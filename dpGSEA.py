'''
dpGSEA 0.9.4

A drug-gene target enrichment technique utilizing a modified GSEA approach. 

dpGSEA enriches on a proto-matrix based on a user-defined cutoff (these matrices need to be built by the user or the user can utilize the ones included with the script). This proto-matrix summarizes drug-gene perturbation profiles in a format iterable by this enrichment. The output includes three different measures along with the "leading-edge" genes as a .csv output file.

- Enrichment score (ES) - this score is interpreted the same way the standard GSEA enrichment score. It reflects the degree to which a complimentary or matching drug gene profile is overrepresented at the top of a ranked list.
- Enrichment score p-value (ES_pvalue) - the statistical significance of the enrichment score for a single drug gene set.
- Target compatibility p-value (TC_pvalue) - a p-value reflecting the quantity and magnitude of statistical significance of differentially expressed genes that 
  match or antagonize a drug profile. This statistical test compares the modulation of the leading edge genes against random modulation.
- Driver Genes aka leading edge genes (driver_genes) - this lists genes that appear in the ranked list at or before the point at which the running sum reaches its maximum deviation from zero. These genes are often interpreted as the genes driving an enrichment or modulation of drug-gene and differential expression analysis.

'''


# Importing the depencies
# This script is dependant on indexing provided by pandas
import pandas as pd

# Numpy has slightly faster list parsing versus native python
import numpy as np

# Random is used to generate null distributions for hypothesis testing
import random as ran

# For cmdline arguments
import argparse

'''

This script contains 3 classes. The first class (inputArgs) parses the arguments and reads the tables using pandas. The second class (tablePreprocessing) merges and ranks the protolist against the differential expression topTable. The last class (dpGSEA) performs the enrichment. The code ladder distributes the computing.

'''



# This class uses the argparse library to set flags for arguments from commandline and reads the tables using pandas
class inputArgs:
	
	# Initialize flags for the specific arguments needed in cmdline
	def __init__(self):
		
		# Add flags needed here
		self.helper = {'tt'  : 'Table (output from limma topTable) of differential expression analysis (comma delimited)',
                       'dr'  : 'User generated drug reference table',
                       'ma'  : 'Find matching drug profiles rather than compliment',
                       'i'   : 'Number of permutations for null',
                       'sd'  : 'Set seed for permutations (default 1)',
                       'o'   : 'Output file names'}
		
		# Initialize argument parser
		self.ap = argparse.ArgumentParser()
		
		# Add flags above to instance of argument parsing
		self.ap.add_argument('-tt',  '--toptable', help = self.helper['tt'], type = str)
		self.ap.add_argument('-dr',  '--drugref',  help = self.helper['dr'], default = 'contact.txt', type = str)
		self.ap.add_argument('-ma',  '--match',    help = self.helper['ma'],       default = False, action='store_true')
		self.ap.add_argument('-i',   '--iterations',      help = self.helper['i'],  default = 1000, type = int)
		self.ap.add_argument('-sd',  '--setseed',  help = self.helper['sd'],       default = 1, type = int)
		self.ap.add_argument('-o',   '--out',      help = self.helper['o'],  type = str)
		
		# Pass all flags to parser instance
		self.args = self.ap.parse_args()
		
	# Read the topTable and renames to geneTable
	def readGeneTable(self):
		
		# Comma delimited file, use write.csv with no arguments in R to generate the topTable
		return pd.read_csv(self.args.toptable, sep = ',')
		
	# Read the proto matrix table and renames it to drugRef
	def readDrugRef(self):
		
		# Comma delimited file with the column names 'drug', 'gene', 'drug_dir'
		return pd.read_csv(self.args.drugref, sep = ',')



# This class parses and adds columns needed for enrichment
class tablePreprocessing:
	
	# The parsing and merging of the tables take place during initiliazation
	def __init__(self, geneTable, drugRef):
		
		# Calculates absolute value of the T-statistic used for ranking and renames some columns
		geneTable['abs_t'] = abs(geneTable['t'])
		geneTable.columns  = ['gene'] + list(geneTable.columns[1:])
		drugRef.columns    = ['drug', 'gene', 'drug_dir']
		
		# Merges the columns on genes based on the drugRef, remove any NAs generated from merge and ranks by abs_t
		rankTable = pd.merge(drugRef, geneTable, on = 'gene', how = 'left')
		rankTable = rankTable[rankTable.t.notna()]
		rankTable = rankTable.sort_values(by = 'abs_t', ascending = False).reset_index()
		
		# Determines the gene direction based on the FC of the DE, 1 is up-regulated while 0 is down-regulated
		rankTable.loc[rankTable.logFC > 0, 'gene_dir'] = 1
		rankTable.loc[rankTable.logFC < 0, 'gene_dir'] = 0
		
		# Remove signal weaker than can be represented by float
		rankTable = rankTable[rankTable.t != 0]
		
		# Assign as ints for faster comp later on
		rankTable.drug_dir = rankTable.drug_dir.astype(int)
		rankTable.gene_dir = rankTable.gene_dir.astype(int)
		rankTable['abs_t'] = rankTable.abs_t.round(6)
		
		# Add to instance
		self.rankTable = rankTable
		
	# Partition the drug list for multi-processing
	def partitionDrugList(self, numChunks):
		
		# Read rankTable from instance
		rankTable = self.rankTable
		
		# Get all unique drugs from the rankTable
		drugList = list(set(rankTable.drug))
		
		# Return table as chunks based on the number of processes
		return list(drugList[i::numChunks] for i in range(numChunks))



# This class performs the enrichment itself, multiple instances of the class is called when performing multi-processing
class dpGSEA:
	
	def __init__(self, rankTable, iterations, seed, matchProfile = False):
		ran.seed(seed)
		self.rankTable    = rankTable
		self.indexLen     = len(self.rankTable)
		self.iterations   = iterations		
		self.matchProfile = matchProfile
		
	def drugList(self):
		return self.rankTable['drug'].unique()
		
	def getDrugIndexes(self, drug):
		rankTable    = self.rankTable
		matchProfile = self.matchProfile
		
		if matchProfile: ind = rankTable[(rankTable.gene_dir == rankTable.drug_dir) & (rankTable.drug == drug)].index
		else: ind = rankTable[(rankTable.gene_dir != rankTable.drug_dir) & (rankTable.drug == drug)].index
		
		if ind.size != 0: return np.asarray(ind)
		else: return None
		
	def getNullIndexes(self, drug):
		iterations  = self.iterations
		try:
			resampleNum = len(self.getDrugIndexes(drug))
			if resampleNum != 0: return np.array([np.random.choice(self.indexLen, resampleNum, replace = False) for _ in range(iterations)])
			else: return None
		except:
			return None
		
	def getMaxDeviations(self, index, getTable = False):
		if index is not None:
			if len(index.shape) == 1:
				# Assigns variable to instance rank table
				rankTable = self.rankTable
				
				# Finds total sum of for brownian bridge
				totalSum  = sum(rankTable.abs_t)
				
				# Calculates the total sum for hits
				hitSum   = sum(rankTable.abs_t[index])
				
				# Negative step for "misses" weighted by the T-statistic
				rankTable['step'] = -1 * rankTable.abs_t/(totalSum - hitSum)
				
				# Calculates the "hit" steps (the comprehension loop will save time on smaller group sizes)
				# rankTable.ix[index, 'step'] = rankTable.ix[index].abs_t/hitSum
				rankTable.loc[index, 'step'] = [rankTable.abs_t[n]/hitSum for n in index]
				
				# Calculates the cumulative sum for the brownian bridge
				rankTable['cumsum'] = np.cumsum(rankTable.step)
				
				# Calculates cumulative sum and finds max deviation and index
				maxDeviation      = max(rankTable['cumsum'])
				maxDeviationIndex = float(rankTable['cumsum'].idxmax())
				
				if getTable:
					return rankTable
					
				else:
					return {'maxDeviation': maxDeviation,
                            'maxDeviationIndex': 1 - (maxDeviationIndex/self.indexLen)}
					
			else:
				# Assigns variable to instance rank table
				rankTable  = self.rankTable
				iterations = self.iterations
				# Iterate through all indexes
				totalSum  = sum(rankTable.abs_t)
				
				maxDeviationList      = np.array([])
				maxDeviationIndexList = np.array([])
				
				for i in index:
					
					# Calculates the total sum for hits
					# hitSum   = sum(rankTable.abs_t[i])
					hitSum   = sum(rankTable.abs_t[n] for n in i)
					
					# Negative step for "misses" weighted by the T-statistic
					rankTable['step'] = -1 * rankTable.abs_t/(totalSum - hitSum)
					
					# Calculates the "hit" steps (the comprehension loop will save time on smaller group sizes)
					# rankTable.ix[i, 'step'] = rankTable.ix[i].abs_t/hitSum
					rankTable.loc[i, 'step'] = [rankTable.abs_t[n]/hitSum for n in i] # faster for shorter <200 lists 
					
					# Calculates cumulative sum and finds max deviation and index
					cumSum            = np.cumsum(rankTable.step)
					maxDeviation      = max(cumSum)
					maxDeviationIndex = float(cumSum.idxmax())
					
					# Adds to the max deviations and index of max deviations arrays
					maxDeviationList      = np.append(maxDeviationList, maxDeviation)
					maxDeviationIndexList = np.append(maxDeviationIndexList, maxDeviationIndex)
					
				maxDeviationListNorm      = maxDeviationList/np.mean(maxDeviationList)
				maxDeviationIndexList     = 1 - (maxDeviationIndexList/self.indexLen)
				maxDeviationIndexListNorm = maxDeviationIndexList/np.mean(maxDeviationIndexList)
					
				return {'maxDeviation'          : maxDeviationList,
                        'maxDeviationNorm'      : maxDeviationListNorm,
                        'maxDeviationIndex'     : maxDeviationIndexList,
                        'maxDeviationIndexNorm' : maxDeviationIndexListNorm}
				
	def getStats(self, drug, drugIndex = None, drugMaxDev = None, nullIndex = None, nullMaxDev = None):
		if drugIndex is None and drugMaxDev is None:
			drugIndex  = self.getDrugIndexes(drug)
			drugMaxDev = self.getMaxDeviations(drugIndex)
			
		if nullIndex is None and nullMaxDev is None:
			nullIndex  = self.getNullIndexes(drug)
			nullMaxDev = self.getMaxDeviations(nullIndex)
		
		iterations = self.iterations + 0.
		rankTable  = self.rankTable
		
		# Sets the enrichment score, enrichment score p, and target compatibility score
		es  = drugMaxDev['maxDeviation']
		nes = es/np.mean(nullMaxDev['maxDeviation'])
		esp = sum(nullMaxDev['maxDeviation'] > drugMaxDev['maxDeviation'])/iterations
		
		tcs  = drugMaxDev['maxDeviationIndex']
		ntcs = tcs/np.mean(nullMaxDev['maxDeviationIndex'])
		tcp  = sum(nullMaxDev['maxDeviationIndex'] > drugMaxDev['maxDeviationIndex'])/iterations
		
		# Finds leading edge genes
		driverGeneIndexes = drugIndex[drugIndex <= (-(drugMaxDev['maxDeviationIndex'] - 1) * self.indexLen) + 0.5]
		genes = list(rankTable.loc[driverGeneIndexes, 'gene'])
		
		# Returns dict of results
		return {'drug'  : drug,
                'ES'    : es,
                'NES'   : nes,
                'ES_p'  : esp,
                'TCS'   : tcs,
                'NTCS'  : ntcs,
                'TCS_p' : tcp,
                'genes' : genes}
		
	def resultsTable(self):
		drugList    = self.drugList()
		drugListLen = len(drugList)
		iterations  = self.iterations
		
		resultsTable = pd.DataFrame(columns = ['drug', 'ES', 'NES', 'ES_p', 'TCS', 'NTCS', 'TCS_p', 'genes'])
		nullDistDict = {}
		nullNESDist  = []
		nullNTCSDist  = []
		
		drugCounter = 1
		for drug in drugList:
			print('DxCL: ' + str(drugCounter) + ' of ' + str(drugListLen) + ', ' + drug)
			drugCounter += 1
			drugIndex  = self.getDrugIndexes(drug)
			
			if drugIndex is not None:
				drugMaxDev = self.getMaxDeviations(drugIndex)
				nullKey    = len(drugIndex)
				
				if nullKey in nullDistDict.keys():
					nullMaxDev = nullDistDict[nullKey]
					
				else:
					nullIndex  = self.getNullIndexes(drug)
					nullMaxDev = self.getMaxDeviations(nullIndex)
					nullDistDict[nullKey] = nullMaxDev
					
				drugStats    = self.getStats(drug, drugIndex = drugIndex, drugMaxDev = drugMaxDev, nullIndex = nullIndex, nullMaxDev = nullMaxDev)
				resultsTable = resultsTable.append(drugStats, ignore_index = True)
				nullNESDist  = nullNESDist + list(nullMaxDev['maxDeviationNorm'])
				nullNTCSDist = nullNTCSDist + list(nullMaxDev['maxDeviationIndexNorm'])
				
		print('Calculating FDRs...')
		
		quantiles = [90, 95]
		for u in quantiles:
			nes_threshold  = np.percentile(nullNESDist, u)
			NES_FDR = resultsTable[resultsTable.NES >= nes_threshold].index
			col_name = 'NES_' + str(u)
			resultsTable[col_name] = 0
			resultsTable.loc[NES_FDR, col_name] = 1
			
			ntcs_threshold = np.percentile(nullNTCSDist, u)
			NTCS_FDR = resultsTable[resultsTable.NTCS >= ntcs_threshold].index
			col_name = 'NTCS_' + str(u)
			resultsTable[col_name] = 0
			resultsTable.loc[NTCS_FDR, col_name] = 1
			
		# resultsTable = resultsTable[['drug', 'ES', 'NES', 'ES_p', 'NES_FDR', 'TCS', 'NTCS', 'TCS_p', 'NTCS_FDR', 'genes']]
		# resultsTable = resultsTable[['drug', 'ES', 'NES', 'ES_p', 'TCS', 'NTCS', 'TCS_p', 'genes']]
		
		return resultsTable
		
# script control ###################################################################################################

if __name__ == "__main__":
	
	# Initialize the argument parsing library
	inp = inputArgs()
	
	# Read the topTable and the proto matrix
	geneTable      = inp.readGeneTable()
	drugRef        = inp.readDrugRef()
	matchProfile   = inp.args.match
	iterations     = inp.args.iterations
	seed           = inp.args.setseed
	outputFileName = inp.args.out
	
	# Perform the processing and merging of the tables to create the rank table
	tp             = tablePreprocessing(geneTable, drugRef)
	rankTable      = tp.rankTable
	
	dt      = dpGSEA(rankTable, iterations = iterations, seed = seed, matchProfile = matchProfile)
	results = dt.resultsTable()
	
	print('Writing results...')
	results.to_csv(path_or_buf = outputFileName, index = False, sep = '\t')


