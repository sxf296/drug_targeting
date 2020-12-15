Version 0.9.3 for python2 or python3

Please contact Mike Fang at sxf296@case.edu for any questions.

# Dependencies

numpy, pandas

# Drug Perturbation GSEA (dpGSEA), previously known as dtGSEA

A drug-gene target enrichment technique utilizing a modified GSEA approach, and uses prior drug-defined gene sets in the form of proto-matrix files: PM(protomatrix), L1000 or CMAP(derived from), FC or P(ranked by), top 20 or 50(number of genes).

# How to use

python dpGSEA.py -h for help

The following flags are listed below:

* -tt TOPTABLE, --toptable TOPTABLE (This refers to the TopTable limma output, but one could use any ranked table as long as there are columns "logFC" and "t")
* -dr DRUGREF, --drugref DRUGREF (This refers to the proto-matrix: PM_L1000_FC20.csv, PM_L1000_FC50.csv, etc.)
* -ma, --match (indicating whether to search for matching profiles)
* -i ITERATIONS (we recommend, at least 1000 iterations for generating significance, default is at 1000)
* -sd SETSEED, --setseed SETSEED
* -o OUT, --out OUT (output file names, results will be tab delimited)

Toy example: python dpGSEA.py -tt CD71pos_nonResvsRes.csv -dr PM_L1000_FC50.csv -i 1000 -o ResVsnonResResults.tsv

# Output file

* drug - specific drug
* ES - enrichment score
* NES - normalized enrichment score
* ES_p - enrichment score p value
* TCS - target compatibility score
* NTCS - normalized target compatibility score
* TCS_p - target compatibility score p value
* genes - leading edge genes
* NTCS_num - FDR cutoff for NTCS, denoted with 1 if drug reaches sig. threshold at default 0.90 and 0.95 confidence levels.
* NES_num - FDR cutoff for NES

Please note, to add different FDR cutoffs, please ctrl+F "quantiles = [90, 95]" within the script and add confidence levels as needed.

# Results

* Enrichment score (ES) - this score is interpreted the same way the standard GSEA enrichment score. It reflects the degree to which a complimentary or matching drug gene profile is overrepresented at the top of a ranked list.

* Enrichment score p-value (ES_pvalue) - the statistical significance of the enrichment score for a single drug gene set.

* Target compatibility p-value (TC_pvalue) - a p-value reflecting the quantity and magnitude of statistical significance of differentially expressed genes that match or antagonize a drug profile. This statistical test compares the modulation of the leading edge genes against random modulation.

* Driver Genes aka leading edge genes (driver_genes) - this lists genes that appear in the ranked list at or before the point at which the running sum reaches its maximum deviation from zero. These genes are often interpreted as the genes driving an enrichment or modulation of drug-gene and differential expression analysis.
