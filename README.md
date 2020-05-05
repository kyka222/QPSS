# QPSS

## Introduction

Quantitative Phenotype Scan Statistic (QPSS) is a tool for scan-statistic approach handling continuous outcomes to identify genomic regions where rare quantitative-phenotype-associated variants cluster. 

## Requirements
Programming language: Python2\
Others: PLINK 1.9, R 2.10 or higher, R packages goft\
\
## How to use
Options:\
--bfile = input PLINK binary fileset (i.e., .bed/bim/fam) (required)\
--chr = chromosome number (required)\
--G-position = start and end positions of G (default: start and end positions in PLINK bfile)\
--pheno = phenotype file with FID (first column), IID (second column)\
--pheno-name = column name for phenotype in phenotype file\
--max-maf - max threshold of MAF to define rare variant (default: 0.05)\
--min-ac = min threshold of MAC to define rare variant (default: 5)\
--extract = list to extract variants\
--exclude = list to exclude variants\
--keep = list to keep individuals\
--remove = list to remove individuals
--W-position or --W = start and end positions of W\
--W-file = list of start and end positions\
--W-fixed or --wf = window size for sliding window approach (default: 2000)\
--W-slide or --ws = slid window size for sliding window approach (default: half of window size specified by --W-fixed)\
--out" = output file name\
--perm or --p = whether p is computed. Add 'gpd' as --perm gpd for GPD approximation\
--threads = number of threads used (default: 10)\
--max-sim = max number of simulation for permutation test (default: 10000)\


Example 1: sliding window approach (W = 2000/S = 1000) for entire chromosome 1 for outcome in a phenotype file to output logLR only\
\
python2 QPSS.v1.py --bfile plink_fileset_name --pheno phenotype_file_name --phenoname outcome_name --chr 1 --out output_name\
\
\
Example 2: Calculate logLR and p-value by permutation test for a specfic window (start is 1234 and end is 5678) on chromosome 1\
\
python2 QPSS.v1.py --bfile plink_fileset_name --pheno phenotype_file_name --phenoname outcome_name --chr 1 --W-position 1234 5678\
\
\
Example 3: Calculate logLR and p-value with GPD approximation for a specfic window (start is 1234 and end is 5678) within a specfic G region (start 123 to end 456789) on chromosome 1 --perm\
\
python2 QPSS.v1.py --bfile plink_fileset_name --pheno phenotype_file_name --phenoname outcome_name --chr 1 --W-position 1234 5678 --G-position 123 456789 --perm gpd\
python2 QPSS.v1.py --bfile plink_fileset_name --pheno phenotype_file_name --phenoname outcome_name --chr 1 --W-position 1234 5678 --G-position 123 456789 --perm gpd\SSSSSS
\