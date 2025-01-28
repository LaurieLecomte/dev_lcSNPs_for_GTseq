#!/bin/bash

# srun -p small -c 1 --mem=50G -J check_sites -o log/check_sites_%j.log /bin/sh 01_scripts/check_sites.sh "02_infos/sites_by_chr/background/sites_all_maf0.01_pctind0.50_maxdepth10_canonical_minmaj.list" &

# VARIABLES
# Files
GENOME="03_genome/genome.fasta"
BAMLIST="02_infos/bam.filelist"
CHR_LIST="02_infos/chrs.txt"
REGION_LIST="02_infos/regions_to_keep.txt"

POP_FILE1="02_infos/pops.txt"
ID_POP="02_infos/ID_POP.txt"

# Directories
SITES_DIR="02_infos/sites_by_chr"
BAM_DIR="04_bam"
SNP_DIR="05_cand_SNPs"
PCA_DIR="06_pca"
ADMIX_DIR="07_admix"
MAF_DIR="08_maf_by_pop"
FST_DIR="09_fst"
SELECT_DIR="10_SNP_selection"

# Values
## SNP calling - data filtering
MIN_MAPQ=10 # min read mapping quality
MIN_Q=20 # minimum base quality score

## SNP calling - sample filtering
MIN_DEPTH=1 # min sequencing depth by sample

## SNP calling - site filtering
MIN_MAF=0.05 # min allele frequency for a given site (across all samples)
MAX_DEPTH_FACTOR=10 # this value * sample number = max combined depth for a given site (across all samples) -> e.g. 10 * Nsamples gives max allowed combined depth
# advice: we usually set it at 2-4 times the expected coverage to remove repeated regions

PERCENT_IND=0.65 # min fraction of samples required to have min depth at given site
#advice: as percentage, avoid going below 50% and also consider the whole number of individuals. (it may makes sense to use 50% with 100 ind/pop, but you may want 90% with 9 ind/pop

## Fst 
WINDOW=25000  #window size for sliding window FST & Thetas
WINDOW_STEP=5000  #window step
NSITES=500000 #to make realSFS goes faster -reduce the number of sites considered
GROUP=pops

## Admixture
K_MIN=2 #min nb of pop to consider for NGS admix
K_MAX=5 #maximum nb of pop to consider for NGS admix

PVAL_THRESHOLD=0.001

# Execs
NGSPARALOG="/project/lbernatchez/users/lalec31/softwares/ngsParalog/ngsParalog"
NGSADMIX="/project/lbernatchez/users/lalec31/softwares/NGSadmix"
REALSFS="/prg/angsd/0.937/misc/realSFS"


BACKGROUND=$1 # list of less stringent, background SNPs

# LOAD REQUIRED MODULES
module load python/3.7
module load R/4.2


# 4. Score SNPs and determine thresholds of complexity, GC contents, MAF sum and number of neighbor background SNPs (from list of background SNPs produced by scrips 02.1, 03.1 and 04.1)
MIN_AFD=0.1
WIN=100
python3 01_scripts/utils/check_sites.py $SELECT_DIR/maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_combined.afds_"$MIN_AFD".ids $BACKGROUND $GENOME $WIN $SELECT_DIR/maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_combined_"$MIN_AFD"_scored_"$WIN".site_check.txt

GOOD_POS=$(less $SELECT_DIR/maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_combined_"$MIN_AFD"_scored_"$WIN".site_check.txt | grep 'True' | wc -l)
WRONG_POS=$(less $SELECT_DIR/maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_combined_"$MIN_AFD"_scored_"$WIN".site_check.txt | grep 'False' | wc -l)
echo "$GOOD_POS sites with right position"
echo "$WRONG_POS sites with wrong position"
