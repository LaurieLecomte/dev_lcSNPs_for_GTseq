#!/bin/bash

# srun -p medium -c 1 --mem=50G -J 09_preselect_SNPs -o log/09_preselect_SNPs_%j.log /bin/sh 01_scripts/09_preselect_SNPs.sh &

# VARIABLES
# Files
GENOME="03_genome/genome.corrected.fasta"
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
GROUP=pop

## Admixture
K_MIN=2 #min nb of pop to consider for NGS admix
K_MAX=5 #maximum nb of pop to consider for NGS admix

PVAL_THRESHOLD=0.001

# Execs
NGSPARALOG="/project/lbernatchez/users/lalec31/softwares/ngsParalog/ngsParalog"
NGSADMIX="/project/lbernatchez/users/lalec31/softwares/NGSadmix"
REALSFS="/prg/angsd/0.937/misc/realSFS"


# LOAD REQUIRED MODULES
module load python/3.7
module load R/4.2

# 1. Combine per-site mafs across all populations into a single file
Rscript 01_scripts/utils/combine_mafs.R 02_infos/sites_all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_canonical_minmaj.list $POP_FILE1 $MAF_DIR/$POP/"$POP"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".mafs

# 2. Calculate pairwise AFDs (Allele Frequency Differences)
python3 01_scripts/utils/01_compute_pairwise_AFDs.py $MAF_DIR/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_combined.mafs $SELECT_DIR/maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_combined.afds

# 3. Prefilter SNPs on AFDs
MIN_AFD=0.1 # decr to 0.1
python3 01_scripts/utils/02_pre_filter_SNPs_on_pairwise_AFDs.py $SELECT_DIR/maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_combined.afds $MIN_AFD $SELECT_DIR/maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_combined.afds_"$MIN_AFD".ids

# 4. Score SNPs and determine thresholds of complexity, GC contents, MAF sum and number of neighbor background SNPs (from list of background SNPs produced by scrips 02.1, 03.1 and 04.1)
WIN=100
python3 01_scripts/utils/03_score_SNPs_for_GTseq.py $SELECT_DIR/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_combined.afds_"$MIN_AFD".ids $SITES_DIR/background/sites_all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_canonical_minmaj.list $GENOME $WIN $SELECT_DIR/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_combined_"$MIN_AFD"_scored_"$WIN"

# 4. Plot metrics to decide thresholds
Rscript 01_scripts/utils/04.1_plot_metrics.R $SELECT_DIR/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_combined_"$MIN_AFD"_scored_"$WIN"
