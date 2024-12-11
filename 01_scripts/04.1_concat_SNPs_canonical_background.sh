#!/bin/bash

# srun -p small -c 1 -J 04_concat_SNPs_canonical_background -o log/04_concat_SNPs_canonical_background_%j.log /bin/sh 01_scripts/04.1_concat_SNPs_canonical_background.sh &

# VARIABLES
# Files
GENOME="03_genome/genome.corrected.fasta"
BAMLIST="02_infos/bam.filelist"
CHR_LIST="02_infos/chrs.txt"
REGION_LIST="02_infos/regions_to_keep.txt"

POP_FILE1="02_infos/pop.txt"
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
VALID_DIR="11_SNP_validation"

# Values # BEWARE: different thresholds than in script 02_call_SNPs_all.sh
## SNP calling - data filtering
MIN_MAPQ=10 # min read mapping quality
MIN_Q=20 # minimum base quality score

## SNP calling - sample filtering
MIN_DEPTH=1 # min sequencing depth by sample

## SNP calling - site filtering
MIN_MAF=0.01   #### 
MAX_DEPTH_FACTOR=10 

PERCENT_IND=0.50  ####

# Execs
NGSPARALOG="/project/lbernatchez/users/lalec31/softwares/ngsParalog/ngsParalog"
NGSADMIX="/project/lbernatchez/users/lalec31/softwares/NGSadmix"
REALSFS="/prg/angsd/0.937/misc/realSFS"

NB_CPU=4

# LOAD REQUIRED MODULES
module load angsd/0.937
module load samtools/1.15
module load R/4.2

# Increase number of opened file limit
ulimit -S -n 2048

# BEAGLE
# 1. Extract header for 1st chr : 
## what is the first chromosome ? 
FIRST_CHR=$(less $CHR_LIST | head -n1) 
## Extract header from beagle for first chr and initialize output file
zless $SNP_DIR/background/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr"$FIRST_CHR"_canon.beagle.gz | head -n1 > $SNP_DIR/background/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_all_chrs_canon.beagle

# 2. Append beagle’s contents for all chromosomes
less $CHR_LIST | while read CHR
do
  # extract the right beagle file for a given chr 
	BEAGLE_FILE=$(ls -1 $SNP_DIR/background/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr*_canon.beagle.gz | grep $CHR) # extract the right beagle file for a given chr 
 echo "appending file $BEAGLE_FILE"
	# Extract all lines except first one and append to ALL_CHR.beagle
  zless $BEAGLE_FILE  | grep -v ^marker >> $SNP_DIR/background/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_all_chrs_canon.beagle
done

# 3. Compress
bgzip -f $SNP_DIR/background/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_all_chrs_canon.beagle

# MAFS
# 1. Extract header for 1st chr : 
## Extract header from maf for first chr and initialize output file. We have already identified the 1st chromosome in previous step.
zless $SNP_DIR/background/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr"$FIRST_CHR"_canon.mafs.gz | head -n1 > $SNP_DIR/background/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_all_chrs_canon.mafs

# 2. Append beagles contents for all chromosomes
less $CHR_LIST | while read CHR
do
  # extract the right beagle file for a given chr 
	MAFS_FILE=$(ls -1 $SNP_DIR/background/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr*_canon.mafs.gz | grep $CHR) # extract the right beagle file for a given chr 
	echo "appending file $MAFS_FILE"
  # Extract all lines except first one and append to ALL_CHR.mafs
  zless $MAFS_FILE  | grep -v ^chromo >> $SNP_DIR/background/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_all_chrs_canon.mafs
done

