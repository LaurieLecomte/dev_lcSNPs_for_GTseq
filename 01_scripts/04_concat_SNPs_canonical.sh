#!/bin/bash

# srun -p small -c 1 -J 04_concat_SNPs_canonical -o log/04_concat_SNPs_canonical_%j.log /bin/sh 01_scripts/04_concat_SNPs_canonical.sh &

# VARIABLES
GENOME="03_genome/genome.fasta"
BAM_DIR="04_bam"
SNP_DIR="05_cand_SNPs"
SITES_DIR="02_infos/sites_by_chr"
BAMLIST="02_infos/bam.filelist"

#REGION=$1

# Filtering
MIN_MAF=0.05 # will keep SNP above this allele frequency (over all individuals)
MIN_DEPTH=1 # will keep positions with at least MIN_DEPTH reads for each individual, This is not necessaily for all individuals, we consider a PERCENT_IND (percentage of individuals over all individuals in step 03, and within each pop at step 07)
#advice: as min depth use a value that is a bit below what you expected. we use 1 for 1X of coverage but if you have 5X it may make sense to put the bar a bit higher to 2 or 3
PERCENT_IND=0.65 #advice: as percentage, avoid going below 50% and also consider the whole number of individuals. (it may makes sense to use 50% with 100 ind/pop, but you may want 90% with 9 ind/pop
MAX_DEPTH_FACTOR=10 #will keep SNP with at least a coverage of this factor multiplied by the number of ind - across all ind. 
# advice: we usually set it at 2-4 times the expected coverage to remove repeated regions
MIN_MAPQ=10
MIN_Q=20

# Window size for Fst a
WINDOW=25000  #window size for sliding window FST & Thetas
WINDOW_STEP=5000  #window step

# NGSadmix
K_MIN=2 #min nb of pop to consider for NGS admix
K_MAX=5 #maximum nb of pop to consider for NGS admix


NB_CPU=4 #change accordingly to the -c argument in srun, 

PVAL_THRESHOLD=0.001
#REGION_NUM="$1" 
#REGION=$(head -n $REGION_NUM 02_info/regions.txt | tail -n 1)

# Execs
NGSPARALOG="/project/lbernatchez/users/lalec31/softwares/ngsParalog/ngsParalog"

CHR_LIST="02_infos/chrs.txt"

# LOAD REQUIRED MODULES
module load angsd/0.931
module load samtools/1.15
module load R/4.2

# Increase number of opened file limit
ulimit -S -n 2048


# BEAGLE
# 1. Extract header for 1st chr : 
## what is the first chromosome ? 
FIRST_CHR=$(less $CHR_LIST | head -n1) 
## Extract header from beagle for first chr and initialize output file
zless $SNP_DIR/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr"$FIRST_CHR"_canon.beagle.gz | head -n1 > $SNP_DIR/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_all_chrs_canon.beagle

# 2. Append beagle’s contents for all chromosomes
less $CHR_LIST | while read CHR
do
  # extract the right beagle file for a given chr 
	BEAGLE_FILE=$(ls -1 $SNP_DIR/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr*_canon.beagle.gz | grep $CHR) # extract the right beagle file for a given chr 
 echo "appending file $BEAGLE_FILE"
	# Extract all lines except first one and append to ALL_CHR.beagle
  zless $BEAGLE_FILE  | grep -v ^marker >> $SNP_DIR/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_all_chrs_canon.beagle
done

# 3. Compress
bgzip -f $SNP_DIR/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_all_chrs_canon.beagle

# MAFS
# 1. Extract header for 1st chr : 
## Extract header from maf for first chr and initialize output file. We have already identified the 1st chromosome in previous step.
zless $SNP_DIR/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr"$FIRST_CHR"_canon.mafs.gz | head -n1 > $SNP_DIR/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_all_chrs_canon.mafs

# 2. Append beagles contents for all chromosomes
less $CHR_LIST | while read CHR
do
  # extract the right beagle file for a given chr 
	MAFS_FILE=$(ls -1 $SNP_DIR/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr*_canon.mafs.gz | grep $CHR) # extract the right beagle file for a given chr 
	echo "appending file $MAFS_FILE"
  # Extract all lines except first one and append to ALL_CHR.mafs
  zless $MAFS_FILE  | grep -v ^chromo >> $SNP_DIR/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_all_chrs_canon.mafs
done

