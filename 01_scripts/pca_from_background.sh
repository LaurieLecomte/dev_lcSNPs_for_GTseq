#!/bin/bash

# srun -p small -c 4 -J pca_from_background -o log/pca_from_background_%j.log /bin/sh 01_scripts/pca_from_background.sh & 

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


##activate a conda environnement in which you have install the necessary library and python version
#conda create --name pcangsd_test python=2.7 
#conda activate pcangsd_test 
#conda install ipython numpy scipy six pandas python-dateutil cython numba

#you may need to edit the name of the environnemnt depending on what you chose
#conda activate pcangsd_test 

# LOAD REQUIRED MODULES
module load pcangsd/1.10
module load R/4.2


#this is the input file for the pca
INPUT="$SNP_DIR/background/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_all_chrs_canon.beagle.gz"

echo "analyse covariance matrix on all individuals"
pcangsd --threads $NB_CPU -b $INPUT -o $PCA_DIR/background/"$(basename -s '.beagle.gz' $INPUT)"

echo "transform covariance matrix into PCA"
COV_MAT=$PCA_DIR/background/"$(basename -s '.beagle.gz' $INPUT)".cov
#Rscript 01_scripts/Rscripts/make_pca_simple.r "$COV_MAT" "$BAMLIST"
Rscript 01_scripts/utils/pca_simple.R "$COV_MAT" "$BAMLIST" "$INPUT" $ID_POP

