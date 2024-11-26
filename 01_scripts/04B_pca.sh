#!/bin/bash

# srun -p small -c 4 -J pca -o log/pca_%j.log /bin/sh 01_scripts/
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


#maybe edit
NB_CPU=4 #change accordingly in SLURM header

##load the adequate python environment that you have created (1st intitialize conda with YOUR PATH to miniconda, except if you already have it in your .bashrc)

##activate a conda environnement in which you have install the necessary library and python version
#conda create --name pcangsd_test python=2.7 
#conda activate pcangsd_test 
#conda install ipython numpy scipy six pandas python-dateutil cython numba

#you may need to edit the name of the environnemnt depending on what you chose
conda activate pcangsd_test 



#this is the input file for the pca
INPUT="$SNP_DIR/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_all_chrs_canon.beagle.gz"

echo "analyse covariance matrix on all individuals"
pcangsd -threads $NB_CPU \
	-beagle $INPUT -o 04_pca/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"

echo "transform covariance matrix into PCA"
COV_MAT=04_pca/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".cov
Rscript 01_scripts/Rscripts/make_pca_simple.r "$COV_MAT" "$BAM_LIST"

