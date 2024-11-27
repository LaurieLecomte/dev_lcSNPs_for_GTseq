#!/bin/bash

# srun -p small -c 4 -J 06_ngs_admix -o log/06_ngs_admix_%j.log /bin/sh 01_scripts/06_ngs_admix.sh &

###this script will work on all individuals using the beagle genotype likelihood and perform an admixture analysis 
#this requires NGSadmix to be installed and its path export in the bashrc
#NGS admiw will explore all number of population between K_MIN and K_MAX as provided in the 01_config.sh

# VARIABLES
GENOME="03_genome/genome.fasta"
BAM_DIR="04_bam"
SNP_DIR="05_cand_SNPs"
SITES_DIR="02_infos/sites_by_chr"
BAMLIST="02_infos/bam.filelist"

PCA_DIR="06_pca"
ADMIX_DIR="07_admix"

CHR_LIST="02_infos/chrs.txt"

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
NGSADMIX="/project/lbernatchez/users/lalec31/softwares/NGSadmix"

# LOAD REQUIRED MODULES
module load angsd/0.937
#module load samtools/1.15
#module load R/4.2

#INPUT="$SNP_DIR/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_all_chrs_canon.beagle.gz"

# Loop from K_MIN to K_MAX
for i in $(seq $K_MIN $K_MAX)
do 
	echo $i
	$NGSADMIX -P $NB_CPU -likes $SNP_DIR/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_all_chrs_canon.beagle.gz \
	-minMaf $MIN_MAF -K $i -o $ADMIX_DIR/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_K"$i"
done
