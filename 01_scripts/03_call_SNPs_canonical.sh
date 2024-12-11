#!/bin/bash

### This script will work on all bamfiles and use ANGSD to call SNPs from the previously create list of canonical SNPs,
### and output genotype likelihood (beagle.gz) and MAF (maf.gz)

# less 02_infos/chrs.txt | parallel -j10 srun -p small -c 4 --mem=20G -J 03_call_SNPs_canonical_{} -o log/03_call_SNPs_canonical_{}_%j.log /bin/sh 01_scripts/03_call_SNPs_canonical.sh {} &
# srun -p small -c 4  --mem=20G -J 03_call_SNPs_canonical_NC_036838.1 -o log/03_call_SNPs_canonical_NC_036838.1_%j.log /bin/sh 01_scripts/03_call_SNPs_canonical.sh "NC_036838.1" &

# VARIABLES
GENOME="03_genome/genome.corrected.fasta"
BAM_DIR="04_bam"
SNP_DIR="05_cand_SNPs"
SITES_DIR="02_infos/sites_by_chr"
BAMLIST="02_infos/bam.filelist"
REGION_FILE="02_infos/regions_to_keep.txt"

CHR=$1

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


# LOAD REQUIRED MODULES
module load angsd/0.937
module load samtools/1.15
module load R/4.2

# Increase number of opened file limit
ulimit -S -n 2048


# 1. Call MAF and likelihoods
echo " Calculate the MAF and GL for all individuals listed in 02_info/bam.filelist" on canonical SNP sites
echo "keep loci with at leat one read for n individuals = $MIN_IND, which is $PERCENT_IND % of total $N_IND individuals"
echo "filter on allele frequency = $MIN_MAF"

## Calculate the MAF and GL
angsd -P $NB_CPU -nQueueSize 50 \
-domaf 1 -GL 2 -doGlf 2 -doMajorMinor 1 \
-ref $GENOME \
-sites $SITES_DIR/sites_all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr"$CHR"_canonical \
-remove_bads 1 -minMapQ $MIN_MAPQ -minQ $MIN_Q -skipTriallelic 1 \
-uniqueOnly 1 -only_proper_pairs 1 \
-r "$CHR" \
-b "$BAMLIST" \
-out $SNP_DIR/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr"$CHR"_canon
#-rf $SITES_DIR/"$(basename -s '.txt' $REGIONS)"_"$CHR".txt \

#main features
#-P nb of threads -nQueueSize maximum waiting in memory (necesary to optimize CPU usage
# -doMaf 1 (allele frequencies)  -dosaf (prior for SFS) -GL (Genotype likelihood 2 GATK method - export GL in beagle format  -doGLF2) 
# -doMajorMinor 1 use the most frequent allele as major
# -anc provide a ancestral sequence = reference in our case -fold 1 (car on utilise la ref comme ancestral
# -rf (file with the region written) work on a defined region : OPTIONAL
# -b (bamlist) input file
# -out  output file

