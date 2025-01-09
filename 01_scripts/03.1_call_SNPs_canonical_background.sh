#!/bin/bash

### This script will work on all bamfiles and use ANGSD to call SNPs from the previously create list of canonical SNPs,
### and output genotype likelihood (beagle.gz) and MAF (maf.gz)

# less 02_infos/chrs.txt | parallel -j10 srun -p small -c 4 --mem=20G -J 03.1_call_SNPs_canonical_background_{} -o log/03.1_call_SNPs_canonical_background_{}_%j.log /bin/sh 01_scripts/03.1_call_SNPs_canonical_background.sh {} &
# srun -p small -c 4  --mem=20G -J 03.1_call_SNPs_canonical_background_NC_036838.1 -o log/03.1_call_SNPs_canonical_background_NC_036838.1_%j.log /bin/sh 01_scripts/03.1_call_SNPs_canonical_background.sh "NC_036838.1" &

# VARIABLES
# Files
GENOME="03_genome/genome.fasta"
BAMLIST="02_infos/bam.filelist"
CHR_LIST="02_infos/chrs.txt"
REGION_LIST="02_infos/regions_to_keep.txt"
REGION_BED="02_infos/regions_to_keep.bed"

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
CHR=$1

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
-sites $SITES_DIR/background/sites_all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr"$CHR"_canonical \
-remove_bads 1 -minMapQ $MIN_MAPQ -minQ $MIN_Q -skipTriallelic 1 \
-uniqueOnly 1 -only_proper_pairs 1 \
-r "$CHR:" \
-b "$BAMLIST" \
-out $SNP_DIR/background/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr"$CHR"_canon

#main features
#-P nb of threads -nQueueSize maximum waiting in memory (necesary to optimize CPU usage
# -doMaf 1 (allele frequencies)  -dosaf (prior for SFS) -GL (Genotype likelihood 2 GATK method - export GL in beagle format  -doGLF2) 
# -doMajorMinor 1 use the most frequent allele as major
# -anc provide a ancestral sequence = reference in our case -fold 1 (car on utilise la ref comme ancestral
# -rf (file with the region written) work on a defined region : OPTIONAL
# -b (bamlist) input file
# -out  output file

