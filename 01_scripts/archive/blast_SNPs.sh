#!/bin/bash

# srun -p small -c 1 --mem=50G -J blast_SNPs -o log/blast_SNPs_%j.log /bin/sh 01_scripts/blast_SNPs.sh &

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

MIN_AFD=0.4
WIN=50
MAX_MAF=2
MAX_NUM_SNP=10
MIN_COMPLEX=60
MIN_GC=0.2
MAX_GC=0.8

TARGET_SUM=12.6
EXP=2
MIN_DIST=10000

# Blast alignment criteria
MIN_IDY=90 #pident
ALN_LEN=160 #length
MAX_HIT=1 #

SNP_LIST="$SELECT_DIR/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_combined_"$MIN_AFD"_scored_"$WIN".good."$TARGET_SUM".tsv"
SNP_SEQS="$SELECT_DIR/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_combined_"$MIN_AFD"_scored_"$WIN""

MULTIHITS="$SELECT_DIR/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_combined_"$MIN_AFD"_scored_"$WIN".good."$TARGET_SUM".multihits"

DIST_SNP=200000

# LOAD REQUIRED MODULES
module load ncbiblast/2.6.0
module load bedtools/2.31.1

#less $SNP_LIST | while read line;
#do
#  CHR=$(echo "$line" | cut -f1)
#  POS=$(echo "$line" | cut -f2)
#  SEQ=$(less $SNP_SEQS | awk -v a="$CHR" -v b="$POS" '$1==a && $2 == b' | cut -f15 | perl -pe "s/[\[\]\ \,\']+//g")
#  echo  '>'"$CHR"_"$POS" >> "$SELECT_DIR/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_combined_"$MIN_AFD"_scored_"$WIN".good."$TARGET_SUM".fasta"
#  echo "$SEQ" >> "$SELECT_DIR/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_combined_"$MIN_AFD"_scored_"$WIN".good."$TARGET_SUM".fasta"
#done

# Make a fasta file for SNP flanking sequences 
#Rscript 01_scripts/utils/flanking_seqs.fasta $SNP_LIST $GENOME $WIN

# Blast flanking seqs against unmasked reference genome
## Fist generate database from unmasked reference
#makeblastdb -in 03_genome/genome.fasta -dbtype nucl 

## Blast
#blastn -query "$SELECT_DIR/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_combined_"$MIN_AFD"_scored_"$WIN".good."$TARGET_SUM".fasta"  -db 03_genome/genome.fasta -out $SELECT_DIR/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_combined_"$MIN_AFD"_scored_"$WIN".good."$TARGET_SUM".blast  -outfmt 6 -perc_identity 0.9 -num_alignments 1 -qcov_hsp_perc 80 

blastn -query "$SELECT_DIR/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_combined_"$MIN_AFD"_scored_"$WIN".good."$TARGET_SUM".fasta" -subject 03_genome/genome.fasta -out $SELECT_DIR/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_combined_"$MIN_AFD"_scored_"$WIN".good."$TARGET_SUM".blast -outfmt 6 #-num_alignments 1 #-perc_identity 0.9 -num_alignments 1 -qcov_hsp_perc 80 

# Filter blast results
## List SNPs that have >1 mapping 
less $SELECT_DIR/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_combined_"$MIN_AFD"_scored_"$WIN".good."$TARGET_SUM".blast | cut -f1 | uniq -d > $MULTIHITS

## Remove multimappings and filter on aln length and identity
less $SELECT_DIR/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_combined_"$MIN_AFD"_scored_"$WIN".good."$TARGET_SUM".blast | grep -vFf $MULTIHITS | awk -v a="$MIN_IDY" -v b="$ALN_LEN" '$3>=a && $4>=b' > $SELECT_DIR/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_combined_"$MIN_AFD"_scored_"$WIN".good."$TARGET_SUM"_filtered.blast


# Trim to one SNP/200 kbp
## Covert blast output to bed 
#less $SELECT_DIR/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_combined_"$MIN_AFD"_scored_"$WIN".good."$TARGET_SUM"_filtered.blast | awk -v OFS="\t" '{ print $2, $9, $9+1, $1}' >  $SELECT_DIR/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_combined_"$MIN_AFD"_scored_"$WIN".good."$TARGET_SUM"_filtered.bed

## Groups SNPs together if they are located less than DIST_SNP bp apart
#bedtools cluster -i $SELECT_DIR/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_combined_"$MIN_AFD"_scored_"$WIN".good."$TARGET_SUM"_filtered.bed -d $DIST_SNP > $SELECT_DIR/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_combined_"$MIN_AFD"_scored_"$WIN".good."$TARGET_SUM"_filtered_clustered.bed 

## Keep only first SNP in each cluster
#less $SELECT_DIR/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_combined_"$MIN_AFD"_scored_"$WIN".good."$TARGET_SUM"_filtered_clustered.bed  | awk '!seen[$5]++' > $SELECT_DIR/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_combined_"$MIN_AFD"_scored_"$WIN".good."$TARGET_SUM"_filtered_trimmed"$(( $DIST_SNP / 1000 ))"kb.bed 