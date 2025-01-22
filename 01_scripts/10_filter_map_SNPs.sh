#!/bin/bash

# srun -p small -c 5 --mem=50G -J 10_filter_map_SNPs -o log/10_filter_map_SNPs_%j.log /bin/sh 01_scripts/10_filter_map_SNPs.sh &

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

MIN_AFD=0.1
WIN=100
MAX_MAF=2
MAX_NUM_SNP=10
MIN_COMPLEX=100
MIN_GC=0.2
MAX_GC=0.7

#TARGET_SUM=$1
EXP=2
MIN_DIST=200000 # min distance between final SNP to trim final list

# Blast alignment criteria
MIN_IDY=90 #pident
ALN_LEN=160 #length
MAX_HIT=1 #

# Execs
NGSPARALOG="/project/lbernatchez/users/lalec31/softwares/ngsParalog/ngsParalog"
NGSADMIX="/project/lbernatchez/users/lalec31/softwares/NGSadmix"
REALSFS="/prg/angsd/0.937/misc/realSFS"

CPU=5

# LOAD REQUIRED MODULES
module load python/3.7
module load R/4.2

module load ncbiblast/2.6.0
module load bedtools/2.31.1


# 1. Filter SNPs based on established thresholds
Rscript 01_scripts/utils/04.2_filter_SNPs.R $SELECT_DIR/maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_combined_"$MIN_AFD"_scored_"$WIN" $MAX_MAF $MAX_NUM_SNP $MIN_COMPLEX $MIN_GC $MAX_GC

# 2. Blast each SNP's flanking sequence against reference
## Build a fasta of flanking sequences
FASTA="$SELECT_DIR/maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_combined_"$MIN_AFD"_scored_"$WIN".good.fasta"
less $SELECT_DIR/maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_combined_"$MIN_AFD"_scored_"$WIN".good.tsv | tail -n+2 | awk '{print ">"$1 "_" $2 "\n" $14}' > $FASTA

## Fist generate database from reference
if [[ ! -f "$GENOME".nsq ]]
  then makeblastdb -in $GENOME -dbtype nucl 
fi

## Run blast
blastn -query $FASTA -subject $GENOME -out $SELECT_DIR/maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_combined_"$MIN_AFD"_scored_"$WIN".good.blast -outfmt 6 -num_threads $CPU

## Filter blast results: get SNPs that have >1 mapping 
MULTIMAP="$SELECT_DIR/maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_combined_"$MIN_AFD"_scored_"$WIN".good.blast.multimap"
less $SELECT_DIR/maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_combined_"$MIN_AFD"_scored_"$WIN".good.blast | cut -f1 | uniq -d > $MULTIMAP

## Remove multimappings and filter on aln length and identity
less $SELECT_DIR/maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_combined_"$MIN_AFD"_scored_"$WIN".good.blast | grep -vFf $MULTIMAP | awk -v a="$MIN_IDY" -v b="$ALN_LEN" '$3>=a && $4>=b' > $SELECT_DIR/maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_combined_"$MIN_AFD"_scored_"$WIN".good.filtered_idy"$MIN_IDY"_len"$ALN_LEN".blast


## Extract remaining SNPs from scored table
Rscript 01_scripts/utils/extract_scored_SNPs_from_blast.R $SELECT_DIR/maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_combined_"$MIN_AFD"_scored_"$WIN".good.filtered_idy"$MIN_IDY"_len"$ALN_LEN".blast $SELECT_DIR/maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_combined_"$MIN_AFD"_scored_"$WIN".good.tsv $SELECT_DIR/maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_combined_"$MIN_AFD"_scored_"$WIN".good_idy"$MIN_IDY"_len"$ALN_LEN".tsv

