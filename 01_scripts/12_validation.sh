#!/bin/bash

# srun -p small -c 1 --mem=50G -J 12_validation -o log/12_validation_%j.log /bin/sh 01_scripts/12_validation.sh 12.6 &

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
VALID_DIR="11_SNP_validation"

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

TARGET_SUM=$1
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

CPU=1

# LOAD REQUIRED MODULES
module load python/3.7
module load R/4.2
module load pcangsd/1.10


# 1. Run pca on final set of SNPs
INPUT="$SELECT_DIR/maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_combined_"$MIN_AFD"_scored_"$WIN".good_idy"$MIN_IDY"_len"$ALN_LEN"."$TARGET_SUM".beagle.gz"
echo "analyse covariance matrix on all individuals"
pcangsd --threads $CPU -b $INPUT -o $VALID_DIR/"$(basename -s '.beagle.gz' $INPUT)"

echo "transform covariance matrix into PCA"
COV_MAT=$VALID_DIR/"$(basename -s '.beagle.gz' $INPUT)".cov
Rscript 01_scripts/utils/pca_simple.R "$COV_MAT" "$BAMLIST" "$INPUT" $ID_POP


# 2. Run NGSadmix on final set of SNPs
# Loop from K_MIN to K_MAX
for i in $(seq $K_MIN $K_MAX)
do 
	echo $i
	$NGSADMIX -P $CPU -likes $INPUT -minMaf $MIN_MAF -K $i -o $VALID_DIR/"$(basename -s '.beagle.gz' $INPUT)"_K"$i"
  Rscript 01_scripts/utils/plot_admix.R $VALID_DIR/"$(basename -s '.beagle.gz' $INPUT)"_K"$i" $ID_POP
done



# 3. Run pca and NGSadmix on subsets of final SNP set

if [[ ! -d $VALID_DIR/subsets ]]
then
  mkdir $VALID_DIR/subsets
fi

# Generate subsets 
Rscript 01_scripts/utils/subsets_from_final_beagle.R $INPUT 10 100 $VALID_DIR/subsets

## bgzip
for file in $(ls -1 $VALID_DIR/subsets/*subset*sites.beagle); do bgzip $file -f ; done

## Loop over these subset files
for file in $(ls -1 $VALID_DIR/subsets/*subset*sites.beagle.gz);
do

  echo $file
  # Run pca
  echo "analyse covariance matrix on all individuals"
  pcangsd --threads $CPU -b $file -o $VALID_DIR/subsets/"$(basename -s '.beagle.gz' $file)"
  
  echo "transform covariance matrix into PCA"
  COV_MAT=$VALID_DIR/subsets/"$(basename -s '.beagle.gz' $file)".cov
  Rscript 01_scripts/utils/pca_simple.R "$COV_MAT" "$BAMLIST" "$file" $ID_POP
  
  # Run NGSadmix
  for i in $(seq $K_MIN $K_MAX); 
  do 
    $NGSADMIX -P $CPU -likes $file -minMaf $MIN_MAF -K $i -o $VALID_DIR/subsets/"$(basename -s '.beagle.gz' $file)"_K"$i"
    Rscript 01_scripts/utils/plot_admix.R $VALID_DIR/subsets/"$(basename -s '.beagle.gz' $file)"_K"$i" $ID_POP
  done
  
  #echo "done for $file"
done
