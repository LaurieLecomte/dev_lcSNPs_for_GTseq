#!/bin/bash

#Produce a list of background SNPs to filter candidate SNPs called using more stringent filters 

# less 02_infos/chrs.txt | parallel -j10 srun -p small -c 4 --mem=20G -J 02.1_call_SNPs_all_background_{} -o log/02.1_call_SNPs_all_background_{}_%j.log /bin/sh 01_scripts/02.1_call_SNPs_background.sh {} &
# srun -p small -c 4 --mem=20G -J 02.1_call_SNPs_all_background_NC_036838.1 -o log/02.1_call_SNPs_all_NC_036838.1_background_%j.log /bin/sh 01_scripts/02.1_call_SNPs_background.sh "NC_036838.1" &

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
CHR=$1

# LOAD REQUIRED MODULES
module load angsd/0.937
module load samtools/1.15
module load R/4.2

# Increase number of opened file limit
ulimit -S -n 2048


# Create a new directory for background SNPs
if [[ ! -d $SNP_DIR/background ]]
then
  mkdir $SNP_DIR/background
fi

if [[ ! -d $SITES_DIR/background ]]
then
  mkdir $SITES_DIR/background
fi


# 0. Prepare variables 
N_IND=$(wc -l $BAMLIST | cut -d " " -f 1)
MIN_IND_FLOAT=$(echo "($N_IND * $PERCENT_IND)" | bc -l)
MIN_IND=${MIN_IND_FLOAT%.*} 
MAX_DEPTH=$(echo "($N_IND * $MAX_DEPTH_FACTOR)" | bc -l)

# 0. Prepare regions
less $REGION_LIST | grep $CHR > $SITES_DIR/"$(basename -s '.txt' $REGION_LIST)"_"$CHR".txt # already done in main script


# 1. Calculate the MAF and HWE
echo "Calculate the MAF and HWE for all individuals listed in 02_info/bam.filelist"
echo "keep loci with at least $MIN_DEPTH read for n individuals = $MIN_IND, which is $PERCENT_IND % of total $N_IND individuals"
echo "filter on allele frequency = $MIN_MAF"

angsd -P $NB_CPU -nQueueSize 50 \
-doMaf 1 -doHWE 1 -GL 2 -doMajorMinor 1 -doCounts 1 \
-remove_bads 1 -minMapQ $MIN_MAPQ -minQ $MIN_Q -skipTriallelic 1 \
-uniqueOnly 1 -only_proper_pairs 1 \
-minInd $MIN_IND -minMaf $MIN_MAF -setMaxDepth $MAX_DEPTH -setMinDepthInd $MIN_DEPTH \
-b $BAMLIST \
-rf $SITES_DIR/"$(basename -s '.txt' $REGION_LIST)"_"$CHR".txt -out $SNP_DIR/background/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr"$CHR"

#main features
#-P nb of threads -nQueueSize maximum waiting in memory (necesary to optimize CPU usage
# -doMaf 1 (allele frequencies) -GL (Genotype likelihood 2 GATK method)
# -doHWE 1 (calculate deviation from Hardy-Weinberg Equilibrium)
# -doMajorMinor 1 use the most frequent allele as major
# -rf (file with the region written) work on a defined region : OPTIONAL
# -b (bamlist) input file
# -out  output file

#main filters
#filter on bam files -remove_bads (remove files with flag above 255) -minMapQ minimum mapquality -minQ (minimum quality of reads?)
#filter on frequency -minInd (minimum number of individuals with at least one read at this locus) we set it to 50%
#filter on allele frequency -minMaf, set to 0.05 
#filter on reads with several hits -uniqueOnly
#filter on pairs of reads not properly mapped -only_proper_pairs 1 (by default in ANGSD)

# 2. Extract SNP which passed the MIN_MAF and PERCENT_IND filters & their Major-minor alleles
echo "from the maf file, extract a list of SNP chr, position, major all, minor all"

## Unzip .maf file
gunzip -f $SNP_DIR/background/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr"$CHR".mafs.gz 

## Filter sites
INFILE=$SNP_DIR/background/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr"$CHR".mafs
OUTFILE_sites=$SITES_DIR/background/sites_all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr"$CHR"
BED_FILE=$SITES_DIR/background/sites_all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr"$CHR".bed

Rscript 01_scripts/utils/make_sites_list_maxdepth_simple.R "$INFILE" "$OUTFILE_sites"

## Index sites
angsd sites index $OUTFILE_sites

## Convert site file to bed for ngsparalog
awk '{print $1"\t"$2-1"\t"$2}' $OUTFILE_sites > $BED_FILE

# 3. Run mpileup and ngsParalog without intermidate files to calculate a likelihood ratio (LR) of mismapping reads covering each site
samtools mpileup -b $BAMLIST -l $BED_FILE -r $CHR -q 0 -Q 0 --ff UNMAP,DUP |
$NGSPARALOG calcLR \
    -infile - \
    -outfile $SNP_DIR/background/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr"$CHR".ngsparalog \
    -minQ 20 -minind $MIN_IND -mincov $MIN_DEPTH -allow_overwrite 1

## Convert ngsparalog output in list of canonical and deviant SNPs based on p-value threshold
Rscript 01_scripts/utils/convert_ngsparalog_to_sitelist.R \
    $SNP_DIR/background/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr"$CHR".ngsparalog \
    $OUTFILE_sites $PVAL_THRESHOLD

## Index
angsd sites index "$OUTFILE_sites"_deviant
angsd sites index "$OUTFILE_sites"_canonical
