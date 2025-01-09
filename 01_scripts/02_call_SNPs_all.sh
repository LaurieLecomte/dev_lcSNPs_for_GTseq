#!/bin/bash

### This script will work on all bamfiles and use ANGSD to call SNPs that pass coverage and MAF filters,
### then calculate the likelihoods that reads wer mismapped at the position of snps using ngsparalog
### to produce list of canonical and deviant SNPs


# less 02_infos/chrs.txt | parallel -j10 srun -p small -c 4 --mem=50G -J 02_call_SNPs_all_{} -o log/02_call_SNPs_all_{}_%j.log /bin/sh 01_scripts/02_call_SNPs_all.sh {} &
# srun -p small -c 4 --mem=50G -J 02_call_SNPs_all_NC_036838.1 -o log/02_call_SNPs_all_NC_036838.1_%j.log /bin/sh 01_scripts/02_call_SNPs_all.sh "NC_036838.1" &

# VARIABLES
GENOME="03_genome/genome.fasta"
BAM_DIR="04_bam"
SNP_DIR="05_cand_SNPs"
SITES_DIR="02_infos/sites_by_chr"
BAMLIST="02_infos/bam.filelist"
REGION_LIST="02_infos/regions_to_keep.txt"
REGION_BED="02_infos/regions_to_keep.bed"

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
module load bedtools/2.31.1

# Increase number of opened file limit
ulimit -S -n 2048



if [[ ! -d $SNP_DIR/all ]]
then
  mkdir $SNP_DIR/all
fi

if [[ ! -d $SITES_DIR/all ]]
then
  mkdir $SITES_DIR/all
fi




# 0. Prepare variables 
N_IND=$(wc -l $BAMLIST | cut -d " " -f 1)
MIN_IND_FLOAT=$(echo "($N_IND * $PERCENT_IND)" | bc -l)
MIN_IND=${MIN_IND_FLOAT%.*} 
MAX_DEPTH=$(echo "($N_IND * $MAX_DEPTH_FACTOR)" | bc -l)

# 0. Prepare regions
#less $REGION_LIST | grep $CHR > $SITES_DIR/"$(basename -s '.txt' $REGION_LIST)"_"$CHR".txt


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
-r "$CHR:" -out $SNP_DIR/all/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr"$CHR"

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
gunzip -f $SNP_DIR/all/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr"$CHR".mafs.gz 

## Filter sites
INFILE=$SNP_DIR/all/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr"$CHR".mafs
OUTFILE_sites=$SITES_DIR/all/sites_all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr"$CHR"
BED_FILE=$SITES_DIR/all/sites_all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr"$CHR".bed
FILT_BED_FILE=$SITES_DIR/all/sites_all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr"$CHR"_filt.bed

Rscript 01_scripts/utils/make_sites_list_maxdepth_simple.R "$INFILE" "$OUTFILE_sites"

## Index sites
angsd sites index $OUTFILE_sites

# 3. Filter out SNPs located in excluded regions
## Convert site file to bed for ngsparalog
awk '{print $1"\t"$2-1"\t"$2}' $OUTFILE_sites > $BED_FILE

## Convert canonical list to bed, then filter out sites located in unwanted regions
bedtools intersect -a $BED_FILE -b $REGION_BED > $FILT_BED_FILE
echo "before filtering for forbidden regions: "$(wc -l $BED_FILE)" sites"
echo "after filtering for forbidden regions: "$(wc -l $FILT_BED_FILE)" sites"
    

# 3. Run mpileup and ngsParalog without intermidate files to calculate a likelihood ratio (LR) of mismapping reads covering each site
#samtools mpileup -b $BAMLIST -l $BED_FILE -r $CHR -q 0 -Q 0 --ff UNMAP,DUP |
samtools mpileup -b $BAMLIST -l $FILT_BED_FILE -r $CHR -q 0 -Q 0 --ff UNMAP,DUP |
$NGSPARALOG calcLR \
    -infile - \
    -outfile $SNP_DIR/all/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr"$CHR".ngsparalog \
    -minQ $MIN_Q -minind $MIN_IND -mincov $MIN_DEPTH -allow_overwrite 1

## Convert ngsparalog output in list of canonical and deviant SNPs based on p-value threshold
Rscript 01_scripts/utils/convert_ngsparalog_to_sitelist.R \
    $SNP_DIR/all/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr"$CHR".ngsparalog \
    $OUTFILE_sites $PVAL_THRESHOLD

## Index
angsd sites index "$OUTFILE_sites"_deviant
angsd sites index "$OUTFILE_sites"_canonical
echo "$(less "$OUTFILE_sites"_canonical | wc -l) sites remaining after removing deviant SNPs"

## Convert to bed for reference if needed
awk '{print $1"\t"$2-1"\t"$2}' "$OUTFILE_sites"_deviant > "$OUTFILE_sites"_deviant.bed
awk '{print $1"\t"$2-1"\t"$2}' "$OUTFILE_sites"_canonical > "$OUTFILE_sites"_canonical.bed
