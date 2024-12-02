#!/bin/bash

# srun -p medium -c 4 --mem=30G -J 07_maf_by_pop -o log/07_maf_by_pop_%j.log /bin/sh 01_scripts/07_maf_by_pop.sh &

# VARIABLES
GENOME="03_genome/genome.fasta"
BAM_DIR="04_bam"
SNP_DIR="05_cand_SNPs"
SITES_DIR="02_infos/sites_by_chr"
BAMLIST="02_infos/bam.filelist"

PCA_DIR="06_pca"
ADMIX_DIR="07_admix"
MAF_DIR="08_maf_by_pop"


CHR_LIST="02_infos/chrs.txt"
ID_POP="02_infos/ID_POP.txt"

POP_FILE1="02_infos/pop_list.txt"

REGION_LIST="02_infos/regions_to_keep.txt"

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



ulimit -S -n 2048


# Before analysis : combine per-chromosome canonical SNP lists
cat $SITES_DIR/sites_all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr*_canonical > 02_infos/sites_all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_canonical
## Index this sites file
angsd sites index 02_infos/sites_all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_canonical

# Do maf for all population listed
cat $POP_FILE1 | while read i
do
  echo $i
  mkdir $MAF_DIR/$i
  
  # Make a list of bam files for a given population #### to correct 
  #less $ID_POP | awk -v val="$i" 'BEGIN{FS="\t"} $2 == val {print}' | cut -f1 > $MAF_DIR/$i/"$i".samples
  #less $BAMLIST | grep -f $MAF_DIR/$i/"$i".samples > $MAF_DIR/$i/"$i"bam.filelist
  
  less $BAMLIST | grep $i > 02_infos/"$i"bam.filelist
  
  N_IND=$(wc -l $MAF_DIR/$i/"$i"bam.filelist | cut -d " " -f 1)
  MIN_IND_FLOAT=$(echo "($N_IND * $PERCENT_IND)"| bc -l)
  MIN_IND=${MIN_IND_FLOAT%.*} 
  
  echo "working on pop $i, $N_IND individuals, will use the sites file provided"
  echo "will filter for sites with at least one read in $MIN_IND individuals, which is $PERCENT_IND of the total"
  
  angsd -P $NB_CPU -nQueueSize 50 \
  -doMaf 1 -GL 2 -doMajorMinor 3 \
  -ref $GENOME \
  -rf $REGION_LIST \
  -remove_bads 1 -minMapQ 30 -minQ 20 -minInd $MIN_IND -setMinDepthInd $MIN_DEPTH \
  -sites 02_infos/sites_all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_canonical \
  -b 02_infos/"$i"bam.filelist -out $MAF_DIR/$i/"$i"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"
  
  gunzip -f $MAF_DIR/$i/"$i"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".mafs.gz
done

