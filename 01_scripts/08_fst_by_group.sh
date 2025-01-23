#!/bin/bash

# srun -p medium -c 4 --mem=50G -J 08_fst_by_group -o log/08_fst_by_group_%j.log /bin/sh 01_scripts/08_fst_by_group.sh &


###this script will 
#1 make a subset of bamlist to ahev equal nb of samples
#2calculate the saf by pop
#3 calculate 2dSFS and FST
#it uses the last version of angsd (unfold for saf and fold for realSFS)

GENOME="03_genome/genome.fasta"
BAM_DIR="04_bam"
SNP_DIR="05_cand_SNPs"
SITES_DIR="02_infos/sites_by_chr"
BAMLIST="02_infos/bam.filelist"

PCA_DIR="06_pca"
ADMIX_DIR="07_admix"
MAF_DIR="08_maf_by_pop"
FST_DIR="09_fst"

CHR_LIST="02_infos/chrs.txt"
ID_POP="02_infos/ID_POP.txt" 

POP_FILE1="02_infos/pops.txt"

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


PVAL_THRESHOLD=0.001
#REGION_NUM="$1" 
#REGION=$(head -n $REGION_NUM 02_info/regions.txt | tail -n 1)

# Execs
NGSPARALOG="/project/lbernatchez/users/lalec31/softwares/ngsParalog/ngsParalog"
NGSADMIX="/project/lbernatchez/users/lalec31/softwares/NGSadmix"
REALSFS="/prg/angsd/0.937/misc/realSFS"
#REALSFS="/prg/angsd/0.931/misc/realSFS"

#maybe edit
NB_CPU=4 #change accordingly in SLURM header
NSITES=500000 #to make realSFS goes faster -reduce the number of sites considered
GROUP=pops #the subgroup on whcih we are making the fst comparison -> it should be a file like GROUP.txt in the folder 02_info
#POP_FILE1=02_info/"$GROUP".txt #choose on which list of pop run the analyses





# LOAD REQUIRED MODULES
#module load angsd/0.937
#module load samtools/1.15
#module load R/4.2
module load angsd/0.937 #only with this version the SFS/FST script runs well (edit in sept 2022)

#make sure you index the sites file with the same version 
ulimit -S -n 2048

# Before analysis and if not done before : combine per-chromosome canonical SNP lists


#make a folder in which write new results
mkdir $FST_DIR/$GROUP


# 1. Subset bamfilelist to have equal number of bam files per pop -> necessary ?

Rscript 01_scripts/utils/subset_random_Nind.r "$GROUP" $FST_DIR

#2 do saf for all population listed
cat $POP_FILE1 | while read i
do
  echo $i

  N_IND=$(wc -l $FST_DIR/$GROUP/"$i"subsetbam.filelist | cut -d " " -f 1) # or N_IND=$(wc -l $MAF_DIR/$i/"$i"bam.filelist | cut -d " " -f 1)
  MIN_IND_FLOAT=$(echo "($N_IND * $PERCENT_IND)"| bc -l)
  MIN_IND=${MIN_IND_FLOAT%.*} 

  echo "working on pop $i, $N_IND individuals, will use the sites file provided"
  echo "will filter for sites with at least one read in $MIN_IND individuals, which is $PERCENT_IND of the total"

  # CHRECK FOR -anc or -ref and -doMajorMinor
  # CORRECT MIN_IND
  angsd -P $NB_CPU \
  -dosaf 5 -GL 2 -doMajorMinor 4 \
  -ref $GENOME \
  -rf $CHR_LIST \
  -remove_bads 1 -minMapQ $MIN_MAPQ -minQ $MIN_Q -minInd $MIN_IND -setMinDepthInd $MIN_DEPTH \
  -sites 02_infos/sites_all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_canonical_sites \
  -b $FST_DIR/$GROUP/"$i"subsetbam.filelist -out $FST_DIR/$GROUP/"$i"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"
done



# 3 Calculate FST
#prepare variables - avoid to modify
num_pops=$(wc -l "$POP_FILE1" | cut -d " " -f 1)

# Estimate pairwise FST for all populations listed

for i in $(seq $num_pops)
do
	pop1=$(cat "$POP_FILE1" | head -"$i" | tail -1)
	for j in $(seq $[ $i + 1 ] $num_pops)
	do
		pop2=$(cat "$POP_FILE1" | head -"$j" | tail -1)
		echo "FST between $pop1 and $pop2"
		echo "$pop1"
		echo "$pop2"
		
		echo "calcualte the 2dsfs priors"
		$REALSFS $FST_DIR/$GROUP/"$pop1"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".saf.idx \
    $FST_DIR/$GROUP/"$pop2"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".saf.idx \
    -P $NB_CPU -maxIter 30 -nSites $NSITES > $FST_DIR/$GROUP/"$pop1"_"$pop2"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"."$NSITES"
    
    file=$FST_DIR/$GROUP/"$pop1"_"$pop2"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"."$NSITES"
    
    Rscript 01_scripts/utils/sum_sites_2dsfs.r "$file"
		
		echo " prepare the fst for easy window analysis etc"
		$REALSFS fst index $FST_DIR/$GROUP/"$pop1"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".saf.idx \
    $FST_DIR/$GROUP/"$pop2"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".saf.idx \
    -sfs $FST_DIR/$GROUP/"$pop1"_"$pop2"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"."$NSITES".2dsfs \
    -P $NB_CPU -fstout $FST_DIR/$GROUP/"$pop1"_"$pop2"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"

		echo "print SFS priori for each position"
		$REALSFS fst print $FST_DIR/$GROUP/"$pop1"_"$pop2"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".fst.idx \
    -P $NB_CPU > $FST_DIR/$GROUP/"$pop1"_"$pop2"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".bypos.sfs
		
		echo "get the global estimate of FST throughout the genome"
		$REALSFS fst stats $FST_DIR/$GROUP/"$pop1"_"$pop2"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".fst.idx \
    -P $NB_CPU > $FST_DIR/$GROUP/"$pop1"_"$pop2"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".fst
		
		echo "calculate FST by slidingwindow, window size=$WINDOW and step=$WINDOW_STEP, as given in 01_config.sh"
		$REALSFS fst stats2 $FST_DIR/$GROUP/"$pop1"_"$pop2"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".fst.idx \
    -win $WINDOW -step $WINDOW_STEP -P $NB_CPU > $FST_DIR/$GROUP/"$pop1"_"$pop2"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".slidingwindow
    
    # Output paiwise Fst in a file
    FST_UNW=$(less $FST_DIR/$GROUP/"$pop1"_"$pop2"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".fst | cut -f1)
    FST_W=$(less $FST_DIR/$GROUP/"$pop1"_"$pop2"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".fst | cut -f2)
    echo -e "$pop1\t$POP2\t$FST_UNW\t$FST_W\n$pop2\t$pop1\t$FST_UNW\t$FST_W" >> $FST_DIR/paiwise_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".txt
	done
done


# Output paiwise Fst in a single file
#for i in $(ls 09_fst/pop/*10.fst); do PAIR=$(echo $i | sed -E 's/.+\/([A-Za-z]+\_[A-Za-z]+)\_maf.+/\1/'); POP1=$(echo ${PAIR%_*} ); POP2=$(echo "${PAIR##*_}"); FST1=$(less $i | cut -f1); FST2=$(less $i | cut -f2); echo -e "$POP1\t$POP2\t$FST1\t$FST2\n$POP2\t$POP1\t$FST1\t$FST2" >> 09_fst/pairwise_fst.txt; done