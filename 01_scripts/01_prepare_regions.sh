#!/bin/bash

# First make bed files of regions to exclude from analysis and put them in 02_infos/regions_to_exclude

## Run RepeatModeler and RepeatMasker to get repeats, and convert the gff to bed:
#module load bedops/2.4.40
#GFF="RepeatMasker/genome.fasta.out.gff"
#EXCL_DIR="02_infos/regions_to_exclude"
#less $GFF | gff2bed | cut -f1-3 > $EXCL_DIR/repeats.bed 

## Extract N positions in the reference fasta
#module load seqkit/2.8.0
#seqkit locate -P -i -r -G -p "[^A^T^C^G]+" --bed 03_genome/genome.fasta | cut -f1-3 > 02_infos/regions_to_exclude/ambiguous_bases.bed

## Add bed files of sex-linked regions to 02_infos/regions_to_exclude

## Run this script
# /bin/sh 01_scripts/01_prepare_regions.sh &


# VARIABLES
FASTA="03_genome/genome.fasta"

CHR_LIST="02_infos/chrs.txt"

EXCL_DIR="02_infos/regions_to_exclude"

# LOAD REQUIRED MODULES
module load bedtools/2.31.1
module load bedops/2.4.40
module load samtools/1.15

# Create a bed from indexed genome, for all chromosomes and contigs
less "$FASTA".fai | awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' > "$FASTA".bed

# Create a bed only for chromosomes
if [[ -f 02_infos/chrs.bed ]]
then
  rm 02_infos/chrs.bed
fi
 
less $CHR_LIST | while read CHR || [ -n "$line" ]; do grep -Fw $CHR "$FASTA".bed >> 02_infos/chrs.bed; done


# Combine bed files of regions to exclude together
cat $(find $EXCL_DIR -type f -name "*.bed") | bedtools sort > 02_infos/excluded_regions_cat.bed

# "Collapse" together overlapping regions 
bedtools merge -i 02_infos/excluded_regions_cat.bed > 02_infos/excluded_intervals.bed

# Get a list of regions to exclude for SNP calling
sort-bed 02_infos/chrs.bed > 02_infos/chrs.sorted.bed
sort-bed 02_infos/excluded_intervals.bed > 02_infos/excluded_intervals.sorted.bed



# Get list of regions to keep 
bedtools subtract -a 02_infos/chrs.sorted.bed -b 02_infos/excluded_intervals.sorted.bed > 02_infos/regions_to_keep.bed
less 02_infos/regions_to_keep.bed | sed -E 's/\t/:/' | sed -E 's/\t/-/' > 02_infos/regions_to_keep.txt




# Mask regions to exclude in the genome fasta
#bedtools maskfasta -mc N -fi $FASTA -bed 02_infos/excluded_intervals.sorted.bed -fo "${FASTA%.*}".corrected.fasta

# Index 
#samtools faidx "${FASTA%.*}".corrected.fasta


#bedops --partition 02_infos/chrs.sorted.bed 02_infos/excluded_intervals.sorted.bed > 02_infos/chrs_excluded_intervals.partition.bed

#bedops --not-element-of 1 02_infos/chrs_excluded_intervals.partition.bed 02_infos/excluded_intervals.sorted.bed | sort -k1,2 > 02_infos/regions_to_keep.bed
## Convert to a list readable by angsd
#less 02_infos/regions_to_keep.bed | sed -E 's/\t/:/' | sed -E 's/\t/-/' > 02_infos/regions_to_keep.txt