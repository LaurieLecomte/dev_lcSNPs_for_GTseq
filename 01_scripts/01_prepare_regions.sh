#!/bin/bash

# /bin/sh 01_scripts/01_prepare_regions.sh &

# VARIABLES
GENOME="03_genome/genome.fasta"

CHR_LIST="02_infos/chrs.txt"

EXCL_DIR="02_infos/regions_to_exclude"

# LOAD REQUIRED MODULES
module load bedtools/2.31.1
module load bedops/2.4.40

# Create a bed from indexed genome, for all chromosomes and contigs
less "$GENOME".fai | awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' > "$GENOME".bed

# Create a bed only for chromosomes
if [[ -f 02_infos/chrs.bed ]]
then
  rm 02_infos/chrs.bed
fi
 
less $CHR_LIST | while read CHR || [ -n "$line" ]; do grep -Fw $CHR "$GENOME".bed >> 02_infos/chrs.bed; done


# Combine bed files of regions to exclude together
cat $(find $EXCL_DIR -type f -name "*.bed") | sort -k1,2 > 02_infos/excluded_regions_cat.bed

# "Collapse" together overlapping regions 
bedtools merge -i 02_infos/excluded_regions_cat.bed > 02_infos/excluded_intervals.bed

# Get a list of regions to keep for SNP calling
sort-bed 02_infos/chrs.bed > 02_infos/chrs.sorted.bed
sort-bed 02_infos/excluded_intervals.bed > 02_infos/excluded_intervals.sorted.bed

bedops --partition 02_infos/chrs.sorted.bed 02_infos/excluded_intervals.sorted.bed > 02_infos/chrs_excluded_intervals.partition.bed

bedops --not-element-of 1 02_infos/chrs_excluded_intervals.partition.bed 02_infos/excluded_intervals.sorted.bed | sort -k1,2 > 02_infos/regions_to_keep.bed

## Convert to a list readable by angsd
less 02_infos/regions_to_keep.bed | sed -E 's/\t/:/' | sed -E 's/\t/-/' > 02_infos/regions_to_keep.txt