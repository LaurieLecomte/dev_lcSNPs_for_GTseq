#!/bin/bash

# parallel -a 02_infos/samples40.txt -j 10 srun -p small -c 1 -J coverage_{} -o log/coverage_{}_%j.log /bin/sh 01_scripts/coverage.sh {} &

# VARIABLES
SAMPLE=$1
BAM="04_bam/"$SAMPLE".bam"
COV_DIR="coverage"
REGIONS="02_infos/chrs.bed"

# LOAD REQUIRED MODULES
module load mosdepth/0.3.6


if [[ ! -d $COV_DIR ]]
then
  mkdir $COV_DIR
fi

mosdepth -n --fast-mode --by 1000 $COV_DIR/"$SAMPLE" $BAM

MEAN_COV=$(zless $COV_DIR/"$SAMPLE".regions.bed.gz | awk '{sum+=$4} END { print sum/NR}')
echo -e "$SAMPLE\t$MEAN_COV" > $COV_DIR/"$SAMPLE".meancov.txt
