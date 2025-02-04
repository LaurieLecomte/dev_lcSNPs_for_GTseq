#!/bin/sh

# Run on Manitou
# srun -p medium -c 10 --mem=50G --time=7-00:00:00 -J RepeatModeler -o log/RepeatModeler_%j.log /bin/sh 01_scripts/RepeatModeler.sh &


# VARIABLES
GENOME_DIR="03_genome"
GENOME="$GENOME_DIR/genome.fasta"

RMOD_DIR="RepeatModeler"
RMAS_DIR="RepeatMasker"
DB_NAME="saal"

CPU=10

# LOAD REQUIRED MODULES 
module load RepeatModeler/2.0.1

# 0. Make directory for RepeatModeler if needed
if [[ ! -d $RMOD_DIR ]]
then
  mkdir $RMOD_DIR
fi

# De novo detection of repeats in whole genome
cd "$RMOD_DIR"
BuildDatabase -name $DB_NAME ../$GENOME
RepeatModeler -pa $CPU -database $DB_NAME

