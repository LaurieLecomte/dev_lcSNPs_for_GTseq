#!/bin/sh

# Detect repeats in the final SaFo assembly using the custom repeat elements library made with RepeatModeler AND repeats for Salmonidae in famdb

# Run on Manitou
# srun -p medium -c 10 --mem=50G --time=7-00:00:00 -J 02.2_final_RepeatMasker -o log/02.2_final_RepeatMasker_%j.log /bin/sh 01_scripts/02.2_final_RepeatMasker.sh &


# VARIABLES
GENOME_DIR="03_genome"
GENOME="$GENOME_DIR/genome.fasta"

RMOD_DIR="RepeatModeler"
RMAS_DIR="RepeatMasker"
DB_NAME="saal"

FAMILY="Salmonidae"
SPECIES="Salvelinus"
CLASS_LIB="$RMOD_DIR/"$DB_NAME"-families.fa" #outputted by RepeatModeler.sh
FAM_LIB="$RMAS_DIR/famdb_$FAMILY_desc.fa" #created in current script

COMB_LIB="$RMAS_DIR/combined_"$DB_NAME"_"$FAMILY"_desc.fa" # we use the custom library produced by the final_NCBI_RepeatMasker_famdb.sh script

CPU=10

# LOAD REQUIRED MODULES 
module load gnu-openmpi/4.1.4
module load exonerate/2.4.0
module load RepeatMasker/4.1.2
module load ncbiblast/2.6.0
module load python/2.7



# 1. Create custom lib using the famdb utility
famdb.py families --format 'fasta_name' --descendants $FAMILY --include-class-in-name > $FAM_LIB

# 2. Combine with custom lib produced by RepeatModeler
## Remove previous 
if [[ -f "${CLASS_LIB%.*}"_renamed.fa ]]
then
  rm "${CLASS_LIB%.*}"_renamed.fa
fi
## First edit sequence header in the RepeatModeler lib
less $CLASS_LIB | while read line; do echo $line | sed -E "s/(^.+)\ \(.+/\1\ \@$SPECIES/" >> "${CLASS_LIB%.*}"_renamed.fa ; done

## Concatenate
cat $RMAS_DIR/famdb_"$FAMILY"_desc.fa "${CLASS_LIB%.*}"_renamed.fa > $COMB_LIB

# 3. Run RepeatMasker on all contigs, using the combined custom library
RepeatMasker -pa $CPU $GENOME -dir $RMAS_DIR -gff -lib $COMB_LIB
