
argv <- commandArgs(T)
SELECT <- argv[1]
GENOME <- argv[2]
WIN <- argv[3]

#SELECT <- "/project/lbernatchez/users/lalec31/projets_labo/Bastien/GTseq_saal_202410/dev_lcSNPs_for_GTseq/10_SNP_selection/all_maf0.05_pctind0.65_maxdepth10_combined_0.4_scored_50.good.12.6.tsv"
#GENOME <- '03_genome/genome.corrected.fasta'
#WIN <- 100

# Import
selected_SNPs <- read.delim(SELECT, header=FALSE, col.names = c('CHR', 'POS'))
# Get position range on reference
ref_range <- GenomicRanges::GRanges(seqnames = selected_SNPs$CHR,
                                    ranges = IRanges::IRanges(start = (selected_SNPs$POS - WIN), end = (selected_SNPs$POS + WIN) )
)
                                    
# Extract sequence from reference sequence
seq_IDs <- paste0('>', selected_SNPs$CHR, '_', selected_SNPs$POS)

ref_seqs <- Rsamtools::scanFa(GENOME, ref_range)

seqs <- (unname(as.character(ref_seqs)))

# Write sequences to output fasta
for (i in 1:length(seqs)){
  cat(seq_IDs[i], file = paste0(unlist(strsplit(SELECT, '.tsv')), '.fasta'), sep="\n", append=TRUE)
  cat(seqs[i], file = paste0(unlist(strsplit(SELECT, '.tsv')), '.fasta'), sep="\n", append=TRUE)
}
 
