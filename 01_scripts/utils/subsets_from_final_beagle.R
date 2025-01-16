# Generate subsets of SNPs from the final SNP set

# 1. Access files ---------------------------------------------------------
argv <- commandArgs(T)
BEAGLE_SELEC <- argv[1]
ITERATIONS <- argv[2]
SUBSET_SIZE <- argv[3]

beagle_selec <- read.delim(BEAGLE_SELEC)


# 2. Make N subsets of final SNPs
for (i in 1:iterations){
  snp_subset <- beagle_selec[sort(sample(1:nrow(beagle_selec), size = 100, replace = FALSE)), ]
  
  write.table(snp_subset, file = paste0(strsplit(BEAGLE_SELEC, split = '.beagle'), '_subset', i, '.beagle'),
              col.names = sapply(X = colnames(remaining_snp), FUN = function(x) { head(unlist(strsplit(x, split = '.', fixed = TRUE)), 1) }),
              row.names = FALSE, quote = FALSE, sep = "\t")
}
