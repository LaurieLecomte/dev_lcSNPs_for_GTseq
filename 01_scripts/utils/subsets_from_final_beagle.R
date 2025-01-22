# Generate subsets of SNPs from the final SNP set

# 1. Access files ---------------------------------------------------------
argv <- commandArgs(T)
BEAGLE_SELEC <- argv[1]
ITERATIONS <- as.numeric(argv[2])
SUBSET_SIZE <- as.numeric(argv[3])
OUT_DIR <- argv[4]

beagle_selec <- read.delim(BEAGLE_SELEC)


# 2. Make N subsets of final SNPs
for (i in 1:ITERATIONS){
  snp_subset <- beagle_selec[sort(sample(1:nrow(beagle_selec), size = SUBSET_SIZE, replace = FALSE)), ]
  
  write.table(snp_subset, file = paste0(OUT_DIR, '/', strsplit(basename(BEAGLE_SELEC), split = '.beagle.gz'), '_subset', i, '_', SUBSET_SIZE, 'sites', '.beagle'),
              col.names = sapply(X = colnames(beagle_selec), FUN = function(x) { head(unlist(strsplit(x, split = '.', fixed = TRUE)), 1) }),
              row.names = FALSE, quote = FALSE, sep = "\t")
}
