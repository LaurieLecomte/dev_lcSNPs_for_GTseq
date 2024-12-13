argv <- commandArgs(T)
SITES <- argv[1]
POP_LIST <- argv[2]
#SUFFIX <- argv[3]

SITES <- "/project/lbernatchez/users/lalec31/projets_labo/Bastien/GTseq_saal_202410/dev_lcSNPs_for_GTseq/02_infos/sites_all_maf0.05_pctind0.65_maxdepth10_canonical_minmaj.list"
POP_FILE <- "02_infos/pop.txt"
SUFFIX <- "all_maf0.05_pctind0.65_maxdepth10"

library(data.table)
options(scipen = 999)

sites <- (fread(SITES, header = TRUE))

pops <- read.table(POP_FILE, col.names = 'pop')



combined <- sites[, c('chromo', 'position')]

for (pop in pops$pop){
  pop_maf <- paste0("08_maf_by_pop/", pop, '/', pop, SUFFIX)
  #assign(x = pop, value = as.data.frame(fread(pop_maf)[, c(1,2,6)]))
  pop_df <- (fread(pop_maf)[, c(1,2,6)])
  colnames(pop_df)[3] <- eval(pop)
  combined <- merge(combined, pop_df, by = c('chromo', 'position'), all = TRUE)
}

colnames(combined) <- c('ChromName', 'pos', colnames(combined)[3:length(colnames(combined))])

# Remove NAs
combined <- na.omit(combined)
write.table(combined, file = paste0('08_maf_by_pop/', SUFFIX, '_combined.mafs'), 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
