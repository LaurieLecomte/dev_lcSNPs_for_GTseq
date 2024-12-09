# 1. Access files ---------------------------------------------------------
argv <- commandArgs(T)
PAIRWISE_FST <- argv[1]

#PAIRWISE_FST <- "/project/lbernatchez/users/lalec31/projets_labo/Bastien/GTseq_saal_202410/dev_lcSNPs_for_GTseq/09_fst/pairwise_fst_both.txt"

fst <- read.table(PAIRWISE_FST, col.names = c('POP1', 'POP2', 'FST_unw', 'FST_w'))

fst <- fst[, c(1,2,4)]

library(reshape2)
fst <- melt(fst)


library(ggplot2)
fst_map <- 
ggplot(data = fst, aes(x=POP1, y=POP2, fill=value)) + 
  geom_tile() +
  geom_text(aes(POP2, POP1, label = round(value, digits = 4)), color = "black", size = 2) + 
  scale_fill_gradient(low = "royalblue", high = "red4")

fst_map
jpeg(file = paste0(strsplit(x = PAIRWISE_FST, split = '.txt')[[1]],
                   "_heatmap.rds"))
saveRDS(fst_map, file = paste0(strsplit(x = PAIRWISE_FST, split = '.txt')[[1]],
                                  "_heatmap.rds"))
