# Extract final selected SNP sites from the beagle file 

# 1. Access files ---------------------------------------------------------
argv <- commandArgs(T)
BEAGLE <- argv[1]
SELEC <- argv[2]
#OUT <- argv[3]

#BEAGLE <- "/project/lbernatchez/users/lalec31/projets_labo/Bastien/lcSNPs_for_GTseq_RIN/05_cand_SNPs/all/all_maf0.05_pctind0.65_maxdepth10_all_chrs_canon.beagle.gz"
beagle <- read.delim(BEAGLE)

#SELEC <- "/project/lbernatchez/users/lalec31/projets_labo/Bastien/lcSNPs_for_GTseq_RIN/10_SNP_selection/maf0.05_pctind0.65_maxdepth10_combined_0.1_scored_100.good_idy90_len160.5.tsv"
selec <- read.delim(SELEC, header=FALSE, col.names = c('CHR', 'POS'))


# 2. Extract SNPs ---------------------------------------------------------
# Paste chr and pos together to match the blast output site name
selec$site <- paste0(selec$CHR, '_', selec$POS)

# Substract SNPs from the scored table
remaining_snp <- subset(beagle, marker %in% selec$site)

# Export
write.table(remaining_snp, 
            #file = paste0(OUT, '.beagle'),
            file = paste0(strsplit(SELEC, split = '.tsv'), '.beagle'),
            col.names = sapply(X = colnames(remaining_snp), FUN = function(x) { head(unlist(strsplit(x, split = '.', fixed = TRUE)), 1) }), 
            row.names = FALSE, quote = FALSE, sep = "\t")


#library(ggplot2)
#library(dplyr)
#remaining_plot <- select(remaining_snp, marker)

#remaining_plot$CHR <- sapply(X = remaining_plot$marker, FUN = function (x) sub("([A-Z0-9_.]+)_([0-9]+)", "\\1", x))
#remaining_plot$POS <- sapply(X = remaining_plot$marker, FUN = function (x) as.numeric(sub("([A-Z0-9_.]+)_([0-9]+)", "\\2", x)))


#chrs <- read.delim("/project/lbernatchez/users/lalec31/projets_labo/Bastien/lcSNPs_for_GTseq_RIN/02_infos/chrs.bed", 
#                   header=FALSE,
#                   col.names = c('CHR', 'START', 'STOP'))
#chrs0 <- chrs

#for (i in 1:nrow(chrs)){
#  if (i > 1){
#    chrs$START[i] <- (chrs$STOP[i - 1]) + 1
#    chrs$STOP[i] <- chrs$START[i] + chrs$STOP[i] 
#  }
#}

#remaining_plot <- merge(x = remaining_plot, y = chrs, by = 'CHR')
#remaining_plot$POS_updated <- remaining_plot$POS + remaining_plot$START
#remaining_plot$STOP_updated <- remaining_plot$STOP + remaining_plot$START

#plot(remaining_plot$POS_updated, remaining_plot$POS_updated)
#abline(v = remaining_plot$STOP)
