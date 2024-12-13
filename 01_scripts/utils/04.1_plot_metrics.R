# Plot scored SNPs metrics
argv <- commandArgs(T)
SCORED <- argv[1]

#SCORED <- "/project/lbernatchez/users/lalec31/projets_labo/Bastien/GTseq_saal_202410/dev_lcSNPs_for_GTseq/08_maf_by_pop/all_maf0.05_pctind0.65_maxdepth10_combined_0.05_scored_50"

library(data.table)

# 1. Load data
d <- fread(SCORED, header = TRUE)

# 2. Plot columns of interest to explore filter thresholds
d <- d[, -c("Position", "Chrom", "Position2", "Chrom2", "Major", "Minor", "Ancestral", "SequenceMod", "Sequence", "NumSamples")]

jpeg(file = paste0(SCORED, ".jpg"), width = 800, height = 500, quality = 90)
p <- plot(d[1:10000, ], pch=19, col="#00000011", cex=0.5)
p
dev.off()

