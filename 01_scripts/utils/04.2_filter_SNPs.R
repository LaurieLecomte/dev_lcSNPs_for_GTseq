# Filter SNPs
argv <- commandArgs(T)
SCORED <- argv[1]

max_SumMAFs <- argv[2]
max_NumSNPs <- argv[3]
min_complexity <- argv[4]
min_GCcontent <- argv[5]
max_GCcontent <- argv[6]

#SCORED <- "/project/lbernatchez/users/lalec31/projets_labo/Bastien/GTseq_saal_202410/dev_lcSNPs_for_GTseq/08_maf_by_pop/all_maf0.05_pctind0.65_maxdepth10_combined_0.05_scored_50"


library(data.table)

# 1. Load data
d <- fread(SCORED, header = TRUE)

# 2. Filter
d_filt <- d[d$SumMAFs <= max_SumMAFs &
             d$NumSNPs <= max_NumSNPs &
             d$Complexity >= min_complexity &
             d$GCcontent >= min_GCcontent &       
             d$GCcontent <= max_GCcontent, ]      



# 3. Plot after filters
d_filt_sub <- d_filt[, -c("Position", "Chrom", "Position2", "Chrom2", "Major", "Minor", "Ancestral", "SequenceMod", "Sequence", "NumSamples")]
jpeg(file = paste0(SCORED, "_filt.jpg"), width = 800, height = 500, quality = 90)

plot(d_filt_sub[1:10000], pch=19, col="#00000004", cex=0.5)
dev.off()

# 4. Write file
write.table(d_filt, paste0(SCORED, ".good.tsv"), quote=F, row.names=F, sep="\t")
