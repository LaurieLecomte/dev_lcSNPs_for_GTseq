# Apply filters to SNPs based on flanking SNPs and sequence properties

argv <- commandArgs(T)
SCORED <- argv[1]

SCORED <- "/project/lbernatchez/users/lalec31/projets_labo/Bastien/GTseq_saal_202410/dev_lcSNPs_for_GTseq/08_maf_by_pop/all_maf0.05_pctind0.65_maxdepth10_combined_0.05_scored_50"

# Cleanup
rm(list=ls())

# Modules
#library(tidyverse)
library(data.table)

# Load data
d = fread(SCORED, header = TRUE)

# Plot columns of interest to explore filter thresholds
dsub = d
dsub = dsub[, -c("Position", "Chrom", "Position2", "Chrom2", "Major", "Minor", "Ancestral", "SequenceMod", "Sequence", "NumSamples")]
plot(dsub[1:10000, ], pch=19, col="#00000011", cex=0.5)

# Set thresholds
max_SumMAFs = 2
max_NumSNPs = 10
#min_complexity = 100
min_complexity = 90
min_GCcontent = 0.2
max_GCcontent = 0.8

# Filter
d_filt = d[d$SumMAFs <= max_SumMAFs &
        d$NumSNPs <= max_NumSNPs &
        d$Complexity >= min_complexity &
        d$GCcontent >= min_GCcontent &       ## tell Eric d$
        d$GCcontent <= max_GCcontent, ]      ## tell Eric d$



# Plot after filters
d_filt_sub = d_filt[, -c("Position", "Chrom", "Position2", "Chrom2", "Major", "Minor", "Ancestral", "SequenceMod", "Sequence", "NumSamples")]
plot(d_filt_sub[1:10000], pch=19, col="#00000004", cex=0.5)

# Write file
write.table(d_filt, paste0(SCORED, ".good.tsv"), quote=F, row.names=F, sep="\t")
