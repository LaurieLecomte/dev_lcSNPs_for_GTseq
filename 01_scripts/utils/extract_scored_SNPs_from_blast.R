# Extract sites that passed blast filters from the scored SNP table 

# 1. Access files ---------------------------------------------------------
argv <- commandArgs(T)
BLAST <- argv[1]
SCORED <- argv[2]
OUT <- argv[3]

blast_snp <- read.delim(BLAST, header=FALSE)
scored_snp <- read.delim(SCORED, header=TRUE)


# 2. Extract SNPs ---------------------------------------------------------
# Paste chr and pos together to match the blast output site name
scored_snp$site <- paste0(scored_snp$Chrom, '_', scored_snp$Position)

# Substract SNPs from the scored table
remaining_snp <- subset(scored_snp, site %in% blast_snp$V1)

# Export
BLAST_STR <- sub("^.+filtered_([idylen0-9_.]+).blast", "\\1", BLAST)
write.table(remaining_snp, 
            file = OUT,
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")