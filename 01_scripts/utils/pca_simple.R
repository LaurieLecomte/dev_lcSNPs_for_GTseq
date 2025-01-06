# Run PCA on covariance matrix

# 1. Access files ---------------------------------------------------------
argv <- commandArgs(T)
COV_MAT <- argv[1]
BAM_LIST <- argv[2]
BEAGLE <- argv[3]
ID_POP <- argv[4]


# 2. Run PCA on cov mat ---------------------------------------------------
print(paste("read cov matrix", COV_MAT))
cov_mat <- as.matrix(read.table(COV_MAT), header = FALSE)
pca <- eigen(cov_mat)

pca.mat <- as.matrix(pca$vectors %*% (diag(pca$values))^0.5)


# Add column names
nPC <- dim(pca$vectors)[2]
col_PC <- vector(length = nPC)
for (i in 1 : nPC) {col_PC[i] <- paste0("PC", i)}
colnames(pca.mat) <- c(col_PC)

# Add row names based on bam list
bam <- read.table(BAM_LIST, sep = "\t", header = FALSE, col.names = 'bam') 
sample_IDs <- sapply(X = bam$bam, FUN = function(x) sub(pattern = '.bam', replacement = '', basename(x)))

#beagle_names <- read.table(BEAGLE, nrows = 1)
#beagle_names <- as.vector(paste(beagle_names[1,]))
#sample_names <- unique(beagle_names[4:length(beagle_names)])

rownames(pca.mat) <- sample_IDs

# 3. Compute PCA stats ----------------------------------------------------
# Calculate varsum (eigen_mats$values[eigen_mats$values>=0]
var1 <- round(pca$values[1]*100/sum(pca$values[pca$values >= 0]), 2)
var2 <- round(pca$values[2]*100/sum(pca$values[pca$values >= 0]), 2)
var3 <- round(pca$values[3]*100/sum(pca$values[pca$values >= 0]), 2)
var4 <- round(pca$values[4]*100/sum(pca$values[pca$values >= 0]), 2)

# Make kmeans for n expected groups on PC1
## Get number of expected clusters = number of pops
ID_pop <- read.table(ID_POP, sep = "\t", header = FALSE) [, 1:2]
colnames(ID_pop) <- c('ID', 'pop')
pop_num <- length(unique(ID_pop$pop))

## Get kmeans
kmeans_res <- kmeans(x = as.matrix(pca.mat[, 1]), c(min(pca.mat[, 1]), median(pca.mat[, 1]), max(pca.mat[, 1])))
#kmeans_res <- kmeans(x = as.matrix(pca.mat[, 1]), pop_num)
k_ss <- round(kmeans_res$betweenss/kmeans_res$totss, pop_num)

# Save 4PCS eigenvalues and k means SS
write.table(pca.mat[, 1:3], paste0(COV_MAT, ".pca"), quote = FALSE)
write.table(c(var1, var2, var3, var3), paste0(COV_MAT, ".eig"), quote = FALSE)


# 4. Plot PCA -------------------------------------------------------------
jpeg(file = paste0(COV_MAT, ".pca.jpg"))
par(mfrow = c(1, 1))
plot(pca.mat[, 1], pca.mat[, 2], pch = 20, ylab = paste("PC2", var2), xlab = paste("PC1", var1), 
     col = kmeans_res$cluster, main = paste("k_SS", k_ss))
#plot(pca.mat[,3], pca.mat[,4], pch=20, ylab=paste("PC4", var4), xlab=paste("PC3",var3))


# Upgrade with ggplot2
library(ggplot2)
library(dplyr)

pca_df <- as.data.frame(pca.mat)

## Import ID pop file for coloring points
#ID_pop <- read.table(ID_POP, sep = "\t", header = FALSE) [, 1:2]
#colnames(ID_pop) <- c('ID', 'pop')

## Assign pop values according to IDs 
POP <- c()
ID <- c()

for (i in rownames(pca_df)) {
  # extract this ID in ID_POP
  sample_df <- subset(ID_pop, ID == i)
  ID[i] <- sample_df$ID
  POP[i] <- sample_df$pop
}

pca_df <- cbind(pca_df, ID, POP)

## Add cluster info
if (all(rownames(pca_df) == names(kmeans_res$cluster) )) {
  pca_df$cluster <- kmeans_res$cluster
}

## Plot
### PC1-2
jpeg(file = paste0(COV_MAT, ".PC1_PC2.jpg"))
par(mfrow = c(1,1))
plot_PC1_PC2 <- ggplot(data = pca_df, aes(x = PC1, y = PC2, label = ID)) +
  geom_point(aes(col = POP, shape = factor(cluster))) +
  geom_text(size = 1, hjust = 0, 
            nudge_x = -0.002, nudge_y = 0.002) +
  labs(x = paste("PC1", var1), y = paste("PC2", var2), 
       title = paste("k_SS" ,k_ss) ) #+
  #scale_color_manual(values = c("red", "blue"))
plot_PC1_PC2

dev.off()

ggsave(plot_PC1_PC2, 
       filename = paste0(COV_MAT, ".PC1_PC2.pdf"),
       device = "pdf", dpi = 320,
       height = 6, width = 5, units = "in")

### PC1-3
#jpeg(file = paste0(COV_MAT, ".PC1_PC3.jpg"))
#plot_PC1_PC3 <- ggplot(data = pca_df, aes(x = PC1, y = PC3, label = ID)) +
#  geom_point(aes(col = POP, shape = factor(cluster))) +
#  geom_text(size = 1.5, hjust = 0, 
#            nudge_x = -0.002, nudge_y = 0.002, check_overlap = TRUE) +
#  labs(x = paste("PC1", var1), y = paste("PC3", var3), 
#       title = paste("k_SS", k_ss) ) 
#+ scale_color_manual(values = c("red", "blue"))
#plot_PC1_PC3
#dev.off()

### PC2-3
#jpeg(file = paste0(COV_MAT, ".PC2_PC3.jpg"))
#plot_PC2_PC3 <- ggplot(data = pca_df, aes(x = PC2, y = PC3, label = ID)) +
#  geom_point(aes(col = POP, shape = factor(cluster))) +
#  geom_text(size = 1.5, hjust = 0, 
#            nudge_x = -0.002, nudge_y = 0.002, check_overlap = TRUE) +
#  labs(x = paste("PC2", var2), y = paste("PC3", var3), 
#       title = paste("k_SS" ,k_ss) ) 
 #+ scale_color_manual(values = c("red", "blue"))
#plot_PC2_PC3
#dev.off()

### PC3-4
jpeg(file = paste0(COV_MAT, ".PC3_PC4.jpg"))
plot_PC3_PC4 <- ggplot(data = pca_df, aes(x = PC3, y = PC4, label = ID)) +
  geom_point(aes(col = POP, shape = factor(cluster))) +
  geom_text(size = 1, hjust = 0, 
            nudge_x = -0.002, nudge_y = 0.002) +
  labs(x = paste("PC3", var3), y = paste("PC4", var4), 
       title = paste("k_SS", k_ss) ) 

ggsave(plot_PC3_PC4, 
       filename = paste0(COV_MAT, ".PC3_PC4.pdf"),
       device = "pdf", dpi = 320,
       height = 6, width = 5, units = "in")

## Store as rds object for easier handling
saveRDS(plot_PC1_PC2, file = paste0(strsplit(x = COV_MAT, split = '.cov')[[1]],
                                    "_PC1_PC2.rds"))
#saveRDS(plot_PC1_PC3, file = paste0(strsplit(x = COV_MAT, split = '.cov')[[1]],
#                                    "_PC1_PC3.rds"))
#saveRDS(plot_PC2_PC3, file = paste0(strsplit(x = COV_MAT, split = '.cov')[[1]],
#                                    "_PC2_PC3.rds"))
saveRDS(plot_PC3_PC4, file = paste0(strsplit(x = COV_MAT, split = '.cov')[[1]],
                                    "_PC3_PC4.rds"))

# 5. Format nicer figures -------------------------------------------------
pca_df$POP <- as.factor(pca_df$POP)

library(ggpubr)
plot_PC1_2 <- 
ggplot(data = pca_df, aes(x = round(PC1, 2), y = round(PC2, 2), label = ID)) +
  geom_point(aes(col = POP, shape = POP)) +
  stat_ellipse(linewidth = 0.5, aes(group = POP, col = POP), show.legend = FALSE) + 
  labs(y = paste0("PC2 (", var2, ' %)'), x = paste0("PC1 (", var1, ' %)'),
       color = 'Population', shape = 'Population') 
#+ scale_color_manual(values = c("red", "blue"))

saveRDS(plot_PC1_2, file = paste0(strsplit(x = COV_MAT, split = '.cov')[[1]],
                                    "_PC1_PC2_formatted.rds"))


plot_PC3_4 <- 
ggplot(data = pca_df, aes(x = round(PC3, 2), y = round(PC4, 2), label = ID)) +
  geom_point(aes(col = POP, shape = POP)) +
  stat_ellipse(linewidth = 0.5, aes(group = POP, col = POP), show.legend = FALSE) + 
  labs(x = paste0("PC3 (", var3, ' %)'), y = paste0("PC4 (", var4, ' %)'),
       color = 'Population', shape = 'Population') 
#+ scale_color_manual(values = c("red", "blue"))

saveRDS(plot_PC3_4, file = paste0(strsplit(x = COV_MAT, split = '.cov')[[1]],
                                    "_PC3_PC4_formatted.rds"))
