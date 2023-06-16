# Make heatmaps and volcano plots from gene expression analysis

# Load packages
library(ggplot2)
library(pheatmap)

# Read in DEseq2 normalized data and calculate row means and log2 ratios
df_norm <- read.delim("norm_allsamples.txt", header = TRUE, row.names = 1)

# Separate out male samples
df_male <- df_norm[c(1:18)]
male_mean <- rowMeans(df_male[1:6])
male_log2 <- log((df_male+1)/(male_mean+1),2)

# Separate out female samples
df_female <- df_norm[c(19:36)]
female_mean <- rowMeans(df_female[1:6])
female_log2 <- log((df_female+1)/(female_mean+1),2)

# Remove outlier samples
# df_log2 <- df_log2[c(1:5,9:18)]

# Get the log2 ratios and FDR from DESeq2 results
femalegenes <- read.delim("DESeq2_high_females.txt", header = TRUE, row.names = 1)
malegenes <- read.delim("DESeq2_high_males.txt", header = TRUE, row.names = 1)
femalegenes <- femalegenes[c(2,6)]
malegenes <- malegenes[c(2,6)]
  
# Select genes of interest from DESeq2 results
sigfemale <- subset(femalegenes, femalegenes$padj < 0.05, select = c(1:2))
sigmale <- subset(malegenes, malegenes$padj < 0.05, select = c(1:2))

# mygenes <- merge(siglow,sighigh, by=0, all=TRUE)
# Retrieve nomalized data for significant genes
female_list <- as.vector(rownames(sigfemale))
femalematrix <- female_log2[rownames(female_log2) %in% female_list, ]
male_list <- as.vector(rownames(sigmale))
malematrix <- male_log2[rownames(male_log2) %in% male_list, ]

# Sort matrix high to low
mymean <- rowMeans(malematrix[13-18])
mymatrix <- malematrix[order(mymean, decreasing = TRUE),]

# Set legend break parameters
brks = seq(-1.1, 1.1, 0.0275)
lgnd_brks = c(-1, -0.5, 0, 0.5, 1)

# Make heatmap
pheatmap(mymatrix, breaks = brks, legend_breaks = lgnd_brks, fontsize = 10, 
         clustering_distance_cols = "euclidean", clustering_method = "complete",
         cluster_row = FALSE, show_rownames = FALSE)

# Create dataframe for volcano plot
malegenes$expression <- as.factor(malegenes$padj < 0.1 & abs(malegenes$log2FoldChange) > 0.585)

# Create Volcano plot
ggplot(data = malegenes, aes(x=log2FoldChange, y=-log10(padj), colour=expression)) +
  geom_point(alpha=0.4, size=1.8) + xlab("Log2 Fold Change") +
  ylab("Adjusted p-value") + 
  theme(axis.title.x = element_text(face = "bold", size = 12),
  axis.text.x = element_text(face = "bold", size = 12)) +
  theme(axis.title.y = element_text(face = "bold", size = 12),
  axis.text.y = element_text(face = "bold", size = 12)) +
  scale_colour_discrete(name = "Fold > 1.5, FDR < 0.1") +
  xlim(c(-3, 3)) + 
  ylim(c(0, 10)) +
  geom_vline(xintercept=c(-0.585,0.585),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = 1,lty=4,col="black",lwd=0.8) +
  theme(legend.title = element_text(face = "bold", size = 12)) +
  theme(legend.text = element_text(size = 10))