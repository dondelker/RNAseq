# Combine, normalize and create log2 ratio data files
# Load libraries
library(DESeq2)
library(pheatmap)

# Read in count files
jcounts <- read.delim("countsFemales.txt", header = TRUE, row.names = 1)

# Remove genes with low reads
jsums <- rowSums(jcounts)
jdf <- subset(jcounts, jsums > 9)

# Create data matrix
m <- as.matrix(jdf)

# Read the sample information and experimental factors
sampleinfo = read.delim("female_info.txt")

# Create separate factor objects for each experiment factor.
Treatment = factor(sampleinfo$Treatment)

# Create a DESeqDataSet object. Make sure that the last factor is the
# desired comparison for differential expression.
dds = DESeqDataSetFromMatrix(countData = m, colData = sampleinfo,
                             design = ~Treatment)

# Estimate size factors and normalize count matrix
dds = estimateSizeFactors(dds)
df_norm <- counts(dds, normalized=TRUE)

# Run DESeq2 and get the results
dds = DESeq(dds)
dds_results = results(dds)

# Save the results to file
write.table(dds_results, file = "DEseq2_all_females.txt", quote = F, 
            sep = "\t", col.names = NA, row.names = TRUE)

# Create log2 ratios heatmaps
df_mean <- rowMeans(df_norm)
df_log2 <- log((df_norm+1)/(df_mean+1),2)

# Set legend break parameters
brks = seq(-1.1, 1.1, 0.02)
lgnd_brks = c(-1, -0.5, 0, 0.5, 1)

# Make heatmap
pheatmap(df_log2, breaks = brks, legend_breaks = lgnd_brks, fontsize = 12, 
         clustering_distance_cols = "correlation", clustering_method = "single",
         show_rownames = FALSE)
