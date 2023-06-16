# Transform DESeq2 normalized data to log2 ratios

# Read in count matrix
df <- read.delim('DESeq2_norm_counts_nodups.txt', header = TRUE, row.names = 1)

# Force dataframe as numeric
df[] <- lapply(df, as.numeric)

# Calculate average expression of all genes
dfmean <- rowMeans(df)

# Calculate log2 ratios
dfratio <- log((df+1)/(dfmean+1),2)

# Add gene names back to data frame
genes <- rownames(df)
mydata <- cbind(genes, dfratio)

# Save transformed matrix
write.table(mydata, "Norm_log2ratios.txt", row.names = FALSE, sep="\t")
