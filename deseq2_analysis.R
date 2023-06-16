# Run DESeq multivariate (multi-factor) analysis

# Load packages.
library(DESeq2)

# Read in count data file and make a matrix of the counts
mydata <- read.delim("Teratogen26_Counts_nodups.txt", header = TRUE, row.names = 1)
m <- as.matrix(mydata) 

# Read the sample information and experimental factors
sampleinfo = read.delim("sample_info.txt")

# Create separate factor objects for each experiment factor.
Vehicle = factor(sampleinfo$Vehicle)
Treatment = factor(sampleinfo$Treatment)

# Create a DESeqDataSet object. Make sure that the last factor is the
# desired comparison for differential expression.
dds = DESeqDataSetFromMatrix(countData = m, colData = sampleinfo, 
                             design = ~Vehicle+Treatment)

# Run DESeq2 and get the results
dds = DESeq(dds)
dds_results = results(dds)

# Save the results to a file
write.table(dds_results, file = "Teratogen26_results.txt", quote = F, 
            sep = "\t", col.names = NA, row.names = TRUE)