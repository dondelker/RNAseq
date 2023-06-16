# Merge DEseq2 results files

# Read in data files
male_low <- read.delim("DEseq2_low_males.txt", header = TRUE, row.names = 1)
male_high <- read.delim("DEseq2_high_males.txt", header = TRUE, row.names = 1)
female_low <- read.delim("DESeq2_low_females.txt", header = TRUE, row.names = 1)
female_high <- read.delim("DESeq2_high_females.txt", header = TRUE, row.names = 1)

# Subset data files and create df list (repeat for females)
mlow <- male_low[,c(2,6)]
mhigh <- male_high[,c(2,6)]

# Merge low and high dose files (repeat for females)
df_m <- merge(mlow, mhigh, by=0, all=TRUE) 
row.names(df_m) = df_m$Row.names
df_m <- df_m[-c(1)]
colnames(df_m) <- c("Mlow_log2fc","Mlow_FDR","Mhigh_log2fc","Mhigh_FDR")

# Merge male and female files
final_df <- merge(df_f, df_m, by=0, all=TRUE)

# Save the results to a file
write.table(final_df, "DEseq2_exp3_results.txt", sep = "\t", col.names = NA)


