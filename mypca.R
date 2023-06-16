# Perform PCA on normalized count data
# Load libraries
library(tidyverse)
library(rgl)

# Read in normalized count table
mydata <- read.delim("Norm_log2ratios.txt", header = TRUE, row.names = 1)

# Perform PCA
mydata.pca <- prcomp(mydata, center = TRUE)
summary(mydata.pca)

# Get eigenvalues and eigenvectors
eigenvalues <- mydata.pca$sdev^2
eigenvectors <- as.data.frame(mydata.pca$rotation)

# Write combined results to file
myresults <- t(rbind(eigenvalues, eigenvectors))
write.table(myresults, "PCA_results.txt", row.names = TRUE, col.names = NA, 
            sep = "\t")

# Make plot of explained variance
cumpro <- cumsum(eigenvalues / sum(eigenvalues))
plot(cumpro[0:8], ylim = c(0, 1), xlab = "PC #", col = "red", pch = 16, 
     cex = 2, cex.main = 2, cex.lab = 2, cex.axis = 2, 
     ylab = "Amount of explained variance", main = "Cumulative variance plot")
abline(v = 3, col="blue", lty=5)
abline(h = 0.76, col="blue", lty=5)
legend("topleft", legend=c("Cut-off @ PC3"),
       col=c("blue"), lty=2, cex=1)

# Make table of genes showing the most variability
mytable <- as.data.frame(mydata.pca$x)
mytable <- mytable %>% 
  as_tibble(rownames = "gene")

# Select gene and PC2 variables
top_genes <- mytable %>% 
  select(gene, PC2) %>%
  # convert to a long format 
  pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
  # for each PC
  group_by(PC) %>% 
  # arrange by descending order of loading
  arrange(desc(abs(loading))) %>% 
  # take the 50 top rows
  slice(1:50) %>% 
  # pull the gene column as a vector
  pull(gene) %>% 
  # ensure only unique genes are retained
  unique()

# Subset the eigenvector table to top 50 genes
top_table <- mytable %>%
  filter (gene %in% top_genes)
write.table(top_table, "mytop50genes.txt", row.names = TRUE, col.names = NA, 
            sep = "\t")

# Create matrix for PCA
my3Dvalues <- eigenvectors[, 1:3]

# Color samples by group
mycol <-rep("black", 18)
mycol[c(1:3)]<-"blue"
mycol[c(4:6)]<-"green"
mycol[c(7:9)]<-"orange"
mycol[c(10:12)]<-"black"
mycol[c(13:15)]<-"red"
mycol[c(16:18)]<-"cyan"

# Make 3D PCA plot and print to file
plot3d(my3Dvalues, xlab="", ylab="", zlab="", size=8, col=mycol)
rgl.snapshot("PCA.png")

# Make fast 2D PCA plot of all genes
library(ggfortify)
autoplot(mydata.pca)