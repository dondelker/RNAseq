# Combine, normalize and create log2 ratio data files
# Load libraries
library(rgl)

# Read in PCA data
mydata <- read.delim("PCA_pup_results.txt", header = TRUE, row.names = 1)

# Create matrix for PCA
mymatrix <- t(mydata[26:37])
row.names(mymatrix) <- gsub("Interim.", "", row.names(mymatrix))
my3Dvalues <- mymatrix[,1:3]

# Color samples by group
mycol <-rep("black", 12)
mycol[c(1:6)]<-"blue"
mycol[c(7:12)]<-"red"

# Make 3D PCA plot and print to file
plot3d(my3Dvalues, xlab="", ylab="", zlab="", type = 's', size=3, col=mycol)
x=my3Dvalues[,1]
y=my3Dvalues[,2]
z=my3Dvalues[,3]
text3d(x=x,y=y,z=z,text=rownames(my3Dvalues), adj=1.4)
rgl.snapshot("interim_PCA3.png")

