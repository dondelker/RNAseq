# Load libraries
library(BioVenn)

# Read in gene lists
mylists <-read.delim("Chemo_Lists.txt", header = TRUE)
mydrug01 <- mylists$Chemo01
mydrug031 <- mylists$Chemo031

# Create two-way Venn
draw.venn(list_x = mydrug01, list_y = mydrug031, list_z = 0, title = "", subtitle = "", xtitle = "", ytitle = "", ztitle = "", nr_s = 0, z_c = "white", nr_c = "clear")