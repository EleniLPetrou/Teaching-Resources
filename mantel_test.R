# The purpose of this script is to conduct a Mantel test on two matrices
# (distance matrix and fst matrix).

# The Mantel test null hypothesis: H0:
# The distances among objects in matrix Dy are not (linearly or monotonically) 
# related to the corresponding distances in Dx (Legendre & Legendre 2012) 
# 
# Alternative hypothesis: H1: Small values of Dy correspond to small values of 
# Dx and large values of Dx to large values of Dy.

# This H0 is tested using a permutation procedure: If H0 is true, then permuting 
#the rows and columns of the matrix should be equally likely to produce a larger 
# or a smaller coefficient.The test statistic is the Pearson product-moment 
# correlation coefficient r. r falls in the range of -1 to +1, where being close
#  to -1 indicates strong negative correlation and +1 indicates strong positive 
#  correlation. An r value of 0 indicates no correlation.


# Load the necessary libraries
library(ggplot2)
library(tidyr)
library(dplyr)
library(vegan)

######################################################################################
MYPATH <- "~/chum_salmon_martin"
setwd(MYPATH)
list.files()

FSTFILE <- "Fst_matrix.txt"  
DISTFILE <- "distance_matrix.txt"

# Read in data
fst_df <- read.delim(FSTFILE)
head(fst_df)

dist_df <- data.table(read.delim(DISTFILE))
head(dist_df)

# Save the data as matrix objects
dist_mat <- as.matrix(dist_df, rownames = "Population")
fst_mat <- as.matrix(fst_df, rownames = "Population")

rownames(dist_mat)
rownames(fst_mat)

# mantel test on the full dataset

mantel_full <- mantel(dist_mat , 
                      fst_mat, 
                      method = "kendall", 
                      permutations = 999, 
                      na.rm = TRUE)

print(mantel_full)

###############################################################################
# Cluster 1 Mantel test: comparisons with FST > 0.20

cluster1_mat <- fst_mat
cluster1_mat[!cluster1_mat > 0.2] <- NA
which(cluster1_mat > 0.2)


mantel_cluster1 <- mantel(cluster1_mat , 
                      dist_mat, 
                      method = "kendall", 
                      permutations = 999, 
                      na.rm = TRUE)

print(mantel_cluster1)

#Mantel statistic r: -0.01914
#Significance: 0.552 
###############################################################################
# Cluster 2 Mantel test: comparisons with 0.1 < FST < 0.20

cluster2_mat <- fst_mat
cluster2_mat[!cluster2_mat > 0.1 & cluster2_mat < 0.2] <- NA
which(cluster2_mat > 0.1 & cluster2_mat < 0.2)

mantel_cluster2 <- mantel(cluster2_mat , 
                          dist_mat, 
                          method = "kendall", 
                          permutations = 999, 
                          na.rm = TRUE)

print(mantel_cluster2)

#Mantel statistic r: -0.03781 
#Significance: 0.749

################################################################################
# Cluster 3 Mantel test: comparisons with FST < 0.1 

cluster3_mat <- fst_mat
cluster3_mat[!cluster3_mat < 0.1] <- NA


mantel_cluster3 <- mantel(cluster3_mat , 
                          dist_mat, 
                          method = "kendall", 
                          permutations = 999, 
                          na.rm = TRUE)

print(mantel_cluster3)

#Mantel statistic r: 0.1559
#Significance: 0.003  
