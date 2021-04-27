# This script was written by Eleni on 20210413
# It calculates Weir and Cockerham pairwise population FST (1984) using the R package hierfstat
# and a genepop file as input

# Install necessary packages
if (!require("tibble")) install.packages("tibble")
if (!require("tidyr")) install.packages("tidyr")
if (!require("dplyr")) install.packages("dplyr")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("hierfstat")) install.packages("hierfstat")
if (!require("adegenet")) install.packages("adegenet")
if (!require("reshape2")) install.packages("reshape2")
if (!require("diveRsity")) install.packages("diveRsity")

# Load the necessary libraries
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(adegenet)
library(hierfstat)
library(reshape2)
library(diveRsity)


############################################################
# Specify the path to the genepop file on your computer
MYPATH <- "D:/sequencing_data/Herring_Coastwide_PopulationStructure/output_stacks_populations/filtered_haplotypesANDsnps_1104indiv_6718loci/data_for_Lorenz_class"

# set the working directory
setwd(MYPATH)
list.files()

# Specify the name of your genepop file
INFILE <- "HerringRADGP.gen"

##############################################################
# PART 1: Read in your data

# Read in the genepop data with the R package adegenet and save it as an adegenet (genind) object. Specify that alleles are coded with 3 characters (ncode =3)
mydata_gen <- read.genepop(INFILE, ncode = 3) 
class(mydata_gen)

# Do a little bit of data wrangling. Save a vector containing the population names and use that to rename the populations in the adegenet object
pop_vec <- gsub("_[0-9][0-9][0-9]", "", mydata_gen@pop) 
head(pop_vec)
mydata_gen@pop <- as.factor(pop_vec)

# Subset the data so you only keep loci on Chromosome 1

LG1_loci <- mydata_gen@loc.fac[1:576]

mydata_gen_LG1 <- mydata_gen[loc = LG1_loci]

# Transform the genepop file into a hierfstat data frame
hierfstat_df <- genind2hierfstat(mydata_gen_LG1)


##############################################################
# PART 2: Using the R package hierfstat, calculate pairwise population FST using the equation from Weir and Cockerham (1984). This takes about 5-10 minutes to run.

fst_mat <- pairwise.WCfst(hierfstat_df, diploid = TRUE)
head(fst_mat)


##############################################################
# PART 3: Do a little data wrangling so you can plot the data

# Round the Fst values to 4 decimal places
fst_mat <- round(fst_mat, 4)

# Specify a custom order for the rows and columns of the matrix. This order is roughly by spawn-timing
my_order <- c("Case07", "Squa14", "Port14", 
              "SmBy15", "PtGb14", "QlBy14", 
              "Gabr15", "ElBy15", "Kwak15",
              "ChPt14", "ChPt16", "Skid14", "Mass16")

fst_mat_ordered <- fst_mat[my_order, my_order]

# Note that a correlation matrix has redundant information (the upper triangle is the same as the lower triangle). We'll use the function below to keep only half of the correlation matrix.

get_upper_tri <- function(Fstmat){
  Fstmat[lower.tri(Fstmat)] <- NA
  return(Fstmat)
}

upper_tri <- get_upper_tri(fst_mat_ordered)
head(upper_tri)

# Use the melt function from the reshape2 package to reformat the matrix so it is a long table
melted_Fstmat <- melt(upper_tri, na.rm = TRUE)
head(melted_Fstmat)



########################################################################
#PART 4: Make a heatmap and visualize the FST values

# identify the maximum and minimum value of the heat map
mymax <-  max(melted_Fstmat$value)
mymin <- min(melted_Fstmat$value)
mymedian <- median(melted_Fstmat$value)


heatmap_plot <- ggplot(data = melted_Fstmat, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = mymedian, limit = c(mymin, mymax), space = "Lab", 
                       name = "Pairwise Fst") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 8, hjust = 1),
        axis.text.y = element_text(angle = 0, vjust = 1, size = 8, hjust = 1)) +
  ylab("Population A") +
  xlab("Population B") +
  coord_fixed()

heatmap_plot

########################################################################
#PART 5: Use the diversity package to estimate statistical significance of pairwise FST values

test_results <- chiCalc(infile = INFILE, outfile = NULL, pairwise = TRUE, mcRep = 5000)

# do alittle bit of data wrangling to rename the populations and make a plottable dataframe
sig_df <- test_results$multilocus_pw %>%
  tidyr::separate(pops, c("pop1", "pop2"), sep = " vs ", remove = FALSE)

sig_df$pop1 <- gsub("*_[0-9][0-9][0-9],", "", sig_df$pop1)
sig_df$pop2 <- gsub("*_[0-9][0-9][0-9],", "", sig_df$pop2)

head(sig_df)

sigmax <-  max(sig_df$p.value)
sigmin <- min(sig_df$p.value)
sigmedian <- median(sig_df$p.value)


sig_plot <- ggplot(data = sig_df, aes(pop1, pop2, fill = p.value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "red", 
                       midpoint = 0.90, limit = c(sigmin, sigmax), space = "Lab", 
                       name = "p value") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 8, hjust = 1),
        axis.text.y = element_text(angle = 0, vjust = 1, size = 8, hjust = 1)) +
  ylab("Population A") +
  xlab("Population B") +
  coord_fixed()

sig_plot


###########################################################################

# save results to file
write.table(upper_tri, file = "results_Fst_matrix.txt", sep = "\t", quote = FALSE)
write.table(sig_df, file = "results_Fst_significance.txt", sep = "\t", quote = FALSE)

# Save the plot
ggsave("pairwise_FST_plot.pdf", heatmap_plot)

