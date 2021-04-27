# Eleni Petrou. 20210426. R version 3.6.1 (2019-07-05)

# Purpose of script: 

# This script takes as input a genepop file. It subsequently allows the user to calculate distributions of FIS, visualize them, and save them to a text document.

############################################################
# Install R packages
if (!require("tibble")) install.packages("tibble")
if (!require("tidyr")) install.packages("tidyr")
if (!require("dplyr")) install.packages("dplyr")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("hierfstat")) install.packages("hierfstat")
if (!require("adegenet")) install.packages("adegenet")


# Load required packages
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(hierfstat)
library(adegenet)

############################################################
# Specify the path to the genepop file on your computer
MYPATH <- "D:/sequencing_data/Herring_Coastwide_PopulationStructure/output_stacks_populations/filtered_haplotypesANDsnps_1104indiv_6718loci/data_for_Lorenz_class"

# set the working directory
setwd(MYPATH)
list.files()

# Specify the name of your genepop file
INFILE <- "HerringRADGP.gen"

##############################################################
# PART 1: Read in your data and calculate Ho, Hs, and Fis

# Read in the genepop data with the R package adegenet and save it as an adegenet object. Specify that alleles are coded with 3 characters (ncode =3)
mydata_gen <- read.genepop(INFILE, ncode = 3) 


# Do a little bit of data wrangling. Save a vector containing the population names and use that to rename the populations in the adegenet object
pop_vec <- gsub("_[0-9][0-9][0-9]", "", mydata_gen@pop) 
head(pop_vec)
mydata_gen@pop <- as.factor(pop_vec)


# Using the R package hierfstat, calculate observed heterozygosities, gene diversities, and FIS for each sampling location
my_stats <- basic.stats(mydata_gen)

# Take a peek at some of the output
# Results over all sampling locations
my_stats$overall

# Take a peek at results reported for seperately for each sampling location
head(my_stats$Hs)
head(my_stats$Ho)
head(my_stats$Fis)


# Do a little bit more data wrangling to get the FIS results into a user-friendly format
# Use the tibble package to save the rownames of this dataframe to a column
fis_df <- as.data.frame(my_stats$Fis) %>%
  tibble::rownames_to_column( var = "locus")
  
head(fis_df)

# Use tidyr library to change the format of the dataframe (from wide to long) so you can plot the data easily with ggplot2. The "gather" function Gather multiple columns and collapses them into key-value pairs. The data is now in "tidy" format
fis_tidy_df <- tidyr::gather(fis_df, "location", "fis", pop_vec)
head(fis_tidy_df)


###########################################
# PART 2: Plot the distribution of FIS in each sampling location using ggplot

fis_plot <- ggplot(fis_tidy_df, aes(x = fis, fill = location)) +
  geom_histogram(binwidth = 0.15) +
  facet_wrap(~location) +
  xlab(expression(italic(F[IS]))) +
  ylab("Number of loci") +
  xlim(-1,1) +
  geom_hline(yintercept = 0, color = "white") +
  theme_bw() 

fis_plot

# Save the FIS data to your working directory as a text file
write.table(fis_tidy_df, file = "HerringRADGP_FIS_results.txt", quote = FALSE, sep = "\t" , row.names = FALSE)

# Save the plot
ggsave("FIS_plot.pdf", fis_plot)

