# Eleni Petrou. 20210426. R version 3.6.1 (2019-07-05)

# Purpose of script: 

# This script takes as input a genepop file and calculates observed vs. expected heterozygosities.
# It also calculates Hardy-Weinberg Equilibrium
############################################################
# Install R packages
if (!require("tibble")) install.packages("tibble")
if (!require("tidyr")) install.packages("tidyr")
if (!require("dplyr")) install.packages("dplyr")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("hierfstat")) install.packages("hierfstat")
if (!require("adegenet")) install.packages("adegenet")
if (!require("pegas")) install.packages("pegas")


# Load required packages
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(hierfstat)
library(adegenet)
library(pegas)

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

##############################################################
# PART 2: Calculate overall observed and expected heterozygosities using adegenet

div <- summary(mydata_gen)
head(div)

# The output is in a funny format, so let's do some data wrangling to save it as a table
div_df <- as.data.frame(cbind(div$Hobs,div$Hexp))

div_df <- div_df %>%
  tibble::rownames_to_column(var = "locus") %>%
  dplyr::rename("Hobs" = V1, "Hexp" = V2)

head(div_df)

div_plot <- ggplot(div_df, aes(x = Hobs, y = Hexp)) +
  geom_point() +
  ggtitle("Expected heterozygosity as a function of observed heterozygosity per locus") +
  theme_bw()

div_plot

##############################################################
# PART 3: Is there a statistically significant difference between observed and expected heterozygosities?
# We will use the Bartlett test of homogeneity of variances to investigate this question
# Our null hypothesis: Hexp = Hobs

bartlett.test(list(div$Hexp, div$Hobs))

##############################################################
# PART 4: Test for Hardy-Weinberg Equilibrium
# We will test each population for deviations from HWE using the function hw.test() from the pegas package.

# overall test- HWE across all sampling locations
HWE_all <- hw.test(mydata_gen, B = 100)

# separate test for each sampling location
hwe.pop <- seppop(mydata_gen) %>% 
  lapply(hw.test, B = 100)

##############################################################
# PART 5: Do some data wrangling to place these HWE results into a nice dataframe

# hwe.pop is a list of matrices
# use do.call function to do r bind to the list of 
# matrices and turn them into a dataframe

hwe_df <- as.data.frame(do.call(rbind, hwe.pop))
hwe_df$locus <- row.names(hwe_df)
head(hwe_df)

# create a temporary vector with the population names
temp_vec <- (names(hwe.pop))
temp_vec

# run a little function to repeat the population names by the number of loci
nloci <- nlevels(mydata_gen$loc.fac)

name_func <- function(x) {
  m <- c()
  for (i in x) {
    y <- (rep(i,nloci))
    m <- c(m, y)
  }
  return(m)
} 

name_vec <- name_func(temp_vec)
name_vec

# append those population names to the dataframe
hwe_df$Location <- name_vec


# pegas has added a funny .pop number after each locus name, to designate that that locus was being evaluated in a specific population. Let's remove this delimiter, since we already have a column in our dataframe designating the population. 
hwe_df$locus_name <- gsub("\\..*","", hwe_df$locus)
head(hwe_df)


# Split the hwe_df into two distinct dataframes, based on the chromosome that a locus is on.
alpha  <- 0.05

final_df <- hwe_df %>%
  dplyr::select(-locus) %>%
  dplyr::rename(chi2_pval = `Pr(chi^2 >)`) %>%
  dplyr::rename(exact_pval = `Pr.exact`) %>%
  tidyr::separate(locus_name, c("chromosome", "position"), sep = "_", remove = FALSE) %>%
  dplyr::mutate(significance = if_else(exact_pval > alpha, 1, 0))



hwe_LG1_df <- final_df %>%
  filter(chromosome == "LG1")

hwe_LG15_df <- final_df %>%
  filter(chromosome == "LG15")


# Chromosome 1 plot  

LG1_plot <- ggplot(data = hwe_LG1_df, aes(exact_pval, fill = Location)) + 
  geom_histogram(binwidth = 0.05) +
  ggtitle("LG1: Distribution of p-values of exact Hardy-Weinberg test in all populations") +
  facet_wrap(~Location) +
  theme_bw() +
  geom_hline(yintercept = 0, color = "white") +
  xlab("p-value of exact Hardy-Weinberg test")

LG1_plot

# Chromosome 15 plot 

LG15_plot <- ggplot(data = hwe_LG15_df, aes(exact_pval, fill = Location)) + 
  geom_histogram(binwidth = 0.05) +
  ggtitle("LG15: Distribution of p-values of exact Hardy-Weinberg test in all populations") +
  facet_wrap(~Location) +
  theme_bw() +
  geom_hline(yintercept = 0, color = "white") +
  xlab("p-value of exact Hardy-Weinberg test")

LG15_plot

# Visualize loci out of HWE on each chromosome
hwe_LG1_sig_df <- filter(hwe_LG1_df, exact_pval < alpha)
hwe_LG15_sig_df <- filter(hwe_LG15_df, exact_pval < alpha)

(hwe_LG1_sig_plot <- ggplot(hwe_LG1_sig_df, aes(x = Location, y = locus_name, fill = exact_pval)) + 
    geom_tile(width = 1, height = 1) +
    ylab("Locus name") +
    theme_bw() +
    theme(axis.text.y = element_text(size = 4)))


(hwe_LG15_sig_plot <- ggplot(hwe_LG15_sig_df, aes(x = Location, y = locus_name, fill = exact_pval)) + 
  geom_tile(width = 1, height = 1) +
  ylab("Locus name") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 4)))


##########################################################
# Save your results to files

# Save the HWE data to your working directory as a text file
write.table(final_df, file = "results_hwe.txt", quote = FALSE, sep = "\t" , row.names = FALSE)

# Save the plots
ggsave("HWE_distro_LG1.pdf", LG1_plot)
ggsave("HWE_distro_LG15.pdf", LG15_plot)
ggsave("HWE_locus_LG1.pdf", hwe_LG1_sig_plot)
ggsave("HWE_locus_LG15.pdf", hwe_LG15_sig_plot)








