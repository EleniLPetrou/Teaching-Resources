# Install necessary packages
if (!require("tibble")) install.packages("tibble")
if (!require("tidyr")) install.packages("tidyr")
if (!require("dplyr")) install.packages("dplyr")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("adegenet")) install.packages("adegenet")
if (!require("viridis")) install.packages("viridis")


# Load the necessary libraries
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(adegenet)
library(viridis)


############################################################
# Specify the path to the genepop file on your computer
MYPATH <- "D:/sequencing_data/Herring_Coastwide_PopulationStructure/output_stacks_populations/filtered_haplotypesANDsnps_1104indiv_6718loci/data_for_Lorenz_class"

# set the working directory
setwd(MYPATH)
list.files()

# Specify the name of your genepop file
INFILE <- "HerringRADGP.gen"
METAFILE <- "Julian_date_metadata.txt"

##############################################################
# PART 1: Read in your data
# read in the sample metadata
meta_df <- read.delim(METAFILE)

# Read in the genepop data with the R package adegenet and save it as an adegenet (genind) object. Specify that alleles are coded with 3 characters (ncode =3)
mydata_gen <- read.genepop(INFILE, ncode = 3) 
class(mydata_gen)

# Do a little bit of data wrangling. Save a vector containing the population names and use that to rename the populations in the adegenet object
pop_vec <- gsub("_[0-9][0-9][0-9]", "", mydata_gen@pop) 
head(pop_vec)
mydata_gen@pop <- as.factor(pop_vec)
levels(mydata_gen@pop)
# OPTIONAL: Subset the data so you only keep loci on Chromosome 1
#LG1_loci <- mydata_gen@loc.fac[1:576]
#mydata_gen <- mydata_gen[loc = LG1_loci]

##############################################################
# PART 2: Conduct a discriminant analysis of principal components
# The DAPC function transforms the data using PCA and then performs a discriminant analysis on the retained principal components.

(mypcas <- length(levels(mydata_gen$loc.fac)))
(mydas <- length(levels(mydata_gen$pop)))

dapc_all <- dapc(mydata_gen, mydata_gen$pop, n.pca = mypcas, n.da = mydas)

# Look at the optim_a score and choose the optimal number of PCS based on the highest value
test_a_score <- optim.a.score(dapc_all)
test_a_score$best

dapc_opt <- dapc(mydata_gen, mydata_gen$pop, n.pca = test_a_score$best, n.da = mydas)
dapc_opt

# Now, save the output of the DAPC as a dataframe. Here I saved all rows, but only the first six discriminant axes. You can access the individual coordinates of a DAPC by typing dapc_all$ind.coord
DAPC_df <- as.data.frame((dapc_opt$ind.coord[ , 1:6]))

# Do a little bit of data wrangling to make a nice dataframe for plotting
temp_df <- DAPC_df %>%
  tibble::rownames_to_column(var = "individual") %>%
  tidyr::separate(individual, c("Population", "id"), sep = "_", remove = FALSE)

plotting_df <- dplyr::left_join(temp_df, meta_df, by = "Population")

##############################################################
# PART 3: Plot the results of the discriminant analysis of principal components

# set the breaks for your color ramp
mybreaks = c(0, 30, 60, 90, 120, 150, 180)
mylabels = c("January", "February", "March", "April", "May", "June", "July")

dapc_plot <- ggplot() +
  geom_point(data = plotting_df, aes(x = LD1, y = LD2, color = julian_date, shape = Region), size = 2, alpha = 0.6) +
  ylab("DA 1") +
  xlab("DA 2") +
  scale_color_viridis(option = "plasma", name = "Sampling date", 
                      breaks = mybreaks,
                      labels = mylabels,
                      begin = 0, end = 0.95) +
  theme_bw() +
  theme(panel.grid = element_blank())

dapc_plot

# Save the plot
ggsave("DAPC_plot.pdf", dapc_plot)
