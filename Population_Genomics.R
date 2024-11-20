#Population Genomics Analyses for Epidendrum fulgens adaptation genomics manuscript#
#J. Mattos, 2024

#Install packages#
install.packages("vcfR")
install.packages("tidyverse")
install.packages("reshape2")
install.packages("RColorBrewer")
install.packages("patchwork")
install.packages("ggrepel")
install.packages("poppr")
install.packages("maps")
install.packages("mapplots")
install.packages("adegenet")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ggtree")

#Load packages#
library(tidyverse)
library(ggplot2)
library(vcfR)
library(adegenet)
library(RColorBrewer)
library(patchwork)
library(ggrepel)
library(poppr)
library(ggtree)
library(reshape2)
library(maps)
library(mapplots)

#Load VCF file and sample info#

vcf_file <- "~/Documents/Jac/Population Genomics/VCF/vcf_pruned_by_plink.vcf"
vcf <- read.vcfR(vcf_file, verbose = FALSE)
sample_info <- read.csv2("~/Documents/Jac/Population Genomics/Tabelas/populations_table_efulgens_POPMAP.csv", header=TRUE)
pop_info <- read.csv2("~/Documents/Jac/Population Genomics/Tabelas/Samples_and_Populations-Habitats.csv", header=TRUE)


#Analysis packages tend to convert files to their own formats so we will use can interpret easily.
#We will use a “genlight” format, which is good for storing variant call data. 
vcf.gl <- vcfR2genlight(vcf)

#In order to use the function properly, like the case of vcfR2genlight, a lot of parameters need to be specified.
pop(vcf.gl) <- sample_info$population
ploidy(vcf.gl) <- 2


vcf.gl
head(vcf.gl@ind.names)
head(vcf.gl@pop)

vcf.pca <- glPca(vcf.gl)
vcf.pca

#convert scores of vcf.pca into a tibble
vcf.pca.scores <- as_tibble(vcf.pca$scores)

#add the population data into a column of vcf.pca.scores tibble
vcf.pca.scores$population <- sample_info$population

# We will also determine the variance each PC contributes the data, which will help us understand potential drivers of patterns in our dataset. Lets plot the eigenvectors to try an understand this a bit more.
barplot(100 * vcf.pca$eig / sum(vcf.pca$eig), col="green")
title(ylab = "Percent of variance explained") 
title(xlab = "Eigenvalues")

#Lets extract the variance associated with the top 4 PCs, so we can use them in our plots.

#first we sum all the eigenvalues
eig.total <- sum(vcf.pca$eig)

#sum the variance
PC1.variance <- formatC(head(vcf.pca$eig)[1]/eig.total * 100)
PC2.variance <- formatC(head(vcf.pca$eig)[2]/eig.total * 100)
PC3.variance <- formatC(head(vcf.pca$eig)[3]/eig.total * 100)
PC4.variance <- formatC(head(vcf.pca$eig)[4]/eig.total * 100)

#Lets check that this has worked
PC1.variance  #10.45 of the variance is explained by the first PC
PC2.variance  #3.46 of the variance is explained by the second PC

#plot PC1 and PC2 with each sample as a point
plot12 <- ggplot(vcf.pca.scores, aes(PC1, PC2)) + geom_point()
plot12

#We’ll add some axis labels, and incorporate the variance information to describe  the relative importance of the spread of the data
#paste0() essentially connects the strings together
xlabel <- paste0("PC1 variance = ",PC1.variance,"%")
ylabel <- paste0("PC2 variance = ", PC2.variance, "%")
plot12 <- plot12 + labs(x = xlabel, y = ylabel)
plot12

#Adding colors and plotting PCA
cols <- colorRampPalette(brewer.pal(8, "Set1"))(17)

#Colorblind-friendly color palette
# The palette with black:
cbp2 <- c("#999933", "#E69F00", "#56B4E9", "#44AA99",
          "#661100", "#332288", "#CC6677", "#CC79A7")

install.packages("ggthemes")
library(ggthemes)

plot12 <- plot12 + 
  geom_point(aes(col = population)) + 
  scale_colour_manual(values=cbp2, name="Population") +
  theme_clean() +
  #stat_ellipse(aes(group = population), linetype = 2)
plot12

#Plotting everything all at once with ggplot2

ggplot(vcf.pca.scores, aes(PC1, PC2))+
  geom_point(aes(col= population))+
  scale_color_manual(values=cbp2, name="Population")+
  theme_clean()+
  labs(x=xlabel, y=ylabel)+
  stat_ellipse(aes(group = population), linetype = 2)
  




#Lets quickly look at PC3/PC4, and compare to the first plot.
plot34 <- ggplot(vcf.pca.scores, aes(PC3, PC4)) + 
  geom_point(aes(col = population)) + 
  labs(x = paste0("PC3 variance = ", PC3.variance,"%"), y = paste0("PC4 variance = ", PC4.variance, "%")) + 
  scale_colour_manual(values = cols) 

plot12 + plot34


#Exploring the genetic data using phylogenetic trees#
#Analysis of pairwise genetic distance using a tree#

#Generating pairwise distances between samples that we will plot in a tree format
tree_data <- aboot(vcf.gl, tree = "upgma", distance = bitwise.dist, sample = 100, showtree = F, cutoff = 50) 

#--- make and plot the tree 
tree_plot <- ggtree(tree_data) + 
  geom_tiplab(size = 2, color = cbp2[pop(vcf.gl)]) + 
  xlim(-0.1, 0.3) + 
  geom_nodelab(size = 2, nudge_x = -0.006, nudge_y = 1) + 
  theme_tree2(legend.position = 'centre')

tree_plot

#______________________________________________________________________________________________________
#Tutorial: https://popgen.nescent.org/StartSNP.html
#Calculating population genomics/genetics basic statistics from SNP data

#Install and load packages
install.packages("hierfstat")
library("adegenet")
library("hierfstat")
library("pegas")

data <- read.table("C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics_RNAseq/Final_VCF/PLINK_LD/vcf_pruned_by_plink.vcf")
dim(data) #matrix of dimensions 75973 x 89

#Transforming the table we imported to a "genind" object, from adegenet#
colnames(data) <- gsub("\\.", "_", colnames(data))




#_______________________________________________________________________________________________________
#### Population structure with LEA ####

#Load the library
library(LEA)

#Convert vcf to geno 

LEA::vcf2geno("C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics/Final_VCF/PLINK_LD/vcf_pruned_by_plink.vcf",
              output.file = "vcf_pruned_by_plink.geno")


#### Model Ancestry Proportions ####

#Now that we have our input file in the expected format, we will estimate individual admixture coefficients using sNMF. This program is implemented as the snmf() function in the R package LEA. In short, this function provides results very similar to programs such as STRUCTURE or ADMIXTURE. Assuming K ancestral populations, the snmf function provides least-squares estimates of ancestry proportions rather than maximum likelihood estimates (Frichot 2014).

#The results allow us to determine what is the best K value (i.e. the most likely number of genetic clusters)

# testing K populations: from K = 1 to K = 10
obj <- LEA::snmf("vcf_pruned_by_plink.geno", K = 1:10, ploidy = 2,
                 entropy = TRUE, CPU = 4, project = "new")



#The snmf function computes an entropy criterion, which assesses the fit of the statistical model to the data using a cross-validation approach. The entropy criterion can help choosing the number of ancestral populations that best explains the genotypic data.

#_Find the best K value 
# plot cross-entropy
plot(obj, col = "blue4", cex = 1.4, pch = 19) # best is 3 here, lowest value
# choose the best LEA run
best = which.min(cross.entropy(obj, K = 3))

#Here the lowest cross-entropy value is clearly at K = 3 (lowest value), suggesting there are three genetic clusters within the dataset. Often, the plot shows a less clear pattern, and choosing the "knee/elbow/inflection" point is a generally good approach.

#The next step is to display a barplot of the ancestry matrix (also called the Q-matrix).

#Plot ancestry proportions across samples 

LEA::barchart(obj, K=3,run=best,border=NA,space=0,
              col=c("orchid4","darkolivegreen4", "darkturquoise"),
              xlab = "Individuals", ylab = "Ancestry proportions (K=3)", main = "Ancestry matrix (K=3)") -> bp_3
axis(1, at = 1:length(bp_3$order), labels = bp_3$order, las=1, cex.axis = .6)



#Plotting with different K values just to check how they look

LEA::barchart(obj, K=2,run=best,border=NA,space=0,
              col=c("orchid4","darkolivegreen4"),
              xlab = "Individuals", ylab = "Ancestry proportions (K=2)", main = "Ancestry matrix (K=2)") -> bp_2
axis(1, at = 1:length(bp_2$order), labels = bp_2$order, las=1, cex.axis = .6)


#Select the best run for K = 4 clusters 
best = which.min(cross.entropy(project, K = 4)) 
my.colors <- c("orchid4", "lightblue", "olivedrab", "darkturquoise") 
LEA::barchart(obj, K = 4, run = best, border = NA, space = 0, col = my.colors, 
              xlab = "Individuals", ylab = "Ancestry proportions (K=4)", 
              main = "Ancestry matrix (K=4)")-> bp_4
axis(1, at = 1:length(bp_4$order), labels = bp_4$order, las=1, cex.axis = .6)




#### Admixture plots for population structure visualization ####

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

### K = 3 ###

# Set parameters
K <- 3  # Specify the value of K
Q_file <- paste0("~/Documents/Jac/Population Genomics/Results/Admixture/pruneddata.", K, ".Q")  # File name with ancestral proportions

# Load ADMIXTURE output
Q_data <- read.table(Q_file, header = FALSE)

# Rename columns to indicate ancestry proportions
colnames(Q_data) <- paste0("Ancestry_", 1:K)

# Add individual labels (e.g., 1, 2, 3...) if desired
Q_data$Individual <- factor(1:nrow(Q_data))

# Reshape the data for plotting
Q_long <- Q_data %>%
  pivot_longer(cols = starts_with("Ancestry_"), 
               names_to = "Ancestry", 
               values_to = "Proportion")


# Plotting according to population 
# Need to use either the .fam or the .ped files that we use as input for Admixture

fam_file <- read.table("~/Documents/Jac/Population Genomics/Results/Admixture/pruneddata.fam", header = FALSE)
individual_ids <- fam_file$V2

#combining the individual IDs from fam file with the results from Admixture in the .Q file 

Q_data$Individual <- individual_ids

#also combining the population of each individual in this file 
#the file sample_file has the same order of individuals and populations as the order in the .Q file now after joining with the .fam file

Q_data$Population <- sample_file$population
colnames(Q_data) <- c("Ancestry_1", "Ancestry_2", "Ancestry_3", "Individual", "Population")


#Changing the order of the populations on the dataframe so it will follow the gradient 
Q_data$Population <- factor(Q_data$Population, levels = c("Ubatuba", "Bertioga", "Cardoso", "Florianopolis", "Torres", "Itapua", "Arambare", "Pelotas"))
Q_data_ordered <- Q_data %>% arrange(Population)


# Explicitly set 'Individual' as a factor (this step helps ensure correct ordering and axis text behavior)
Q_data_ordered$Individual <- factor(Q_data_ordered$Individual, levels = unique(Q_data_ordered$Individual))

#Checking the levels of the factor columns (for the individual orders)
levels(Q_data_ordered$Individual)

#Reshaping the data again for plotting 
Q_long <- Q_data_ordered %>% pivot_longer(cols = starts_with("Ancestry_"), names_to = "Ancestry", values_to = "Proportion")

#Individual data are now as "levels" of a factor

#colour palette
install.packages("wesanderson")
library(wesanderson)
names(wes_palettes)
wes_palettes


# Plot the data
ggplot(Q_long, aes(x = Individual, y = Proportion, fill = Ancestry)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Ancestry_1" = "#9986A5", "Ancestry_2" = "#ABDDDE", "Ancestry_3" = "#446455"))
  labs(x = "Individual", y = "Ancestry Proportion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        plot.margin = margin(t = 10, r = 10, b = 30, l = 10))




### K = 3 ###

# Set parameters
K <- 3  # Specify the value of K
Q_file <- paste0("~/Documents/Jac/Population Genomics/Results/Admixture/pruneddata.", K, ".Q")  # File name with ancestral proportions

# Load ADMIXTURE output
Q_data <- read.table(Q_file, header = FALSE)

# Rename columns to indicate ancestry proportions
colnames(Q_data) <- paste0("Ancestry_", 1:K)

# Add individual labels (e.g., 1, 2, 3...) if desired
Q_data$Individual <- factor(1:nrow(Q_data))

# Reshape the data for plotting
Q_long <- Q_data %>%
  pivot_longer(cols = starts_with("Ancestry_"), 
               names_to = "Ancestry", 
               values_to = "Proportion")


# Plotting according to population 
# Need to use either the .fam or the .ped files that we use as input for Admixture

fam_file <- read.table("~/Documents/Jac/Population Genomics/Results/Admixture/pruneddata.fam", header = FALSE)
individual_ids <- fam_file$V2

#combining the individual IDs from fam file with the results from Admixture in the .Q file 

Q_data$Individual <- individual_ids

#also combining the population of each individual in this file 
#the file sample_file has the same order of individuals and populations as the order in the .Q file now after joining with the .fam file

Q_data$Population <- sample_file$population
colnames(Q_data) <- c("Ancestry_1", "Ancestry_2", "Ancestry_3", "Individual", "Population")


#Changing the order of the populations on the dataframe so it will follow the gradient 
Q_data$Population <- factor(Q_data$Population, levels = c("Ubatuba", "Bertioga", "Cardoso", "Florianopolis", "Torres", "Itapua", "Arambare", "Pelotas"))
Q_data_ordered <- Q_data %>% arrange(Population)


# Explicitly set 'Individual' as a factor (this step helps ensure correct ordering and axis text behavior)
Q_data_ordered$Individual <- factor(Q_data_ordered$Individual, levels = unique(Q_data_ordered$Individual))

#Checking the levels of the factor columns (for the individual orders)
levels(Q_data_ordered$Individual)

#Reshaping the data again for plotting 
Q_long <- Q_data_ordered %>% pivot_longer(cols = starts_with("Ancestry_"), names_to = "Ancestry", values_to = "Proportion")

#Individual data are now as "levels" of a factor

#colour palette
install.packages("wesanderson")
library(wesanderson)
names(wes_palettes)
wes_palettes


# Plot the data
ggplot(Q_long, aes(x = Population, y = Proportion, fill = Ancestry)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Ancestry_1" = "#9986A5", "Ancestry_2" = "#ABDDDE", "Ancestry_3" = "#446455"))
  labs(x = "Individual", y = "Ancestry Proportion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        plot.margin = margin(t = 10, r = 10, b = 30, l = 10))
  


### K = 2 ###
  
# Set parameters
K <- 2  # Specify the value of K
Q_file <- paste0("~/Documents/Jac/Population Genomics/Results/Admixture/pruneddata.", K, ".Q")  # File name with ancestral proportions
  
# Load ADMIXTURE output
Q_data_2 <- read.table("Q_file", header = FALSE)
  
# Rename columns to indicate ancestry proportions
colnames(Q_data_2) <- paste0("Ancestry_", 1:K)
  
# Plotting according to population 
# Need to use either the .fam or the .ped files that we use as input for Admixture

fam_file <- read.table("~/Documents/Jac/Population Genomics/Results/Admixture/pruneddata.fam", header = FALSE)
individual_ids <- fam_file$V2

#combining the individual IDs from fam file with the results from Admixture in the .Q file 

Q_data_2$Individual <- individual_ids
  
#also combining the population of each individual in this file 
#the file sample_file has the same order of individuals and populations as the order in the .Q file now after joining with the .fam file

Q_data_2$Population <- sample_file$population
colnames(Q_data_2) <- c("Ancestry_1", "Ancestry_2", "Individual", "Population")


#Changing the order of the populations on the dataframe so it will follow the gradient 
Q_data_2$Population <- factor(Q_data_2$Population, levels = c("Ubatuba", "Bertioga", "Cardoso", "Florianopolis", "Torres", "Itapua", "Arambare", "Pelotas"))
Q_data_2_ordered <- Q_data_2 %>% arrange(Population)


# Explicitly set 'Individual' as a factor (this step helps ensure correct ordering and axis text behavior)
Q_data_2_ordered$Individual <- factor(Q_data_2_ordered$Individual, levels = unique(Q_data_2_ordered$Individual))

#Reshaping the data again for plotting 
Q_long_2 <- Q_data_2_ordered %>% pivot_longer(cols = starts_with("Ancestry_"), names_to = "Ancestry", values_to = "Proportion")


# Plot the data
ggplot(Q_long_2, aes(x = Individual, y = Proportion, fill = Ancestry)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Ancestry_1" = "#9986A5", "Ancestry_2" = "#446455"))
labs(x = "Individual", y = "Ancestry Proportion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        plot.margin = margin(t = 10, r = 10, b = 30, l = 10))


### K = 4 ###

# Set parameters
K <- 4  # Specify the value of K
#Q_file <- paste0("~/Documents/Jac/Population Genomics/Results/Admixture/pruneddata.", K, ".Q")  # File name with ancestral proportions

# Load ADMIXTURE output
Q_data_4 <- read.table("~/Documents/Jac/Population Genomics/Results/Admixture/pruneddata.4.Q", header = FALSE)

# Rename columns to indicate ancestry proportions
colnames(Q_data_4) <- paste0("Ancestry_", 1:K)

# Plotting according to population 
# Need to use either the .fam or the .ped files that we use as input for Admixture

fam_file <- read.table("~/Documents/Jac/Population Genomics/Results/Admixture/pruneddata.fam", header = FALSE)
individual_ids <- fam_file$V2

#combining the individual IDs from fam file with the results from Admixture in the .Q file 

Q_data_4$Individual <- individual_ids

#also combining the population of each individual in this file 
#the file sample_file has the same order of individuals and populations as the order in the .Q file now after joining with the .fam file

Q_data_4$Population <- sample_file$population
colnames(Q_data_4) <- c("Ancestry_1", "Ancestry_2", "Ancestry_3", "Ancestry_4", "Individual", "Population")


#Changing the order of the populations on the dataframe so it will follow the gradient 
Q_data_4$Population <- factor(Q_data_4$Population, levels = c("Ubatuba", "Bertioga", "Cardoso", "Florianopolis", "Torres", "Itapua", "Arambare", "Pelotas"))
Q_data_4_ordered <- Q_data_4 %>% arrange(Population)


# Explicitly set 'Individual' as a factor (this step helps ensure correct ordering and axis text behavior)
Q_data_4_ordered$Individual <- factor(Q_data_4_ordered$Individual, levels = unique(Q_data_4_ordered$Individual))

#Reshaping the data again for plotting 
Q_long_4 <- Q_data_4_ordered %>% pivot_longer(cols = starts_with("Ancestry_"), names_to = "Ancestry", values_to = "Proportion")


# Plot the data
ggplot(Q_long_4, aes(x = Population, y = Proportion, fill = Ancestry)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Ancestry_1" = "#9986A5", "Ancestry_2" = "#446455", "Ancestry_3" = "#ABDDDE", "Ancestry_4" = "#446499"))
labs(x = "Individual", y = "Ancestry Proportion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        plot.margin = margin(t = 10, r = 10, b = 30, l = 10))




#----------------------------------------------------#


#### Population differentiation using Fst statistics ####

#Load the library
library(hierfstat)
library(adegenet)
library(vcfR)


#loading vcf file with vcfR package
vcf_file <- "C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics/Final_VCF/PLINK_LD/vcf_pruned_by_plink.vcf"
vcf <- read.vcfR(vcf_file, verbose = FALSE)

#loading population info data
pop_info <- read.csv("C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics/Sample_info/trinity_efulgens_table.csv", header=TRUE) 
head(pop_info)
ind <- as.character(pop_info$sample)
pop <- as.character(pop_info$population)


#Transforming the vcf we imported to a "genind" object, from adegenet#
my_geneIND <- vcfR2genind(vcf)
class(my_geneIND)
my_geneIND


#Transforming the genIND object into a hierfstat object
mydata <- genind2hierfstat(my_geneIND, pop = pop) 
class(mydata)


#Getting basic stats from hierfstat object
basic.stats(mydata) # Fst following Nei (1987) on genind object
wc(mydata) # Weir and Cockerham's estimate

#Pairwise Fst
genet.dist(mydata, method = "WC84")




#----------------------------------------------------#

#### Fst visualization with Manhattan plots ####
#It needs the file .fst generated by VCFtools on the server
#This is a measure of "per-site" Fst


# load package
library(ggplot2)
library(dplyr)
library(tidyverse)
library(stats)

# load FST file
Fst <- read.table("populations_efulgens_plinkpruned_windows100kb_Fst.windowed.weir.fst", header = TRUE)
# calculate the middle point of each bin
Fst$midPos <- (Fst$BIN_START + Fst$BIN_END) / 2
# explore the object
head(Fst)

#selecting only the chromosomes/ 12 largest scaffolds 
Fst_chromosomes <- Fst[1:7399,]


# plot
ggplot(data = Fst_chromosomes, aes(x = midPos / 1000000, y = WEIGHTED_FST, col = CHROM)) +
  geom_point() + 
  geom_smooth() + 
  theme_classic() +
  facet_grid(cols = vars(CHROM), scales = "free_x", space = "free_x") +
  labs(x = "position (in MB)")


#### Identifying outliers of differentiation from the Fst file#### 


# identify the 95% and 99% percentile
quantile(Fst_chromosomes$WEIGHTED_FST, c(0.975, 0.995), na.rm = T)

# identify the 95% percentile
my_threshold <- quantile(Fst_chromosomes$WEIGHTED_FST, 0.975, na.rm = T)


# make an outlier column in the data.frame
Fst_chromosomes <- Fst_chromosomes %>% mutate(outlier = ifelse(WEIGHTED_FST > my_threshold, "outlier", "background"))


#how many outlier SNPs do we have compared to the background ones?
Fst_chromosomes %>% group_by(outlier) %>% tally()


#plotting the outliers
ggplot(Fst_chromosomes, aes(midPos, WEIGHTED_FST, colour = outlier)) + geom_point()

ggplot(data = Fst_chromosomes, aes(x = midPos / 1000000, y = WEIGHTED_FST, col = outlier)) +
  geom_point() +  
  theme_classic() +
  facet_grid(cols = vars(CHROM), scales = "free_x", space = "free_x") +
  labs(x = "position (in MB)")




#----------------------------------------------------#


#### Basic Population Genetic Statistics from SNPs data ####

# Here we'll calculate 
# 1. Genetic Diversity 
# 2. Test Hardy Weinberg Equilibrium 
# 3. Fis and global Fst 

# load packages 

library("adegenet")
library("hierfstat")
library("pegas")
library("vcfR")

# Creating a GenIND object and converting to a hierfstat object 


#loading vcf file with vcfR package
vcf_file <- "C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics/Final_VCF/PLINK_LD/vcf_pruned_by_plink.vcf"
vcf <- read.vcfR(vcf_file, verbose = FALSE)

#loading population info data
pop_info <- read.csv("C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics/Sample_info/trinity_efulgens_table.csv", header=TRUE) 
head(pop_info)
ind <- as.character(pop_info$sample)
pop <- as.character(pop_info$population)


#Transforming the vcf we imported to a "genind" object, from adegenet#
my_geneIND <- vcfR2genind(vcf)
class(my_geneIND)
my_geneIND


#Transforming the genIND object into a hierfstat object
mydata <- genind2hierfstat(my_geneIND, pop = pop) 
class(mydata)

#### Genetic diversity - Observed and Expected Heterozygosity per locus #### 

#Using the genIND object to create a summary with Adegenet
div <- adegenet::summary(my_geneIND)
names(div)


#Plotting Observed Heterozygosity 

plot(div$Hobs, xlab="Loci number", ylab="Observed Heterozygosity", 
     main="Observed heterozygosity per locus")


#Plotting Expected Heterozygosity 

plot(div$Hobs,div$Hexp, xlab="Hobs", ylab="Hexp", 
     main="Expected heterozygosity as a function of observed heterozygosity per locus")


##Bartlett test of homogeneity of variances

bartlett.test(list(div$Hexp, div$Hobs)) #H0: Hexp = Hobs

#Bartlett K-squared = 312.42, df = 1, p-value < 2.2e-16
#We see a difference in observed and expected heterozygosity 


#Basic stats (Ho, Hs, Fis and global Fst)

basicstats <- basic.stats(mydata, diploid = TRUE, digits = 2) 
names(basicstats)

#Compute confidence interval for FIs
boot.ppfis(mydata)

#$fis.ci
#ll     hl
#1 0.0807 0.0882
#2 0.1547 0.1627
#3 0.1015 0.1078
#4 0.0797 0.0864
#5 0.0753 0.0815
#6 0.1082 0.1153
#7 0.1361 0.1433
#8 0.1155 0.1226


#### Testing for Hardy-Weinberg Equilibrium ####

#We used the function hw.test() from the pegas package.

hw.test(my_geneIND, B=1000)

#With this, we get for each locus a test of significance of the null hypothesis: H0 - the locus is in HW equilibrium in the population.




#----------------------------------------------------#


#### Per population Basic Summary Stats ####

#Getting basic stats from hierfstat object
basic.stats(mydata) # Fst following Nei (1987) on genind object
wc(mydata) # Weir and Cockerham's estimate

#A table –with np (number of populations) columns and nl (number of loci) rows– of genotype counts
n.ind.samp(mydata)

#A list containing allele frequencies. Each element of the list is one locus. For each locus, Populations are in columns and alleles in rows
pop.freq(mydata)

#	A table –with np (number of populations) columns and nl (number of loci) rows– of observed heterozygosities
Ho(mydata)

#A table –with np (number of populations) columns and nl (number of loci) rows– of observed gene diversities
Hs(mydata)

#	A table –with np (number of populations) columns and nl (number of loci) rows–of observed Fis
Fis(mydata)




#----------------------------------------------------#


#### Genetic Diversity metrics from VCFtools ####

#Nucleotide Diversity#

nucleotide_diversity <- read.table("efulgens_10kb.windowed.pi", header = TRUE)

#making a histogram 

hist(nucleotide_diversity$PI, br=20)

#making a boxplot 

boxplot(nucleotide_diversity$PI, ylab="nucleotide diversity")


#Tajima's D#

tajimas <- read.table("efulgens.tajima.Tajima.D", header = TRUE)
hist(tajimas$TajimaD)
boxplot(tajimas$TajimaD)

#Plotting tajima's D along a chromosome#
#Subseting a chromosome 

tajima.chr1 <- subset(tajimas, CHROM == "scaffold_1")

plot(tajima.chr1$BIN_START,tajima.chr1$TajimaD, xlab="Position", ylab="Tajima's D")


#----------------------------------------------------#


#### Pairwise Fst - COASTAL x INLAND populations ####

# load package
library(ggplot2)
library(dplyr)
library(tidyverse)
library(stats)

# load FST file
Fst_coastal_inland <- read.table("Fst_efulgens_coastal_vs_inland.windowed.weir.fst", header = TRUE)
# calculate the middle point of each bin
Fst_coastal_inland$midPos <- (Fst$BIN_START + Fst$BIN_END) / 2
# explore the object
head(Fst_coastal_inland)

#selecting only the chromosomes/ 12 largest scaffolds 
Fst_chrm_coastal_inland <- Fst[1:7399,]


# plot
ggplot(data = Fst_chrm_coastal_inland, aes(x = midPos / 1000000, y = WEIGHTED_FST, col = CHROM)) +
  geom_point() + 
  geom_smooth() + 
  theme_classic() +
  facet_grid(cols = vars(CHROM), scales = "free_x", space = "free_x") +
  labs(x = "position (in MB)")






#----------------------------------------------------#


#### Pairwise Fst with STAMPP and Heatmaps ####
#Physalia Courses - Adaptation Genomics#
#https://github.com/MafaldaSFerreira/physalia_adaptation_course-2024/blob/main/02_day2/Tutorial_day2_IBD_optional.md


# load packages
library(dplyr)
library(magrittr)
library(tibble)
library(gplots)
library(RColorBrewer)
library(corrplot)

#creating a function to be able to plot a matrix into a heatmap#

makeSymm <- function(m, position) {
  # add symmetrical triangle matrix (upper or lower)
  if (position == "upper") {
    m[upper.tri(m)] <- t(m)[upper.tri(m)]
    return(m)
  }
  if (position == "lower") {
    m[lower.tri(m)] <- t(m)[lower.tri(m)]
    return(m)
  }
}


# load the FST matrix for all SNPs
#C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics/Population_Genomics/Fst/StAMPP
fst.mat <- read.table("~/Documents/Jac/Population Genomics/Tabelas/StAMPP/efulgens_fst_matrix.txt")

# use the given function to fill the upper diagonal of the matrix
fst.all.mat <- fst.mat %>%
  as.matrix(.) %>%
  makeSymm(., "upper")

fst.all.mat[is.na(fst.all.mat)] <- 0 # replace NAs by 0 (NAs unaccepted for the heatmap function)
fst.all.mat[1:5, 1:5] # check the fst_matrix

# visualise values
corrplot(fst.all.mat, is.corr = FALSE, method = "number", addgrid.col = FALSE, diag = FALSE, type = "lower", number.digits = 3, number.cex = 0.7)


# visualize pairwise FST with a heatmap plot
gplots::heatmap.2(fst.all.mat, trace = "none",
                  #col = colorRampPalette(brewer.pal(9, "Reds"))(15),
                  #col = colorRampPalette(brewer.pal(9, "Greens")),
                  col = colorRampPalette(brewer.pal(9, "PuBuGn")),
                  key.xlab = "FST")




#----------------------------------------------------#

#### Isolation by Distance (IBD) Analyses ####

# if not done already
install.packages("geosphere")

# load packages
library(reshape2)
library(dplyr)
library(magrittr)
library(tibble)
library(ggplot2)
library(geosphere)

# import information about populations
pop_coord <- read.csv("~/Documents/Jac/Population Genomics/Tabelas/Coordinates_Pops.csv", header=TRUE) 
head(pop_coord)

pop_info <- read.csv2("~/Documents/Jac/Population Genomics/Tabelas/Samples_and_Populations-Habitats.csv", header=TRUE) 
pop <- as.character(pop_info$population)


# calculate geographic (euclidian) distances between all pairs of populations
distance <- geosphere::distm(pop_coord[, c(3, 4)], fun = distGeo) %>%
  as.matrix(.)

#distance <- geosphere::distm(info_pop[, c(4, 3)], fun = distGeo) %>%
#  as.matrix(.)  # correct lat and long order?

# change it from meters to km
distance <- distance / 1000

# set the colnames and rownames of the distance matrix
dimnames(distance) <- list(pop_info$population, pop_info$population)
distance

# prepare datasets
# linearize the distance matrix
dist.melt <- reshape2::melt(distance) %>%
  set_colnames(., c("pop1", "pop2", "distance"))
head(dist.melt)
dist.melt

# linearize the fst matrix
fst.melt <- reshape2::melt(fst.all.mat) %>%
  set_colnames(., c("pop1", "pop2", "FST"))

# join the distance and fst
IBD.df <- left_join(dist.melt, fst.melt, by = c("pop1", "pop2")) %>%
  filter(., distance > 0)
head(IBD.df)

# test association with FST
cor.test(log(IBD.df$distance), IBD.df$FST / (1 - IBD.df$FST))

# plot IBD
ggplot(IBD.df) + aes(x = log(distance), y = FST / (1 - FST)) +
  geom_point(color="black", shape=18) +
  geom_smooth(method = "lm", formula = y~x, color="#56B4E9") +
  theme_light() +
  labs(
    x="Distance (log)",
    y="Fst/(1-Fst)",
    title="Correlation: Isolation by Distance"
  ) 





#----------------------------------------------------#

#### Isolation by Environment (IBE) Analysis with Mantel Test ####

library(vegan)
library(geosphere)


# import information about populations
pop_coord <- read.csv("~/Documents/Jac/Population Genomics/Tabelas/Coordinates_Pops.csv", header=TRUE) 
head(pop_coord)

pop_info <- read.csv2("~/Documents/Jac/Population Genomics/Tabelas/Samples_and_Populations-Habitats.csv", header=TRUE) 
pop <- as.character(pop_info$population)


## Genetic distances ## 
#These will be from the Fst matrix from previous section
gen_dist <- fst.all.mat


## Geographic distances ##
# calculate geographic (euclidian) distances between all pairs of populations
distance <- geosphere::distm(pop_coord[, c(3, 4)], fun = distGeo) %>%
  as.matrix(.)

# change it from meters to km
distance <- distance / 1000

# set the colnames and rownames of the distance matrix
dimnames(distance) <- list(pop_info$population, pop_info$population)
distance

geo_dist <- distance

## Environmental distances ##
pops_allBIOCLIM <- read.csv2("~/Documents/Jac/Population Genomics/Tabelas/only_pops_data_ALL_BIOCLIM.csv", header = TRUE, sep=",")
str(pops_allBIOCLIM)

# selection of environmental predictors based on VIF:
bioclim <- pops_allBIOCLIM[, 4:22]

#predictors <- all_bioclim %>% select(bio1, bio2, bio8, bio14)
predictors_mantel <- pops_allBIOCLIM[, c("bio1", "bio2", "bio8", "bio14")]


# Convert factor or character columns to numeric if necessary
predictors_mantel[] <- lapply(predictors_mantel, function(x) {
  if(is.factor(x) || is.character(x)) {
    return(as.numeric(as.character(x)))  # Convert factors/characters to numeric
  } else {
    return(x)
  }
})
#scale the variables so they are all comparable
env_bioclim_scale <- scale(predictors_mantel)

# Calculate the environmental distance matrix using Euclidean distance
env_dist <- dist(env_bioclim_scale, method = "euclidean")
print(as.matrix(env_dist))
env_dist <- as.matrix(env_dist)


## Perform the partial Mantel test ## 
partial_mantel_result <- mantel.partial(gen_dist, env_dist, geo_dist, method = "pearson", permutations = 999)
print(partial_mantel_result)

#Mantel’s r (partial): This is the correlation between genetic distance and environmental distance, controlling for geographic distance.
#p-value: Indicates whether the correlation is statistically significant.


# Plotting the Mantel Test:)

# Convert distance matrices to vectors
genetic_vector <- as.vector(gen_dist)  # Flatten the genetic distance matrix
environmental_vector <- as.vector(env_dist)  # Flatten the environmental distance matrix

# Check if both vectors have the same length
length(genetic_vector) == length(environmental_vector)


# Scatter plot of genetic vs. environmental distances
plot(genetic_vector, environmental_vector, 
     xlab = "Genetic Distance", 
     ylab = "Environmental Distance", 
     main = "Mantel Test: Genetic vs. Environmental Distances", 
     pch = 19, col = "#927BC9")

# Optionally, add a regression line to visualize the trend
abline(lm(environmental_vector ~ genetic_vector), col = "#56B4E9")


#Plotting with ggplot
#Need to transform the data into a dataframe in order to use geom_point() for the scatterplot.

data_mantel_ggplot <- data.frame(
  Genetic_Distance = genetic_vector,
  Environmental_Distance = environmental_vector
)

ggplot(data_mantel_ggplot, aes(x= Genetic_Distance, y= Environmental_Distance)) +
  geom_point(color="black", shape=18)+
  geom_smooth(method = "lm", formula = y~x, color="#56B4E9") +
  labs(
    x="Genetic Distance",
    y="Environmental Distance",
    title="Mantel Test: Isolation by Environment"
  ) +
  theme_light()



# Add Mantel result as text on the plot
mantel_r <- round(partial_mantel_result$statistic, 3)
mantel_p <- round(partial_mantel_result$signif, 3)

text(x = max(genetic_vector), y = min(environmental_vector), 
     labels = paste("r =", mantel_r, "\np =", mantel_p), 
     pos = 4, col = "#784")



## Plotting a Partial Mantel Test

# Convert distance matrices to vectors
genetic_vector <- as.vector(gen_dist)
environmental_vector <- as.vector(env_dist)
geographic_vector <- as.vector(geo_dist)

# Use linear models to extract residuals
# Fit linear model: environmental ~ geographic
env_geo_lm <- lm(environmental_vector ~ geographic_vector)

# Get the residuals from the environmental vs. geographic model
env_residuals <- resid(env_geo_lm)

# Now fit genetic ~ geographic to control for geographic distance
gen_geo_lm <- lm(genetic_vector ~ geographic_vector)

# Get the residuals from the genetic vs. geographic model
gen_residuals <- resid(gen_geo_lm)

# Now, you can plot the residuals (genetic residuals vs. environmental residuals)
plot(gen_residuals, env_residuals, 
     xlab = "Genetic Distance Residuals", 
     ylab = "Environmental Distance Residuals", 
     main = "Partial Mantel Test: Residuals Plot", 
     pch = 19, col = "#927BC9")

# Add a regression line
abline(lm(env_residuals ~ gen_residuals), col = "#56B4E9")

# Extract Mantel's r and p-value
mantel_r <- round(partial_mantel_result$statistic, 3)  # Mantel's r value
mantel_p <- round(partial_mantel_result$signif, 3)     # p-value

# Add text to the plot
text(x = max(gen_residuals), y = min(env_residuals), 
     labels = paste("r =", mantel_r, "\np =", mantel_p), 
     pos = 4, col = "black")




#-----------------------------------------------------#

#### Investigate outliers of differentiation - With OUTFLANK ####

#OutFLANK is an R package that implements the method developed by Whitlock and Lotterhos using likelihood on a trimmed distribution of FST values to infer the distribution of FST for neutral markers. This distribution is then used to assign q-values to each locus to detect outliers that may be due to spatially heterogeneous selection.#

#Whitlock, M. C., and K. J. Lotterhos. 2015. Reliable detection of loci responsible for local adaptation: Inference of a neutral model through trimming the distribution of FST. The American Naturalist. 186:S24–S36.


# load packages
library(OutFLANK)
library(vcfR)
library(ggplot2)

# use the library vcfR to convert the VCF into the OutFLANK format
vcf_file <- "C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics/Final_VCF/PLINK_LD/vcf_pruned_by_plink.vcf"
obj.vcfR <- read.vcfR(vcf_file, verbose = FALSE)

# extract information about SNP id and position
position <- getPOS(obj.vcfR) # positions in bp
chromosome <- getCHROM(obj.vcfR) # chromosome information
id_snp <- getID(obj.vcfR) # ID of the SNP

# gather this info in a dataframe
chr_pos <- as.data.frame(cbind(id_snp, chromosome, position)) # save info about id, chr, position
str(chr_pos) # explore the column types

# R is sometimes not good at categorizing columns, and here we had a problem that bp position was converted to a character and we need it as a number 
# use this command to transform this column into numeric
chr_pos$position <- as.numeric(as.character(chr_pos$position)) 


# we expect that it will be useful for subsequent analysis to have a file with snp id and position, so let's save this data frame as a text file in our directory
write.table(chr_pos, "C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics/Population_Genomics/Fst_Outliers/Outflank/table_snpID_positions.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# extract and format the genotype matrix
geno <- extract.gt(obj.vcfR) # character matrix containing the genotypes

# create an empty matrix, (9 stands for missing data)
G <- matrix(9, nrow = nrow(geno), ncol = ncol(geno))

# that we fill with genotypes
G[geno %in% c("0/0", "0|0")] <- 0
G[geno %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G[geno %in% c("1/1", "1|1")] <- 2

# an overview of our data and its first 10 rows/10 columns
table(as.vector(G))
dim(G)
G[1:10,1:10]

# as it will be useful later, we suggest to save this object as a text file
write.table(G, "C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics/Population_Genomics/Fst_Outliers/Outflank/geno_matrix.txt", sep = "\t", col.names = FALSE, row.names = FALSE)

# import population info
pop_info <- read.csv2("C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics/Population_Genomics/Samples_Info/Samples_and_Populations-Habitats.csv", header=TRUE) 
pop_vector <- as.character(pop_info$population)

#We will now use OutFLANK to calculate FST for each locus. It needs the information about populations. For OutFLANK we will keep only the pop column. Then we will calculate a FST value for each SNP

# FST matrix with OutFLANK
my_fst <- MakeDiploidFSTMat(t(G), locusNames = id_snp, popNames = pop_vector)





#-----------------------------------------------------#

#### Investigate outliers of differentiation - With BayPASS ####

# load package
library(ggplot2)
library(OutFLANK)
library(vcfR)
library(ggplot2)

# use the library vcfR (the one with all SNPs, and not only the pruned snps)
vcf_file <- "C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics/Final_VCF/epidendrum-final-variants-removed.maf_miss.vcf"
obj.vcfR_allSNPs <- read.vcfR(vcf_file, verbose = FALSE)

# extract information about SNP id and position
position <- getPOS(obj.vcfR_allSNPs) # positions in bp
chromosome <- getCHROM(obj.vcfR_allSNPs) # chromosome information
id_snp <- getID(obj.vcfR_allSNPs) # ID of the SNP

# gather this info in a dataframe
chr_pos_allSNPs <- as.data.frame(cbind(id_snp, chromosome, position)) # save info about id, chr, position
str(chr_pos_allSNPs) # explore the column types

# load XtX values
xtx_allsnps <- read.table("allsnps_epidendrum.controlled.output_summary_pi_xtx.out", header = TRUE)
head(xtx_allsnps)
# we will mostly work with M_XtX

# load position info about the SNPs
#SNP_pos <- read.table("02_data/SNP_pos.txt", header = TRUE)
SNP_pos <- chr_pos_allSNPs

# should be same number of rows
dim(xtx_allsnps)
dim(SNP_pos)


xtx_pos <- cbind(SNP_pos, xtx_allsnps)

ggplot(xtx_pos, aes(x = position, y = M_XtX, colour = chromosome)) + 
  geom_point() +
  theme_classic() +
  facet_grid(cols = vars(chromosome), scales = "free_x", space = "free_x")

#Now we realised that we really need to know at which value we put the threshold:

# load XtX values from simulatd data
xtx_simu <- read.table("simulate_controlled.output_summary_pi_xtx.out", header = TRUE)
head(xtx_simu)

# calculate the threshold
threshold_fdr0.01 = quantile(xtx_simu$M_XtX, probs = 0.99)
threshold_fdr0.05 = quantile(xtx_simu$M_XtX, probs = 0.95)


#selecting only the 12 chromosomes
xtxpos_12chrs <- dplyr::filter(xtx_pos, chromosome %in% c("scaffold_1", "scaffold_2", "scaffold_3", "scaffold_4", "scaffold_5", "scaffold_6", "scaffold_7", "scaffold_8", "scaffold_9", "scaffold_10", "scaffold_11", "scaffold_12"))


# add it on the plot
ggplot(xtxpos_12chrs, aes(x = position, y = M_XtX, colour = chromosome)) + 
  geom_point() +
  theme_classic() +
  facet_grid(cols = vars(chromosome), scales = "free_x", space = "free_x") +
  geom_hline(aes(yintercept = threshold_fdr0.05), linetype = "dotted", linewidth = 1, col = "red", show.legend = FALSE) +
  geom_hline(aes(yintercept = threshold_fdr0.01), linetype = "dotted", linewidth = 1, show.legend = FALSE)

# output outliers
xtx_pos[xtx_pos$M_XtX >= threshold_fdr0.05, ]





#-----------------------------------------------------#

#### Genotype-Environment Association Analyses - With RDA (Redundancy Analysis)####

install.packages("LEA")

library(vegan)    # Used to run PCA & RDA
library(lfmm)     # Used to run LFMM
library(LEA)
library(vcfR)
library(ggplot2)

# use the library vcfR to convert the VCF into the OutFLANK format
vcf_file <- "~/Documents/Jac/Population Genomics/Tabelas/vcf_pruned_by_plink.vcf"
obj.vcfR <- read.vcfR(vcf_file, verbose = FALSE)


# extract information about SNP id and position
position <- getPOS(obj.vcfR) # positions in bp
chromosome <- getCHROM(obj.vcfR) # chromosome information
id_snp <- getID(obj.vcfR) # ID of the SNP

# gather this info in a dataframe
chr_pos <- as.data.frame(cbind(id_snp, chromosome, position)) # save info about id, chr, position
str(chr_pos) # explore the column types

# R is sometimes not good at categorizing columns, and here we had a problem that bp position was converted to a character and we need it as a number 
# use this command to transform this column into numeric
chr_pos$position <- as.numeric(as.character(chr_pos$position)) 

# extract and format the genotype matrix
geno <- extract.gt(obj.vcfR) # character matrix containing the genotypes
#SNP position table
SNP_pos <- chr_pos


geno# transpose data and give meaningful colnames
gen <- t(geno)
colnames(gen) <- paste(SNP_pos$chromosome, SNP_pos$position, sep = "_")
gen [1:10, 1:10]
dim(gen)


# Converting values on "geno" matrix to 0, 1, and 2; instead of 0/0, 0/1, 1/0, 1/1. 
#0/0 = 0 = homozygote 
#1/0 and 0/1 = 1 = heterozygote 
#1/1 = 2 = homozygote for the alternative allele 

library(dplyr)

#transforming the genotype matrix to a dataframe
genotype_df <- as.data.frame(gen, stringsAsFactors = FALSE)

#replacing specific values using dplyr::recode()

genotype_df[] <- lapply(genotype_df, function(col) {recode(col, "0/0" = 0, "0/1" = 1, "1/0" = 1, "1/1" = 2)})
print(genotype_df)


# replace 9 by NA - 9 is already NAs in my matrix!
gen[which(gen == "9")] <- NA

# evaluate % of missing
sum(is.na(genotype_df))/(dim(genotype_df)[1]*dim(genotype_df)[2]) # 4% of missing data
sum(is.na(genotype_df))

# impute missing with the most common geno
gen.imp <- apply(genotype_df, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))

# we can't have NAs from now on
sum(is.na(gen.imp)) 
gen.imp

str(gen.imp)

    
#Let's look now at our environmental/phenotype matrix.
install.packages("psych")
library(psych)

env <- read.csv2("~/Documents/Jac/Population Genomics/Tabelas/pop_ind_lat_long_env_data_important_env_variables.csv", header = TRUE)
head(env)

# Make individual names characters (not factors)
env$individual_ID <- as.character(env$individual_ID)

# Confirm that genotypes and environmental data are in the same order
identical(rownames(gen.imp), env[,1]) 

env[,1]
row.names(gen.imp)

pairs.panels(env[,5:7], scale=T)


#______________________________________________________#

#Run RDA on environmental variable and test it

## RDA FOR TEMPERATURE (mean annual temperature) ##

# load package
library(vegan)

# run rda
temp.rda <- vegan::rda(gen.imp ~ env$temperature, scale = TRUE)
temp.rda

#Looking at the fraction of variance explained by the 1st axis of the RDA
RsquareAdj(temp.rda)


#Test the significance of the model using permutation tests.
temp.signif.full <- anova.cca(temp.rda, parallel = getOption("mc.cores")) # default is permutation = 999
temp.signif.full

#Analyse the RDA output
#We can plot the RDA. We’ll start with simple triplots from vegan. Here we’ll use scaling=3 (also known as “symmetrical scaling”) for the ordination plots. This scales the SNP and individual scores by the square root of the eigenvalues so that we can easily visualize them in the sample plot. Here, the SNPs are in red (in the center of each plot), and the individuals are colour-coded by population. The blue vectors are the environmental predictors. The relative arrangement of these items in the ordination space reflects their relationship with the ordination axes, which are linear combinations of the predictor variables.

#jpeg("C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics/Population_Genomics/Gene-Environment-Association/RDA/RDA-temperature")

plot(temp.rda, scaling = 3) 
points(temp.rda, display = "sites", pch = 20, cex = 1.3, col = as.factor(pop_info$population), scaling = 3)



#______________________________________________________#

#### RDA for all three environmental variables ####

# load package
library(vegan)

# run rda
env.rda <- vegan::rda(gen.imp ~ ., data=env[,5:7], scale = TRUE)
env.rda

#Looking at the fraction of variance explained by the 1st axis of the RDA
RsquareAdj(env.rda)

#$adj.r.squared
#[1] 0.09908134

#Our constrained ordination explains very little about the variation (9%); this low explanatory power is not surprising given that we expect that most of the SNPs in our dataset will not show a relationship with the environmental predictors (e.g., most SNPs will be neutral).

#The eigenvalues for the constrained axes reflect the variance explained by each canonical axis:

summary(eigenvals(env.rda, model = "constrained"))
screeplot(env.rda)


#Test the significance of the model using permutation tests.
#Now let’s check our RDA model for significance using formal tests. We can assess both the full model and each constrained axis using F-statistics (Legendre et al, 2010). The null hypothesis is that no linear relationship exists between the SNP data and the environmental predictors. See ?anova.cca for more details and options.

env.signif.full <- anova.cca(env.rda, parallel = getOption("mc.cores")) # default is permutation = 999
env.signif.full

#Analyse the RDA output
#We can plot the RDA. We’ll start with simple triplots from vegan. Here we’ll use scaling=3 (also known as “symmetrical scaling”) for the ordination plots. This scales the SNP and individual scores by the square root of the eigenvalues so that we can easily visualize them in the sample plot. Here, the SNPs are in red (in the center of each plot), and the individuals are colour-coded by population. The blue vectors are the environmental predictors. The relative arrangement of these items in the ordination space reflects their relationship with the ordination axes, which are linear combinations of the predictor variables.

## extract % explained by the first 2 axes
perc <- round(100*(summary(env.rda)$cont$importance[2, 1:2]), 2)
perc

#associating populations with colors

populations <- env$population
pop_colors <- c("Ubatuba" = "#CC79A7", "Bertioga" = "darkolivegreen2", "Cardoso" = "#927BC9", "Floripa" = "#F0E442", "Torres" = "#009E73", "Itapua" = "#56B4E9", "Arambare" = "#E69F00", "Pelotas" = "#000000")

#Map the population labels to their corresponding colors
point_colors <- pop_colors[populations]
point_colors


plot(env.rda, type="n", scaling = 3, xlab = paste0("RDA1 (", perc[1], "%)"), ylab = paste0("RDA2 (", perc[2], "%)") ) 
points(env.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3)  # the SNPs
points(env.rda, display = "sites", pch = 20, cex = 1.7, col=point_colors, scaling = 3)   # the individuals
text(env.rda, scaling=3, display="bp", col="#0868ac", cex=0.9)                           # the predictors
legend("topleft", legend = unique(populations), bty = "n", col = "gray32", pch = 21, cex = 0.9, pt.bg = pop_colors[unique(populations)])



#checking data to color populations

length(populations) #80
nrow(sc_si) #80

unique(populations)


#______________________________________________________#

#simple plots
plot(env.rda, type="n", scaling = 3)
points(env.rda, display = "sites", pch = 20, cex = 1.3, col = point_colors, scaling = 3)


#simple triplot
ordiplot(env.rda, scaling = 1, type = "text")


#custom triplot, step by step

## extract scores - these are coordinates in the RDA space
sc_si <- scores(env.rda, display="sites", choices=c(1,2), scaling=1)
sc_sp <- scores(env.rda, display="species", choices=c(1,2), scaling=1)
sc_bp <- scores(env.rda, display="bp", choices=c(1, 2), scaling=1)

## extract % explained by the first 2 axes
perc <- round(100*(summary(env.rda)$cont$importance[2, 1:2]), 2)
perc

#### Custom RDA triplot, step by step ####

# Set up a blank plot with scaling, axes, and labels
plot(env.rda,
     scaling = 1, # set scaling type 
     type = "none", # this excludes the plotting of any points from the results
     frame = FALSE,
     # set axis limits
     xlim = c(-4,4), 
     ylim = c(-4,4),
     # label the plot (title, and axes)
     main = "Triplot RDA - scaling 1",
     xlab = paste0("RDA1 (", perc[1], "%)"), 
     ylab = paste0("RDA2 (", perc[2], "%)") 
)
# add points for site scores
points(sc_si, 
       pch = 21, # set shape (here, circle with a fill colour)
       col = "black", # outline colour
       bg = point_colors, # fill colour
       cex = 1.2) # size

# add points for species scores
#points(sc_sp, 
#       pch = 22, # set shape (here, square with a fill colour)
#       col = "black",
#       bg = "#D0B663", 
#       cex = 1.2)

# add text labels for species abbreviations
text(sc_si + c(0.03, 0.09), # adjust text coordinates to avoid overlap with points 
     labels = rownames(sc_si), 
     col = "grey40", 
     font = 2, # bold
     cex = 0.6)

# the predictors
# add arrows for effects of the explanatory variables
#arrows(0,0, # start them from (0,0)
#       sc_bp[,1], sc_bp[,2], # end them at the score value
#       col = "#33a02c", 
#       lwd = 3)

# add environmental predictors
text(env.rda, scaling=1, display="bp", col="#33a02c", cex=0.9)                           # the predictors


# add text labels for arrows
text(x = sc_bp[,1] +3, # adjust text coordinate to avoid overlap with arrow tip
     y = sc_bp[,2] +3, 
     labels = rownames(sc_bp), 
     col = "#33a02c", 
     cex = 0.7, 
     font = 2)

legend("topleft", legend = unique(populations), bty = "n", col = "gray32", pch = 21, cex = 0.9, pt.bg = pop_colors[unique(populations)])


#______________________________________________________#

#### RDA with PC scores from environmental PCAs ####

env <- read.csv2("C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics/Population_Genomics/Samples_Info/Environmental_Data/pop_ind_lat_long_env_data_PCs.csv", header = TRUE)

PCs_data <- env[,5:8]

# Run all steps above first, and then continue from here:

# run rda
env.rda_PCs <- vegan::rda(gen.imp ~ ., data=env[,5:8], scale = TRUE)
env.rda_PCs

#Looking at the fraction of variance explained by the 1st axis of the RDA
RsquareAdj(env.rda_PCs)

#$adj.r.squared
#[1] 0.1124486

#Our constrained ordination explains very little about the variation (11%); this low explanatory power is not surprising given that we expect that most of the SNPs in our dataset will not show a relationship with the environmental predictors (e.g., most SNPs will be neutral).

#The eigenvalues for the constrained axes reflect the variance explained by each canonical axis:

summary(eigenvals(env.rda_PCs, model = "constrained"))
screeplot(env.rda_PCs)


#Test the significance of the model using permutation tests.
#Now let’s check our RDA model for significance using formal tests. We can assess both the full model and each constrained axis using F-statistics (Legendre et al, 2010). The null hypothesis is that no linear relationship exists between the SNP data and the environmental predictors. See ?anova.cca for more details and options.

env.signif.full <- anova.cca(env.rda_PCs, parallel = getOption("mc.cores")) # default is permutation = 999
env.signif.full

#Analyse the RDA output
#We can plot the RDA. We’ll start with simple triplots from vegan. Here we’ll use scaling=3 (also known as “symmetrical scaling”) for the ordination plots. This scales the SNP and individual scores by the square root of the eigenvalues so that we can easily visualize them in the sample plot. Here, the SNPs are in red (in the center of each plot), and the individuals are colour-coded by population. The blue vectors are the environmental predictors. The relative arrangement of these items in the ordination space reflects their relationship with the ordination axes, which are linear combinations of the predictor variables.

#associating populations with colors

populations <- env$population
pop_colors <- c("Ubatuba" = "#CC79A7", "Bertioga" = "darkolivegreen2", "Cardoso" = "#927BC9", "Floripa" = "#F0E442", "Torres" = "#009E73", "Itapua" = "#56B4E9", "Arambare" = "#E69F00", "Pelotas" = "#000000")

#Map the population labels to their corresponding colors
point_colors <- pop_colors[populations]
point_colors

## extract % explained by the first 2 axes
perc_PCs <- round(100*(summary(env.rda_PCs)$cont$importance[2, 1:2]), 2)
perc_PCs


plot(env.rda_PCs, type="n", scaling = 3, xlab = paste0("RDA1 (", perc_PCs[1], "%)"), ylab = paste0("RDA2 (", perc_PCs[2], "%)")) 
points(env.rda_PCs, display="species", pch=20, cex=0.7, col="gray32", scaling=3)  # the SNPs
points(env.rda_PCs, display = "sites", pch = 20, cex = 1.7, col=point_colors, scaling = 3)   # the individuals
text(env.rda_PCs, scaling=3, display="bp", col="#0868ac", cex=0.9)                           # the predictors
legend("bottomright", legend = unique(populations), bty = "n", col = "gray32", pch = 21, cex = 0.9, pt.bg = pop_colors[unique(populations)])


#custom triplot, step by step

## extract scores - these are coordinates in the RDA space
sc_si_PCs <- scores(env.rda_PCs, display="sites", choices=c(1,2), scaling=1)
sc_sp_PCs <- scores(env.rda_PCs, display="species", choices=c(1,2), scaling=1)
sc_bp_PCs <- scores(env.rda_PCs, display="bp", choices=c(1, 2), scaling=1)



#______________________________________________________#

#### RDA for ALL BIOCLIM variables ####


#BIOCLIM <- read.csv("C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics/Population_Genomics/Environmental_Variables/extracted_ALL_BIOCLIM_10min.csv", header = TRUE)

env_allBIOCLIM <- read.csv2("~/Documents/Jac/Population Genomics/Tabelas/pop_ind_lat_long_env_data_ALL_BIOCLIM.csv", header = TRUE)

all_bioclim <- env_allBIOCLIM[,5:23]
str(all_bioclim)

#investigating the environmental variables 
library(psych)    # Used to investigate correlations among predictors

#acommodating the window size
x11()
pairs.panels(all_bioclim, scale=T)

#RDA is a regression-based method, and so can be subject to problems when using highly correlated predictors (Dormann et al., 2013). 
#Generally, the |r| > 0.7 “rule of thumb” is a good guideline for removing correlated predictors. 
#We will also check for multicollinearity using Variance Inflation Factors (VIF), below.
#Variable reduction should be guided by an ecological interpretation of the relevance of possible predictors. 
#Here, we use the function pairs.panels to visualize correlations among our predictors. 
#Correlation coefficients are in the upper right diagonal, with their size scaled to their |r|. The lower left shows scatter plots, while the diagonal shows histograms of the data. See ?pairs.panels for more information.


#baseando no pairs.panels, escolhi sete variaveis que eram menos correlacionadas com todas as outras (R>0.7)

pred <- subset(all_bioclim, select=c(bio2, bio5, bio7, bio9, bio14, bio17, bio18))
pred

#conferindo a correlacao dessas variaveis que escolhi usando o pairs.panels de novo

pairs.panels(pred, scale=T)


# load package
library(vegan)

# run rda
env_allBIOCLIM.rda <- vegan::rda(gen.imp ~ ., data=all_bioclim, scale = TRUE)
env_allBIOCLIM.rda

env_pred.rda <- vegan::rda(gen.imp ~ ., data=pred, scale = TRUE)


#Looking at the fraction of variance explained by the 1st axis of the RDA
RsquareAdj(env_pred.rda)

#$adj.r.squared
#[1] 0.1374833

#Our constrained ordination explains very little about the variation (13%); this low explanatory power is not surprising given that we expect that most of the SNPs in our dataset will not show a relationship with the environmental predictors (e.g., most SNPs will be neutral).

#The eigenvalues for the constrained axes reflect the variance explained by each canonical axis:

summary(eigenvals(env_pred.rda, model = "constrained"))
screeplot(env_pred.rda)


#Test the significance of the model using permutation tests.
#Now let’s check our RDA model for significance using formal tests. We can assess both the full model and each constrained axis using F-statistics (Legendre et al, 2010). The null hypothesis is that no linear relationship exists between the SNP data and the environmental predictors. See ?anova.cca for more details and options.

env.signif.full <- anova.cca(env_pred.rda, parallel = getOption("mc.cores")) # default is permutation = 999
env.signif.full

#Analyse the RDA output
#We can plot the RDA. We’ll start with simple triplots from vegan. Here we’ll use scaling=3 (also known as “symmetrical scaling”) for the ordination plots. This scales the SNP and individual scores by the square root of the eigenvalues so that we can easily visualize them in the sample plot. Here, the SNPs are in red (in the center of each plot), and the individuals are colour-coded by population. The blue vectors are the environmental predictors. The relative arrangement of these items in the ordination space reflects their relationship with the ordination axes, which are linear combinations of the predictor variables.

## extract % explained by the first 2 axes
perc <- round(100*(summary(env_pred.rda)$cont$importance[2, 1:2]), 2)
perc

#associating populations with colors

populations <- env_allBIOCLIM$population
pop_colors <- c("Ubatuba" = "#CC79A7", "Bertioga" = "darkolivegreen2", "Cardoso" = "#927BC9", "Floripa" = "#F0E442", "Torres" = "#009E73", "Itapua" = "#56B4E9", "Arambare" = "#E69F00", "Pelotas" = "#000000")

#Map the population labels to their corresponding colors
point_colors <- pop_colors[populations]
point_colors


plot(env_pred.rda, type="n", scaling = 3, xlab = paste0("RDA1 (", perc[1], "%)"), ylab = paste0("RDA2 (", perc[2], "%)") ) 
points(env_pred.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3)  # the SNPs
points(env_pred.rda, display = "sites", pch = 20, cex = 1.7, col=point_colors, scaling = 3)   # the individuals
text(env_pred.rda, scaling=3, display="bp", col="#0868ac", cex=0.9)                           # the predictors
legend("topleft", legend = unique(populations), bty = "n", col = "gray32", pch = 21, cex = 0.9, pt.bg = pop_colors[unique(populations)])


#simple plots
plot(env_allBIOCLIM.rda, type="n", scaling = 3)
points(env_allBIOCLIM.rda, display = "sites", pch = 20, cex = 1.3, col = point_colors, scaling = 3)
ordiplot(env_allBIOCLIM.rda, scaling = 1, type = "text")


#custom triplot, step by step
## extract scores - these are coordinates in the RDA space
sc_si <- scores(env_allBIOCLIM.rda, display="sites", choices=c(1,2), scaling=1)
sc_sp <- scores(env_allBIOCLIM.rda, display="species", choices=c(1,2), scaling=1)
sc_bp <- scores(env_allBIOCLIM.rda, display="bp", choices=c(1, 2), scaling=1)

## extract % explained by the first 2 axes
perc <- round(100*(summary(env_allBIOCLIM.rda)$cont$importance[2, 1:2]), 2)
perc

## Custom RDA triplot, step by step ##

# Set up a blank plot with scaling, axes, and labels
plot(env_allBIOCLIM.rda,
     scaling = 1, # set scaling type 
     type = "none", # this excludes the plotting of any points from the results
     frame = FALSE,
     # set axis limits
     xlim = c(-4,4), 
     ylim = c(-4,4),
     # label the plot (title, and axes)
     main = "Triplot RDA - scaling 1",
     xlab = paste0("RDA1 (", perc[1], "%)"), 
     ylab = paste0("RDA2 (", perc[2], "%)") 
)
# add points for site scores
points(sc_si, 
       pch = 21, # set shape (here, circle with a fill colour)
       col = "black", # outline colour
       bg = point_colors, # fill colour
       cex = 1.2) # size

# add points for species scores
#points(sc_sp, 
#       pch = 22, # set shape (here, square with a fill colour)
#       col = "black",
#       bg = "#D0B663", 
#       cex = 1.2)

# add text labels for species abbreviations
text(sc_si + c(0.03, 0.09), # adjust text coordinates to avoid overlap with points 
     labels = rownames(sc_si), 
     col = "grey40", 
     font = 2, # bold
     cex = 0.6)

# the predictors
# add arrows for effects of the explanatory variables
#arrows(0,0, # start them from (0,0)
#       sc_bp[,1], sc_bp[,2], # end them at the score value
#       col = "#33a02c", 
#       lwd = 3)

# add environmental predictors
text(env_allBIOCLIM.rda, scaling=1, display="bp", col="#33a02c", cex=0.9)                           # the predictors


# add text labels for arrows
text(x = sc_bp[,1] +3, # adjust text coordinate to avoid overlap with arrow tip
     y = sc_bp[,2] +3, 
     labels = rownames(sc_bp), 
     col = "#33a02c", 
     cex = 0.7, 
     font = 2)

legend("topleft", legend = unique(populations), bty = "n", col = "gray32", pch = 21, cex = 0.9, pt.bg = pop_colors[unique(populations)])


#______________________________________________________#

#### RDA for final set of predictors based on VIF ####

# load package
library(vegan)
library(dplyr)

env_allBIOCLIM <- read.csv2("~/Documents/Jac/Population Genomics/Tabelas/pop_ind_lat_long_env_data_ALL_BIOCLIM.csv", header = TRUE)

all_bioclim <- env_allBIOCLIM[,5:23]
str(all_bioclim)

# selection of environmental predictors based on VIF:
bioclim <- env_allBIOCLIM[, 5:23]
bioclim <- as.matrix(bioclim)

corr <- vifcor(bioclim, th=0.7)
#kept bio1, bio2, bio8 and bio14

#predictors <- all_bioclim %>% select(bio1, bio2, bio8, bio14)
predictors <- all_bioclim[, c("bio1", "bio2", "bio8", "bio14")]


# run rda
predictors.rda <- vegan::rda(gen.imp ~ ., data=predictors, scale = TRUE)
predictors.rda 


#Looking at the fraction of variance explained by the 1st axis of the RDA
RsquareAdj(predictors.rda)

#$adj.r.squared
#[1] 0.1065415

#Our constrained ordination explains very little about the variation (13%); this low explanatory power is not surprising given that we expect that most of the SNPs in our dataset will not show a relationship with the environmental predictors (e.g., most SNPs will be neutral).

#The eigenvalues for the constrained axes reflect the variance explained by each canonical axis:

summary(eigenvals(predictors.rda, model = "constrained"))
screeplot(predictors.rda)


#Test the significance of the model using permutation tests.
#Now let’s check our RDA model for significance using formal tests. We can assess both the full model and each constrained axis using F-statistics (Legendre et al, 2010). The null hypothesis is that no linear relationship exists between the SNP data and the environmental predictors. See ?anova.cca for more details and options.

env.signif.full <- anova.cca(predictors.rda, parallel = getOption("mc.cores")) # default is permutation = 999
env.signif.full

#Analyse the RDA output
#We can plot the RDA. We’ll start with simple triplots from vegan. Here we’ll use scaling=3 (also known as “symmetrical scaling”) for the ordination plots. This scales the SNP and individual scores by the square root of the eigenvalues so that we can easily visualize them in the sample plot. Here, the SNPs are in red (in the center of each plot), and the individuals are colour-coded by population. The blue vectors are the environmental predictors. The relative arrangement of these items in the ordination space reflects their relationship with the ordination axes, which are linear combinations of the predictor variables.

## extract % explained by the first 2 axes
perc <- round(100*(summary(predictors.rda)$cont$importance[2, 1:2]), 2)
perc

perc_3_4 <- round(100*(summary(predictors.rda)$cont$importance[2, 3:4]), 2)
perc_3_4


#associating populations with colors

populations <- env_allBIOCLIM$population
pop_colors <- c("Ubatuba" = "#CC79A7", "Bertioga" = "darkolivegreen2", "Cardoso" = "#927BC9", "Florianopolis" = "#F0E442", "Torres" = "#009E73", "Itapua" = "#56B4E9", "Arambare" = "#E69F00", "Pelotas" = "#000000")

#Map the population labels to their corresponding colors
point_colors <- pop_colors[populations]
point_colors

#plotting RDA1 and RDA2

plot(predictors.rda, type="n", scaling = 3, xlab = paste0("RDA1 (", perc[1], "%)"), ylab = paste0("RDA2 (", perc[2], "%)") ) 
points(predictors.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3)  # the SNPs
points(predictors.rda, display = "sites", pch = 20, cex = 1.7, col=point_colors, scaling = 3)   # the individuals
text(predictors.rda, scaling=3, display="bp", col="#0868ac", cex=0.9)                           # the predictors
legend("topleft", legend = unique(populations), bty = "n", col = "gray32", pch = 21, cex = 0.9, pt.bg = pop_colors[unique(populations)])


#plotting RDA3 and RDA4 to check collinearity among predictors 
plot(predictors.rda, type="n", choices= c(3,4), scaling = 3, xlab = paste0("RDA3 (", perc_3_4[1], "%)"), ylab = paste0("RDA4 (", perc_3_4[2], "%)") ) 
points(predictors.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3, choices= c(3,4))  # the SNPs
points(predictors.rda, display = "sites", pch = 20, cex = 1.7, col=point_colors, scaling = 3, choices= c(3,4))   # the individuals
text(predictors.rda, scaling=3, display="bp", col="#0868ac", cex=0.9, choices= c(3,4))                           # the predictors
legend("topleft", legend = unique(populations), bty = "n", col = "gray32", pch = 21, cex = 0.9, pt.bg = pop_colors[unique(populations)])


#plotting with ggplot2
#first we need to extract the values of the scores of the rda (species, sites and predictors)
species_scores <- as.data.frame(scores(predictors.rda, display = "species", scaling = 3))
sites_scores <- as.data.frame(scores(predictors.rda, display = "sites", scaling = 3))
predictors_scores <- as.data.frame(scores(predictors.rda, display = "bp", scaling = 3))

# Add relevant columns
species_scores$type <- "species"
sites_scores$type <- "sites"
sites_scores$population <- populations
sites_scores$point_color <- point_colors
predictors_scores$type <- "predictors"

# Scale explained variance for axes labels
xlab <- paste0("RDA1 (", perc[1], "%)")
ylab <- paste0("RDA2 (", perc[2], "%)")

#Re-scaling the predictor values so the arrows don't look too small on the RDA plot
#We multiply the RDA1 and RDA2 values for predictors by a scaling factor (e.g., arrow_scale). Adjust this factor to balance arrow size.

arrow_scale <- 8

predictors_scores$RDA1 <- predictors_scores$RDA1 * arrow_scale
predictors_scores$RDA2 <- predictors_scores$RDA2 * arrow_scale


ggplot() +
  # Add SNP points (species)
  geom_point(data = species_scores, aes(x = RDA1, y = RDA2), 
             color = "gray32", size = 0.7) +
  # Add individual points (sites)
  geom_point(data = sites_scores, aes(x = RDA1, y = RDA2, fill = point_color), 
             shape = 21, size = 1.7) +
  # Add predictors (arrows)
  geom_segment(data = predictors_scores, 
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "black") +
  geom_text(data = predictors_scores, 
            aes(x = RDA1, y = RDA2, label = rownames(predictors_scores)), 
            color = "black", size = 4) +
  # Customize the legend
  scale_fill_identity(guide = "legend", 
                      labels = unique(populations), 
                      breaks = unique(point_colors)) +
  theme_light() +
  labs(x = xlab, y = ylab, fill = "Population") +
  theme(legend.position = "top")




#simple plots
plot(env_allBIOCLIM.rda, type="n", scaling = 3)
points(env_allBIOCLIM.rda, display = "sites", pch = 20, cex = 1.3, col = point_colors, scaling = 3)
ordiplot(env_allBIOCLIM.rda, scaling = 1, type = "text")


#custom triplot, step by step
## extract scores - these are coordinates in the RDA space
sc_si <- scores(env_allBIOCLIM.rda, display="sites", choices=c(1,2), scaling=1)
sc_sp <- scores(env_allBIOCLIM.rda, display="species", choices=c(1,2), scaling=1)
sc_bp <- scores(env_allBIOCLIM.rda, display="bp", choices=c(1, 2), scaling=1)

## extract % explained by the first 2 axes
perc <- round(100*(summary(env_allBIOCLIM.rda)$cont$importance[2, 1:2]), 2)
perc

## Custom RDA triplot, step by step ##

# Set up a blank plot with scaling, axes, and labels
plot(env_allBIOCLIM.rda,
     scaling = 1, # set scaling type 
     type = "none", # this excludes the plotting of any points from the results
     frame = FALSE,
     # set axis limits
     xlim = c(-4,4), 
     ylim = c(-4,4),
     # label the plot (title, and axes)
     main = "Triplot RDA - scaling 1",
     xlab = paste0("RDA1 (", perc[1], "%)"), 
     ylab = paste0("RDA2 (", perc[2], "%)") 
)
# add points for site scores
points(sc_si, 
       pch = 21, # set shape (here, circle with a fill colour)
       col = "black", # outline colour
       bg = point_colors, # fill colour
       cex = 1.2) # size

# add points for species scores
#points(sc_sp, 
#       pch = 22, # set shape (here, square with a fill colour)
#       col = "black",
#       bg = "#D0B663", 
#       cex = 1.2)

# add text labels for species abbreviations
text(sc_si + c(0.03, 0.09), # adjust text coordinates to avoid overlap with points 
     labels = rownames(sc_si), 
     col = "grey40", 
     font = 2, # bold
     cex = 0.6)

# the predictors
# add arrows for effects of the explanatory variables
#arrows(0,0, # start them from (0,0)
#       sc_bp[,1], sc_bp[,2], # end them at the score value
#       col = "#33a02c", 
#       lwd = 3)

# add environmental predictors
text(env_allBIOCLIM.rda, scaling=1, display="bp", col="#33a02c", cex=0.9)                           # the predictors


# add text labels for arrows
text(x = sc_bp[,1] +3, # adjust text coordinate to avoid overlap with arrow tip
     y = sc_bp[,2] +3, 
     labels = rownames(sc_bp), 
     col = "#33a02c", 
     cex = 0.7, 
     font = 2)

legend("topleft", legend = unique(populations), bty = "n", col = "gray32", pch = 21, cex = 0.9, pt.bg = pop_colors[unique(populations)])







#______________________________________________________#


#### Identifying candidate SNPs for local adaptation ####
##Proceeding using the RDA with the three environmental predictors selected on the niche modelling chapter##


#We’ll use the loadings of the SNPs in the ordination space to determine which SNPs are candidates for local adaptation. The SNP loadings are stored as species in the RDA object. We’ll extract the SNP loadings from the three significant constrained axes

load.rda <- scores(predictors.rda, choices=c(1:3), display="species")  # Species scores for the first three constrained axes

#If we look at histograms of the loadings on each RDA axis, we can see their (relatively normal) distributions. SNPs loading at the center of the distribution are not showing a relationship with the environmental predictors; those loading in the tails are, and are more likely to be under selection as a function of those predictors (or some other predictor correlated with them).

hist(load.rda[,1], main="Loadings on RDA1")
hist(load.rda[,2], main="Loadings on RDA2")
hist(load.rda[,3], main="Loadings on RDA3") 


#I’ve written a simple function to identify SNPs that load in the tails of these distributions. We’ll start with a 3 standard deviation cutoff (two-tailed p-value = 0.0027). As with all cutoffs, this can be modified to reflect the goals of the analysis and our tolerance for true positives vs. false positives. For example, if you needed to be very conservative and only identify those loci under very strong selection (i.e., minimize false positive rates), you could increase the number of standard deviations to 3.5 (two-tailed p-value = 0.0005). This would also increase the false negative rate. If you were less concerned with false positives, and more concerned with identifying as many potential candidate loci as possible (including those that may be under weaker selection), you might choose a 2.5 standard deviation cutoff (two-tailed p-value = 0.012).

#I define the function here as outliers, where x is the vector of loadings and z is the number of standard deviations to use

outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}


#Now let’s apply it to each significant constrained axis

cand1 <- outliers(load.rda[,1],3) 
cand2 <- outliers(load.rda[,2],3) 
cand3 <- outliers(load.rda[,3],3) 

length(cand1)
length(cand2)
length(cand3)

ncand <- length(cand1) + length(cand2) + length(cand3)
ncand # 278 candidate SNPs for local adaptation

#Next, we’ll organize our results by making one data frame with the axis, SNP name, loading, & correlation with each predictor:
  
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))

colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")

cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

#Let’s add in the correlations of each candidate SNP with the environmental predictors:
  
foo <- matrix(nrow=(ncand), ncol=4)  # 4 columns for 4 predictors
colnames(foo) <- c("BIO1", "BIO2", "BIO8", "BIO14")


for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- gen.imp[,nam]
  foo[i,] <- apply(predictors,2,function(x) cor(x,snp.gen))
}

cand <- cbind.data.frame(cand,foo)  
head(cand)

#Now we have a data frame of the candidate SNPs and their correlation with our environmental predictors!

#Investigating the candidate SNPs
#We’ll start off by looking for duplicate detection. These are SNPs that are identified as candidates on more than one RDA axis.

length(cand$snp[duplicated(cand$snp)])  # 0 duplicate detections
cand <- cand[!duplicated(cand$snp),] # remove duplicate detections

#Next, we’ll see which of the predictors each candidate SNP is most strongly correlated with:
  
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,8] <- names(which.max(abs(bar[4:6]))) # gives the variable
  cand[i,9] <- max(abs(bar[4:6]))              # gives the correlation
}


colnames(cand)[8] <- "predictor"
colnames(cand)[9] <- "correlation"

table(cand$predictor) 


#BIO1 BIO2 BIO8 
#90  121   67 

#90 SNPs associated with BIO1
#121 SNPs associated with BIO2
#67 SNPs associated with BIO8

#Note that, in some cases, correlations may be strong for multiple variables (depending on collinearity among predictors). 
#It may be useful to consider how candidate SNPs are correlated with multiple predictors. 
#We could, for example, look at the cand object and investigate correlations with predictors other than the 
#predictor with the highest correlation coefficient. 
#However, for this analysis we will focus on the strongest correlations of each SNP with one predictor.

#export table 

write.table(cand, "~/Documents/Jac/Population Genomics/Results/RDA/Candidate SNPs/candidate_SNPs.txt", row.names = FALSE, quote = FALSE, sep = "\t")
write.csv(cand, "~/Documents/Jac/Population Genomics/Results/RDA/Candidate SNPs/candidate_SNPs.csv", row.names = FALSE)


#### Plotting the strongly associated SNPs ####
#Let’s look at RDA plots again, but this time focus in on the SNPs in the ordination space. 
#We’ll color code the SNPs based on the predictor variable that they are most strongly correlated with. 

sel <- cand$snp
env <- cand$predictor

env[env=="BIO1"] <- '#1f78b4'
env[env=="BIO2"] <- '#6a3d9a'
env[env=="BIO8"] <- '#33a02c'

# color by predictor:
col.pred <- rownames(predictors.rda$CCA$v) # pull the SNP names
predictors.rda$CCA$v


# Set a default color for all SNPs
col.pred[] <- '#f1eef6'  # Light gray for non-candidate SNPs

# Color candidate SNPs based on their predictor
for (i in 1:length(sel)) {
  foo <- match(sel[i], col.pred)
  if (!is.na(foo)) {
    col.pred[foo] <- env[i]
  }
}

# Ensure all values in col.pred are valid colors
col.pred <- ifelse(col.pred %in% c('#1f78b4', '#6a3d9a', '#33a02c', '#f1eef6'), 
                   col.pred, '#f1eef6')

empty <- col.pred
empty[empty == "#f1eef6"] <- rgb(0,1,0, alpha=0)  # transparent
empty.outline <- ifelse(empty == "#00FF0000", "#00FF0000", "gray32")
#col.pred[grep("chr",col.pred)] <- '#f1eef6' # non-candidate SNPs
bg <- c('#1f78b4','#6a3d9a','#33a02c')


# axes 1 & 2
plot(predictors.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1))
points(predictors.rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3)
points(predictors.rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3)
text(predictors.rda, scaling=3, display="bp", col="#0868ac", cex=1)
legend("bottomright", legend=c("BIO4", "BIO13", "BIO18"), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)


#Debugging

print(table(col.pred))
print(table(empty))
print(table(empty.outline))

print(head(sel))
print(table(env))

for (i in 1:length(sel)) {
  foo <- match(sel[i], rownames(predictors.rda$CCA$v))
  if (!is.na(foo)) {
    col.pred[foo] <- env[i]
    print(paste("Matched SNP:", sel[i], "to color:", env[i]))
  } else {
    print(paste("Failed to match SNP:", sel[i]))
  }
}

print(table(col.pred))

plot(predictors.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1))
points(predictors.rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3)
points(predictors.rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3)
text(predictors.rda, scaling=3, display="bp", col="#0868ac", cex=1)
legend("bottomright", legend=c("BIO1", "BIO2", "BIO8"), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)




#______________________________________________________#

#### GEA with BayPass ####

# load package
library(ggplot2)
library(vcfR)

# use the library vcfR (the one with all SNPs, and not only the pruned snps)
vcf_file <- "~/Documents/Jac/Population Genomics/VCF/epidendrum-final-variants-removed.maf_miss.vcf"
obj.vcfR_allSNPs <- read.vcfR(vcf_file, verbose = FALSE)

# extract information about SNP id and position
position <- getPOS(obj.vcfR_allSNPs) # positions in bp
chromosome <- getCHROM(obj.vcfR_allSNPs) # chromosome information
id_snp <- getID(obj.vcfR_allSNPs) # ID of the SNP

# gather this info in a dataframe
chr_pos_allSNPs <- as.data.frame(cbind(id_snp, chromosome, position)) # save info about id, chr, position
str(chr_pos_allSNPs) # explore the column types

#Like we did before we can plot the XtX but we will mostly be interested in the BF value (Bayesian factor of association with an environmental variable). 
#You will find it in the files ending with betai_reg.out in the column BF.dB. 
#"the Bayes Factor (column BF(dB)) in dB units (i.e., 10 × log10(BF)) measuring the support of the association of each SNP with each population covariable and the corresponding regression coefficients βi (column Beta_is)" 
#BF can be informative in itself, good candidate are usually above 20. Following the rule of Jeffrey, we can consider BF as meaning

#< 3 --> nothing
#3 to 10 --> weak support
#10 to 20 --> interesting
#more than 20 --> strong support

# load Bayesian Factors values from the new table (from real SNPs dataset)
BF_allsnps <- read.table("~/Documents/Jac/Population Genomics/Results/BayPass/GEA/allsnps_environment_controlled_new.output_summary_betai_reg.out", header = TRUE)
head(BF_allsnps)
# we will be working with the BF (the Bayesian Factor) values now

# load position info about the SNPs
#SNP_pos <- read.table("02_data/SNP_pos.txt", header = TRUE)
SNP_pos <- chr_pos_allSNPs

# should be same number of rows
dim(BF_allsnps)
dim(SNP_pos)


BF_pos <- cbind(SNP_pos, BF_allsnps)

ggplot(BF_pos, aes(x = position, y = BF_allsnps$BF.dB., colour = chromosome)) + 
  geom_point() +
  theme_classic() +
  facet_grid(cols = vars(chromosome), scales = "free_x", space = "free_x")


#Now we realised that we really need to know at which value we put the threshold:

# load BF (bayes factor) values from simulated data
BF_simu <- read.table("~/Documents/Jac/Population Genomics/Results/BayPass/GEA/simulated_enviroment_controlled_new.output_summary_betai_reg.out", header = TRUE)
head(BF_simu)

# calculate the threshold
threshold_fdr0.01 = quantile(BF_simu$BF.dB, probs = 0.99)
threshold_fdr0.05 = quantile(BF_simu$BF.dB, probs = 0.95)
threshold_fdr0.001 = quantile(BF_simu$BF.dB, probs = 0.999)
threshold_BF <- BF_simu %>% filter(BF_simu$BF.dB. > 10)


#selecting only the 12 chromosomes
BFpos_12chrs <- dplyr::filter(BF_pos, chromosome %in% c("scaffold_1", "scaffold_2", "scaffold_3", "scaffold_4", "scaffold_5", "scaffold_6", "scaffold_7", "scaffold_8", "scaffold_9", "scaffold_10", "scaffold_11", "scaffold_12"))


# add it on the plot
ggplot(BFpos_12chrs, aes(x = position, y = BF.dB., colour = chromosome)) + 
  geom_point() +
  theme_classic() +
  facet_grid(cols = vars(chromosome), scales = "free_x", space = "free_x") +
  geom_hline(aes(yintercept = threshold_fdr0.01), linetype = "dotted", linewidth = 1, col = "red", show.legend = FALSE) +
  geom_hline(aes(yintercept = threshold_fdr0.001), linetype = "dotted", linewidth = 1, show.legend = FALSE) 


# output outliers
outliers <- BF_pos[BF_pos$BF.dB. >= threshold_fdr0.001, ]

#### Exporting the outliers that are associated with the environmental variables ####

#We can export the list of outlier SNPs for subsequent analysis. Here is the code for the controlled models.

# load bf values
bf_allsnps <- read.table("~/Documents/Jac/Population Genomics/Tabelas/BayPass/allsnps_environment_controlled.output_summary_betai_reg.out", header = TRUE)
bf_pos <- cbind(SNP_pos, bf_allsnps)

# load bf values from simulatd data
bf_simu <- read.table("~/Documents/Jac/Population Genomics/Tabelas/BayPass/simulated_enviroment_controlled.output_summary_betai_reg.out", header = TRUE)

# calculate the threshold from simulations (or you can use BF = 10)
threshold_fdr0.01 = quantile(bf_simu$BF.dB, probs = 0.99)
threshold_fdr0.05 = quantile(bf_simu$BF.dB, probs = 0.95)
threshold_BF <- bf_simu %>% filter(bf_simu$BF.dB. > 10)


outliers <- bf_pos[bf_pos$BF.dB >= threshold_fdr0.001, ]
outliers
write.table(outliers, "~/Documents/Jac/Population Genomics/Results/BayPass/GEA/outliers_bp.txt", row.names = FALSE, quote = FALSE, sep = "\t")
write.csv(outliers, "~/Documents/Jac/Population Genomics/Results/BayPass/GEA/outliers_bp.csv", row.names = FALSE)

#______________________________________________________#

#### Identify outliers detected by the two different methods (RDA and BayPass) ####

#It is often recommended to keep only outliers (or SNPs with strong GEA) detected by more than one method to avoid false positive. This will nevertheless also reduce the power of the analysis... a matter of choice? Here we will run a few R command to keep the intersection of our BayPass and RDA outliers. 
#It is worth noting that since RDA and BayPass works differently, it may not be surprising to have a limited overlap.

# load package
library(dplyr)

# load outliers tables
outlier_temp_rda <- read.table("~/Documents/Jac/Population Genomics/Results/RDA/Candidate SNPs/candidate_SNPs.txt", header = TRUE)
head(outlier_temp_rda)
nRDA <- dim(outlier_temp_rda)[1]
nRDA # how many outliers?

outlier_temp_bp <- read.table("~/Documents/Jac/Population Genomics/Results/BayPass/GEA/outliers_bp.txt", header = TRUE)
head(outlier_temp_bp)
outlier_temp_bp <- outlier_temp_bp[, c(2, 3, 10)] # we keep snp id, chr, pos and BF
dim(outlier_temp_bp)
nBP <- dim(outlier_temp_bp)[1]
nBP # how many outliers?

#concatenating column chr and position to make a new column named SNP_id just as we have on the RDA table
outlier_temp_bp$snp <- paste(outlier_temp_bp$chromosome, outlier_temp_bp$position, sep = "_")


# join outliers keeping positions present in both the 1st or the 2nd database
outlier_temp_fulljoin <- dplyr::full_join(outlier_temp_rda, outlier_temp_bp, by = c("snp" = "snp"))
head(outlier_temp_fulljoin)
nALL <- dim(outlier_temp_fulljoin)[1]
nALL # how many in total?

# join outliers keeping positions present in either the 1st or the 2nd database (only the intersection between them)
outlier_temp_innerjoin <- dplyr::inner_join(outlier_temp_rda, outlier_temp_bp, by = c("snp" = "snp"))
head(outlier_temp_innerjoin)
dim(outlier_temp_innerjoin)
nboth <- dim(outlier_temp_innerjoin)[1]
nboth # how many joint outliers?

# visualize
library(ggVennDiagram)
ggVennDiagram(list(RDA = 1:nRDA, BP = (nRDA + 1 - nboth):(nRDA - nboth + nBP))) + 
                scale_fill_gradient(high = "#9986A5", low = "#ABDDDE")




#________________________________________________________#

#### Preparing files for identifying the intersections between genes and annotated SNPs ####

# load file
outlier_temp_rda <- read.table("~/Documents/Jac/Population Genomics/Results/RDA/Candidate SNPs/candidate_SNPs.txt", header = TRUE)

# have a quick look
head(outlier_temp_rda)

# what's its dimension?
dim(outlier_temp_rda)

# which size around the SNP
window <- 10000
print(window)

# add a vector with start position
print(head(outlier_temp_rda$position))
print(length(outlier_temp_rda$position))

outlier_temp_rda$start <- outlier_temp_rda$position - (window / 2)

# start position can't be negative! replace negative by 0
outlier_temp_rda$start[outlier_temp_rda$start < 0] <- 0

# add a vector with stop position
outlier_temp_rda$stop <- outlier_temp_rda$position + (window / 2)

# have a look
head(outlier_temp_rda)

# which columns should we keep?
outlier_temp_rda_bed <- outlier_temp_rda[, c(2, 5, 6, 1)]

# save your file
write.table(outlier_temp_rda_bed, "05_bed/outlier_temp_rda.bed", row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)

























