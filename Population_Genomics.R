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

vcf_file <- "C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics_RNAseq/Final_VCF/PLINK_LD/vcf_pruned_by_plink.vcf"
vcf <- read.vcfR(vcf_file, verbose = FALSE)
sample_info <- read.csv("C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics_RNAseq/Sample_info/trinity_efulgens_table.csv", header=TRUE)


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
cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

plot12 <- plot12 + geom_point(aes(col = population)) + scale_colour_manual(values=cbp2) 
plot12

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


#### Basic Population Genetic Statistics from SNPS data ####

# Here we'll calculate 
# 1. Genetic Diversity 
# 2. Test Hardy Weinberg Equilibrium 
# 3. Fis and global Fst 

# load packages 

library("adegenet")
library("hierfstat")
library("pegas")
library("vcfR")

# Creating a GenIND object and onverting to a hierfstat object 


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
fst.mat <- read.table("C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics/Population_Genomics/Fst/StAMPP/efulgens_fst_matrix.txt")

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
                  col = colorRampPalette(brewer.pal(9, "Reds"))(15),
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
pop_coord <- read.csv("C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics/Population_Genomics/Samples_Info/Coordinates_Pops.csv", header=TRUE) 
head(pop_coord)

pop_info <- read.csv2("C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics/Population_Genomics/Samples_Info/Samples_and_Populations-Habitats.csv", header=TRUE) 
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
  geom_point() +
  geom_smooth(method = "lm", formula = y~x) +
  theme_bw()


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

library(vegan)    # Used to run PCA & RDA
library(lfmm)     # Used to run LFMM
library(LEA)
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

env <- read.csv2("C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics/Population_Genomics/Samples_Info/Environmental_Data/pop_ind_lat_long_env_data_important_env_variables.csv", header = TRUE)
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


#### Identifying candidate SNPs for local adaptation ####

#We’ll use the loadings of the SNPs in the ordination space to determine which SNPs are candidates for local adaptation. The SNP loadings are stored as species in the RDA object. We’ll extract the SNP loadings from the three significant constrained axes

load.rda <- scores(env.rda, choices=c(1:3), display="species")  # Species scores for the first three constrained axes

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

cand1 <- outliers(load.rda[,1],3) # 3
cand2 <- outliers(load.rda[,2],3) # 290
cand3 <- outliers(load.rda[,3],3) # 154

ncand <- length(cand1) + length(cand2) + length(cand3)
ncand # 447 candidate SNPs for local adaptation




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

env_allBIOCLIM <- read.csv2("C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics/Population_Genomics/Samples_Info/Environmental_Data/pop_ind_lat_long_env_data_ALL_BIOCLIM.csv", header = TRUE)

all_bioclim <- env_allBIOCLIM[,5:23]
all_bioclim$bio1

# load package
library(vegan)

# run rda
env_allBIOCLIM.rda <- vegan::rda(gen.imp ~ ., data=all_bioclim, scale = TRUE)
env_allBIOCLIM.rda

#Looking at the fraction of variance explained by the 1st axis of the RDA
RsquareAdj(env_allBIOCLIM.rda)

#$adj.r.squared
#[1] 0.09908134

#Our constrained ordination explains very little about the variation (9%); this low explanatory power is not surprising given that we expect that most of the SNPs in our dataset will not show a relationship with the environmental predictors (e.g., most SNPs will be neutral).

#The eigenvalues for the constrained axes reflect the variance explained by each canonical axis:

summary(eigenvals(env_allBIOCLIM.rda, model = "constrained"))
screeplot(env_allBIOCLIM.rda)


#Test the significance of the model using permutation tests.
#Now let’s check our RDA model for significance using formal tests. We can assess both the full model and each constrained axis using F-statistics (Legendre et al, 2010). The null hypothesis is that no linear relationship exists between the SNP data and the environmental predictors. See ?anova.cca for more details and options.

env.signif.full <- anova.cca(env_allBIOCLIM.rda, parallel = getOption("mc.cores")) # default is permutation = 999
env.signif.full

#Analyse the RDA output
#We can plot the RDA. We’ll start with simple triplots from vegan. Here we’ll use scaling=3 (also known as “symmetrical scaling”) for the ordination plots. This scales the SNP and individual scores by the square root of the eigenvalues so that we can easily visualize them in the sample plot. Here, the SNPs are in red (in the center of each plot), and the individuals are colour-coded by population. The blue vectors are the environmental predictors. The relative arrangement of these items in the ordination space reflects their relationship with the ordination axes, which are linear combinations of the predictor variables.

## extract % explained by the first 2 axes
perc <- round(100*(summary(env_allBIOCLIM.rda)$cont$importance[2, 1:2]), 2)
perc

#associating populations with colors

populations <- env$population
pop_colors <- c("Ubatuba" = "#CC79A7", "Bertioga" = "darkolivegreen2", "Cardoso" = "#927BC9", "Floripa" = "#F0E442", "Torres" = "#009E73", "Itapua" = "#56B4E9", "Arambare" = "#E69F00", "Pelotas" = "#000000")

#Map the population labels to their corresponding colors
point_colors <- pop_colors[populations]
point_colors


plot(env_allBIOCLIM.rda, type="n", scaling = 3, xlab = paste0("RDA1 (", perc[1], "%)"), ylab = paste0("RDA2 (", perc[2], "%)") ) 
points(env_allBIOCLIM.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3)  # the SNPs
points(env_allBIOCLIM.rda, display = "sites", pch = 20, cex = 1.7, col=point_colors, scaling = 3)   # the individuals
text(env_allBIOCLIM.rda, scaling=3, display="bp", col="#0868ac", cex=0.9)                           # the predictors
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

#### Custom RDA triplot, step by step ####

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

#### GEA with BayPass ####

# load package
library(ggplot2)
library(vcfR)

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

#Like we did before we can plot the XtX but we will mostly be interested in the BF value (Bayesian factor of association with an environmental variable). You will find it in the files ending with betai_reg.out in the column BF.dB. "the Bayes Factor (column BF(dB)) in dB units (i.e., 10 × log10(BF)) measuring the support of the association of each SNP with each population covariable and the corresponding regression coefficients βi (column Beta_is)" BF can be informative in itself, good candidate are usually above 20. Following the rule of Jeffrey, we can consider BF as meaning

#< 3 --> nothing
#3 to 10 --> weak support
#10 to 20 --> interesting
#more than 20 --> strong support

# load Bayesian Factors values from the new table (from real SNPs dataset)
BF_allsnps <- read.table("C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics/Population_Genomics/Gene-Environment-Association/BayPass/allsnps_environment_controlled.output_summary_betai_reg.out", header = TRUE)
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
BF_simu <- read.table("C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics/Population_Genomics/Gene-Environment-Association/BayPass/simulated_enviroment_controlled.output_summary_betai_reg.out", header = TRUE)
head(BF_simu)

# calculate the threshold
threshold_fdr0.01 = quantile(BF_simu$BF.dB, probs = 0.99)
threshold_fdr0.05 = quantile(BF_simu$BF.dB, probs = 0.95)


#selecting only the 12 chromosomes
BFpos_12chrs <- dplyr::filter(BF_pos, chromosome %in% c("scaffold_1", "scaffold_2", "scaffold_3", "scaffold_4", "scaffold_5", "scaffold_6", "scaffold_7", "scaffold_8", "scaffold_9", "scaffold_10", "scaffold_11", "scaffold_12"))


# add it on the plot
ggplot(BFpos_12chrs, aes(x = position, y = BF.dB., colour = chromosome)) + 
  geom_point() +
  theme_classic() +
  facet_grid(cols = vars(chromosome), scales = "free_x", space = "free_x") +
  geom_hline(aes(yintercept = threshold_fdr0.05), linetype = "dotted", linewidth = 1, col = "red", show.legend = FALSE) +
  geom_hline(aes(yintercept = threshold_fdr0.01), linetype = "dotted", linewidth = 1, show.legend = FALSE)

# output outliers
BF_pos[BF_pos$BF.dB. >= threshold_fdr0.05, ]

#### Exporting the outliers that are associated with the environmental variables ####

#We can export the list of outlier SNPs for subsequent analysis. Here is the code for the controlled models.

# load bf values
bf_allsnps <- read.table("C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics/Population_Genomics/Gene-Environment-Association/BayPass/allsnps_environment_controlled.output_summary_betai_reg.out", header = TRUE)
bf_pos <- cbind(SNP_pos, bf_allsnps)

# load bf values from simulatd data
bf_simu <- read.table("C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics/Population_Genomics/Gene-Environment-Association/BayPass/simulated_enviroment_controlled.output_summary_betai_reg.out", header = TRUE)

# calculate the threshold from simulations (or you can use BF = 10)
threshold_fdr0.01 = quantile(bf_simu$BF.dB, probs = 0.99)
threshold_fdr0.05 = quantile(bf_simu$BF.dB, probs = 0.95)

outliers <- bf_pos[bf_pos$BF.dB >= threshold_fdr0.05, ]
outliers
write.table(outliers, "C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics/Population_Genomics/Gene-Environment-Association/BayPass/outlier_temp_bp.txt", row.names = FALSE, quote = FALSE, sep = "\t")










