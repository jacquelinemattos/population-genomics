#Population Genetics Summary Statistics for Epidendrum fulgens adaptation genomics manuscript#
#J. Mattos, 2024

#Install packages#
install.packages(c("vcfR", "poppr", "pegas", "adegenet"))
install.packages("hierfstat")
library(vcfR)
library(poppr)
library(pegas)
library(adegenet)
library(hierfstat)


# Creating a GenIND object and converting to a hierfstat object 

#loading vcf file with vcfR package
vcf_file <- "~/Documents/Jac/Population Genomics/Tabelas/vcf_pruned_by_plink.vcf"
vcf <- read.vcfR(vcf_file, verbose = FALSE)

#loading population info data
#all populations/individuals
pop_info <- read.csv2("~/Documents/Jac/Population Genomics/Tabelas/pop_ind_lat_long_env_data_ALL_BIOCLIM.csv", header=TRUE) 
head(pop_info)
ind <- as.character(pop_info$individual_ID)
pop <- as.character(pop_info$population)



#separating the vcf into populations, to be able to get the summary statistics per population
pop_map <- read.csv2("~/Documents/Jac/Population Genomics/Tabelas/populations_table_efulgens_POPMAP.csv")
individuals <- colnames(vcf@gt)[-1] 
individuals

# Get individuals for the current population
ubatuba <- pop_map$individual_ID[pop_map$population == "Ubatuba"]
bertioga <- pop_map$individual_ID[pop_map$population == "Bertioga"]
cardoso <- pop_map$individual_ID[pop_map$population == "Cardoso"]
floripa <- pop_map$individual_ID[pop_map$population == "Florianopolis"]
torres <- pop_map$individual_ID[pop_map$population == "Torres"]
arambare <- pop_map$individual_ID[pop_map$population == "Arambare"]
pelotas <- pop_map$individual_ID[pop_map$population == "Pelotas"]
itapua <- pop_map$individual_ID[pop_map$population == "Itapua"]
  

#Transforming the vcf we imported to a "genind" object, from adegenet#
my_geneIND <- vcfR2genind(vcf)
class(my_geneIND)
my_geneIND

# Check population assignments
my_geneIND@pop

# Assign populations to the genind object based on the individual names
my_geneIND@pop <- factor(pop_map$population[match(indNames(my_geneIND), pop_map$individual_ID)])


# Now you can proceed with splitting
genind_list <- seppop(my_geneIND)

genind_Ubatuba <- genind_list[["Ubatuba"]]
genind_Bertioga <- genind_list[["Bertioga"]]  
genind_Cardoso <- genind_list[["Cardoso"]]
genind_Floripa <- genind_list[["Floripa"]]
genind_Torres <- genind_list[["Torres"]]
genind_Arambare <- genind_list[["Arambare"]]
genind_Pelotas <- genind_list[["Pelotas"]]
genind_Itapua <- genind_list[["Itapua"]]


#Transforming the genIND object into a hierfstat object
mydata <- genind2hierfstat(my_geneIND, pop = pop) 
class(mydata)


#creating a genIND for each population (just in case)
Ubatuba <- genind2hierfstat(genind_Ubatuba)
class(Ubatuba)
Ubatuba$pop <- c("Ubatuba", "Ubatuba", "Ubatuba","Ubatuba", "Ubatuba", "Ubatuba", "Ubatuba", "Ubatuba", "Ubatuba", "Ubatuba")
str(Ubatuba)
as.numeric(Ubatuba$pop)

#selecting each population from hierfstat object
library(dplyr)

Ubatuba_Hierfstat <- mydata %>% filter(pop == "Ubatuba")



#Using the genIND object to create a summary with Adegenet
div <- adegenet::summary(my_geneIND)
names(div)


#div_Ubatuba <- adegenet::summary(genind_Ubatuba)


############################################
#### Per population Basic Summary Stats ####
############################################


#Basic stats (Ho, Hs, Fis) for all populations :)

basicstats <- basic.stats(mydata, diploid = TRUE, digits = 2) 

names(basicstats)
Ho <- basicstats$Ho
Hs <- basicstats$Hs
Fis <- basicstats$Fis
Overall <- basicstats$overall

#getting the mean values for each parameter and each population

Ho <- as.data.frame(Ho)
Hs <- as.data.frame(Hs)
Fis <- as.data.frame(Fis)

mean_values_Ho <- Ho %>% summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))
mean_values_Hs <- Hs %>% summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))
mean_values_Fis <- Fis %>% summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))

mean_summary_stats <- rbind(mean_values_Ho, mean_values_Hs, mean_values_Fis)
row.names(mean_summary_stats) <- c("Ho", "Hs", "Fis")

#Saving this dataframe into a .csv file

write.csv(mean_summary_stats, file="~/Documents/Jac/Population Genomics/Results/Mean_Summary_Stats.csv", row.names = TRUE)


#Getting basic stats from hierfstat object
basic.stats(mydata) # Fst following Nei (1987) on genind object
wc(mydata) # Weir and Cockerham's estimate

#A list containing allele frequencies. Each element of the list is one locus. For each locus, Populations are in columns and alleles in rows
pop_freq <- pop.freq(mydata)
pop_freq


#	A table –with np (number of populations) columns and nl (number of loci) rows– of observed heterozygosities
obs_het <- Ho(mydata)
obs_het


#A table –with np (number of populations) columns and nl (number of loci) rows– of observed gene diversities
gene_diversity <- Hs(mydata)
gene_diversity


#### Testing for Hardy-Weinberg Equilibrium ####
#We used the function hw.test() from the pegas package.

hw.test(my_geneIND, B=1000)

#With this, we get for each locus a test of significance of the null hypothesis: H0 - the locus is in HW equilibrium in the population.


#### Creating boxplots for visualization ####

# The palette with black:
cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

boxplot(Ho, main = "Observed Heterozygosity", 
        ylab = "Values", xlab = "Populations", 
        col = cbp2)

boxplot(Hs, main = "Expected Heterozygosity", 
        ylab = "Values", xlab = "Populations", 
        col = cbp2)

boxplot(Fis, main = "Inbreeding Coefficient", 
        ylab = "Values", xlab = "Populations", 
        col = cbp2)





