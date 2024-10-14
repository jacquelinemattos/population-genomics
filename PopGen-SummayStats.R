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



#### VCFTools Per Sample Calculations ####

het_vcftools <- read.csv2("~/Documents/Jac/Population Genomics/Results/Mean Summary Stats/VCFTools/efulgens_not_LDpruned.het.het", header = TRUE, sep = "\t")

#Encontrando heterozigosidade 
#Numero de sites - Homozygous sites / Numero de sites
#Column4 - Column2 / Column4

het_vcftools$heterozygosity <- ((het_vcftools$N_SITES) - (het_vcftools$O.HOM.)) / (het_vcftools$N_SITES) 

#Separando cada populacao em um dataframe 

het_Itapua <- het_vcftools[1:10,]
het_Torres <- het_vcftools[11:20,]
het_Ubatuba <- het_vcftools[21:30,]
het_Bertioga <- het_vcftools[31:40,]
het_Cardoso <- het_vcftools[41:50,]
het_Floripa <- het_vcftools[51:60,]
het_Arambare <- het_vcftools[61:70,]
het_Pelotas <- het_vcftools[71:80,]


#Tirando a media dos valores para cada populacao - Obs Het; Exp Het; Fis

#Transformando colunas em valores numericos
str(het_Ubatuba)
het_Ubatuba$heterozygosity <- as.numeric(het_Ubatuba$heterozygosity)
het_Ubatuba$F <- as.numeric(het_Ubatuba$F)
mean_values_Ubatuba <- het_Ubatuba %>% summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))

str(het_Bertioga)
het_Bertioga$heterozygosity <- as.numeric(het_Bertioga$heterozygosity)
het_Bertioga$F <- as.numeric(het_Bertioga$F)
mean_values_Bertioga <- het_Bertioga %>% summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))

str(het_Cardoso)
het_Cardoso$heterozygosity <- as.numeric(het_Cardoso$heterozygosity)
het_Cardoso$F <- as.numeric(het_Cardoso$F)
mean_values_Cardoso <- het_Cardoso %>% summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))

str(het_Floripa)
het_Floripa$heterozygosity <- as.numeric(het_Floripa$heterozygosity)
het_Floripa$F <- as.numeric(het_Floripa$F)
mean_values_Floripa <- het_Floripa %>% summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))

str(het_Torres)
het_Torres$heterozygosity <- as.numeric(het_Torres$heterozygosity)
het_Torres$F <- as.numeric(het_Torres$F)
mean_values_Torres <- het_Torres %>% summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))

str(het_Arambare)
het_Arambare$heterozygosity <- as.numeric(het_Arambare$heterozygosity)
het_Arambare$F <- as.numeric(het_Arambare$F)
mean_values_Arambare <- het_Arambare %>% summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))

str(het_Pelotas)
het_Pelotas$heterozygosity <- as.numeric(het_Pelotas$heterozygosity)
het_Pelotas$F <- as.numeric(het_Pelotas$F)
mean_values_Pelotas <- het_Pelotas %>% summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))

str(het_Itapua)
het_Itapua$heterozygosity <- as.numeric(het_Itapua$heterozygosity)
het_Itapua$F <- as.numeric(het_Itapua$F)
mean_values_Itapua <- het_Itapua %>% summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))


#Juntando tudo em um dataframe 

mean_summary_stats_VCFTools <- rbind(mean_values_Ubatuba, mean_values_Bertioga, mean_values_Cardoso, mean_values_Floripa, mean_values_Torres, mean_values_Arambare, mean_values_Pelotas, mean_values_Itapua)
row.names(mean_summary_stats_VCFTools) <- c("Ubatuba", "Bertioga", "Cardoso", "Florianopolis", "Torres", "Arambare", "Pelotas", "Itapua")

mean_summary_stats_VCFTools$O.HOM. <- NULL
mean_summary_stats_VCFTools$N_SITES <- NULL

colnames(mean_summary_stats_VCFTools) <- c("Fis", "Obs Heterozygosity")
write.csv(mean_summary_stats_VCFTools, "~/Documents/Jac/Population Genomics/Results/Mean Summary Stats/VCFTools/mean_summary_stats_VCFTools.csv", row.names = TRUE)



#### Creating boxplots for visualization ####

#Heterozygosity
heterozygosity_all_pops <- cbind(het_Ubatuba$heterozygosity, het_Bertioga$heterozygosity, het_Cardoso$heterozygosity, het_Floripa$heterozygosity, 
                                 het_Torres$heterozygosity, het_Arambare$heterozygosity, het_Pelotas$heterozygosity, het_Itapua$heterozygosity)
heterozygosity_all_pops <- as.data.frame(heterozygosity_all_pops)

colnames(heterozygosity_all_pops) <- c("Ubatuba", "Bertioga", "Cardoso", "Florianopolis", "Torres", "Arambare", "Pelotas", "Itapua")

cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

boxplot(heterozygosity_all_pops, main = "Observed Heterozygosity", 
        ylab = "Values", xlab = "Populations", 
        col = cbp2)



#Inbreeding coefficient
Fis_all_pops <- cbind(het_Ubatuba$F , het_Bertioga$F, het_Cardoso$F, het_Floripa$F, 
                      het_Torres$F, het_Arambare$F, het_Pelotas$F, het_Itapua$F)
Fis_all_pops <- as.data.frame(Fis_all_pops)

colnames(Fis_all_pops) <- c("Ubatuba", "Bertioga", "Cardoso", "Florianopolis", "Torres", "Arambare", "Pelotas", "Itapua")

boxplot(Fis_all_pops, main = "Inbreeding Coefficient", 
        ylab = "Values", xlab = "Populations", 
        col = cbp2)









