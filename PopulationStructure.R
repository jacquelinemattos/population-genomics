### Population Genomics Analyses for Epidendrum fulgens ###
### Jacqueline Mattos - 2024 ###


#### Population structure with LEA ####

# load the library
library(LEA)

# convert vcf to geno 
# note that we are working in the directory 02-day2 telling the program that the path of the populations.snps.vcf is inside the folder populations_2lin_random, and the output should go inside the folder LEA, with the name population_2lin_random.geno.
LEA::vcf2geno("populations_2lin_random/populations.snps.vcf",
              output.file = "LEA/population_2lin_random.geno")


#### Model Ancestry Proportions ####

#Now that we have our input file in the expected format, we will estimate individual admixture coefficients using sNMF. This program is implemented as the snmf() function in the R package LEA. In short, this function provides results very similar to programs such as STRUCTURE or ADMIXTURE. Assuming K ancestral populations, the snmf function provides least-squares estimates of ancestry proportions rather than maximum likelihood estimates (Frichot 2014).

#The results allow us to determine what is the best K value (i.e. the most likely number of genetic clusters)

# testing K populations: from K = 1 to K = 10
obj <- LEA::snmf("LEA/population_2lin_random.geno", K = 1:10, ploidy = 2,
                 entropy = TRUE, CPU = 4, project = "new")



#The snmf function computes an entropy criterion, which assesses the fit of the statistical model to the data using a cross-validation approach. The entropy criterion can help choosing the number of ancestral populations that best explains the genotypic data.

#_Find the best K value 
# plot cross-entropy
plot(obj, col = "blue4", cex = 1.4, pch = 19) # best is 2 here, lowest value
# choose the best LEA run
best = which.min(cross.entropy(obj, K = 2))

#Here the lowest cross-entropy value is clearly at K = 2 (lowest value), suggesting there are two genetic clusters within the dataset. Often, the plot shows a less clear pattern, and choosing the "knee/elbow/inflection" point is a generally good approach.

#The next step is to display a barplot of the ancestry matrix (also called the Q-matrix).

#Plot ancestry proportions across samples 

LEA::barchart(obj, K=2,run=best,border=NA,space=0,
              col=c("red","yellow"),
              xlab = "Individuals", ylab = "Ancestry proportions (K=2)", main = "Capelin lineages")


