#####Jacqueline Mattos, Setembro2024####
###Código para Redundancy Analysis###


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

#jpeg("C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics/Population_Genomics/Gene-Environment-Association/RDA/RDA-temperature")

populations <- c("Ubatuba", "Bertioga", "Cardoso", "Floripa", "Torres", "Itapua", "Arambare", "Pelotas")
pop_colors <- c("Ubatuba" = "#D0B663", "Bertioga" = "#1f78b4", "Cardoso" = "#ffff33", "Floripa" = "#a6cee3", "Torres" = "#33a02c", "Itapua" = "#7CACBA", "Arambare" = "#7BC9A2", "Pelotas" = "#927BC9")


plot(env.rda, type="n", scaling = 3) 
points(env.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3)  # the SNPs
points(env.rda, display = "sites", pch = 20, cex = 1.7, col=pop_colors[populations], scaling = 3)   # the individuals
text(env.rda, scaling=3, display="bp", col="#0868ac", cex=0.9)                           # the predictors
legend("topleft", legend=pop_info$population, bty="n", col="gray32", pch=21, cex=0.9, pt.bg=pop_colors)
