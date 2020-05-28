#Set working directory
setwd("/srv/scratch/z3516573/gitrepo/ldeli_growth")

#Load libraries
library(dplyr)
library(magrittr)
library(MCMCglmm)

##Read in data
# Dataset
data_DA <- read.csv("data/growth/processed/analysis/Ldeli_quangen_growth_DA.csv", stringsAsFactors = F)
str(data_DA)
dim(data_DA)

#Format date and treatment
data_DA %<>% mutate(treatment = as.factor(treatment),
                    liz_id = as.factor(liz_id),
                    dam_id = as.factor(dam_id), 
                    treatent = as.factor(treatment))

#Ginverse matrix
Ga_inverse <- read.csv("output/G/complete_preddata_InverseGa.csv", row.names = 1) %>% as.matrix()
Ga_inverse <- as(Ga_inverse, "CsparseMatrix")

# Set up priors for model and configure model
prior.3 <- list(G = list(G1 = list(V = diag(3), nu = 0.0002),
                         G2 = list(V = diag(3), nu = 0.0002)),
                R = list(V = 1, nu = 0.002))

nitt = 153000
burnin = 3000
thin = 100

#Run 1 chain of mod 1
mod.3.2.DA <- MCMCglmm(lnMass ~ z_hatch_since_start + z_lnEggMass + 
                       treatment + z_days_since_hatch +  z_days_since_hatch_I2,
                       random = ~us(1+z_days_since_hatch + z_days_since_hatch_I2):liz_id + us(1+z_days_since_hatch + z_days_since_hatch_I2):dam_id, 
                       family = "gaussian",
                       ginverse = list(liz_id = Ga_inverse),
                       prior = prior.3, 
                       nitt = nitt,
                       burnin = burnin,
                       thin = thin,
                       data = data_DA, 
                       verbose = T)

saveRDS(mod.3.2.DA, "output/rds/mod.3.2.DA")

