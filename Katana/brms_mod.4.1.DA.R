#Set working directory
setwd("/srv/scratch/z3516573/gitrepo/ldeli_growth")

#Load libraries
library(dplyr)
library(magrittr)
library(brms)

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

#G matrix
G_VCV <- read.csv("output/G/Ga_SNPready.csv", row.names = 1) %>% as.matrix()

#The model
brm_4.1 <- brm(lnMass ~ 1 +
                 (1 + z_days_since_hatch + z_days_since_hatch_I2 | liz_id) + 
                 (1 + z_days_since_hatch + z_days_since_hatch_I2 | dam_id) + 
                 (1 + z_days_since_hatch + z_days_since_hatch_I2 | id),
               family = gaussian(),
               cov_ranef = list(liz_id = G_VCV),
               data = data_DA, 
               chains = 4, cores = 4, iter = 4000, warmup = 1500, thin = 5,
               control = list(adapt_delta = 0.98))

saveRDS(brm_4.1, "output/rds/brm_4.1")