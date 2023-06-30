#Set working directory
setwd("/srv/scratch/z3516573/gitrepo/ldeli_growth")

#Load libraries
library(dplyr)
library(magrittr)
library(brms)

##Read in data
# Dataset
cold_DA <- read.csv("data/growth/processed/analysis/COLD_quangen_growth_DA.csv", stringsAsFactors = F)
str(cold_DA)
dim(cold_DA)

#Format date and treatment
cold_DA %<>% mutate(treatment = as.factor(treatment),
                    liz_id = as.factor(liz_id),
                    dam_id = as.factor(dam_id), 
                    treatent = as.factor(treatment))

#G matrix
G_VCV <- read.csv("output/G/Ga_SNPready.csv", row.names = 1) %>% as.matrix()

cold_brm_5.5_mod <- bf(
  lnMass ~ 1 +
    (1 + z_days_since_hatch + z_days_since_hatch_I2 | liz_id) + 
    (1 + z_days_since_hatch + z_days_since_hatch_I2 | dam_id) + 
    (1 | id), 
  sigma ~  z_days_since_hatch 
)

cold_brm_5.5 <- brm(cold_brm_5.5_mod,
                   family = gaussian(),
                   cov_ranef = list(liz_id = G_VCV),
                   data = hot_DA, 
                   chains = 4, cores = 4, iter = 4000, warmup = 1500, thin = 5,
                   control = list(adapt_delta = 0.98))

saveRDS(cold_brm_5.5, "output/rds/cold_brm_5.5")