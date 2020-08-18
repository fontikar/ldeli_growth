#Set working directory
setwd("/srv/scratch/z3516573/gitrepo/ldeli_growth")

#Load libraries
library(dplyr)
library(magrittr)
library(brms)

##Read in data
# Dataset
hot_DA <- read.csv("data/growth/processed/analysis/HOT_quangen_growth_DA.csv", stringsAsFactors = F)
str(hot_DA)
dim(hot_DA)

#Format date and treatment
hot_DA %<>% mutate(treatment = as.factor(treatment),
                   liz_id = as.factor(liz_id),
                   dam_id = as.factor(dam_id), 
                   treatent = as.factor(treatment))

#G matrix
G_VCV <- read.csv("output/G/Ga_SNPready.csv", row.names = 1) %>% as.matrix()

hot_brm_5.5_mod <- bf(
  lnMass ~ 1 +
    (1 + z_days_since_hatch + z_days_since_hatch_I2 | liz_id) + 
    (1 + z_days_since_hatch + z_days_since_hatch_I2 | dam_id) + 
    (1 | id), 
  sigma ~  z_days_since_hatch 
)

hot_brm_5.5 <- brm(hot_brm_5.5_mod,
               family = gaussian(),
               cov_ranef = list(liz_id = G_VCV),
               data = hot_DA, 
               chains = 4, cores = 4, iter = 4000, warmup = 1500, thin = 5,
               control = list(adapt_delta = 0.98))

saveRDS(hot_brm_5.5, "output/rds/hot_brm_5.5")