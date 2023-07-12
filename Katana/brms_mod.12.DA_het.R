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
                    dam_id = as.factor(dam_id))

#G matrix
G_VCV <- read.csv("output/G/Ga_SNPready.csv", row.names = 1) %>% as.matrix()

# Set some prirors
priors <- c(prior(normal(0, 10), "Intercept"),
            prior(student_t(3, 0, 10), class = "sd"))

#The model
brm_12_het_mod <- 
  bf(lnMass ~ 1 + (1 + z_days_since_hatch + z_days_since_hatch_I2 | gr(F1_Genotype, cov = G_VCV)) + 
                  (1 + z_days_since_hatch + z_days_since_hatch_I2 | dam_id), 
      sigma ~ z_days_since_hatch)

brm_12_het <- brm(brm_12_het_mod,
               family = gaussian(),
               prior = priors,
               data2 = list(G_VCV = G_VCV),
               data = data_DA, 
               chains = 4, cores = 4, iter = 6000, warmup = 1000, thin = 10,
               control = list(adapt_delta = 0.98), save_pars = save_pars(all = TRUE))

saveRDS(brm_12_het, "output/rds/brm_12_het")