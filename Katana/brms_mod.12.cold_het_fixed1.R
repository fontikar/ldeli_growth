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

priors <- c(prior(normal(0, 5), "Intercept"),
            prior(normal(0, 10), "b"),
            prior(student_t(3, 0, 10), class = "sd"))

#The model

mods_cold <- bf(lnMass ~ z_days_since_hatch + z_days_since_hatch_I2 + 
                 (1 + z_days_since_hatch + z_days_since_hatch_I2 | gr(F1_Genotype, cov = G_VCV)) + 
                 (1 + z_days_since_hatch + z_days_since_hatch_I2 | dam_id),
                 sigma ~ z_days_since_hatch)

brm_12_cold_het_fixed1 <- brm(mods_cold,
               family = gaussian(),
               prior = priors,
               data2 = list(G_VCV = G_VCV),
               data = cold_DA, 
               chains = 4, cores = 4, iter = 6000, warmup = 1000, thin = 10,
               control = list(adapt_delta = 0.98, max_treedepth=12), save_pars = save_pars(all = TRUE))

saveRDS(brm_12_cold_het_fixed1, "output/rds/brm_12_cold_het_fixed1")
