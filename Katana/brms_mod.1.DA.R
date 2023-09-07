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

priors <- c(prior(normal(0, 10), "Intercept"),
            prior(student_t(3, 0, 10), class = "sd"),
            prior(student_t(3, 0, 10), class = "sigma"))

#The model
brm_1 <- brm(lnMass ~ 1 +
                 (1 | gr(F1_Genotype, cov = G_VCV)) + 
                 (1 | dam_id) + 
                 (1 | id),
               family = gaussian(),
               data2 = list(G_VCV = G_VCV),
               prior = priors,
               data = data_DA, 
               chains = 4, cores = 4, iter = 6000, warmup = 1000, thin = 10,
               control = list(adapt_delta = 0.98), save_pars = save_pars(all = TRUE))

add_criterion(brm_1, c("loo", "waic"), moment_match = TRUE)

saveRDS(brm_1, "output/rds/brm_1")

# Test model to make sure conversions back to observed scale for Vr are done correctly
#The model

priors <- c(prior(normal(0, 10), "Intercept"),
            prior(student_t(3, 0, 10), class = "sd"))

mod_sigma <- bf(lnMass ~ 1 +
                 (1  | gr(F1_Genotype, cov = G_VCV)) + 
                 (1 | dam_id),
                 sigma ~ 1)

brm_1_sigma <- brm(mod_sigma, 
               family = gaussian(),
               data2 = list(G_VCV = G_VCV),
               prior = priors,
               data = data_DA, 
               chains = 4, cores = 4, iter = 6000, warmup = 1000, thin = 10,
               control = list(adapt_delta = 0.98), save_pars = save_pars(all = TRUE))

saveRDS(brm_1_sigma, "output/rds/brm_1_sigma")

# Test 2

mod_sigma2 <- bf(lnMass ~ 1 +
                 (1  | gr(F1_Genotype, cov = G_VCV)) + 
                 (1 | dam_id),
                 sigma ~ z_days_since_hatch)

brm_1_sigma2 <- brm(mod_sigma2, 
               family = gaussian(),
               data2 = list(G_VCV = G_VCV),
               prior = priors,
               data = data_DA, 
               chains = 4, cores = 4, iter = 6000, warmup = 1000, thin = 10,
               control = list(adapt_delta = 0.98), save_pars = save_pars(all = TRUE))
