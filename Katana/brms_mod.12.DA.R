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
            prior(student_t(3, 0, 10), class = "sd"),
            prior(student_t(3, 0, 10), class = "sigma"))

#The model
brm_12 <- brm(lnMass ~ 1 +
                 (1 + z_days_since_hatch + z_days_since_hatch_I2 | gr(F1_Genotype, cov = G_VCV)) + 
                  (1 + z_days_since_hatch + z_days_since_hatch_I2 | dam_id),
               family = gaussian(),
               prior = priors,
               data2 = list(G_VCV = G_VCV),
               data = data_DA, 
               chains = 4, cores = 4, iter = 6000, warmup = 1000, thin = 10,
               control = list(adapt_delta = 0.98), save_pars = save_pars(all = TRUE))

saveRDS(brm_12, "output/rds/brm_12")


# We can also compare the above model to a Von Bertalanffy growth model
growth_formula <- bf( mass ~ Linf * (1 - exp(-K * (z_days_since_hatch - t0))), 
                      Linf + K + t0 ~ 1 + (1|gr(F1_Genotype, cov = G_VCV)) + (1|dam_id), nl = T)


growth_prior <- c(set_prior("normal(0, 0.5)", nlpar = "Linf", lb = 0), 
                  set_prior("normal(0,0.005)", nlpar = "K", lb = 0), 
                  set_prior("normal(0,50)", nlpar = "t0"), 
                  set_prior("normal(0,10)", class = "sigma"), 
                  set_prior("normal(0,1)", class = "sd", group = "F1_Genotype", nlpar = "Linf"), 
                  set_prior("normal(0,1e-10)", class = "sd", group = "F1_Genotype", nlpar = "K"), 
                  set_prior("normal(0,1)", class = "sd", group = "dam_id", nlpar = "Linf"), 
                  set_prior("normal(0,1e-10)", class = "sd", group = "dam_id", nlpar = "K")) 

brm_12_Von_Bert <- brm(growth_formula, family = gaussian(), 
               data = data_DA, 
               data2 = list(G_VCV = G_VCV),
               prior = growth_prior, 
               sample_prior = T, cores = 4, chains = 4, warmup = 4000, iter = 10000, thin = 10, 
               control = list(adapt_delta = 0.98))

saveRDS(brm_12_Von_Bert, file = "output/rds/brm_12_Von_Bert")