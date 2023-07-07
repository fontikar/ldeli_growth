
#Load libraries
library(dplyr)
library(magrittr)
library(brms)
library(parallel)
library(MASS)

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

# Set some prirors
priors <- c(prior(normal(0, 10), "Intercept"),
            prior(student_t(3, 0, 10), class = "sd"),
            prior(student_t(3, 0, 10), class = "sigma"))

# Reviewer 2 Suggestion to Drop M

brm_8_hot <- brm(lnMass ~ 1 +
                 (1 + z_days_since_hatch + z_days_since_hatch_I2 | gr(F1_Genotype, cov = G_VCV)) +  
                 (1  | id),
               family = gaussian(),
               data2 = list(G_VCV = G_VCV),
			   prior = priors,
               data = hot_DA, 
               chains = 4, cores = 4, iter = 6000, warmup = 1000, thin = 10,
               control = list(adapt_delta = 0.98), save_pars = save_pars(all = TRUE))

add_criterion(brm_8_hot, c("waic", "loo"), moment_match = TRUE)

saveRDS(brm_8_hot, "output/rds/brm_8_hot")

# MCMCglmm

priorB<-list(G = list(G1 = list(V  = diag(3), 
	                               nu = 3, 
	                         alpha.mu = rep(0,3), 
	                          alpha.V = diag(25^2,3,3)),
                        G2 = list(V = diag(1), 
	                               nu = 1, 
	                         alpha.mu = rep(0,1), 
	                          alpha.V = diag(25^2,1,1))),
	             R = list(V = 1, nu = 0.002)) 

invA <- as(ginv(G_VCV), "dgCMatrix")
dimnames(invA) <- list(rownames(G_VCV), rownames(G_VCV))

mod_tets <- mclapply(1:3, function(i) MCMCglmm::MCMCglmm(
		      lnMass ~ 1, 
		      random = ~us(1 + z_days_since_hatch + z_days_since_hatch_I2):F1_Genotype + id,  
		      ginverse= list(F1_Genotype = invA),
		      data = data_DA, family = "gaussian", 
		      nitt = 550000, burnin = 50000, thin = 5, prior = priorB,
		      verbose = T))
