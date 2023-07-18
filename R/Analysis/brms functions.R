
#Function to get model prediction from brms_5.het
func_growth_predictions <- function(day, predat, posterior){
  
  func <- function(liz_id = predat$liz_id[i]){
    # finding dam_id that corresponds to liz_id
    dam_id <- predat$dam_id[predat$liz_id == liz_id]
    treat <- predat$treatment[predat$liz_id == liz_id]
    
    pred <- posterior[, "b_Intercept"] + #Intercept
      posterior[, "b_z_days_since_hatch"] * ztran_DsH(day) + #Linear day effect 
      posterior[, "b_z_days_since_hatch_I2"] * ((ztran_DsH(day))^2) + #Curve day effect 
      posterior[, paste0("r_liz_id[",liz_id, ",Intercept]")] + #Lizard intercept
      posterior[, paste0("r_liz_id[",liz_id, ",z_days_since_hatch]")] + #Lizard slope
      posterior[, paste0("r_liz_id[",liz_id, ",z_days_since_hatch_I2]")] + #Lizard curve z_days_since_hatch_I2
      posterior[, paste0("r_dam_id[",dam_id, ",Intercept]")] + #Dam intercept
      posterior[, paste0("r_dam_id[",dam_id, ",z_days_since_hatch]")] + #Dam slope
      posterior[, paste0("r_dam_id[",dam_id, ",z_days_since_hatch_I2]")]  #Dam curve
    
    if(treat == 29){
      pred <- pred + 
        posterior[, "b_treatment29"] * 1 +
        posterior[, "b_treatment29:z_days_since_hatch"] * 1 + 
        posterior[, "b_treatment29:z_days_since_hatch_I2"] * 1
    }
    
    original_scale_pred <- posterior_summary(pred) %>% exp() #Calculate posterior mean and upper and lower
    
    df <- data.frame(liz_id = liz_id,
                     dam_id = dam_id,
                     treatment = treat,
                     day = day,
                     mean_mass = original_scale_pred[,1],
                     lower = original_scale_pred[,3],
                     upper = original_scale_pred[,4])
    return(df)
  } #Extract posterior mean and upper and lower
  df <- bind_rows(lapply(predat$liz_id, function(x) func(x))) %>% mutate(day = day)
  return(df)
}

#Transforming day to 'z scores' and back tranforming 'z scores of day' back to days
ztran_DsH <- function(day) {
  (day - mean(data$days_since_hatch, na.rm = T)) / sd(data$days_since_hatch, na.rm = T)
}

backztran_DSH <- function(z_day_since_hatch){
  (z_day_since_hatch * sd(data$days_since_hatch, na.rm = T)) + mean(data$days_since_hatch, na.rm = T)
}

#Functions to calculate variance components and h2 and m2 form brms models over x e.g Age
brms_Vcomp <- function(model, x, group_var){
  
  if(group_var == "F1_Genotype"){
    #Strings to search for relevant (co)variance components
    G_vars <- paste0("sd_",group_var)
    G_cors <- paste0("cor_",group_var)
    
    #Extract the relevant sd/(co)variance components 
    SD <- posterior_samples(model, G_vars)
    COR <- posterior_samples(model, G_cors) 
    V <- posterior_samples(model, G_vars)^2
    #Convert correlation to covariance #Cor(1,2) * (SD1 X SD2) = Cov(1,2)
    COV <- cbind(COR[1] * (SD[1] * SD[2]), 
                 COR[2] * (SD[1] * SD[3]),
                 COR[3] * (SD[2] * SD[3]))
    names(COV) <- str_replace(names(COV), "cor", "cov")
    
    # Now, add everything together while accounting for covariances and their respective powers
    Var_comp <- V[1] + (x^2)*V[2] + (x^4)*V[3] +  #The SD of the intercept and linear and quadratic slope
      2*x*COV[1] +    # Covariance of intercept and linear slope
      2*(x^2)*COV[2] + # Covariance of intercept and quadratic slope
      2*(x^3)*COV[3] # Covariance of linear and quadratic slope
    
    df <- data.frame(z_day = x,
                     day = backztran_DSH(x),
                     group_id = group_var,
                     Estimate = posterior_summary(Var_comp)[1],
                     Lower =  posterior_summary(Var_comp)[3],
                     Upper =  posterior_summary(Var_comp)[4])
  }
  
  if(group_var == "dam_id"){
    #Strings to search for relevant (co)variance components
    M_vars <- paste0("sd_",group_var)
    M_cors <- paste0("cor_",group_var)
    
    #Extract the relevant sd/(co)variance components 
    SD <- posterior_samples(model, M_vars)
    COR <- posterior_samples(model, M_cors)
    V <- posterior_samples(model, M_vars)^2
    
    #Convert correlation to covariance #Cor(1,2) * (SD1 X SD2) = Cov(1,2)
    COV <- cbind(COR[1] * (SD[1] * SD[2]), 
                 COR[2] * (SD[1] * SD[3]),
                 COR[3] * (SD[2] * SD[3]))
    names(COV) <- str_replace(names(COV), "cor", "cov")
    
    # Now, add everything together while accounting for covariances and their respective powers
    Var_comp <- V[1] + (x^2)*V[2] + (x^4)*V[3] +  #The SD of the intercept and linear and quadratic slope
      2*x*COV[1] +    # Covariance of intercept and linear slope
      2*(x^2)*COV[2] + # Covariance of intercept and quadratic slope
      2*(x^3)*COV[3] # Covariance of linear and quadratic slope
    
    df <- data.frame(z_day = x,
                     day = backztran_DSH(x),
                     group_id = group_var,
                     Estimate = posterior_summary(Var_comp)[1],
                     Lower =  posterior_summary(Var_comp)[3],
                     Upper =  posterior_summary(Var_comp)[4])
  }
  
  if(group_var == "id"){
    #Strings to search for relevant (co)variance components
    PE_vars <- paste0("sd_",group_var)
    
    #Extract the relevant sd/(co)variance components 
    SD <- posterior_samples(model, PE_vars)
    
    #Squaring SD to get the variance
    Var_comp <- (SD)^2 
    
    df <- data.frame(z_day = x,
                     day = backztran_DSH(x),
                     group_id = group_var,
                     Estimate = posterior_summary(Var_comp)[1],
                     Lower =  posterior_summary(Var_comp)[3],
                     Upper =  posterior_summary(Var_comp)[4])
  }
  
  if(group_var == "sigma"){
    # To get sigma when using a het model we need to extract the fixed effects because we are modelling log(SD)
      SD <- posterior_samples(model, pars = "b_") # Extract the sigma of intercept, linear slope
    
    # Now, compute age-specific SD. Remember that we are modelling log(SD), so we need to exp() at the end to get to the real SD. ## SHIN CHECK.
      SD_comp <- exp(SD[,"b_sigma_Intercept"] + ((x)*SD[,"b_sigma_z_days_since_hatch"])) #The SD of the intercept and linear slope. Remmebr for het models they are modelled as log(SD)
    
    #Squaring SD to get the variance
      Var_comp <- (SD_comp)^2 
    
    df <- data.frame(z_day = x,
                     day = backztran_DSH(x),
                     group_id = group_var,
                     Estimate = posterior_summary(Var_comp)[1],
                     Lower =  posterior_summary(Var_comp)[3],
                     Upper =  posterior_summary(Var_comp)[4])
  }
  
  if(group_var == "total"){
    ##Among ID variance
    #Strings to search for relevant (co)variance components
    G_vars <- paste0("sd_F1_Genotype")
    G_cors <- paste0("cor_F1_Genotype")
    
    #Extract the relevant sd/(co)variance components 
    SD_G <- posterior_samples(model, G_vars)
    COR_G <- posterior_samples(model, G_cors) 
    V_G <- posterior_samples(model, G_vars)^2
    #Convert correlation to covariance #Cor(1,2) * (SD1 X SD2) = Cov(1,2)
    COV_G <- cbind(COR_G[1] * (SD_G[1] * SD_G[2]), 
                   COR_G[2] * (SD_G[1] * SD_G[3]),
                   COR_G[3] * (SD_G[2] * SD_G[3]))
    names(COV_G) <- str_replace(names(COR_G), "cor", "cov")
    
    # Now, add everything together while accounting for covariances and their respective powers
    Vg <- V_G[1] + (x^2)*V_G[2] + (x^4)*V_G[3] +  #The SD of the intercept and linear and quadratic slope
      2*x*COV_G[1] +    # Covariance of intercept and linear slope
      2*(x^2)*COV_G[2] + # Covariance of intercept and quadratic slope
      2*(x^3)*COV_G[3] # Covariance of linear and quadratic slope

    ##Among Dam variance
    #Strings to search for relevant (co)variance components
    M_vars <- paste0("sd_dam_id")
    M_cors <- paste0("cor_dam_id")
    
    #Extract the relevant sd/(co)variance components 
    SD_M <- posterior_samples(model, M_vars)
    COR_M <- posterior_samples(model, M_cors) 
    V_M <- posterior_samples(model, M_vars)^2
    #Convert correlation to covariance #Cor(1,2) * (SD1 X SD2) = Cov(1,2)
    COV_M <- cbind(COR_M[1] * (SD_M[1] * SD_M[2]), 
                   COR_M[2] * (SD_M[1] * SD_M[3]),
                   COR_M[3] * (SD_M[2] * SD_M[3]))
    names(COV_M) <- str_replace(names(COV_M), "cor", "cov")
    
    # Now, add everything together while accounting for covariances and their respective powers
    Vm <- V_M[1] + (x^2)*V_M[2] + (x^4)*V_M[3] +  #The SD of the intercept and linear and quadratic slope
      2*x*COV_M[1] +    # Covariance of intercept and linear slope
      2*(x^2)*COV_M[2] + # Covariance of intercept and quadratic slope
      2*(x^3)*COV_M[3] # Covariance of linear and quadratic slope

       #Residuals

       # To get sigma when using a het model we need to extract the fixed effects because we are modelling log(SD)
      SD_e <- posterior_samples(model, pars = "b_") # Extract the sigma of intercept, linear slope
    
    # Now, compute age-specific SD. Remember that we are modelling log(SD), so we need to exp() at the end to get to the real SD. ## SHIN CHECK.
    SD_comp_e <- exp(SD_e[,"b_sigma_Intercept"] + ((x)*SD_e[,"b_sigma_z_days_since_hatch"])) #The SD of the intercept and linear and quadratic slope
    
    #Squaring SD to get the variance
    Vresid <- (SD_comp_e)^2 
    
    #Calculate total phenotypic variance
    VtotalP <- Vg + Vm +  Vresid
    
    df <- data.frame(z_day = x,
                     day = backztran_DSH(x),
                     group_id = "Total_variance",
                     Estimate = posterior_summary(VtotalP)[1],
                     Lower =  posterior_summary(VtotalP)[3],
                     Upper =  posterior_summary(VtotalP)[4])
    return(df)
  }
  
  
  return(df)
} 

# brms_Vcomp(model = hot_brm_5.5,
#            x = ztran_DsH(seq(0, 500, 30))[10],
#            group_var = "sigma") #it works
# 
#lapply(z_days, function(x) brms_Vcomp(model = hot_brm_5.5, x = x, group_var = "liz_id")) %>% bind_rows() #it works

  
brms_m2 <- function(model, x){
  ##Among ID variance
  #Strings to search for relevant (co)variance components
  G_vars <- paste0("sd_liz_id")
  G_cors <- paste0("cor_liz_id")
  
  #Extract the relevant sd/(co)variance components 
  SD_G <- posterior_samples(model, G_vars)
  COR_G <- posterior_samples(model, G_cors) 
  V_G <- posterior_samples(model, G_vars)^2
  #Convert correlation to covariance #Cor(1,2) * (SD1 X SD2) = Cov(1,2)
  COV_G <- cbind(COR_G[1] * (SD_G[1] * SD_G[2]), 
                 COR_G[2] * (SD_G[1] * SD_G[3]),
                 COR_G[3] * (SD_G[2] * SD_G[3]))
  names(COV_G) <- str_replace(names(COR_G), "cor", "cov")
  
  # Now, add everything together while accounting for covariances and their respective powers
  Vg <- V_G[1] + (x^2)*V_G[2] + (x^4)*V_G[3] +  #The SD of the intercept and linear and quadratic slope
    2*x*COV_G[1] +    # Covariance of intercept and linear slope
    2*(x^2)*COV_G[2] + # Covariance of intercept and quadratic slope
    2*(x^3)*COV_G[3] # Covariance of linear and quadratic slope
  
  ##Among Dam variance
  #Strings to search for relevant (co)variance components
  M_vars <- paste0("sd_dam_id")
  M_cors <- paste0("cor_dam_id")
  
  #Extract the relevant sd/(co)variance components 
  SD_M <- posterior_samples(model, M_vars)
  COR_M <- posterior_samples(model, M_cors) 
  V_M <- posterior_samples(model, M_vars)^2
  #Convert correlation to covariance #Cor(1,2) * (SD1 X SD2) = Cov(1,2)
  COV_M <- cbind(COR_M[1] * (SD_M[1] * SD_M[2]), 
                 COR_M[2] * (SD_M[1] * SD_M[3]),
                 COR_M[3] * (SD_M[2] * SD_M[3]))
  names(COV_M) <- str_replace(names(COV_M), "cor", "cov")
  
  # Now, add everything together while accounting for covariances and their respective powers
  Vm <- V_M[1] + (x^2)*V_M[2] + (x^4)*V_M[3] +  #The SD of the intercept and linear and quadratic slope
    2*x*COV_M[1] +    # Covariance of intercept and linear slope
    2*(x^2)*COV_M[2] + # Covariance of intercept and quadratic slope
    2*(x^3)*COV_M[3] # Covariance of linear and quadratic slope
  
  #Permanent Environment variance
  #Strings to search for relevant (co)variance components
  PE_vars <- paste0("sd_id")
  
  #Extract the relevant sd/(co)variance components 
  SD_PE <- posterior_samples(model, PE_vars)
  
  #Squaring SD to get the variance
  Vpe <- (SD_PE)^2 
  
  #Residuals
  SD_e <- posterior_samples(model, "sigma") # Extract the variance of intercept, linear slope
  
  # Now, add everything together
  SD_comp_e <- SD_e[1] + ((x)*SD_e[2]) #The SD of the intercept and linear and quadratic slope
  
  #Squaring SD to get the variance
  Vresid <- (SD_comp_e)^2 
  
  #Calculate total phenotypic variance
  VtotalP <- Vg + Vm + Vpe + Vresid
  
  #Calculate m2 (proportion of variance explained by maternal effects)
  m2 <- Vm / VtotalP
  
  df <- data.frame(z_day = x,
                   day = backztran_DSH(x),
                   group_id = "m2",
                   Estimate = posterior_summary(m2)[1],
                   Lower =  posterior_summary(m2)[3],
                   Upper =  posterior_summary(m2)[4])
  return(df)
}

# lapply(z_days, function(x) brms_m2(model = hot_brm_5.5, x = x)) %>% bind_rows() #it works
# 
# brms_m2(model = cold_brm_5.5,
#              x = ztran_DsH(seq(0, 500, 30))[10]) #it works

brms_h2 <- function(model, x){
  ##Among ID variance
  #Strings to search for relevant (co)variance components
  G_vars <- paste0("sd_liz_id")
  G_cors <- paste0("cor_liz_id")
  
  #Extract the relevant sd/(co)variance components 
  SD_G <- posterior_samples(model, G_vars)
  COR_G <- posterior_samples(model, G_cors) 
  V_G <- posterior_samples(model, G_vars)^2
  #Convert correlation to covariance #Cor(1,2) * (SD1 X SD2) = Cov(1,2)
  COV_G <- cbind(COR_G[1] * (SD_G[1] * SD_G[2]), 
                 COR_G[2] * (SD_G[1] * SD_G[3]),
                 COR_G[3] * (SD_G[2] * SD_G[3]))
  names(COV_G) <- str_replace(names(COR_G), "cor", "cov")
  
  # Now, add everything together while accounting for covariances and their respective powers
  Vg <- V_G[1] + (x^2)*V_G[2] + (x^4)*V_G[3] +  #The SD of the intercept and linear and quadratic slope
    2*x*COV_G[1] +    # Covariance of intercept and linear slope
    2*(x^2)*COV_G[2] + # Covariance of intercept and quadratic slope
    2*(x^3)*COV_G[3] # Covariance of linear and quadratic slope
  
  ##Among Dam variance
  #Strings to search for relevant (co)variance components
  M_vars <- paste0("sd_dam_id")
  M_cors <- paste0("cor_dam_id")
  
  #Extract the relevant sd/(co)variance components 
  SD_M <- posterior_samples(model, M_vars)
  COR_M <- posterior_samples(model, M_cors) 
  V_M <- posterior_samples(model, M_vars)^2
  #Convert correlation to covariance #Cor(1,2) * (SD1 X SD2) = Cov(1,2)
  COV_M <- cbind(COR_M[1] * (SD_M[1] * SD_M[2]), 
                 COR_M[2] * (SD_M[1] * SD_M[3]),
                 COR_M[3] * (SD_M[2] * SD_M[3]))
  names(COV_M) <- str_replace(names(COV_M), "cor", "cov")
  
  # Now, add everything together while accounting for covariances and their respective powers
  Vm <- V_M[1] + (x^2)*V_M[2] + (x^4)*V_M[3] +  #The SD of the intercept and linear and quadratic slope
    2*x*COV_M[1] +    # Covariance of intercept and linear slope
    2*(x^2)*COV_M[2] + # Covariance of intercept and quadratic slope
    2*(x^3)*COV_M[3] # Covariance of linear and quadratic slope
  
  #Permanent Environment variance
  #Strings to search for relevant (co)variance components
  PE_vars <- paste0("sd_id")
  
  #Extract the relevant sd/(co)variance components 
  SD_PE <- posterior_samples(model, PE_vars)
  
  #Squaring SD to get the variance
  Vpe <- (SD_PE)^2 
  
  #Residuals
  SD_e <- posterior_samples(model, "sigma") # Extract the variance of intercept, linear slope
  
  # Now, add everything together
  SD_comp_e <- SD_e[1] + ((x)*SD_e[2]) #The SD of the intercept and linear and quadratic slope
  
  #Squaring SD to get the variance
  Vresid <- (SD_comp_e)^2 
  
  #Calculate total phenotypic variance
  VtotalP <- Vg + Vm + Vpe + Vresid
  
  #Calculate heritability
  h2 <- Vg / VtotalP
  
  df <- data.frame(z_day = x,
                   day = backztran_DSH(x),
                   group_id = "h2",
                   Estimate = posterior_summary(h2)[1],
                   Lower =  posterior_summary(h2)[3],
                   Upper =  posterior_summary(h2)[4])
  return(df)
}

# brms_h2(model = brm_5.4,
#         x = ztran_DsH(seq(0, 500, 30))[10]) #it works
# 
# lapply(z_days, function(x) brms_h2(model = hot_brm_5.5, x = x)) %>% bind_rows() #it works

get_CV_brms <- function(x, model, group_var){
  #Extract SOL
  model_post <- posterior_samples(model)
  
  #Convert regular day to Z days 
  z <- ztran_DsH(x)
  
  #Make predictions for mean mass at a given age
  pred <- model_post[, "b_Intercept"] + #Intercept
    model_post[, "b_z_days_since_hatch"] *  z + #Linear day effect 
    model_post[, "b_z_days_since_hatch_I2"] * (z^2) #Curve day effect 
  
  #Get variance at a given age
  if(group_var == "F1_Genotype"){
    #Strings to search for relevant (co)variance components
    G_vars <- paste0("sd_",group_var)
    G_cors <- paste0("cor_",group_var)
    
    #Extract the relevant sd/(co)variance components 
      SD <- posterior_samples(model, G_vars)
     COR <- posterior_samples(model, G_cors) 
       V <- posterior_samples(model, G_vars)^2
    
    #Convert correlation to covariance #Cor(1,2) * (SD1 X SD2) = Cov(1,2)
    COV <- cbind(COR[1] * (SD[1] * SD[2]), 
                 COR[2] * (SD[1] * SD[3]),
                 COR[3] * (SD[2] * SD[3]))
    names(COV) <- str_replace(names(COV), "cor", "cov")
   
    # Now, add everything together while accounting for covariances and their powers
    V_g <- V[1] + (z^2)*V[2] + (z^4)*V[3] +  #The variances of the intercept and linear and quadratic slope
      2*z*COV[1] +    # Covariance of intercept and linear slope
      2*(z^2)*COV[2] + # Covariance of intercept and quadratic slope
      2*(z^3)*COV[3] # Covariance of linear and quadratic slope
    
    CV_g <- (100 * (V_g^0.5)) / exp(pred)
    
    #Arrange predictions neatly
    df <- data.frame(   z_day = z,
                          day = x,
                     group_id = group_var,
                     Estimate = posterior_summary(CV_g)[1],
                        Lower =  posterior_summary(CV_g)[3],
                        Upper =  posterior_summary(CV_g)[4])
  }
  
  if(group_var == "dam_id"){
    #Strings to search for relevant (co)variance components
    M_vars <- paste0("sd_",group_var)
    M_cors <- paste0("cor_",group_var)
    
    #Extract the relevant sd/(co)variance components 
      SD <- posterior_samples(model, M_vars)
     COR <- posterior_samples(model, M_cors)
       V <- posterior_samples(model, M_vars)^2
    
    #Convert correlation to covariance #Cor(1,2) * (SD1 X SD2) = Cov(1,2)
    COV <- cbind(COR[1] * (SD[1] * SD[2]), 
                 COR[2] * (SD[1] * SD[3]),
                 COR[3] * (SD[2] * SD[3]))
    names(COV) <- str_replace(names(COV), "cor", "cov")
    
    # Now, add everything together while accounting for covariances and their powers
    V_m <- V[1] + (z^2)*V[2] + (z^4)*V[3] +  #The variances of the intercept and linear and quadratic slope
      2*z*COV[1] +    # Covariance of intercept and linear slope
      2*(z^2)*COV[2] + # Covariance of intercept and quadratic slope
      2*(z^3)*COV[3] # Covariance of linear and quadratic slope
    
    CV_m <- (100 * (V_m^0.5)) / exp(pred)
    
    #Arrange predictions neatly
    df <- data.frame(   z_day = z,
                          day = x,
                     group_id = group_var,
                     Estimate = posterior_summary(CV_m)[1],
                        Lower =  posterior_summary(CV_m)[3],
                        Upper =  posterior_summary(CV_m)[4])
  }
  
  if(group_var == "id"){
    Vpe_vars <- paste0("sd_",group_var)
    Vpe <- (posterior_samples(model, Vpe_vars))^2
    
    CV_pe <- (100 * (Vpe^0.5)) / exp(pred)
    
    #Arrange predictions neatly
    df <- data.frame(z_day = z,
                     day = x,
                     group_id = group_var,
                     Estimate = posterior_summary(CV_pe)[1],
                     Lower =  posterior_summary(CV_pe)[3],
                     Upper =  posterior_summary(CV_pe)[4])
  }
  
  if(group_var == "sigma"){
    
    SD <- posterior_samples(model, pars = "b") # Extract the variance of intercept, linear slope
    
    # Now, add everything together
    SD_comp <- exp(SD[,"b_sigma_Intercept"] + ((x)*SD[,"b_sigma_z_days_since_hatch"])) #The SD of the intercept and linear and quadratic slope
    
    #Squaring SD to get the variance
    Vs <- (SD_comp)^2 
    
    CV_e <- (100 * (Vs^0.5)) / exp(pred)
    
    #Arrange predictions neatly
    df <- data.frame(z_day = z,
                     day = x,
                     group_id = group_var,
                     Estimate = posterior_summary(CV_e)[1],
                     Lower =  posterior_summary(CV_e)[3],
                     Upper =  posterior_summary(CV_e)[4])
  }
  
  if(group_var == "total"){
    ##Among ID variance
    #Strings to search for relevant (co)variance components
    G_vars <- paste0("sd_liz_id")
    G_cors <- paste0("cor_liz_id")
    
    #Extract the relevant sd/(co)variance components 
    SD_G <- posterior_samples(model, G_vars)
    COR_G <- posterior_samples(model, G_cors) 
    V_G <- posterior_samples(model, G_vars)^2
    #Convert correlation to covariance #Cor(1,2) * (SD1 X SD2) = Cov(1,2)
    COV_G <- cbind(COR_G[1] * (SD_G[1] * SD_G[2]), 
                   COR_G[2] * (SD_G[1] * SD_G[3]),
                   COR_G[3] * (SD_G[2] * SD_G[3]))
    names(COV_G) <- str_replace(names(COR_G), "cor", "cov")
    
    # Now, add everything together while accounting for covariances and their respective powers
    Vg <- V_G[1] + (x^2)*V_G[2] + (x^4)*V_G[3] +  #The SD of the intercept and linear and quadratic slope
      2*x*COV_G[1] +    # Covariance of intercept and linear slope
      2*(x^2)*COV_G[2] + # Covariance of intercept and quadratic slope
      2*(x^3)*COV_G[3] # Covariance of linear and quadratic slope
    
    ##Among Dam variance
    #Strings to search for relevant (co)variance components
    M_vars <- paste0("sd_dam_id")
    M_cors <- paste0("cor_dam_id")
    
    #Extract the relevant sd/(co)variance components 
    SD_M <- posterior_samples(model, M_vars)
    COR_M <- posterior_samples(model, M_cors) 
    V_M <- posterior_samples(model, M_vars)^2
    #Convert correlation to covariance #Cor(1,2) * (SD1 X SD2) = Cov(1,2)
    COV_M <- cbind(COR_M[1] * (SD_M[1] * SD_M[2]), 
                   COR_M[2] * (SD_M[1] * SD_M[3]),
                   COR_M[3] * (SD_M[2] * SD_M[3]))
    names(COV_M) <- str_replace(names(COV_M), "cor", "cov")
    
    # Now, add everything together while accounting for covariances and their respective powers
    Vm <- V_M[1] + (x^2)*V_M[2] + (x^4)*V_M[3] +  #The SD of the intercept and linear and quadratic slope
      2*x*COV_M[1] +    # Covariance of intercept and linear slope
      2*(x^2)*COV_M[2] + # Covariance of intercept and quadratic slope
      2*(x^3)*COV_M[3] # Covariance of linear and quadratic slope
    
    #Residuals
    SD_e <- posterior_samples(model, pars = "b") # Extract the variance of intercept, linear slope
    
    # Now, add everything together
    SD_comp_e <- exp(SD[,"b_sigma_Intercept"] + ((x)*SD[,"b_sigma_z_days_since_hatch"])) #The SD of the intercept and linear and quadratic slope
    
    #Squaring SD to get the variance
    Vresid <- (SD_comp_e)^2 
    
    #Calculate total phenotypic variance
    VtotalP <- Vg + Vm + Vresid
    
    CVtot <- (100 * (VtotalP^0.5)) / exp(pred)
    
    #Arrange predictions neatly
    df <- data.frame(   z_day = z,
                          day = x,
                     group_id = group_var,
                     Estimate = posterior_summary(CVtot)[1],
                        Lower =  posterior_summary(CVtot)[3],
                        Upper =  posterior_summary(CVtot)[4])
    
  }
  return(df)
}

# get_CV_brms(model = brm_5.4,
#         x = ztran_DsH(seq(0, 500, 30))[10], group_var = "total") #it works
# 
# lapply(z_days, function(x) get_CV_brms(model = hot_brm_5.4, x = x, group_var = "total")) %>% bind_rows() #it works

get_CV_X2 <- function(x, model, group_var){
  #Extract posterior
  model_post <- posterior_samples(model)
  
  #Convert regular day to Z days 
  z <- ztran_DsH(x)
  #Make predictions for mean mass at a given age
  pred <- model_post[, "b_Intercept"] + #Intercept
    model_post[, "b_z_days_since_hatch"] *  z + #Linear day effect 
    model_post[, "b_z_days_since_hatch_I2"] * (z^2) #Curve day effect 
  
  #Get variance at a given age
  #Among ID
    #Strings to search for relevant (co)variance components
    G_vars <- paste0("sd_","F1_Genotype")
    G_cors <- paste0("cor_","F1_Genotype")
    
    #Extract the relevant sd/(co)variance components 
    G_SD <- posterior_samples(model, G_vars)
    G_COR <- posterior_samples(model, G_cors) 
    G_V <- posterior_samples(model, G_vars)^2
    #Convert correlation to covariance #Cor(1,2) * (SD1 X SD2) = Cov(1,2)
    G_COV <- cbind(G_COR[1] * (G_SD[1] * G_SD[2]), 
                   G_COR[2] * (G_SD[1] * G_SD[3]),
                   G_COR[3] * (G_SD[2] * G_SD[3]))
    names(G_COV) <- str_replace(names(G_COR), "cor", "cov")
    
    # Now, add everything together while accounting for covariances and their powers
    V_g <- G_V[1] + (z^2)*G_V[2] + (z^4)*G_V[3] +  #The variances of the intercept and linear and quadratic slope
      2*z*G_COV[1] +    # Covariance of intercept and linear slope
      2*(z^2)*G_COV[2] + # Covariance of intercept and quadratic slope
      2*(z^3)*G_COV[3] # Covariance of linear and quadratic slope
    
    CV_g <- (100 * (V_g^0.5)) / exp(pred)
    
    #Strings to search for relevant (co)variance components
    M_vars <- paste0("sd_","dam_id")
    M_cors <- paste0("cor_","dam_id")
    
    #Extract the relevant sd/(co)variance components 
     M_SD <- posterior_samples(model, M_vars)
    M_COR <- posterior_samples(model, M_cors)
      M_V <- posterior_samples(model, M_vars)^2
    
    #Convert correlation to covariance #Cor(1,2) * (SD1 X SD2) = Cov(1,2)
    M_COV <- cbind(M_COR[1] * (M_SD[1] * M_SD[2]), 
                   M_COR[2] * (M_SD[1] * M_SD[3]),
                   M_COR[3] * (M_SD[2] * M_SD[3]))
    names(M_COV) <- str_replace(names(M_COR), "cor", "cov")
    
    # Now, add everything together while accounting for covariances and their powers
    V_m <- M_V[1] + (z^2)*M_V[2] + (z^4)*M_V[3] +  #The variances of the intercept and linear and quadratic slope
           2*z*M_COV[1] +    # Covariance of intercept and linear slope
           2*(z^2)*M_COV[2] + # Covariance of intercept and quadratic slope
           2*(z^3)*M_COV[3] # Covariance of linear and quadratic slope
    
    CV_m <- (100 * (V_m^0.5)) / exp(pred)
    
    # Sigma / Residual Variance
    Sigma_SD <- posterior_samples(model, pars = "b") # Extract the variance of intercept, linear slope
    
    # Now, add everything together
    Sigma_comp <- exp(SD[,"b_sigma_Intercept"] + ((x)*SD[,"b_sigma_z_days_since_hatch"])) #The SD of the intercept and linear and quadratic slope
    
    #Squaring SD to get the variance
    Vs <- (Sigma_comp)^2 
    
    CV_e <- (100 * (Vs^0.5)) / exp(pred)
    
  #Calculate heritability
    if(group_var == "h2"){
      CV_X2 = (CV_g / (CV_g + CV_m + CV_e))
      
      df <- data.frame(   z_day = z,
                            day = backztran_DSH(z),
                       group_id = group_var,
                       Mean_age = posterior_summary(pred)[1],
                       Estimate = posterior_summary(CV_X2)[1],
                          Lower =  posterior_summary(CV_X2)[3],
                          Upper =  posterior_summary(CV_X2)[4])
      
    }
    if(group_var == "m2"){
      CV_X2 = (CV_m / (CV_g + CV_m + CV_e))
      
      df <- data.frame(   z_day = z,
                            day = backztran_DSH(z),
                       group_id = group_var,
                       Mean_age = posterior_summary(pred)[1],
                       Estimate = posterior_summary(CV_X2)[1],
                          Lower =  posterior_summary(CV_X2)[3],
                          Upper =  posterior_summary(CV_X2)[4])
    }
    return(df)
}



                  