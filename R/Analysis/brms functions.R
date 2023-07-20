###########################
# Data Creation Functions
###########################

#' @title ztran_DsH
#' @description Takes a variable (in this case age) and z-transforms it to be centered on 0 and have SD units
#' @param day The variable in days
#' @return Returns z-transformed age (in SD units and cenetred on 0)
ztran_DsH <- function(day) {
  (day - mean(data$days_since_hatch, na.rm = T)) / sd(data$days_since_hatch, na.rm = T)
}

#' @title backztran_DSH
#' @description Takes a z-transformed variable (in this case age) and converts it from the z-scale to the raw scale (in days)
#' @param z_day_since_hatch The z-score data 
#' @return Returns age in days
backztran_DSH <- function(z_day_since_hatch){
  (z_day_since_hatch * sd(data$days_since_hatch, na.rm = T)) + mean(data$days_since_hatch, na.rm = T)
}

#' @title extract_V
#' @description Extracts the variances and correlations for intercepts and slopes at a given random effect level from a `brms` model object and converts SDs to variances and cors to covariances.
#' @param model The model that explicitly models 'sigma' (resiudal) across z-scaled age. 
#' @param level A character string that specifies what the grouping variable used to estimate the random effect was called (e.g., "F1_Genotype" is the random effect level for estimating additive genetic variance)
#' @return Returns a list with two dataframes. V is the dataframe of variance compoenents for the random intercept, linear slope and quadratic slope and COV is the covariance between intercept, linear and quadratic slopes
extract_V <- function(model, level = "F1_Genotype"){
        
        #Strings to search for relevant (co)variance components
            vars <- paste0("sd_", level)
            cors <- paste0("cor_", level)

        #Extract the relevant sd/(co)variance components
              SD <- posterior_samples(model, vars)
             COR <- posterior_samples(model, cors) 
               V <- SD^2

        #Convert correlation to covariance #Cor(1,2) * (SD1 X SD2) = Cov(1,2)
                   COV <-  cbind(COR[1] * (SD[1] * SD[2]), 
                                 COR[2] * (SD[1] * SD[3]),
                                 COR[3] * (SD[2] * SD[3]))
            
        #Change column names to reflect the values in the dataframe
            names(COV) <- str_replace(names(COR), "cor", "cov")
              names(V) <- str_replace(names(V), "sd", "v")
        
        # Return a list of V and COV which can be used to calculate variance across age
            return(list(V = V, COV = COV))
}

#' @title generate_data
#' @description Creates a dataframe of day, variable type of the grouping ID (random effect) and the estaimte along with 95% credible interval
#' @param z_age The z-scaled age variable
#' @param type A character string that specifies what the variable is (i.e., h2, m2 or some variance component - e.g., dam_id)
#' @return Returns a dataframe containing the z-scaled age, backtransformed age (in days), the type and the mean estimate and upper and lower 95% credible interval for the variable of interest.
generate_data <- function(x, z_age, type){
        data.frame(   z_day = z_age,
                        day = backztran_DSH(z_age),
                   group_id = type,
                   Estimate = posterior_summary(x)[1],
                      Lower = posterior_summary(x)[3],
                      Upper = posterior_summary(x)[4])
}

###########################
# Model Checking Function
###########################

#' @title brms_model_check
#' @description Checks a 'brms' model by plotting a histogram of residuals and a scatterplot of the observed and predcited log mass
#' @param model The 'brms' model object that models log (mass (grams))
#' @param main Title of the plot
#' @param xlab Label of the x-axis which defaults to 'Residuals'
#' @return Returns z-transformed age (in SD units and cenetred on 0)
brms_model_check <- function(model, main = NULL, xlab = "Residuals"){
  
  # Histogram of residuals - assumed normal - pretty good to me
      resid <-  model$data$lnMass - predict(model, summary = TRUE)[,"Estimate"]
 
  # Look at the residuals - should be normally distributed
     hist(resid, main = main, xlab = xlab)
  
  # We already know roughly from R2 that model does good job predciting observed response, but lets have a look. Little bit of over/underpredciting but nothing serious
      plot(model$data$lnMass ~ predict(model, summary = TRUE)[,"Estimate"], ylab = "Observed lnMass", xlab = "Predicted lnMass", main = main)
      abline(0,1, col = "red")
}

###########################
# Calculation Functions
###########################

#' @title calc_V_across_age
#' @description Calculates the change in variance across age using the variance and covariance in intercept, linear slope and quadratic slope.
#' @param z_age The z-scaled age variable
#' @param V The dataframe containing the estimated variances in intercept, linear slope and quadratic slope
#' @param COV The dataframe containing the covariance between intercept, linear and quadratic slope
#' @return Returns a vector with the estimated variance at a given age
calc_V_across_age <- function(z_age, V, COV){

  # Calculating expected changes in variance across a continuous variable
    V <- V[1] + (z_age^2) * V[2] + (z_age^4) * V[3] +  # The SD of the intercept and linear and quadratic slope
           2 * z_age * COV[1]                       +  # Covariance of intercept and linear slope
           2 * (z_age^2) * COV[2]                   +  # Covariance of intercept and quadratic slope
           2 * (z_age^3) * COV[3]                      # Covariance of linear and quadratic slope
    
    return(V)
}

#' @title create_h_m2
#' @description Calculates heritability (h2) and maternal effects (m2) proportions across age.
#' @param z_age The z-scaled age variable
#' @param G The estimate of additive genetic variance
#' @param M The estimate of maternal effect variance
#' @param E The estimated environmental/residual variance
#' @return Returns a dataframe with the estimated h2 and m2 at a given age along with upper and lower credible intervals
create_h_m2 <- function(z_age, G, M, E, type = c("h2", "m2")){
   type <-  match.arg(type)
   
   if(type == "h2"){
        h2 <- G / (G + M + E)
        df <- generate_data(h2, z_age, type = type)
   }

   if(type == "m2"){
        m2 <- M / (G + M + E)
        df <- generate_data(m2, z_age, type = type)
   }
    
  return(df)
}

#' @title get_Vr_across_age
#' @description Calculates the environmental/residual variance across age.
#' @param model The model that explicitly models 'sigma' (resiudal) across z-scaled age. 
#' @param z_age The z-scaled age variable
#' @return Returns a vector with the estimated environmental/resiudal variance at a given age
get_Vr_across_age <- function(model, z_age){
  # Extract the sigma of intercept, linear slope 
       log_SD_e <- posterior_samples(model, pars = "b_") 
  
  # Calculate the SD across age intercept and linear slope. exp() because SD is on log scale
           SD_e <- exp(log_SD_e[,"b_sigma_Intercept"] + (z_age * log_SD_e[,"b_sigma_z_days_since_hatch"])) 
  
  #Squaring SD to get the variance
             Ve <- (SD_e)^2 
  
  # Return Vr
             return(Ve)
}

###########################
## Fonti's Core Functions
###########################

#Functions to calculate variance components and h2 and m2 form brms models over x e.g Age
brms_Vcomp <- function(model, x, group_var){
  
  if(group_var == "F1_Genotype"){
               G <- extract_V(model, level = group_var)
    G_across_age <- calc_V_across_age(x, G[["V"]], G[["COV"]])
              df <- generate_data(G_across_age, x, type = group_var)
  }
  
  if(group_var == "dam_id"){
               M <- extract_V(model, level = group_var)
    M_across_age <- calc_V_across_age(x, G[["V"]], G[["COV"]])
              df <- generate_data(M_across_age, x, type = group_var)
  }
 
  if(group_var == "sigma"){
               E <-  get_Vr_across_age(model, x)
              df <-  generate_data(E, x, type = group_var)
  }
  
  if(group_var == "total"){   
    #Calculate total phenotypic variance
    VtotalP <- VG + VM + VE
    df <-  generate_data(VtotalP, x, type = group_var)
  }
  return(df)
} 

brms_m2 <- function(model, z_age, G = "F1_Genotype", M = "dam_id") {
    # G - Additive genetic variance
                 G <- extract_V(model, level = G)
      G_across_age <- calc_V_across_age(z_age, G[["V"]], G[["COV"]])

    # M - Maternal effect variance
                 M <- extract_V(model, level = M)
      M_across_age <- calc_V_across_age(z_age, M[["V"]], M[["COV"]])

    # R - Environmental variance
      E_across_age <- get_Vr_across_age(model, z_age)

    # Create the data frame
                df <- create_h_m2(z_age, G_across_age, M_across_age, E_across_age, type = "m2")
                return(df)
}               

brms_h2 <- function(model, z_age, G = "F1_Genotype", M = "dam_id") {
    # G - Additive genetic variance
                 G <- extract_V(model, level = G)
      G_across_age <- calc_V_across_age(z_age, G[["V"]], G[["COV"]])

    # M - Maternal effect variance
                 M <- extract_V(model, level = M)
      M_across_age <- calc_V_across_age(z_age, M[["V"]], M[["COV"]])

    # R - Environmental variance
      E_across_age <- get_Vr_across_age(model, z_age)

    # Create the data frame
                df <- create_h_m2(z_age, G_across_age, M_across_age, E_across_age, type = "h2")
                return(df)
}               

func_growth_predictions <- function(day, predat, posterior){
  
  func <- function(liz_id = predat$F1_Genotype[i]){
    # finding dam_id that corresponds to liz_id
    dam_id <- predat$dam_id[predat$F1_Genotype == liz_id]
    treat <- predat$treatment[predat$F1_Genotype == liz_id]
    
    pred <- posterior[, "b_Intercept"] + #Intercept
      posterior[, "b_z_days_since_hatch"] * ztran_DsH(day) + #Linear day effect 
      posterior[, "b_z_days_since_hatch_I2"] * ((ztran_DsH(day))^2) + #Curve day effect 
      posterior[, paste0("r_F1_Genotype[",liz_id, ",Intercept]")] + #Lizard intercept
      posterior[, paste0("r_F1_Genotype[",liz_id, ",z_days_since_hatch]")] + #Lizard slope
      posterior[, paste0("r_F1_Genotype[",liz_id, ",z_days_since_hatch_I2]")] + #Lizard curve z_days_since_hatch_I2
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
  df <- bind_rows(lapply(predat$F1_Genotype, function(x) func(x))) %>% mutate(day = day)
  return(df)
}

############################
## Functions Below Still need to be simplified. 
###########################
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