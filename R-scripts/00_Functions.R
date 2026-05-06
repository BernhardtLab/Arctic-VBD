## Lilian Chan, University of Guelph
## Arctic vector-borne disease transmission suitability model
##
## Purpose: Functions to load for use in other scripts
## 
## Code adapted from https://github.com/JoeyBernhardt/anopheles-rate-summation/blob/master/R-scripts/working-versions-code/00_RSProjectFunctions.R

##	Table of contents:
##  0. Set-up workspace
##	1. Functions to process TPC model output
##			A. Function to calculate TPC posterior summary statistics across a temperature gradient
##			B. Function to extract full posterior distributions for 3 mean-defining TPC parameters & calculate Tbreadth
##			C. Function to calculate Topt
##			D. Function to calculate Tmin, Max, and Tbreadth for derived TPCs 
##			E. Wrapper function to calculate summary data for and extract parameter posteriors from JAGS fitted TPCs
##			F. Wrapper function to calculate summary data for derived TPCs
##  2. Function to process trait data for plotting
##  3. Function to calculate relative suitability
##  4. Functions for Sensitivity Analysis
##			A. Function for derivative of Briere thermal response
##			B. Function for derivative of quadratic thermal response
##			C. Function for sensitivity analysis - partial derivatives




# 0. Set-up workspace ----------------------------------------------------------

library(tidyverse)



# 1. Functions to process TPC model output -------------------------------------

## A. Function to calculate TPC posterior summary statistics across a temperature gradient ----

calcPostQuants <- function(TPC_predictions, trait_name, temp_gradient) {
  # TPC_predictions: 15000 rows (MCMC iterations) x 451 cols (temp_gradient, 0-45C at 0.1C interval)
  # output: 8 cols (temperature, lowerCI, upperCI, lowerQ, upperQ, mean, median, trait_name) x 451 cols
  
  # Reassign column names to the temperature gradient
  colnames(TPC_predictions) <- temp_gradient
  
  output <- TPC_predictions %>% 
    mutate(iteration = rownames(.)) %>% # add column with iteration number 
    pivot_longer(!iteration, names_to = "temperature", values_to = "trait_value") %>% # convert to long format (3 columns: iteration, temp, & trait value)
    group_by(temperature) %>%
    summarise(lowerCI = quantile(trait_value, probs = 0.025),
              upperCI = quantile(trait_value, probs = 0.975),
              lowerQ = quantile(trait_value, probs = 0.25),
              upperQ = quantile(trait_value, probs = 0.75),
              mean = mean(trait_value),
              median = median(trait_value)) %>% 
    mutate(temperature = as.numeric(temperature)) %>% # make temperature numeric
    arrange(temperature) %>% # re-order rows by ascending temperature (since grouping made it categorical, it is ordered alphabetical)
    mutate(trait = trait_name) # add column with trait name
  
  return(output) # return output
  
}

## B. Function to extract full posterior distributions for 3 mean-defining TPC parameters (q, Tmin, Tmax) & calculate Tbreadth ----

getTPCParamFullPosts <- function (TPC_model, trait_name) {
  # TPC_model is the output from JAGS fitted TPCs
  
  # Extract the Tmin, Tmax, q, and Tbreadth (Tmax - Tmin) of each iteration
  output <- data.frame(iteration = seq(1,length(TPC_model$BUGSoutput$sims.list$cf.T0),1),
                       cf.T0 = TPC_model$BUGSoutput$sims.list$cf.T0[,1],
                       cf.Tm = TPC_model$BUGSoutput$sims.list$cf.Tm[,1],
                       cf.q = TPC_model$BUGSoutput$sims.list$cf.q[,1],
                       Tbreadth = (TPC_model$BUGSoutput$sims.list$cf.Tm[,1] - TPC_model$BUGSoutput$sims.list$cf.T0[,1]),
                       trait = trait_name)
  
  return(output)
}


## C. Function to calculate Topt ----
calcToptQuants <- function(TPC_predictions, trait_name, temp_gradient) {
  # TPC_predictions: 15000 rows (MCMC iterations) x 451 cols (temp_gradient, 0-45C at 0.1C interval)
  
  # Reassign column names to the temperature gradient
  colnames(TPC_predictions) <- temp_gradient
  
  output <- TPC_predictions %>% 
    mutate(iteration = rownames(.)) %>% # add column with iteration number 
    pivot_longer(!iteration, names_to = "temperature", values_to = "trait_value") %>% # convert to long format (3 columns: iteration, temp, & trait value)
    group_by(iteration) %>% 
    slice_max(order_by = trait_value, n = 1, with_ties = FALSE) %>% # for each iteration, select row with highest value for the trait
    ungroup() %>% 
    mutate(temperature = as.numeric(temperature)) %>% # make temperature numeric
    summarise(mean = mean(temperature),
              sd = sd(temperature),
              lowerCI = quantile(temperature, probs = 0.025),
              upperCI = quantile(temperature, probs = 0.975),
              lowerQ = quantile(temperature, probs = 0.25),
              upperQ = quantile(temperature, probs = 0.75),
              median = median(temperature)) %>% 
    mutate(term = "Topt") %>% # add column with calculation type
    mutate(trait = trait_name) # add column with trait name
  
  return(output) # return output
  
}


## D. Function to calculate posteriors for Tmin, Tmax, and Tbreadth for derived TPCs ----

# Use this function for TPCs that are not from JAGS fitted TPCs (e.g. output from the suitability model)

calcDerivedTPCParamPosteriors <- function(TPC_predictions, temp_gradient) {
  # TPC_predictions: 15000 rows (MCMC iterations) x 451 cols (temp_gradient, 0-45C at 0.1C interval)
  
  # Total number of MCMC iterations = 15000
  output_length <- nrow(TPC_predictions)
  # Create a dataframe to store the output
  output.df <- data.frame("cf.T0" = numeric(nrow(TPC_predictions)), 
                          "cf.Tm" = numeric(nrow(TPC_predictions)),
                          "Tbreadth" = numeric(nrow(TPC_predictions)))
  
  # Get length of temp gradient to compare to final index of Tmax
  length_temp_gradient <- length(temp_gradient)
  
  # Loop through each row of TPC_predictions (MCMC step of input traits)
  for (i in 1:output_length) {
    
    ## Create vector of list of indices where S > 0
    index.list <- which(TPC_predictions[i,] > 0)
    length.index.list <- length(index.list) # get the length of index.list 
    
    # Store Tmin
    ifelse(
      # If S > 0 at the lowest temp (should be 0.0ºC, or index 1), then 0ºC will be T0
      index.list[1] == 1, output.df$cf.T0[i] <- temp_gradient[1], 
      # Otherwise, store the corresponding temp of the index before the first value of index.list
      output.df$cf.T0[i] <- temp_gradient[index.list[1] - 1])
    
    
    # Store Tm 
    ifelse(
      # If S > 0 at the highest temp (should be 45.0ºC, or index 451), then 45.0ºC will be Tmax
      index.list[length.index.list] == length_temp_gradient, output.df$cf.Tm[i] <- temp_gradient[index.list[length.index.list]], 
      # Otherwise, store the corresponding temp of the index after the last value of index.list
      output.df$cf.Tm[i] <- temp_gradient[index.list[length.index.list] + 1])
    
    
    # Calculate Tbreadth from Tmin and Tmax
    output.df$Tbreadth[i] <- output.df$cf.Tm[i] - output.df$cf.T0[i]
  }
  
  output.df # return
  
}


## E. Wrapper function to calculate summary data for and extract parameter posteriors from JAGS fitted TPCs ----

## Output: a list of 3 items
##   i. TPC_pred_summary: Posterior summary of TPC predictions across temperatures (output from calcPostQuants)
##  ii. TPC_param_summary_all: Summary statistics of TPC parameters
## iii. TPC_param_full_posts: Full posterior distributions for q, Tmin, Tmax, and Tbreadth (output from getTPCParamFullPosts)

extractTPC <- function(TPC_model, trait_name, temp_gradient) {
  
  # Extract predicted trait values over the temperature gradient
  TPC_predictions <- as.data.frame(TPC_model$BUGSoutput$sims.list$z.trait.mu.pred)
  
  # Calculate TPC posterior summary statistics (means & quantiles)
  TPC_pred_summary <- calcPostQuants(TPC_predictions, trait_name, temp_gradient)
  
  # Extract full posterior distribution for 3 mean-defining TPC parameters (q, Tmin, Tmax) + calculate Tbreadth for each iteration
  TPC_param_full_posts <- getTPCParamFullPosts(TPC_model, trait_name)
  
  # Calculate Tbreadth summary statistics (mean, sd, & quantiles)
  Tbreadth_summary <- data.frame(term = "Tbreadth",
                                 mean = mean(TPC_param_full_posts$Tbreadth),
                                 sd = sd(TPC_param_full_posts$Tbreadth),
                                 lowerCI = quantile(TPC_param_full_posts$Tbreadth, 0.025)[[1]],
                                 median =  quantile(TPC_param_full_posts$Tbreadth, 0.5)[[1]],
                                 upperCI = quantile(TPC_param_full_posts$Tbreadth, 0.975)[[1]],
                                 trait = trait_name)
  
  # Calculate Topt for each iteration and calculate summary statistics (mean, sd, & quantiles)
  Topt_summary <- calcToptQuants(TPC_predictions, trait_name, temp_gradient)
  
  # Pull out parameter summary from the fitted model 
  TPC_param_summary <- as.data.frame(TPC_model$BUGSoutput$summary[1:5,]) %>%
    rownames_to_column(var = "term") %>%
    rename(lowerCI = `2.5%`, median = `50%`, upperCI = `97.5%`) %>% # Rename columns so they are easier to reference & can merge with Topt quantiles
    mutate(trait = trait_name)
  
  # Remove unwanted columns (25% quantile, 75% quantile, Rhat, and n.eff)
  TPC_param_summary <- TPC_param_summary %>%
    dplyr::select(term, mean, sd, lowerCI, median, upperCI, trait)
  
  # Add Topt and Tbreadth to parameters summary data frame
  TPC_param_summary_all <- bind_rows(TPC_param_summary, Topt_summary, Tbreadth_summary)
  
  # Bundle output in a list
  output_list <- list(TPC_pred_summary, TPC_param_summary_all, TPC_param_full_posts)
  
  return(output_list) # return output
  
}

## Use the following function if the JAGS fitted TPC includes random effects
extractTPC_raneff <- function(TPC_model, trait_name, temp_gradient) {
  
  # Extract predicted trait values over the temperature gradient
  TPC_predictions <- as.data.frame(TPC_model$BUGSoutput$sims.list$z.trait.mu.pred.pop) # extract the global-level TPC
  
  # Calculate TPC posterior summary statistics (means & quantiles)
  TPC_pred_summary <- calcPostQuants(TPC_predictions, trait_name, temp_gradient)
  
  # Extract full posterior distribution for 3 mean-defining TPC parameters + calculate Tbreadth for each iteration
  TPC_param_full_posts <- getTPCParamFullPosts(TPC_model, trait_name)
  
  # Calculate Tbreadth summary statistics (mean, sd, & quantiles)
  Tbreadth_summary <- data.frame(term = "Tbreadth",
                                 mean = mean(TPC_param_full_posts$Tbreadth),
                                 sd = sd(TPC_param_full_posts$Tbreadth),
                                 lowerCI = quantile(TPC_param_full_posts$Tbreadth, 0.025)[[1]],
                                 median =  quantile(TPC_param_full_posts$Tbreadth, 0.5)[[1]],
                                 upperCI = quantile(TPC_param_full_posts$Tbreadth, 0.975)[[1]],
                                 trait = trait_name)
  
  # Calculate Topt for each iteration and calculate summary statistics (mean, sd, & quantiles)
  Topt_summary <- calcToptQuants(TPC_predictions, trait_name, temp_gradient)
  
  # Pull out parameter summary from the fitted model 
  TPC_param_summary <- as.data.frame(TPC_model$BUGSoutput$summary[1:8,]) %>%
    rownames_to_column(var = "term") %>%
    rename(lowerCI = `2.5%`, median = `50%`, upperCI = `97.5%`) %>% # Rename columns so they are easier to reference & can merge with Topt quantiles
    mutate(trait = trait_name)
  
  # Remove unwanted columns (25% quantile, 75% quantile, Rhat, and n.eff)
  TPC_param_summary <- TPC_param_summary %>%
    dplyr::select(term, mean, sd, lowerCI, median, upperCI, trait)
  
  # Add Topt and Tbreadth to parameters summary data frame
  TPC_param_summary_all <- bind_rows(TPC_param_summary, Topt_summary, Tbreadth_summary)
  
  # Bundle output in a list
  output_list <- list(TPC_pred_summary, TPC_param_summary_all, TPC_param_full_posts)
  
  return(output_list) # return output
  
}


## F. Wrapper function to calculate summary data for derived TPCs ----
extractDerivedTPC <- function(TPC_predictions, trait_name, temp_gradient) {
  # TPC_predictions: 15000 rows (MCMC iterations) x 451 cols (temp_gradient, 0-45C at 0.1C interval)
  
  # Calculate Tmin, Tmax, and Tbreadth posteriors
  TPC_param_full_posts <- calcDerivedTPCParamPosteriors(TPC_predictions, temp_gradient)
  
  # Calculate Tmin, Tmax, and Tbreadth summary statistics (mean, sd, & quantiles)
  Tmin_summary <- data.frame(term = "cf.T0",
                             mean = mean(TPC_param_full_posts$cf.T0),
                             sd = sd(TPC_param_full_posts$cf.T0),
                             lowerCI = quantile(TPC_param_full_posts$cf.T0, 0.025)[[1]],
                             lowerQ = quantile(TPC_param_full_posts$cf.T0, 0.25)[[1]],
                             median =  quantile(TPC_param_full_posts$cf.T0, 0.5)[[1]],
                             upperQ = quantile(TPC_param_full_posts$cf.T0, 0.75)[[1]],
                             upperCI = quantile(TPC_param_full_posts$cf.T0, 0.975)[[1]],
                             trait = trait_name)
  
  Tmax_summary <- data.frame(term = "cf.Tm",
                             mean = mean(TPC_param_full_posts$cf.Tm),
                             sd = sd(TPC_param_full_posts$cf.Tm),
                             lowerCI = quantile(TPC_param_full_posts$cf.Tm, 0.025)[[1]],
                             lowerQ = quantile(TPC_param_full_posts$cf.Tm, 0.25)[[1]],
                             median =  quantile(TPC_param_full_posts$cf.Tm, 0.5)[[1]],
                             upperQ = quantile(TPC_param_full_posts$cf.Tm, 0.75)[[1]],
                             upperCI = quantile(TPC_param_full_posts$cf.Tm, 0.975)[[1]],
                             trait = trait_name)
  
  Tbreadth_summary <- data.frame(term = "Tbreadth",
                                 mean = mean(TPC_param_full_posts$Tbreadth),
                                 sd = sd(TPC_param_full_posts$Tbreadth),
                                 lowerCI = quantile(TPC_param_full_posts$Tbreadth, 0.025)[[1]],
                                 lowerQ = quantile(TPC_param_full_posts$Tbreadth, 0.25)[[1]],
                                 median =  quantile(TPC_param_full_posts$Tbreadth, 0.5)[[1]],
                                 upperQ = quantile(TPC_param_full_posts$Tbreadth, 0.75)[[1]],
                                 upperCI = quantile(TPC_param_full_posts$Tbreadth, 0.975)[[1]],
                                 trait = trait_name)
  
  # Calculate Topt for each iteration and calculate summary statistics (mean, sd, & quantiles)
  Topt_summary <- calcToptQuants(TPC_predictions, trait_name, temp_gradient)
  
  # Add Topt and Tbreadth to parameters summary data frame
  TPC_param_summary_all <- bind_rows(Tmin_summary, Tmax_summary, Topt_summary, Tbreadth_summary)
  
  return(TPC_param_summary_all) # return output
  
}


# 2. Function to process trait data for plotting -------------------------------
processTraitData <- function (data_input, trait_name) {
  
  output <- data_input %>%
    group_by(temp) %>% 
    summarise(n = n(),
              mean = mean(trait),
              std_error = ifelse(n >= 3, sd(trait)/sqrt(n()), NA)) %>% 
    mutate(trait = trait_name)
  
  return(output)
}


# 3. Function to calculate relative suitability --------------------------------

# Creating a small constant to keep denominators from being zero
ec <- 0.000001

S = function(a, bc, lf, PDR, EFGC, EV, pLA, MDR){
  (a^3 * bc * exp(-(1/(lf+ec))*(1/(PDR+ec))) * EFGC * EV * pLA * MDR * lf^3)^0.5
}


# 4. Functions for Sensitivity Analysis ----------------------------------------

## A. Function for derivative of Briere thermal response ----

d_briere = function(T, T0, Tm, q) {
  
  b <- c()
  
  for (i in 1:length(T)) {
    if (T[i]>T0 && T[i]<Tm) ## When trait value > 0
      {b[i] <- (q*(-5*(T[i]^2) + 3*T[i]*T0 + 4*T[i]*Tm - 2*T0*Tm)/(2*sqrt(Tm-T[i])))}
    else {b[i] <- 0}
  }
  
  b # return output
  
}

## B. Function for derivative of quadratic thermal response ----
d_quad = function(T, T0, Tm, q){
  
  b <- c()
  
  for (i in 1:length(T)){
    if (T[i]>T0 && T[i]<Tm) # When trait value > 0
      {b[i] <- -1*q*(2*T[i] - T0 - Tm)}
    else {b[i] <- 0}
  }
  
  b # return output
  
}


## C. Function for sensitivity analysis - partial derivatives ----

# Arguments: 
#      x_preds = full posterior distribution for TPC predicted values at each temp interval. 
#                15000 rows (MCMC iterations) x 451 cols (temp_gradient, 0-45C at 0.1C interval);
#			 m_param = full posterior distribution for q, Tmin, and Tmax;
#			 m_x = the mean value for each trait over the temperature gradient

SensitivityAnalysis_pd = function(# predicted trait values:
                                  a_preds, bc_preds, lf_preds, PDR_preds, 
                                  EFGC_preds, EV_preds, pLA_preds, MDR_preds,
                                  # Full posterior distribution of TPC parameters:
                                  a_param, bc_param, lf_param, PDR_param, 
                                  EFGC_param, EV_param, pLA_param, MDR_param,
                                  # Trait means at each temperature
                                  m_a, m_bc, m_lf, m_PDR, m_EFGC, m_EV, 
                                  m_pLA, m_MDR) {
  
  # Create matrices to hold results
  dS.da <- dS.dbc <- dS.dlf <- dS.dPDR <- dS.dEFGC <- dS.dEV <- dS.dpLA <- dS.dMDR <- dS.dT <- matrix(NA, nMCMC, N.Temp.xs)
  
  
  # Calculate dy/dt and dS/dy for each MCMC step across the temp gradient
  for(i in 1:nMCMC){ # loop through MCMC steps
    
    # Calculate derivative of all traits with respect to temp (dy/dt) across temp gradient (for a single MCMC step)
    # The sims.list refers to the lists of fitted TPC parameters (T0, Tm, and q)
    da.dT <- d_briere(Temp.xs,
                      a_param$cf.T0[i], # T0
                      a_param$cf.Tm[i], # Tm
                      a_param$cf.q[i]) # q
    
    dbc.dT <- d_quad(Temp.xs,
                     bc_param$cf.T0[i], # T0
                     bc_param$cf.Tm[i], # Tm
                     bc_param$cf.q[i]) # q
    
    dlf.dT <- d_quad(Temp.xs,
                       lf_param$cf.T0[i], # T0
                       lf_param$cf.Tm[i], # Tm
                       lf_param$cf.q[i]) # q
    
    dPDR.dT <- d_briere(Temp.xs,
                        PDR_param$cf.T0[i], # T0
                        PDR_param$cf.Tm[i], # Tm
                        PDR_param$cf.q[i]) # q
    
    dEFGC.dT <- d_quad(Temp.xs,
                         EFGC_param$cf.T0[i], # T0
                         EFGC_param$cf.Tm[i], # Tm
                         EFGC_param$cf.q[i]) # q
    
    dEV.dT <- d_quad(Temp.xs,
                     EV_param$cf.T0[i], # T0
                     EV_param$cf.Tm[i], # Tm
                     EV_param$cf.q[i]) # q
    
    dpLA.dT <- d_quad(Temp.xs,
                      pLA_param$cf.T0[i], # T0
                      pLA_param$cf.Tm[i], # Tm
                      pLA_param$cf.q[i]) # q
    
    dMDR.dT <- d_briere(Temp.xs,
                        MDR_param$cf.T0[i], # T0
                        MDR_param$cf.Tm[i], # Tm
                        MDR_param$cf.q[i]) # q
    
    # Calculate sensitivity (dS/dy * dy/dt) across temp gradient (for a single MCMC step)
    
    # See Mathematica notebook from Shocket et al. 2018 eLife for dR0/dy derivative calculations
    
    dS.da[i, ] <- 3/2 * S(a_preds[i, ], m_bc, m_lf, m_PDR, m_EFGC, m_EV, m_pLA, m_MDR)/(a_preds[i, ]+ec) * da.dT
    dS.dbc[i, ] <- 1/2 * (S(m_a, bc_preds[i, ], m_lf, m_PDR, m_EFGC, m_EV, m_pLA, m_MDR)/(bc_preds[i, ]+ec) * dbc.dT)
    dS.dlf[i, ] <- 1/2 * (S(m_a, m_bc, lf_preds[i, ], m_PDR, m_EFGC, m_EV, m_pLA, m_MDR) * 
                            (1 + 3*lf_preds[i, ]*m_PDR) / ((lf_preds[i, ] + ec)^2 * (m_PDR + ec)) * dlf.dT)
    dS.dPDR[i, ] <- 1/2 * (S(m_a, m_bc, m_lf, PDR_preds[i, ], m_EFGC, m_EV, m_pLA, m_MDR)/((m_lf + ec)*(PDR_preds[i, ]+ec)^2) * dPDR.dT)
    dS.dEFGC[i, ] <- 1/2 * (S(m_a, m_bc, m_lf, m_PDR, EFGC_preds[i, ], m_EV, m_pLA, m_MDR)/(EFGC_preds[i, ]+ec) * dEFGC.dT)
    dS.dEV[i, ] <- 1/2 * (S(m_a, m_bc, m_lf, m_PDR, m_EFGC, EV_preds[i, ], m_pLA, m_MDR)/(EV_preds[i, ]+ec) * dEV.dT)
    dS.dpLA[i, ] <- 1/2 * (S(m_a, m_bc, m_lf, m_PDR, m_EFGC, m_EV, pLA_preds[i, ], m_MDR)/(pLA_preds[i, ]+ec) * dpLA.dT)
    dS.dMDR[i, ] <- 1/2 * (S(m_a, m_bc, m_lf, m_PDR, m_EFGC, m_EV, m_pLA, MDR_preds[i, ])/(MDR_preds[i, ]+ec) * dMDR.dT)
    
    dS.dT[i, ] <-  dS.da[i, ] + dS.dbc[i, ] + dS.dlf[i, ] + dS.dPDR[i, ] + dS.dEFGC[i, ] + dS.dEV[i, ] + dS.dpLA[i, ] + dS.dMDR[i, ]
    
  } # end MCMC loop
  
  # Collect output in a list and return it
  SA_list_out <- list(dS.da, dS.dbc, dS.dlf, dS.dPDR, dS.dEFGC, dS.dEV, dS.dpLA, dS.dMDR, dS.dT)
  SA_list_out
  
} # end function
