## Lilian Chan, University of Guelph
## Arctic vector-borne disease transmission suitability model
##
## Functions to load for use in other scripts
## Code adapted from https://github.com/JoeyBernhardt/anopheles-rate-summation/blob/master/R-scripts/working-versions-code/00_RSProjectFunctions.R

##	Table of contents:
##
##	1. Functions to process TPC model output
##			A. Function to calculate TPC posterior summary statistics across a temperature gradient
##			B. Function to extract full posterior distributions for 3 mean-defining TPC parameters & calculate Tbreadth
##			C. Function to calculate Topt
##			D. Function to calculate Tmin, Max, and Tbreadth for derived TPCs 
##			E. Wrapper function to calculate summary data for and extract parameter posteriors from JAGS fitted TPCs
##			F. Wrapper function to calculate summary data for derived TPCs 


##  5. Function to calculate relative suitability
##  6. Functions for Sensitivity Analysis
##			A. Function for derivative of Briere thermal response
##			B. Function for derivative of quadratic thermal response
##			C. Function for sensitivity analysis #1 - partial derivatives




##########
###### 0. Set-up workspace ----
##########

library(tidyverse)


########################################### 1. Functions to process TPC model output

###### A. Function to calculate TPC posterior summary statistics across a temperature gradient
calcPostQuants = function(input, grad.xs) {
  
  # Get length of gradient
  N.grad.xs <- length(grad.xs)
  
  # Create output dataframe
  output.df <- data.frame("mean" = numeric(N.Temp.xs), "median" = numeric(N.Temp.xs), 
                          "lowerCI" = numeric(N.Temp.xs), "upperCI" = numeric(N.Temp.xs), 
                          "lowerQuartile" = numeric(N.Temp.xs), "upperQuartile" = numeric(N.Temp.xs), temp = grad.xs)
  
  # Calculate mean & quantiles
  for(i in 1:N.grad.xs){
    output.df$mean[i] <- mean(input[ ,i])
    output.df$median[i] <- quantile(input[ ,i], 0.5, na.rm = TRUE)
    output.df$lowerCI[i] <- quantile(input[ ,i], 0.025, na.rm = TRUE)
    output.df$upperCI[i] <- quantile(input[ ,i], 0.975, na.rm = TRUE)
    output.df$lowerQuartile[i] <- quantile(input[ ,i], 0.25, na.rm = TRUE)
    output.df$upperQuartile[i] <- quantile(input[ ,i], 0.75, na.rm = TRUE)
  }
  
  output.df # return output
  
}


# Creating a small constant to keep denominators from being zero
ec <- 0.000001

# Define S(T) with bc as one value
S = function(a, bc, lf, PDR, B, EV, pLA, MDR){
  (a^2 * bc * exp(-(1/(lf+ec))*(1/(PDR+ec))) * B * EV * pLA * MDR * lf^2)^0.5
}

########################################### 6. Functions for Sensitivity Analysis

###### A. Function for derivative of Briere thermal response
d_briere = function(T, T0, Tm, q) {
  
  b <- c()
  
  for (i in 1:length(T)) {
    if (T[i]>T0 && T[i]<Tm) ## When trait value > 0
      {b[i] <- (q*(-5*(T[i]^2) + 3*T[i]*T0 + 4*T[i]*Tm - 2*T0*Tm)/(2*sqrt(Tm-T[i])))}
    else {b[i] <- 0}
  }
  
  b # return output
  
}

###### B. Function for derivative of quadratic thermal response
d_quad = function(T, T0, Tm, q){
  
  b <- c()
  
  for (i in 1:length(T)){
    if (T[i]>T0 && T[i]<Tm) # When trait value > 0
      {b[i] <- -1*q*(2*T[i] - T0 - Tm)}
    else {b[i] <- 0}
  }
  
  b # return output
  
}


###### C. Function for sensitivity analysis #1 - partial derivatives

# Arguments: mod_x_pred = the JAGS model for each trait (for the fitted TPC parameters - T0, Tm, and q);
#			 m_x = the mean value for each trait over the temperature gradient

SensitivityAnalysis_pd = function(mod_a, mod_bc, mod_lf, mod_PDR, mod_B, mod_EV, mod_pLA, mod_MDR,
                                   m_a, m_bc, m_lf, m_PDR, m_B, m_EV, m_pLA, m_MDR) {
  
  # Create matrices to hold results
  dS.da <- dS.dbc <- dS.dlf <- dS.dPDR <- dS.dB <- dS.dEV <- dS.dpLA <- dS.dMDR <- dS.dT <- matrix(NA, nMCMC, N.Temp.xs)
  
  # Extract predicted trait values
  mod_a_preds <- mod_a$BUGSoutput$sims.list$z.trait.mu.pred.pop ## Only get the population-level fit for a
  mod_bc_preds <- mod_bc$BUGSoutput$sims.list$z.trait.mu.pred
  mod_lf_preds <- mod_lf$BUGSoutput$sims.list$z.trait.mu.pred.pop ## Only get the population-level fit for lf
  mod_PDR_preds <- mod_PDR$BUGSoutput$sims.list$z.trait.mu.pred
  mod_B_preds <- mod_B$BUGSoutput$sims.list$z.trait.mu.pred
  mod_EV_preds <- mod_EV$BUGSoutput$sims.list$z.trait.mu.pred
  mod_pLA_preds <- mod_pLA$BUGSoutput$sims.list$z.trait.mu.pred
  mod_MDR_preds <- mod_MDR$BUGSoutput$sims.list$z.trait.mu.pred
  
  # Calculate dy/dt and dS/dy for each MCMC step across the temp gradient
  for(i in 1:nMCMC){ # loop through MCMC steps
    
    # Calculate derivative of all traits with respect to temp (dy/dt) across temp gradient (for a single MCMC step)
    # The sims.list refers to the lists of fitted TPC parameters (T0, Tm, and q)
    da.dT <- d_briere(Temp.xs, 
                      mod_a$BUGSoutput$sims.list[[1]][i], # T0
                      mod_a$BUGSoutput$sims.list[[2]][i], # Tm
                      mod_a$BUGSoutput$sims.list[[3]][i]) # q
    
    dbc.dT <- d_quad(Temp.xs, 
                     mod_bc$BUGSoutput$sims.list[[1]][i], # T0
                     mod_bc$BUGSoutput$sims.list[[2]][i], # Tm
                     mod_bc$BUGSoutput$sims.list[[3]][i]) # q
    
    dlf.dT <- d_briere(Temp.xs, 
                       mod_lf$BUGSoutput$sims.list[[1]][i], # T0
                       mod_lf$BUGSoutput$sims.list[[2]][i], # Tm
                       mod_lf$BUGSoutput$sims.list[[3]][i]) # q
    
    dPDR.dT <- d_briere(Temp.xs,
                        mod_PDR$BUGSoutput$sims.list[[1]][i], # T0
                        mod_PDR$BUGSoutput$sims.list[[2]][i], # Tm
                        mod_PDR$BUGSoutput$sims.list[[3]][i]) # q
    
    dB.dT <- d_briere(Temp.xs,
                      mod_B$BUGSoutput$sims.list[[1]][i], # T0
                      mod_B$BUGSoutput$sims.list[[2]][i], # Tm
                      mod_B$BUGSoutput$sims.list[[3]][i]) # q
    
    dEV.dT <- d_quad(Temp.xs,
                      mod_EV$BUGSoutput$sims.list[[1]][i], # T0
                      mod_EV$BUGSoutput$sims.list[[2]][i], # Tm
                      mod_EV$BUGSoutput$sims.list[[3]][i]) # q
    
    dpLA.dT <- d_quad(Temp.xs,
                      mod_pLA$BUGSoutput$sims.list[[1]][i], # T0
                      mod_pLA$BUGSoutput$sims.list[[2]][i], # Tm
                      mod_pLA$BUGSoutput$sims.list[[3]][i]) # q
    
    dMDR.dT <- d_briere(Temp.xs,
                        mod_MDR$BUGSoutput$sims.list[[1]][i], # T0
                        mod_MDR$BUGSoutput$sims.list[[2]][i], # Tm
                        mod_MDR$BUGSoutput$sims.list[[3]][i]) # q
    
    # Calculate sensitivity (dS/dy * dy/dt) across temp gradient (for a single MCMC step)

    # See Mathematica notebook from Shocket et al. 2018 eLife for dR0/dy derivative calculations
    
    dS.da[i, ] <- S(mod_a_preds[i, ], m_bc, m_lf, m_PDR, m_B, m_EV, m_pLA, m_MDR)/(mod_a_preds[i, ]+ec) * da.dT
    dS.dbc[i, ] <- 1/2 * (S(m_a, mod_bc_preds[i, ], m_lf, m_PDR, m_B, m_EV, m_pLA, m_MDR)/(mod_bc_preds[i, ]+ec) * dbc.dT)
    dS.dlf[i, ] <- 1/2 * (S(m_a, m_bc, mod_lf_preds[i, ], m_PDR, m_B, m_EV, m_pLA, m_MDR) * 
                             (1 + 2*mod_lf_preds[i, ]*m_PDR) / ((mod_lf_preds[i, ] + ec)^2 * (m_PDR + ec)) * dlf.dT)
    dS.dPDR[i, ] <- 1/2 * (S(m_a, m_bc, m_lf, mod_PDR_preds[i, ], m_B, m_EV, m_pLA, m_MDR)/((m_lf + ec)*(mod_PDR_preds[i, ]+ec)^2) * dPDR.dT)
    dS.dB[i, ] <- 1/2 * (S(m_a, m_bc, m_lf, m_PDR, mod_B_preds[i, ], m_EV, m_pLA, m_MDR)/(mod_B_preds[i, ]+ec) * dB.dT)
    dS.dEV[i, ] <- 1/2 * (S(m_a, m_bc, m_lf, m_PDR, m_B, mod_EV_preds[i, ], m_pLA, m_MDR)/(mod_EV_preds[i, ]+ec) * dEV.dT)
    dS.dpLA[i, ] <- 1/2 * (S(m_a, m_bc, m_lf, m_PDR, m_B, m_EV, mod_pLA_preds[i, ], m_MDR)/(mod_pLA_preds[i, ]+ec) * dpLA.dT)
    dS.dMDR[i, ] <- 1/2 * (S(m_a, m_bc, m_lf, m_PDR, m_B, m_EV, m_pLA, mod_MDR_preds[i, ])/(mod_MDR_preds[i, ]+ec) * dMDR.dT)
    
    dS.dT[i, ] <-  dS.da[i, ] + dS.dbc[i, ] + dS.dlf[i, ] + dS.dPDR[i, ] + dS.dB[i, ] + dS.dEV[i, ] + dS.dpLA[i, ] + dS.dMDR[i, ]
    
  } # end MCMC loop
  
  # Collect output in a list and return it
  SA_list_out <- list(dS.da, dS.dbc, dS.dlf, dS.dPDR, dS.dB, dS.dEV, dS.dpLA, dS.dMDR, dS.dT)
  SA_list_out
  
} # end function


## No lifetime egg production ----
SensitivityAnalysis_pd_noB = function(mod_a, mod_bc, mod_lf, mod_PDR, mod_EV, mod_pLA, mod_MDR,
                                  m_a, m_bc, m_lf, m_PDR, m_B, m_EV, m_pLA, m_MDR) {
  
  # Create matrices to hold results
  dS.da <- dS.dbc <- dS.dlf <- dS.dPDR <- dS.dEV <- dS.dpLA <- dS.dMDR <- dS.dT <- matrix(NA, nMCMC, N.Temp.xs)
  
  # Extract predicted trait values
  mod_a_preds <- mod_a$BUGSoutput$sims.list$z.trait.mu.pred.pop ## Only get the population-level fit for a
  mod_bc_preds <- mod_bc$BUGSoutput$sims.list$z.trait.mu.pred
  mod_lf_preds <- mod_lf$BUGSoutput$sims.list$z.trait.mu.pred.pop ## Only get the population-level fit for lf
  mod_PDR_preds <- mod_PDR$BUGSoutput$sims.list$z.trait.mu.pred
  mod_EV_preds <- mod_EV$BUGSoutput$sims.list$z.trait.mu.pred
  mod_pLA_preds <- mod_pLA$BUGSoutput$sims.list$z.trait.mu.pred
  mod_MDR_preds <- mod_MDR$BUGSoutput$sims.list$z.trait.mu.pred
  
  # Calculate dy/dt and dS/dy for each MCMC step across the temp gradient
  for(i in 1:nMCMC){ # loop through MCMC steps
    
    # Calculate derivative of all traits with respect to temp (dy/dt) across temp gradient (for a single MCMC step)
    # The sims.list refers to the lists of fitted TPC parameters (T0, Tm, and q)
    da.dT <- d_briere(Temp.xs, 
                      mod_a$BUGSoutput$sims.list[[1]][i], # T0
                      mod_a$BUGSoutput$sims.list[[2]][i], # Tm
                      mod_a$BUGSoutput$sims.list[[3]][i]) # q
    
    dbc.dT <- d_quad(Temp.xs, 
                     mod_bc$BUGSoutput$sims.list[[1]][i], # T0
                     mod_bc$BUGSoutput$sims.list[[2]][i], # Tm
                     mod_bc$BUGSoutput$sims.list[[3]][i]) # q
    
    dlf.dT <- d_briere(Temp.xs, 
                       mod_lf$BUGSoutput$sims.list[[1]][i], # T0
                       mod_lf$BUGSoutput$sims.list[[2]][i], # Tm
                       mod_lf$BUGSoutput$sims.list[[3]][i]) # q
    
    dPDR.dT <- d_briere(Temp.xs,
                        mod_PDR$BUGSoutput$sims.list[[1]][i], # T0
                        mod_PDR$BUGSoutput$sims.list[[2]][i], # Tm
                        mod_PDR$BUGSoutput$sims.list[[3]][i]) # q
    
    
    dEV.dT <- d_quad(Temp.xs,
                     mod_EV$BUGSoutput$sims.list[[1]][i], # T0
                     mod_EV$BUGSoutput$sims.list[[2]][i], # Tm
                     mod_EV$BUGSoutput$sims.list[[3]][i]) # q
    
    dpLA.dT <- d_quad(Temp.xs,
                      mod_pLA$BUGSoutput$sims.list[[1]][i], # T0
                      mod_pLA$BUGSoutput$sims.list[[2]][i], # Tm
                      mod_pLA$BUGSoutput$sims.list[[3]][i]) # q
    
    dMDR.dT <- d_briere(Temp.xs,
                        mod_MDR$BUGSoutput$sims.list[[1]][i], # T0
                        mod_MDR$BUGSoutput$sims.list[[2]][i], # Tm
                        mod_MDR$BUGSoutput$sims.list[[3]][i]) # q
    
    # Calculate sensitivity (dS/dy * dy/dt) across temp gradient (for a single MCMC step)
    
    # See Mathematica notebook from Shocket et al. 2018 eLife for dR0/dy derivative calculations
    
    dS.da[i, ] <- S(mod_a_preds[i, ], m_bc, m_lf, m_PDR, m_B, m_EV, m_pLA, m_MDR)/(mod_a_preds[i, ]+ec) * da.dT
    dS.dbc[i, ] <- 1/2 * (S(m_a, mod_bc_preds[i, ], m_lf, m_PDR, m_B, m_EV, m_pLA, m_MDR)/(mod_bc_preds[i, ]+ec) * dbc.dT)
    dS.dlf[i, ] <- 1/2 * (S(m_a, m_bc, mod_lf_preds[i, ], m_PDR, m_B, m_EV, m_pLA, m_MDR) * 
                            (1 + 2*mod_lf_preds[i, ]*m_PDR) / ((mod_lf_preds[i, ] + ec)^2 * (m_PDR + ec)) * dlf.dT)
    dS.dPDR[i, ] <- 1/2 * (S(m_a, m_bc, m_lf, mod_PDR_preds[i, ], m_B, m_EV, m_pLA, m_MDR)/((m_lf + ec)*(mod_PDR_preds[i, ]+ec)^2) * dPDR.dT)
    dS.dEV[i, ] <- 1/2 * (S(m_a, m_bc, m_lf, m_PDR, m_B, mod_EV_preds[i, ], m_pLA, m_MDR)/(mod_EV_preds[i, ]+ec) * dEV.dT)
    dS.dpLA[i, ] <- 1/2 * (S(m_a, m_bc, m_lf, m_PDR, m_B, m_EV, mod_pLA_preds[i, ], m_MDR)/(mod_pLA_preds[i, ]+ec) * dpLA.dT)
    dS.dMDR[i, ] <- 1/2 * (S(m_a, m_bc, m_lf, m_PDR, m_B, m_EV, m_pLA, mod_MDR_preds[i, ])/(mod_MDR_preds[i, ]+ec) * dMDR.dT)
    
    dS.dT[i, ] <-  dS.da[i, ] + dS.dbc[i, ] + dS.dlf[i, ] + dS.dPDR[i, ] + dS.dEV[i, ] + dS.dpLA[i, ] + dS.dMDR[i, ]
    
  } # end MCMC loop
  
  # Collect output in a list and return it
  SA_list_out <- list(dS.da, dS.dbc, dS.dlf, dS.dPDR, dS.dEV, dS.dpLA, dS.dMDR, dS.dT)
  SA_list_out
  
} # end function

