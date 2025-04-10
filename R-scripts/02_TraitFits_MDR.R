## Lilian Chan, University of Guelph
## Arctic vector-borne disease transmission suitability model
##
## Purpose: use Bayesian inference (JAGS) to fit TPCs for mosquito development 
## rate (MDR) for Aedes nigripes
##     1) with uniform priors; and 
##     2) with data-informed priors from Aedes sierrensis data
## 
## Table of content:
##    0. Set-up workspace
##
##    1. MCMC settings for all models
##
##    2. Fitting TPC (Briere)
##        A. Fit MDR thermal responses with uniform priors (Ae. nigripes)
##        B. Fit MDR thermal responses for priors (Ae. sierrensis)
##        C. Fit gamma distributions to MDR prior thermal responses
##        D. Fit MDR thermal responses with data-informed priors (Ae. nigripes)
##
##    3. Fitting TPC (Quadratic) (do this later - LC 2025.04.09)
##        A. Fit MDR thermal responses with uniform priors (Ae. nigripes)
##        B. Fit MDR thermal responses for priors (Ae. sierrensis)
##        C. Fit gamma distributions to MDR prior thermal responses
##        D. Fit MDR thermal responses with data-informed priors (Ae. nigripes)
##
##    4. Plotting



##########
###### 0. Set-up workspace ----
##########

library(tidyverse)
library(readxl)
library(janitor)
library(R2jags)
library(mcmcplots) # Diagnostic plots for fits
library(rTPC)
library(nls.multstart)
library(broom)

setwd("~/Documents/UofG/Arctic-VBD")

# Load data
data <- read_csv("data/data-processed/TraitData_MDR.csv")
unique(data$species)

# Subset data
data.MDR.nigripes <- subset(data, species == "nigripes")
data.MDR.sierrensis <- subset(data, species == "sierrensis")

# Plot the data
data %>% ggplot() +
  geom_point(aes(x = temp, y = trait, color = species), position = "jitter") +
  theme_bw()


##########
###### 1. model settings for all models ----
##########

##### inits Function
inits<-function(){list(
  cf.q = 0.01,
  cf.Tm = 35,
  cf.T0 = 5,
  cf.sigma = rlnorm(1))}

##### Parameters to Estimate
parameters <- c("cf.q", "cf.T0", "cf.Tm","cf.sigma", "z.trait.mu.pred")

##### MCMC Settings
# Number of posterior dist elements = [(ni - nb) / nt ] * nc = [ (25000 - 5000) / 8 ] * 3 = 7500
ni <- 110000 # number of iterations in each chain
nb <- 10000 # number of 'burn in' iterations to discard
nt <- 100 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

##### Temp sequence for derived quantity calculations
# For actual fits
# Temp.xs <- seq(1, 45, 0.1)
# N.Temp.xs <-length(Temp.xs)

# For priors - fewer temps for derived calculations makes it go faster
Temp.xs <- seq(5, 45, 0.5)
N.Temp.xs <-length(Temp.xs)

##########
###### 2A. Fit MDR thermal responses with uniform priors (Ae. nigripes) ----
##########

##### Set data
data <- data.MDR.nigripes

##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS -----
MDR.nigripes.bri.uni <- jags(data = jag.data, 
                     inits = inits, 
                     parameters.to.save = parameters, 
                     model.file = "R-scripts/briere_T.txt",
                     n.thin = nt, 
                     n.chains = nc, 
                     n.burnin = nb, 
                     n.iter = ni, 
                     DIC = T, 
                     working.directory = getwd()
)

## Save the model as Rdata 
#save(MDR.nigripes.bri.uni, file = "R-scripts/R2jags-objects/MDR.nigripes.bri.uni.Rdata")

# Read the .Rdata
#load("R-scripts/R2jags-objects/MDR.nigripes.bri.uni.Rdata")


## Diagnostics ----
##### Examine output
MDR.nigripes.bri.uni$BUGSoutput$summary[1:5,]
mcmcplot(MDR.nigripes.bri.uni)

