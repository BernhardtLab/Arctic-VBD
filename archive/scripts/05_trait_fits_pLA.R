## Lilian Chan, University of Guelph
## Arctic vector-borne disease transmission suitability model
##
## Purpose: Fit thermal performance curves (TPCs) for larval-to-adult
## survival (pLA) using Bayesian inference (JAGS). Arctic species models are fit 
## using data-informed priors derived from non-Arctic species.
##
## 
## Table of content:
##    0. Set-up workspace
##
##    1. MCMC settings for all models
##
##    2. Fitting TPC (Briere)
##        A. Fit non-Arctic TPC to generate priors
##        B. Fit gamma distributions to non-Arctic TPC parameters
##        C. Fit Arctic TPC using data-informed priors
##
##    3. Fitting TPC (Quadratic)
##        A. Fit non-Arctic TPC to generate priors
##        B. Fit gamma distributions to non-Arctic TPC parameters
##        C. Fit Arctic TPC using data-informed priors
##
##    4. Compare model fit between Briere and Quadratic models
##    5. Process and save model output for visualization
##
##
## Inputs:
## data-processed/TraitData_pLA.csv - 
##     Synthesized published trait data for pLA
##
## Outputs: 
## R-scripts/R2jags-objects/best-fitting-mods/pLA.arctic.mod.Rdata - 
##     Best-fitting TPC models for Arctic species
##
## R-scripts/R2jags-objects/best-fitting-mods/pLA.nonarctic.mod.Rdata -
##     Best-fitting TPC models for non-Arctic species
##
## data-processed/pLA/pLA.arctic.predictions.summary.csv -
##     Posterior summary of TPC predictions for Arctic species across temperatures
##
## data-processed/pLA/pLA.arctic.params.summary.csv -
##     Summary statistics of TPC parameters (Arctic TPC)
##
## data-processed/pLA/pLA.arctic.params.fullposts.csv -
##     Full posterior distributions for TPC parameters (Arctic TPC)
##
## data-processed/pLA/pLA.nonarctic.predictions.summary.csv -
##     Posterior summary of TPC predictions for non-Arctic species
##
## data-processed/pLA/pLA.nonarctic.params.summary.csv -
##     Summary statistics of TPC parameters (non-Arctic TPC)
##
## data-processed/pLA/pLA.nonarctic.params.fullposts.csv -
##     Full posterior distributions for TPC parameters (non-Arctic TPC)



# 0. Set-up workspace -----------------------------------------------------

library(tidyverse)
library(readxl)
library(janitor)
library(R2jags)
library(mcmcplots) # Diagnostic plots for fits
library(MASS)
library(cowplot)


# Load functions
source("R-scripts/00_Functions.R")


# Load data
data.all <- read_csv("data-processed/TraitData_pLA.csv")
unique(data.all$species)


# Subset data
## Arctic species
data.pLA.arctic <- subset(data.all, type == "Arctic")

## Non-Arctic species
data.pLA.nonarctic <- subset(data.all, type == "non-Arctic")


# Plot the raw data
plot.data.pLA <- data.all %>% 
  ggplot(aes(x = temp, y = trait)) +
  geom_point(aes(colour = species)) +
  labs(y = "Larval survival (%)", x = expression(paste("Temperature (", degree, "C)"))) +
  scale_colour_discrete(name = "Species", labels = c("Ae. flavescens",
                                                     "Ae. nigromaculis",
                                                     "Ae. sollicitans",
                                                     "Ae. triseriatus",
                                                     "Ae. vexans"
  )) +
  facet_grid(rows = vars(type)) +
  theme_bw()

plot.data.pLA



# 1. MCMC Settings for all models ----------------------------------------------

# Number of posterior dist elements = [(ni - nb) / nt] * nc = [(450000 - 50000) / 100] * 3 = 12000
ni <- 450000 # number of iterations in each chain
nb <- 50000 # number of 'burn in' iterations to discard
nt <- 100 # thinning rate - jags saves every nt iterations in each chain
nc <- 3 # number of chains

set.seed(123) # for reproducibility


# 2. Fitting TPC (Briere) ------------------------------------------------------

## 2A. Fit non-Arctic TPC to generate priors -----------------------------------

##### Temp sequence for derived quantity calculations
# For priors - fewer temps for derived calculations makes it go faster
Temp.xs <- seq(0, 45, 0.5)
N.Temp.xs <-length(Temp.xs)


##### Set data
data <- data.pLA.nonarctic
# Although this dataset contains data from multiple species or multiple studies
# of the same species, there are less than 5 species-study combination. 
# Thus we will not incorporate random effects 

## Set priors
prior <- data.frame(q = c(0, 1),
                    T0 = c(0, 20),
                    Tm = c(25, 45))

##### inits Function
inits <- function(){list(
  cf.q = 0.01,
  cf.Tm = 35,
  cf.T0 = 5,
  cf.sigma = rlnorm(1))}


##### Parameters to Estimate
parameters <- c("cf.q", "cf.T0", "cf.Tm","cf.sigma", "z.trait.mu.pred")


##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp


##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, 
                 N.Temp.xs = N.Temp.xs, prior = prior)

##### Run JAGS
pLA.nonarctic.bri.uni <- jags(data = jag.data,
                              inits = inits,
                              parameters.to.save = parameters,
                              model.file = "R-scripts/briereprob.txt",
                              n.thin = nt,
                              n.chains = nc,
                              n.burnin = nb,
                              n.iter = ni,
                              DIC = T,
                              working.directory = getwd()
)


## Save the model as Rdata 
save(pLA.nonarctic.bri.uni, file = "R-scripts/R2jags-objects/all-mods/pLA.nonarctic.bri.uni.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/all-mods/pLA.nonarctic.bri.uni.Rdata")


## Diagnostics
##### Examine output
pLA.nonarctic.bri.uni$BUGSoutput$summary[1:5,]
# mcmcplot(pLA.nonarctic.bri.uni)

# Extract the DIC for future model comparisons
pLA.nonarctic.bri.uni$BUGSoutput$DIC


## Plot data + fit
df.pLA.nonarctic.bri.uni <- data.frame(pLA.nonarctic.bri.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)

head(df.pLA.nonarctic.bri.uni)


##### Plot
plot.pLA.nonarctic.bri.uni <- df.pLA.nonarctic.bri.uni %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "grey", alpha = 0.5) +
  geom_line(aes(y = X50.), color = "#868686FF", linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait), size = 2) +
  
  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Survival probability"
       ) +
  theme_bw()


plot.pLA.nonarctic.bri.uni

ggsave("figures/pLA.nonarctic.bri.uni.png", plot.pLA.nonarctic.bri.uni,
       width = 10.3, height = 5.6)



## 2B. Fit gamma distributions to non-Arctic TPC parameters --------------------

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
pLA.bri.prior.cf.dists <- data.frame(q = as.vector(pLA.nonarctic.bri.uni$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(pLA.nonarctic.bri.uni$BUGSoutput$sims.list$cf.T0),
                                     Tm = as.vector(pLA.nonarctic.bri.uni$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
pLA.bri.prior.gamma.fits = apply(pLA.bri.prior.cf.dists, 2, 
                                 function(df) fitdistr(df, "gamma")$estimate)


save(pLA.bri.prior.gamma.fits, file = "R-scripts/R2jags-objects/priors/pLA.bri.priors.Rsave")



## 2C. Fit Arctic TPC using data-informed priors -------------------------------

# load("R-scripts/R2jags-objects/priors/pLA.bri.priors.Rsave")


##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)


##### Set data
data <- data.pLA.arctic
# Although this dataset contains data from multiple species or multiple studies
# of the same species, there are less than 5 species-study combination. 
# Thus we will not incorporate random effects 

hypers <- pLA.bri.prior.gamma.fits * 0.5
hypers[,2] <- pLA.bri.prior.gamma.fits[,2]

q <- data.frame(hist = rgamma(100000, shape = hypers[1,1], rate = hypers[2,1]))
T0 <-data.frame(hist = rgamma(100000, shape = hypers[1,2], rate = hypers[2,2]))
Tm <- data.frame(hist = rgamma(100000, shape = hypers[1,3], rate = hypers[2,3]))

q <- ggplot() +
  geom_histogram(data = q, aes(x = hist, y = ..density..), 
                 bins = 50, fill = "grey", color = "black", alpha = 0.6) +
  stat_function(fun = dgamma, args = list(hypers[1,1], rate = hypers[2,1]), 
                color = "black", linewidth = 1.2) +
  labs(x = "q",
       y = "Density") +
  theme_bw()

T0 <- ggplot() +
  geom_histogram(data = T0, aes(x = hist, y = ..density..), 
                 bins = 50, fill = "grey", color = "black", alpha = 0.6) +
  stat_function(fun = dgamma, args = list(hypers[1,2], rate = hypers[2,2]), 
                color = "black", linewidth = 1.2) +
  labs(x = "Tmin",
       y = "Density") +
  theme_bw()

Tm <- ggplot() +
  geom_histogram(data = Tm, aes(x = hist, y = ..density..), 
                 bins = 50, fill = "grey", color = "black", alpha = 0.6) +
  stat_function(fun = dgamma, args = list(hypers[1,3], rate = hypers[2,3]), 
                color = "black", linewidth = 1.2) +
  labs(x = "Tmax",
       y = "Density") +
  theme_bw()

priors <- plot_grid(q, T0, Tm, ncol = 3)
priors
ggsave("figures/pLA.bri.priors.png", priors, width = 10.3, height = 5.6)

##### inits Function
inits<-function(){list(
  cf.q = 0.01,
  cf.Tm = 35,
  cf.T0 = 5,
  cf.sigma = rlnorm(1))}

##### Parameters to Estimate
parameters <- c("cf.q", "cf.T0", "cf.Tm","cf.sigma", "z.trait.mu.pred")


##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, 
                 N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
pLA.arctic.bri.inf <- jags(data = jag.data,
                           inits = inits,
                           parameters.to.save = parameters,
                           model.file = "R-scripts/briereprob_inf.txt",
                           n.thin = nt,
                           n.chains = nc,
                           n.burnin = nb,
                           n.iter = ni,
                           DIC = T,
                           working.directory = getwd()
)

## Save the model as Rdata 
save(pLA.arctic.bri.inf, file = "R-scripts/R2jags-objects/all-mods/pLA.arctic.bri.inf.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/all-mods/pLA.arctic.bri.inf.Rdata")


## Diagnostics
##### Examine output
pLA.arctic.bri.inf$BUGSoutput$summary[1:5,]
# mcmcplot(pLA.arctic.bri.inf)

# Extract the DIC for future model comparisons
pLA.arctic.bri.inf$BUGSoutput$DIC

## Plot data + fit
df.pLA.arctic.bri.inf <- data.frame(pLA.arctic.bri.inf$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)

head(df.pLA.arctic.bri.inf)

##### Plot
plot.pLA.arctic.bri.inf <- df.pLA.arctic.bri.inf %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "#4363d8", alpha = 0.5) +
  geom_line(aes(y = X50.), color = "blue", linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait), size = 2) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Survival probability"
  ) +
  theme_bw()

plot.pLA.arctic.bri.inf

ggsave("figures/pLA.arctic.bri.inf.png", plot.pLA.arctic.bri.inf,
       width = 10.3, height = 5.6)




# 3. Fitting TPC (quadratic) ---------------------------------------------------

## 3A. Fit non-Arctic TPC to generate priors -----------------------------------


##### Temp sequence for derived quantity calculations
# For priors - fewer temps for derived calculations makes it go faster
Temp.xs <- seq(0, 45, 0.5)
N.Temp.xs <-length(Temp.xs)


##### Set data
data <- data.pLA.nonarctic


## Set priors
prior <- data.frame(q = c(0, 1),
                    T0 = c(0, 20),
                    Tm = c(25, 45))

##### inits Function
inits <- function(){list(
  cf.q = 0.01,
  cf.Tm = 35,
  cf.T0 = 5,
  cf.sigma = rlnorm(1))}


##### Parameters to Estimate
parameters <- c("cf.q", "cf.T0", "cf.Tm","cf.sigma", "z.trait.mu.pred")


##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp


##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, 
                 N.Temp.xs = N.Temp.xs, prior = prior)

##### Run JAGS
pLA.nonarctic.quad.uni <- jags(data = jag.data,
                              inits = inits,
                              parameters.to.save = parameters,
                              model.file = "R-scripts/quadprob.txt",
                              n.thin = nt,
                              n.chains = nc,
                              n.burnin = nb,
                              n.iter = ni,
                              DIC = T,
                              working.directory = getwd()
)


## Save the model as Rdata 
save(pLA.nonarctic.quad.uni, file = "R-scripts/R2jags-objects/all-mods/pLA.nonarctic.quad.uni.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/all-mods/pLA.nonarctic.quad.uni.Rdata")


## Diagnostics
##### Examine output
pLA.nonarctic.quad.uni$BUGSoutput$summary[1:5,]
# mcmcplot(pLA.nonarctic.quad.uni)

# Extract the DIC for future model comparisons
pLA.nonarctic.quad.uni$BUGSoutput$DIC


## Plot data + fit
df.pLA.nonarctic.quad.uni <- data.frame(pLA.nonarctic.quad.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)

head(df.pLA.nonarctic.quad.uni)


##### Plot
plot.pLA.nonarctic.quad.uni <- df.pLA.nonarctic.quad.uni %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "grey", alpha = 0.5) +
  geom_line(aes(y = X50.), color = "#868686FF", linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait), size = 2) +
  
  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Survival probability"
  ) +
  theme_bw()


plot.pLA.nonarctic.quad.uni

ggsave("figures/pLA.nonarctic.quad.uni.png", plot.pLA.nonarctic.quad.uni,
       width = 10.3, height = 5.6)



## 3B. Fit gamma distributions to non-Arctic TPC parameters --------------------

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
pLA.quad.prior.cf.dists <- data.frame(q = as.vector(pLA.nonarctic.quad.uni$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(pLA.nonarctic.quad.uni$BUGSoutput$sims.list$cf.T0),
                                     Tm = as.vector(pLA.nonarctic.quad.uni$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
pLA.quad.prior.gamma.fits = apply(pLA.quad.prior.cf.dists, 2, 
                                 function(df) fitdistr(df, "gamma")$estimate)


save(pLA.quad.prior.gamma.fits, file = "R-scripts/R2jags-objects/priors/pLA.quad.priors.Rsave")


## 3C. Fit Arctic TPC using data-informed priors -------------------------------

# load("R-scripts/R2jags-objects/priors/pLA.quad.priors.Rsave")


##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)


##### Set data
data <- data.pLA.arctic
hypers <- pLA.quad.prior.gamma.fits * 0.5

q <- data.frame(hist = rgamma(100000, shape = hypers[1,1], rate = hypers[2,1]))
T0 <-data.frame(hist = rgamma(100000, shape = hypers[1,2], rate = hypers[2,2]))
Tm <- data.frame(hist = rgamma(100000, shape = hypers[1,3], rate = hypers[2,3]))

q <- ggplot() +
  geom_histogram(data = q, aes(x = hist, y = ..density..), 
                 bins = 50, fill = "grey", color = "black", alpha = 0.6) +
  stat_function(fun = dgamma, args = list(hypers[1,1], rate = hypers[2,1]), 
                color = "black", linewidth = 1.2) +
  labs(x = "q",
       y = "Density") +
  theme_bw()

T0 <- ggplot() +
  geom_histogram(data = T0, aes(x = hist, y = ..density..), 
                 bins = 50, fill = "grey", color = "black", alpha = 0.6) +
  stat_function(fun = dgamma, args = list(hypers[1,2], rate = hypers[2,2]), 
                color = "black", linewidth = 1.2) +
  labs(x = "Tmin",
       y = "Density") +
  theme_bw()

Tm <- ggplot() +
  geom_histogram(data = Tm, aes(x = hist, y = ..density..), 
                 bins = 50, fill = "grey", color = "black", alpha = 0.6) +
  stat_function(fun = dgamma, args = list(hypers[1,3], rate = hypers[2,3]), 
                color = "black", linewidth = 1.2) +
  labs(x = "Tmax",
       y = "Density") +
  theme_bw()

priors <- plot_grid(q, T0, Tm, ncol = 3)
priors
ggsave("figures/pLA.quad.priors.png", priors, width = 10.3, height = 5.6)

##### inits Function
inits<-function(){list(
  cf.q = 0.01,
  cf.Tm = 35,
  cf.T0 = 5,
  cf.sigma = rlnorm(1))}

##### Parameters to Estimate
parameters <- c("cf.q", "cf.T0", "cf.Tm","cf.sigma", "z.trait.mu.pred")


##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, 
                 N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS 
pLA.arctic.quad.inf <- jags(data = jag.data,
                           inits = inits,
                           parameters.to.save = parameters,
                           model.file = "R-scripts/quadprob_inf.txt",
                           n.thin = nt,
                           n.chains = nc,
                           n.burnin = nb,
                           n.iter = ni,
                           DIC = T,
                           working.directory = getwd()
)

## Save the model as Rdata 
save(pLA.arctic.quad.inf, file = "R-scripts/R2jags-objects/all-mods/pLA.arctic.quad.inf.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/all-mods/pLA.arctic.quad.inf.Rdata")


## Diagnostics 
##### Examine output
pLA.arctic.quad.inf$BUGSoutput$summary[1:5,]
# mcmcplot(pLA.arctic.quad.inf)

# Extract the DIC for future model comparisons
pLA.arctic.quad.inf$BUGSoutput$DIC

## Plot data + fit
df.pLA.arctic.quad.inf <- data.frame(pLA.arctic.quad.inf$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)

head(df.pLA.arctic.quad.inf)

##### Plot
plot.pLA.arctic.quad.inf <- df.pLA.arctic.quad.inf %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "#4363d8", alpha = 0.5) +
  geom_line(aes(y = X50.), color = "blue", linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait), size = 2) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Survival probability"
  ) +
  theme_bw()

plot.pLA.arctic.quad.inf

ggsave("figures/pLA.arctic.quad.inf.png", plot.pLA.arctic.quad.inf,
       width = 10.3, height = 5.6)




# 4. Compare model fit between Briere and Quadratic models ---------------------

##### Find best fitting model #####
# Add an identifying column in each model output dataframe
df.pLA.arctic.bri.inf <- df.pLA.arctic.bri.inf %>% 
  mutate(type = "briere")

df.pLA.arctic.quad.inf <- df.pLA.arctic.quad.inf %>% 
  mutate(type = "quadratic")

# Combine the three dataframes
df.all <- bind_rows(df.pLA.arctic.bri.inf, df.pLA.arctic.quad.inf)


##### Plot
plot.all <- df.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = X50., color = type), linewidth = 1) +
  geom_point(data = data.pLA.arctic, aes(x = temp, y = trait), size = 2) +
  
  # Customize the axes and labels
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Survival probability"
  ) +
  
  # Customize the colours
  ## ribbon
  scale_fill_manual(values = c("briere" = "#4363d8", 
                               "quadratic" = "pink")) +
  ## line
  scale_color_manual(values = c("briere" = "blue", 
                                "quadratic" = "red")) +
  theme_bw()

plot.all

ggsave("figures/pLA.bri.quad.png", plot.all, width = 10.3, height = 5.6)


## DIC
pLA.arctic.bri.inf$BUGSoutput$DIC 
pLA.arctic.quad.inf$BUGSoutput$DIC # This is the best fitting TPC


##### Plot Arctic vs. non-Arctic TPCs for the best fitting TPC #####
df.pLA.nonarctic.quad.uni <- df.pLA.nonarctic.quad.uni %>% 
  mutate(type = "non-Arctic")

df.pLA.arctic.quad.inf <- df.pLA.arctic.quad.inf %>% 
  mutate(type = "Arctic")

df.arctic.nonarctic <- bind_rows(df.pLA.nonarctic.quad.uni, df.pLA.arctic.quad.inf)

plot.arctic.nonarctic <- df.arctic.nonarctic %>% 
  ggplot(aes(x = temp)) +
  geom_point(data = data.all, aes(x = temp, y = trait, colour = type), size = 2) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = X50., color = type), linewidth = 1) +
  
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Survival probability"
  ) +
  
  # Customize the colours
  ## ribbon
  scale_fill_manual(values = c("Arctic" = "#4363d8", 
                               "non-Arctic" = "grey")) +
  ## line
  scale_color_manual(values = c("Arctic" = "blue", 
                                "non-Arctic" = "azure4")) +
  theme_bw()

plot.arctic.nonarctic

ggsave("figures/pLA.arctic.nonarctic.png", plot.arctic.nonarctic, width = 10.3, height = 5.6)


# Save best-fitting TPC in a separate folder
pLA.arctic.mod <- pLA.arctic.quad.inf
pLA.nonarctic.mod <- pLA.nonarctic.quad.uni


## Save the model as Rdata 
save(pLA.arctic.mod, file = "R-scripts/R2jags-objects/best-fitting-mods/pLA.arctic.mod.Rdata")
save(pLA.nonarctic.mod, file = "R-scripts/R2jags-objects/best-fitting-mods/pLA.nonarctic.mod.Rdata")


# 5. Process and save model output for visualization ---------------------------

## Analyze TPC model
# We will create 3 files: 
# a. predictions.summary: showing the mean, median, and 95% credible interval of
#      the predicted trait value at each temp from 0 to 45ºC at a 0.1ºC
# b. params.summary: showing the showing the mean, median, and 95% credible 
#      interval of TPC parameters, Topt, and Tbreadth
# c. params.fullposts: showing the TPC parameter of each MCMC iteration


##### Arctic #####
Temp.xs <- seq(0, 45, 0.1)
pLA.TPC.analysis <- extractTPC(pLA.arctic.quad.inf, "pLA", Temp.xs)
pLA.arctic.predictions.summary <- pLA.TPC.analysis[[1]]
pLA.arctic.params.summary <- pLA.TPC.analysis[[2]]
pLA.arctic.params.fullposts <- pLA.TPC.analysis[[3]]

write_csv(pLA.arctic.predictions.summary, "data-processed/pLA/pLA.arctic.predictions.summary.csv")
write_csv(pLA.arctic.params.summary, "data-processed/pLA/pLA.arctic.params.summary.csv")
write_csv(pLA.arctic.params.fullposts, "data-processed/pLA/pLA.arctic.params.fullposts.csv")

##### non-Arctic #####
Temp.xs <- seq(0, 45, 0.5)
pLA.TPC.analysis <- extractTPC(pLA.nonarctic.quad.uni, "pLA", Temp.xs)
pLA.nonarctic.predictions.summary <- pLA.TPC.analysis[[1]]
pLA.nonarctic.params.summary <- pLA.TPC.analysis[[2]]
pLA.nonarctic.params.fullposts <- pLA.TPC.analysis[[3]]

write_csv(pLA.nonarctic.predictions.summary, "data-processed/pLA/pLA.nonarctic.predictions.summary.csv")
write_csv(pLA.nonarctic.params.summary, "data-processed/pLA/pLA.nonarctic.params.summary.csv")
write_csv(pLA.nonarctic.params.fullposts, "data-processed/pLA/pLA.nonarctic.params.fullposts.csv")
