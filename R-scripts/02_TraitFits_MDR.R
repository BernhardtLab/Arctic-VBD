## Lilian Chan, University of Guelph
## Arctic vector-borne disease transmission suitability model
##
## Purpose: use Bayesian inference (JAGS) to fit TPCs for mosquito development 
## rate (MDR) for Aedes nigripes (Culler et al. 2015)
##     1) with uniform priors; and 
##     2) with data-informed priors from Aedes sierrensis data (Couper et al. 2024)
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
##        E. Plot all three TPCs in the same graph (for comparison)
##
##    3. Fitting TPC (Quadratic)
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
library(MASS)
library(ggsci)

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
# Number of posterior dist elements = [(ni - nb) / nt ] * nc = [ (110000 - 10000) / 100 ] * 5 = 5000
ni <- 110000 # number of iterations in each chain
nb <- 10000 # number of 'burn in' iterations to discard
nt <- 100 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

##### Temp sequence for derived quantity calculations
# For actual fits
# Temp.xs <- seq(0, 45, 0.1)
# N.Temp.xs <-length(Temp.xs)

# For priors - fewer temps for derived calculations makes it go faster
Temp.xs <- seq(0, 45, 0.5)
N.Temp.xs <-length(Temp.xs)

##########
###### 2A. Fit MDR thermal responses with uniform priors (Ae. nigripes): Briere ----
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
# MDR.nigripes.bri.uni <- jags(data = jag.data, 
#                      inits = inits, 
#                      parameters.to.save = parameters, 
#                      model.file = "R-scripts/briere_T.txt",
#                      n.thin = nt, 
#                      n.chains = nc, 
#                      n.burnin = nb, 
#                      n.iter = ni, 
#                      DIC = T, 
#                      working.directory = getwd()
# )

## Save the model as Rdata 
#save(MDR.nigripes.bri.uni, file = "R-scripts/R2jags-objects/MDR.nigripes.bri.uni.Rdata")

# Read the .Rdata
load("R-scripts/R2jags-objects/MDR.nigripes.bri.uni.Rdata")


## Diagnostics ----
##### Examine output
MDR.nigripes.bri.uni$BUGSoutput$summary[1:5,]
mcmcplot(MDR.nigripes.bri.uni)

# Extract the DIC for future model comparisons
MDR.nigripes.bri.uni$BUGSoutput$DIC

## Plot data + fit ----
df.MDR.nigripes.bri.uni <- data.frame(MDR.nigripes.bri.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.MDR.nigripes.bri.uni)

##### Plot
plot.MDR.nigripes.bri.uni <- df.MDR.nigripes.bri.uni %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "#4363d8", alpha = 0.5) +
  geom_line(aes(y = mean), color = "blue", linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Development rate (days-1)"
  ) +
  theme_bw()

plot.MDR.nigripes.bri.uni

# ggsave("figures/MDR.nigripes.bri.uni.png", plot.MDR.nigripes.bri.uni, 
#        width = 10.3, height = 5.6)

##########
###### 2B. Fit MDR thermal responses for priors (Ae. sierrensis): Briere ----
##########

##### Set data
data <- data.MDR.sierrensis

##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS -----

# This code took an hour to run!
# MDR.sierrensis.bri.uni <- jags(data = jag.data,
#                                inits = inits,
#                                parameters.to.save = parameters,
#                                model.file = "R-scripts/briere_T.txt",
#                                n.thin = nt,
#                                n.chains = nc,
#                                n.burnin = nb,
#                                n.iter = ni,
#                                DIC = T,
#                                working.directory = getwd()
# )

## Save the model as Rdata 
#save(MDR.sierrensis.bri.uni, file = "R-scripts/R2jags-objects/MDR.sierrensis.bri.uni.Rdata")

# Read the .Rdata
load("R-scripts/R2jags-objects/MDR.sierrensis.bri.uni.Rdata")


## Diagnostics ----
##### Examine output
MDR.sierrensis.bri.uni$BUGSoutput$summary[1:5,]
mcmcplot(MDR.sierrensis.bri.uni)

# Extract the DIC for future model comparisons
MDR.sierrensis.bri.uni$BUGSoutput$DIC

## Plot data + fit ----
df.MDR.sierrensis.bri.uni <- data.frame(MDR.sierrensis.bri.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.MDR.sierrensis.bri.uni)

##### Plot
plot.df.MDR.sierrensis.bri.uni <- df.MDR.sierrensis.bri.uni %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "grey", alpha = 0.5) +
  geom_line(aes(y = mean), color = "#868686FF", linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait), size = 2 
             , position = "jitter"
             ) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Development rate (days-1)"
  ) +
  theme_bw()

plot.df.MDR.sierrensis.bri.uni

# ggsave("figures/MDR.sierrensis.bri.uni.png", plot.df.MDR.sierrensis.bri.uni, 
#        width = 10.3, height = 5.6)


##########
###### 2C. Fit gamma distributions to MDR prior thermal responses: Briere ----
##########

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
MDR.nigripes.prior.cf.dists <- data.frame(q = as.vector(MDR.sierrensis.bri.uni$BUGSoutput$sims.list$cf.q),
                                          T0 = as.vector(MDR.sierrensis.bri.uni$BUGSoutput$sims.list$cf.T0),
                                          Tm = as.vector(MDR.sierrensis.bri.uni$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
MDR.nigripes.prior.gamma.fits = apply(MDR.nigripes.prior.cf.dists, 2, 
                                      function(df) fitdistr(df, "gamma")$estimate)


MDR.hypers <- MDR.nigripes.prior.gamma.fits
save(MDR.hypers, file = "R-scripts/R2jags-objects/MDRhypers.bri.Rsave")


##########
###### 2D. Fit MDR thermal responses with data-informed priors (Ae. nigripes): Briere ----
##########

load("R-scripts/R2jags-objects/MDRhypers.Rsave")
MDR.nigripes.prior.gamma.fits <- MDR.hypers

##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)

##### Set data
data <- data.MDR.nigripes
hypers <- MDR.nigripes.prior.gamma.fits * 0.1

##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, 
                 N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS -----
# MDR.nigripes.bri.inf <- jags(data = jag.data, 
#                              inits = inits, 
#                              parameters.to.save = parameters, 
#                              model.file = "R-scripts/briere_inf.txt",
#                              n.thin = nt, 
#                              n.chains = nc, 
#                              n.burnin = nb, 
#                              n.iter = ni, 
#                              DIC = T, 
#                              working.directory = getwd()
# )

## Save the model as Rdata 
#save(MDR.nigripes.bri.uni, file = "R-scripts/R2jags-objects/MDR.nigripes.bri.inf.Rdata")

# Read the .Rdata
load("R-scripts/R2jags-objects/MDR.nigripes.bri.inf.Rdata")


## Diagnostics ----
##### Examine output
MDR.nigripes.bri.inf$BUGSoutput$summary[1:5,]
mcmcplot(MDR.nigripes.bri.inf)

# Extract the DIC for future model comparisons
MDR.nigripes.bri.inf$BUGSoutput$DIC

## Plot data + fit ----
df.MDR.nigripes.bri.inf <- data.frame(MDR.nigripes.bri.inf$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.MDR.nigripes.bri.inf)

##### Plot
plot.MDR.nigripes.bri.inf <- df.MDR.nigripes.bri.inf %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "pink", alpha = 0.5) +
  geom_line(aes(y = mean), color = "red", linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Development rate (days-1)"
  ) +
  theme_bw()

plot.MDR.nigripes.bri.inf

# ggsave("figures/MDR.nigripes.bri.inf.png", plot.MDR.nigripes.bri.inf, 
#        width = 10.3, height = 5.6)


##########
###### 2E. Plot all three TPCs in the same graph (for comparison): Briere ----
##########

# Add an identifying column in each model output dataframe
df.MDR.nigripes.bri.uni <- df.MDR.nigripes.bri.uni %>% 
  mutate(type = "Ae. nigripes uniform")

df.MDR.sierrensis.bri.uni <- df.MDR.sierrensis.bri.uni %>% 
  mutate(type = "Ae. sierrensis uniform")

df.MDR.nigripes.bri.inf <- df.MDR.nigripes.bri.inf %>% 
  mutate(type = "Ae. nigripes informative")

# Combine the three dataframes
df.all <- rbind(df.MDR.nigripes.bri.uni, df.MDR.sierrensis.bri.uni, df.MDR.nigripes.bri.inf)

# Plot
plot.all <- df.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = mean, color = type), linewidth = 1) +
  geom_point(data = data.MDR.nigripes, aes(x = temp, y = trait), size = 2) +
  #geom_point(data = data.MDR.sierrensis, aes(x = temp, y = trait), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Development rate (days-1)"
  ) +
  # Customize the colours
  ## ribbon
  scale_fill_manual(values = c("Ae. nigripes uniform" = "#4363d8", 
                             "Ae. sierrensis uniform" = "grey",
                             "Ae. nigripes informative" = "pink")) +
  ## line
  scale_color_manual(values = c("Ae. nigripes uniform" = "blue", 
                               "Ae. sierrensis uniform" = "#868686FF",
                               "Ae. nigripes informative" = "red")) +
  theme_bw()

plot.all

#ggsave("figures/MDR.all.png", plot.all, 
#        width = 10.3, height = 5.6)



##########
###### 3A. Fit MDR thermal responses with uniform priors (Ae. nigripes): Quadratic ----
##########

##### Set data
data <- data.MDR.nigripes

##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# ##### Run JAGS -----
# MDR.nigripes.quad.uni <- jags(data = jag.data, 
#                               inits = inits, 
#                               parameters.to.save = parameters, 
#                               model.file = "R-scripts/quad_T.txt",
#                               n.chains = nc, 
#                               n.burnin = nb, 
#                               n.iter = ni, 
#                               DIC = T, 
#                               working.directory = getwd()
# )

## Save the model as Rdata 
#save(MDR.nigripes.quad.uni, file = "R-scripts/R2jags-objects/MDR.nigripes.quad.uni.Rdata")

# Read the .Rdata
load("R-scripts/R2jags-objects/MDR.nigripes.quad.uni.Rdata")


## Diagnostics ----
##### Examine output
MDR.nigripes.quad.uni$BUGSoutput$summary[1:5,]
mcmcplot(MDR.nigripes.quad.uni)

# Extract the DIC for future model comparisons
MDR.nigripes.quad.uni$BUGSoutput$DIC

## Plot data + fit ----
df.MDR.nigripes.quad.uni <- data.frame(MDR.nigripes.quad.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.MDR.nigripes.quad.uni)

##### Plot
plot.MDR.nigripes.quad.uni <- df.MDR.nigripes.quad.uni %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "#4363d8", alpha = 0.5) +
  geom_line(aes(y = mean), color = "blue", linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Development rate (days-1)"
  ) +
  theme_bw()

plot.MDR.nigripes.quad.uni

# ggsave("figures/MDR.nigripes.quad.uni.png", plot.MDR.nigripes.quad.uni, 
#        width = 10.3, height = 5.6)

##########
###### 3B. Fit MDR thermal responses for priors (Ae. sierrensis): Quadratic ----
##########

##### Set data
data <- data.MDR.sierrensis

##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS -----

# This code took an hour to run!
# MDR.sierrensis.quad.uni <- jags(data = jag.data, 
#                              inits = inits, 
#                              parameters.to.save = parameters, 
#                              model.file = "R-scripts/quad_T.txt",
#                              n.thin = nt, 
#                              n.chains = nc, 
#                              n.burnin = nb, 
#                              n.iter = ni, 
#                              DIC = T, 
#                              working.directory = getwd()
# )

## Save the model as Rdata 
#save(MDR.sierrensis.quad.uni, file = "R-scripts/R2jags-objects/MDR.sierrensis.quad.uni.Rdata")

# Read the .Rdata
load("R-scripts/R2jags-objects/MDR.sierrensis.quad.uni.Rdata")


## Diagnostics ----
##### Examine output
MDR.sierrensis.quad.uni$BUGSoutput$summary[1:5,]
mcmcplot(MDR.sierrensis.quad.uni)

# Extract the DIC for future model comparisons
MDR.sierrensis.quad.uni$BUGSoutput$DIC

## Plot data + fit ----
df.MDR.sierrensis.quad.uni <- data.frame(MDR.sierrensis.quad.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.MDR.sierrensis.quad.uni)

##### Plot
plot.df.MDR.sierrensis.quad.uni <- df.MDR.sierrensis.quad.uni %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "grey", alpha = 0.5) +
  geom_line(aes(y = mean), color = "#868686FF", linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait), size = 2 
             , position = "jitter"
  ) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Development rate (days-1)"
  ) +
  theme_bw()

plot.df.MDR.sierrensis.quad.uni

# ggsave("figures/MDR.sierrensis.quad.uni.png", plot.df.MDR.sierrensis.quad.uni, 
#        width = 10.3, height = 5.6)


##########
###### 3C. Fit gamma distributions to MDR prior thermal responses: Quadratic ----
##########

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
MDR.nigripes.prior.cf.dists <- data.frame(q = as.vector(MDR.sierrensis.quad.uni$BUGSoutput$sims.list$cf.q),
                                          T0 = as.vector(MDR.sierrensis.quad.uni$BUGSoutput$sims.list$cf.T0),
                                          Tm = as.vector(MDR.sierrensis.quad.uni$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
MDR.nigripes.prior.gamma.fits = apply(MDR.nigripes.prior.cf.dists, 2, 
                                      function(df) fitdistr(df, "gamma")$estimate)


MDR.hypers <- MDR.nigripes.prior.gamma.fits
save(MDR.hypers, file = "R-scripts/R2jags-objects/MDRhypers.quad.Rsave")


##########
###### 3D. Fit MDR thermal responses with data-informed priors (Ae. nigripes): Quadratic ----
##########

load("R-scripts/R2jags-objects/MDRhypers.quad.Rsave")
MDR.nigripes.prior.gamma.fits <- MDR.hypers

##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)

##### Set data
data <- data.MDR.nigripes
hypers <- MDR.nigripes.prior.gamma.fits * 0.1

##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, 
                 N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS -----
# MDR.nigripes.quad.inf <- jags(data = jag.data, 
#                              inits = inits, 
#                              parameters.to.save = parameters, 
#                              model.file = "R-scripts/quad_inf.txt",
#                              n.thin = nt, 
#                              n.chains = nc, 
#                              n.burnin = nb, 
#                              n.iter = ni, 
#                              DIC = T, 
#                              working.directory = getwd()
# )

## Save the model as Rdata 
#save(MDR.nigripes.quad.uni, file = "R-scripts/R2jags-objects/MDR.nigripes.quad.inf.Rdata")

# Read the .Rdata
load("R-scripts/R2jags-objects/MDR.nigripes.quad.inf.Rdata")


## Diagnostics ----
##### Examine output
MDR.nigripes.quad.inf$BUGSoutput$summary[1:5,]
mcmcplot(MDR.nigripes.quad.inf)

# Extract the DIC for future model comparisons
MDR.nigripes.quad.inf$BUGSoutput$DIC

## Plot data + fit ----
df.MDR.nigripes.quad.inf <- data.frame(MDR.nigripes.quad.inf$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.MDR.nigripes.quad.inf)

##### Plot
plot.MDR.nigripes.quad.inf <- df.MDR.nigripes.quad.inf %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "pink", alpha = 0.5) +
  geom_line(aes(y = mean), color = "red", linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Development rate (days-1)"
  ) +
  theme_bw()

plot.MDR.nigripes.quad.inf

# ggsave("figures/MDR.nigripes.quad.inf.png", plot.MDR.nigripes.quad.inf, 
#        width = 10.3, height = 5.6)


##########
###### 3E. Plot all three TPCs in the same graph (for comparison) ----
##########

# Add an identifying column in each model output dataframe
df.MDR.nigripes.quad.uni <- df.MDR.nigripes.quad.uni %>% 
  mutate(type = "Ae. nigripes uniform")

df.MDR.sierrensis.quad.uni <- df.MDR.sierrensis.quad.uni %>% 
  mutate(type = "Ae. sierrensis uniform")

df.MDR.nigripes.quad.inf <- df.MDR.nigripes.quad.inf %>% 
  mutate(type = "Ae. nigripes informative")

# Combine the three dataframes
df.all <- rbind(df.MDR.nigripes.quad.uni, df.MDR.sierrensis.quad.uni, df.MDR.nigripes.quad.inf)

# Plot
plot.all <- df.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = mean, color = type), linewidth = 1) +
  geom_point(data = data.MDR.nigripes, aes(x = temp, y = trait), size = 2) +
  #geom_point(data = data.MDR.sierrensis, aes(x = temp, y = trait), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Development rate (days-1)"
  ) +
  # Customize the colours
  ## ribbon
  scale_fill_manual(values = c("Ae. nigripes uniform" = "#4363d8", 
                               "Ae. sierrensis uniform" = "grey",
                               "Ae. nigripes informative" = "pink")) +
  ## line
  scale_color_manual(values = c("Ae. nigripes uniform" = "blue", 
                                "Ae. sierrensis uniform" = "#868686FF",
                                "Ae. nigripes informative" = "red")) +
  theme_bw()

plot.all

#ggsave("figures/MDR.all.png", plot.all, 
#        width = 10.3, height = 5.6)