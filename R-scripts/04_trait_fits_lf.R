## Lilian Chan, University of Guelph
## Arctic vector-borne disease transmission suitability model
##
## Purpose: use Bayesian inference (JAGS) to fit TPCs for mosquito adult lifespan (lf) 
## for Aedes vexans (Costello and Brust 1971) 
##     1) with uniform priors; and 
##     2) with data-informed priors from Aedes sierrensis data (Couper et al. 2024)
##
## Table of content:
##    0. Set-up workspace
##
##    1. MCMC settings for all models
##
##    2. Fitting TPC (Briere)
##        A. Fit lf thermal responses with uniform priors (Ae. vexans)
##        B. Fit lf thermal responses for priors (Ae. sierrensis)
##        C. Fit gamma distributions to lf prior thermal responses
##        D. Fit lf thermal responses with data-informed priors (Ae. vexans)
##        E. Plot all three TPCs in the same graph (for comparison)
##
##    3. Fitting TPC (Quadratic)
##        A. Fit lf thermal responses with uniform priors (Ae. vexans)
##        B. Fit lf thermal responses for priors (Ae. sierrensis)
##        C. Fit gamma distributions to lf prior thermal responses
##        D. Fit lf thermal responses with data-informed priors (Ae. vexans)


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


# Load data
data <- read_csv("data/data-processed/TraitData_lf.csv")
unique(data$species)

# Subset data
data.lf.vexans <- subset(data, species == "vexans")
colnames(data.lf.vexans)
  
data.lf.sierrensis <- subset(data, species == "sierrensis")


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
ni <- 25000 # number of iterations in each chain
nb <- 5000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 3 # number of chains


##########
###### 2A. Fit lf thermal responses with uniform priors (Ae. vexans): Briere ----
##########

##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)


##### Set data
data <- data.lf.vexans

##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS -----
lf.vexans.bri.uni <- jags(data = jag.data,
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
#save(lf.vexans.bri.uni, file = "R-scripts/R2jags-objects/lf.vexans.bri.uni.Rdata")

# Read the .Rdata
load("R-scripts/R2jags-objects/lf.vexans.bri.uni.Rdata")


## Diagnostics ----
##### Examine output
lf.vexans.bri.uni$BUGSoutput$summary[1:5,]
mcmcplot(lf.vexans.bri.uni) # Tmin doesn't look good

# Extract the DIC for future model comparisons
lf.vexans.bri.uni$BUGSoutput$DIC

## Plot data + fit ----
df.lf.vexans.bri.uni <- data.frame(lf.vexans.bri.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.lf.vexans.bri.uni)

##### Plot
plot.lf.vexans.bri.uni <- df.lf.vexans.bri.uni %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "#4363d8", alpha = 0.5) +
  geom_line(aes(y = mean), color = "blue", linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Adult lifespan (days)"
  ) +
  theme_bw()

plot.lf.vexans.bri.uni

# ggsave("figures/lf.vexans.bri.uni.png", plot.lf.vexans.bri.uni, width = 10.3, height = 5.6)


##########
###### 2B. Fit lf thermal responses for priors (Ae. sierrensis): Briere ----
##########

##### Temp sequence for derived quantity calculations
# For priors - fewer temps for derived calculations makes it go faster
Temp.xs <- seq(0, 45, 0.5)
N.Temp.xs <-length(Temp.xs)

##### Set data
data <- data.lf.sierrensis

##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS -----

lf.sierrensis.bri.uni <- jags(data = jag.data,
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
# save(lf.sierrensis.bri.uni, file = "R-scripts/R2jags-objects/lf.sierrensis.bri.uni.Rdata")

# Read the .Rdata
load("R-scripts/R2jags-objects/lf.sierrensis.bri.uni.Rdata")


## Diagnostics ----
##### Examine output
lf.sierrensis.bri.uni$BUGSoutput$summary[1:5,]
mcmcplot(lf.sierrensis.bri.uni)

# Extract the DIC for future model comparisons
lf.sierrensis.bri.uni$BUGSoutput$DIC

## Plot data + fit ----
df.lf.sierrensis.bri.uni <- data.frame(lf.sierrensis.bri.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.lf.sierrensis.bri.uni)

##### Plot
plot.lf.sierrensis.bri.uni <- df.lf.sierrensis.bri.uni %>% 
  ggplot(aes(x = temp)) +
  geom_point(data = data, aes(x = temp, y = trait), size = 2
             , position = "jitter"
  ) +
  stat_summary(data = data, aes(x = temp, y = trait, colour = "red"),
               fun = mean, geom = "point") +
  stat_summary(data = data, aes(x = temp, y = trait, colour = "red"),
               fun.data = "mean_se", geom = "errorbar") +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "grey", alpha = 0.5) +
  geom_line(aes(y = mean), color = "#868686FF", linewidth = 1) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Adult lifespan (days)"
  ) +
  theme_bw()

plot.lf.sierrensis.bri.uni

# ggsave("figures/lf.sierrensis.bri.uni.png", plot.lf.sierrensis.bri.uni, width = 10.3, height = 5.6)

##########
###### 2C. Fit gamma distributions to lf prior thermal responses: Briere ----
##########

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
lf.vexans.prior.cf.dists <- data.frame(q = as.vector(lf.sierrensis.bri.uni$BUGSoutput$sims.list$cf.q),
                                          T0 = as.vector(lf.sierrensis.bri.uni$BUGSoutput$sims.list$cf.T0),
                                          Tm = as.vector(lf.sierrensis.bri.uni$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
lf.vexans.prior.gamma.fits <- apply(lf.vexans.prior.cf.dists, 2, 
                                      function(df) fitdistr(df, "gamma")$estimate)


lf.hypers <- lf.vexans.prior.gamma.fits
#save(lf.hypers, file = "R-scripts/R2jags-objects/lfhypers.bri.Rsave")


##########
###### 2D. Fit lf thermal responses with data-informed priors (Ae. vexans): Briere ----
##########

load("R-scripts/R2jags-objects/lfhypers.bri.Rsave")
lf.vexans.prior.gamma.fits <- lf.hypers

##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)

##### Set data
data <- data.lf.vexans
hypers <- lf.vexans.prior.gamma.fits * 0.1

##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, 
                 N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS -----
lf.vexans.bri.inf <- jags(data = jag.data,
                             inits = inits,
                             parameters.to.save = parameters,
                             model.file = "R-scripts/briere_inf.txt",
                             n.thin = nt,
                             n.chains = nc,
                             n.burnin = nb,
                             n.iter = ni,
                             DIC = T,
                             working.directory = getwd()
)

## Save the model as Rdata 
# save(lf.vexans.bri.inf, file = "R-scripts/R2jags-objects/lf.vexans.bri.inf.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/lf.vexans.bri.inf.Rdata")


## Diagnostics ----
##### Examine output
lf.vexans.bri.inf$BUGSoutput$summary[1:5,]
mcmcplot(lf.vexans.bri.inf)

# Extract the DIC for future model comparisons
lf.vexans.bri.inf$BUGSoutput$DIC

## Plot data + fit ----
df.lf.vexans.bri.inf <- data.frame(lf.vexans.bri.inf$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.lf.vexans.bri.inf)

##### Plot
plot.lf.vexans.bri.inf <- df.lf.vexans.bri.inf %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "pink", alpha = 0.5) +
  geom_line(aes(y = mean), color = "red", linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Adult lifespan (days)"
  ) +
  theme_bw()

plot.lf.vexans.bri.inf

# ggsave("figures/lf.vexans.bri.inf.png", plot.lf.vexans.bri.inf, 
#        width = 10.3, height = 5.6)


##########
###### 2E. Plot all three TPCs in the same graph (for comparison): Briere ----
##########

# Add an identifying column in each model output dataframe
df.lf.vexans.bri.uni <- df.lf.vexans.bri.uni %>% 
  mutate(type = "Ae. vexans uniform")

df.lf.sierrensis.bri.uni <- df.lf.sierrensis.bri.uni %>% 
  mutate(type = "Ae. sierrensis uniform")

df.lf.vexans.bri.inf <- df.lf.vexans.bri.inf %>% 
  mutate(type = "Ae. vexans informative")

# Combine the three dataframes
df.all <- rbind(df.lf.vexans.bri.uni, df.lf.sierrensis.bri.uni, df.lf.vexans.bri.inf)

# Plot
plot.all <- df.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = mean, color = type), linewidth = 1) +
  geom_point(data = data.lf.vexans, aes(x = temp, y = trait), size = 2) +
  #geom_point(data = data.lf.sierrensis, aes(x = temp, y = trait), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Adult lifespan (days)"
  ) +
  # Customize the colours
  ## ribbon
  scale_fill_manual(values = c("Ae. vexans uniform" = "#4363d8", 
                               "Ae. sierrensis uniform" = "grey",
                               "Ae. vexans informative" = "pink")) +
  ## line
  scale_color_manual(values = c("Ae. vexans uniform" = "blue", 
                                "Ae. sierrensis uniform" = "#868686FF",
                                "Ae. vexans informative" = "red")) +
  theme_bw()

plot.all

#ggsave("figures/lf.all.bri.png", plot.all, 
#        width = 10.3, height = 5.6)



##########
###### 3A. Fit lf thermal responses with uniform priors (Ae. vexans): Quadratic ----
##########

##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)

##### Set data
data <- data.lf.vexans

##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# ##### Run JAGS -----
lf.vexans.quad.uni <- jags(data = jag.data,
                              inits = inits,
                              parameters.to.save = parameters,
                              model.file = "R-scripts/quad_T.txt",
                              n.chains = nc,
                              n.burnin = nb,
                              n.iter = ni,
                              DIC = T,
                              working.directory = getwd()
)

## Save the model as Rdata 
#save(lf.vexans.quad.uni, file = "R-scripts/R2jags-objects/lf.vexans.quad.uni.Rdata")

# Read the .Rdata
load("R-scripts/R2jags-objects/lf.vexans.quad.uni.Rdata")


## Diagnostics ----
##### Examine output
lf.vexans.quad.uni$BUGSoutput$summary[1:5,]
mcmcplot(lf.vexans.quad.uni)

# Extract the DIC for future model comparisons
lf.vexans.quad.uni$BUGSoutput$DIC

## Plot data + fit ----
df.lf.vexans.quad.uni <- data.frame(lf.vexans.quad.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.lf.vexans.quad.uni)

##### Plot
plot.lf.vexans.quad.uni <- df.lf.vexans.quad.uni %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "#4363d8", alpha = 0.5) +
  geom_line(aes(y = mean), color = "blue", linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Adult lifespan (days)"
  ) +
  theme_bw()

plot.lf.vexans.quad.uni

# ggsave("figures/lf.vexans.quad.uni.png", plot.lf.vexans.quad.uni, 
#        width = 10.3, height = 5.6)


##########
###### 3B. Fit lf thermal responses for priors (Ae. sierrensis): Quadratic ----
##########

##### Temp sequence for derived quantity calculations
# For priors - fewer temps for derived calculations makes it go faster
Temp.xs <- seq(0, 45, 0.5)
N.Temp.xs <-length(Temp.xs)


##### Set data
data <- data.lf.sierrensis

##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS -----

# This code took an hour to run!
# lf.sierrensis.quad.uni <- jags(data = jag.data, 
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
# save(lf.sierrensis.quad.uni, file = "R-scripts/R2jags-objects/lf.sierrensis.quad.uni.Rdata")

# Read the .Rdata
load("R-scripts/R2jags-objects/lf.sierrensis.quad.uni.Rdata")


## Diagnostics ----
##### Examine output
lf.sierrensis.quad.uni$BUGSoutput$summary[1:5,]
mcmcplot(lf.sierrensis.quad.uni)

# Extract the DIC for future model comparisons
lf.sierrensis.quad.uni$BUGSoutput$DIC

## Plot data + fit ----
df.lf.sierrensis.quad.uni <- data.frame(lf.sierrensis.quad.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.lf.sierrensis.quad.uni)

##### Plot
plot.df.lf.sierrensis.quad.uni <- df.lf.sierrensis.quad.uni %>% 
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
    y = "Adult lifespan (days)"
  ) +
  theme_bw()

plot.df.lf.sierrensis.quad.uni

# ggsave("figures/lf.sierrensis.quad.uni.png", plot.df.lf.sierrensis.quad.uni, 
#        width = 10.3, height = 5.6)


##########
###### 3C. Fit gamma distributions to lf prior thermal responses: Quadratic ----
##########

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
lf.vexans.prior.cf.dists <- data.frame(q = as.vector(lf.sierrensis.quad.uni$BUGSoutput$sims.list$cf.q),
                                       T0 = as.vector(lf.sierrensis.quad.uni$BUGSoutput$sims.list$cf.T0),
                                       Tm = as.vector(lf.sierrensis.quad.uni$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
lf.vexans.prior.gamma.fits = apply(lf.vexans.prior.cf.dists, 2, 
                                   function(df) fitdistr(df, "gamma")$estimate)


lf.hypers <- lf.vexans.prior.gamma.fits
save(lf.hypers, file = "R-scripts/R2jags-objects/lfhypers.quad.Rsave")


##########
###### 3D. Fit lf thermal responses with data-informed priors (Ae. vexans): Quadratic ----
##########

load("R-scripts/R2jags-objects/lfhypers.quad.Rsave")
lf.vexans.prior.gamma.fits <- lf.hypers

##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)

##### Set data
data <- data.lf.vexans
hypers <- lf.vexans.prior.gamma.fits * 0.1

##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, 
                 N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS -----
lf.vexans.quad.inf <- jags(data = jag.data,
                              inits = inits,
                              parameters.to.save = parameters,
                              model.file = "R-scripts/quad_inf.txt",
                              n.thin = nt,
                              n.chains = nc,
                              n.burnin = nb,
                              n.iter = ni,
                              DIC = T,
                              working.directory = getwd()
)

## Save the model as Rdata 
#save(lf.vexans.quad.inf, file = "R-scripts/R2jags-objects/lf.vexans.quad.inf.Rdata")

# Read the .Rdata
load("R-scripts/R2jags-objects/lf.vexans.quad.inf.Rdata")


## Diagnostics ----
##### Examine output
lf.vexans.quad.inf$BUGSoutput$summary[1:5,]
mcmcplot(lf.vexans.quad.inf)

# Extract the DIC for future model comparisons
lf.vexans.quad.inf$BUGSoutput$DIC

## Plot data + fit ----
df.lf.vexans.quad.inf <- data.frame(lf.vexans.quad.inf$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.lf.vexans.quad.inf)

##### Plot
plot.lf.vexans.quad.inf <- df.lf.vexans.quad.inf %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "pink", alpha = 0.5) +
  geom_line(aes(y = mean), color = "red", linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Adult lifespan (days)"
  ) +
  theme_bw()

plot.lf.vexans.quad.inf

# ggsave("figures/lf.vexans.quad.inf.png", plot.lf.vexans.quad.inf, 
#        width = 10.3, height = 5.6)


##########
###### 3E. Plot all three TPCs in the same graph (for comparison) ----
##########

# Add an identifying column in each model output dataframe
df.lf.vexans.quad.uni <- df.lf.vexans.quad.uni %>% 
  mutate(type = "Ae. vexans uniform")

df.lf.sierrensis.quad.uni <- df.lf.sierrensis.quad.uni %>% 
  mutate(type = "Ae. sierrensis uniform")

df.lf.vexans.quad.inf <- df.lf.vexans.quad.inf %>% 
  mutate(type = "Ae. vexans informative")

# Combine the three dataframes
df.all <- rbind(df.lf.vexans.quad.uni, df.lf.sierrensis.quad.uni, df.lf.vexans.quad.inf)

##### Plot
plot.all <- df.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = mean, color = type), linewidth = 1) +
  geom_point(data = data.lf.vexans, aes(x = temp, y = trait), size = 2) +
  #geom_point(data = data.lf.sierrensis, aes(x = temp, y = trait), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Adult lifespan (days)"
  ) +
  # Customize the colours
  ## ribbon
  scale_fill_manual(values = c("Ae. vexans uniform" = "#4363d8", 
                               "Ae. sierrensis uniform" = "grey",
                               "Ae. vexans informative" = "pink")) +
  ## line
  scale_color_manual(values = c("Ae. vexans uniform" = "blue", 
                                "Ae. sierrensis uniform" = "#868686FF",
                                "Ae. vexans informative" = "red")) +
  theme_bw()

plot.all

# ggsave("figures/lf.all.quad.png", plot.all,
#        width = 10.3, height = 5.6)

##### Plot all vexans TPCs for comparison ----
# Add an identifying column in each model output dataframe
df.lf.vexans.bri.uni <- df.lf.vexans.bri.uni %>% 
  mutate(type = "Briere (uni)")

df.lf.vexans.bri.inf <- df.lf.vexans.bri.inf %>% 
  mutate(type = "Briere (inf)")

df.lf.vexans.quad.uni <- df.lf.vexans.quad.uni %>% 
  mutate(type = "Quadratic (uni)")

df.lf.vexans.quad.inf <- df.lf.vexans.quad.inf %>% 
  mutate(type = "Quadratic (inf)")

# Combine the three dataframes
df.all <- rbind(df.lf.vexans.bri.uni, df.lf.vexans.bri.inf, 
                df.lf.vexans.quad.uni, df.lf.vexans.quad.inf)

##### Plot
plot.all <- df.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = mean, color = type), linewidth = 1) +
  geom_point(data = data.lf.vexans, aes(x = temp, y = trait), size = 2) +
  #geom_point(data = data.lf.sierrensis, aes(x = temp, y = trait), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Adult lifespan (days)"
  ) +
  # Customize the colours
  scale_fill_jco() +
  scale_color_jco() +
  theme_bw()

plot.all

# ggsave("figures/lf.all.vexans.png", plot.all,
#        width = 10.3, height = 5.6)