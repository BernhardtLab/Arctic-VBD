## Lilian Chan, University of Guelph
## Arctic vector-borne disease transmission suitability model
##
## Purpose: use Bayesian inference (JAGS) to fit TPCs for mosquito development 
## rate (MDR) for Aedes arctic (Culler et al. 2015)
##     1) with uniform priors; and 
##     2) with data-informed priors from Aedes nonarctic data (Couper et al. 2024)
## 
## Table of content:
##    0. Set-up workspace
##
##    1. MCMC settings for all models
##
##    2. Fitting TPC (Briere)
##        A. Fit MDR thermal responses with uniform priors (Arctic species)
##        B. Fit MDR thermal responses for priors (non-Arctic)
##        C. Fit gamma distributions to MDR prior thermal responses
##        D. Fit MDR thermal responses with data-informed priors (Arctic)
##        E. Plot all three TPCs in the same graph (for comparison)
##
##    3. Fitting TPC (Quadratic)
##        A. Fit MDR thermal responses with uniform priors (Arctic)
##        B. Fit MDR thermal responses for priors (non-Arctic)
##        C. Fit gamma distributions to MDR prior thermal responses
##        D. Fit MDR thermal responses with data-informed priors (Arctic)
##
##    4. Process and save model output for plotting



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

# Load functions
source("R-scripts/00_Functions.R")

# Load data
data <- read_csv("data-processed/TraitData_MDR.csv")
unique(data$species)

## Convert development time (1/MDR) to development rate (MDR)
data <- data %>% 
  mutate(trait = ifelse(trait_name == "1/MDR", 1/trait, trait)) %>% 
  mutate(trait_name = "MDR") 

# Subset data
data.MDR.arctic <- subset(data, species == "nigripes")
data.MDR.nonarctic <- subset(data, species == "sierrensis")


## Plot raw data
plot.data.MDR <- data %>% 
  ggplot(aes(x = temp, y = trait, colour = species)) +
  
  ## Since the Ae. sierrensis has many data, I will just plot the mean±SE
  #geom_point(data = ~filter(.x, type == "Arctic")) +
  # stat_summary(fun = mean, geom = "point") +
  # stat_summary(fun.data = "mean_se", geom = "errorbar") +
  geom_point(position = "jitter") +
  labs(y = "Mosquito development rate (days)", x = "Temperature ºC") +
  scale_colour_discrete(name = "Species", labels = c("Ae. nigripes",
                                                     "Ae. sierrensis"
  )) +
  facet_grid(rows = vars(type), scales = "free_y") +
  theme_bw()

plot.data.MDR

# ggsave("figures/raw_data/plot.data.MDR.png", plot.data.MDR, width = 9.83, height = 6.17)


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
# Number of posterior dist elements = [(ni - nb) / nt] * nc = [(45000 - 5000) / 8] * 3 = 15000
ni <- 45000 # number of iterations in each chain
nb <- 5000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 3 # number of chains



##########
###### 2A. Fit MDR thermal responses with uniform priors (Arctic species): Briere ----
##########

##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)


##### Set data
data <- data.MDR.arctic

##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp

##### Set priors
prior <- data.frame(q = c(0, 1),
                    T0 = c(0, 20),
                    Tm = c(20, 45)
)

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, 
                 N.Temp.xs = N.Temp.xs, prior = prior)

##### Run JAGS -----
MDR.arctic.bri.uni <- jags(data = jag.data,
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
# save(MDR.arctic.bri.uni, file = "R-scripts/R2jags-objects/MDR.arctic.bri.uni.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/MDR.arctic.bri.uni.Rdata")


## Diagnostics ----
##### Examine output
MDR.arctic.bri.uni$BUGSoutput$summary[1:5,]
mcmcplot(MDR.arctic.bri.uni)

# Extract the DIC for future model comparisons
MDR.arctic.bri.uni$BUGSoutput$DIC

## Plot data + fit ----
df.MDR.arctic.bri.uni <- data.frame(MDR.arctic.bri.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.MDR.arctic.bri.uni)

##### Plot
plot.MDR.arctic.bri.uni <- df.MDR.arctic.bri.uni %>% 
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

plot.MDR.arctic.bri.uni

# ggsave("figures/MDR.arctic.bri.uni.png", plot.MDR.arctic.bri.uni, 
#        width = 10.3, height = 5.6)

##########
###### 2B. Fit MDR thermal responses for priors (non-Arctic): Briere ----
##########

##### Temp sequence for derived quantity calculations
# For priors - fewer temps for derived calculations makes it go faster
Temp.xs <- seq(0, 45, 0.5)
N.Temp.xs <-length(Temp.xs)

##### Set data
data <- data.MDR.nonarctic

##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp

## Set priors
prior <- data.frame(q = c(0, 1),
                    T0 = c(0, 20),
                    Tm = c(20, 45)
)


##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, 
                 N.Temp.xs = N.Temp.xs, prior = prior)

##### Run JAGS -----

# This code took half an hour to run!
MDR.nonarctic.bri.uni <- jags(data = jag.data,
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
# save(MDR.nonarctic.bri.uni, file = "R-scripts/R2jags-objects/MDR.nonarctic.bri.uni.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/MDR.nonarctic.bri.uni.Rdata")


## Diagnostics ----
##### Examine output
MDR.nonarctic.bri.uni$BUGSoutput$summary[1:5,]
mcmcplot(MDR.nonarctic.bri.uni)

# Extract the DIC for future model comparisons
MDR.nonarctic.bri.uni$BUGSoutput$DIC

## Plot data + fit ----
df.MDR.nonarctic.bri.uni <- data.frame(MDR.nonarctic.bri.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.MDR.nonarctic.bri.uni)

##### Plot
plot.df.MDR.nonarctic.bri.uni <- df.MDR.nonarctic.bri.uni %>% 
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

plot.df.MDR.nonarctic.bri.uni

# ggsave("figures/MDR.nonarctic.bri.uni.png", plot.df.MDR.nonarctic.bri.uni, 
#        width = 10.3, height = 5.6)


##########
###### 2C. Fit gamma distributions to MDR prior thermal responses: Briere ----
##########

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
MDR.arctic.prior.cf.dists <- data.frame(q = as.vector(MDR.nonarctic.bri.uni$BUGSoutput$sims.list$cf.q),
                                          T0 = as.vector(MDR.nonarctic.bri.uni$BUGSoutput$sims.list$cf.T0),
                                          Tm = as.vector(MDR.nonarctic.bri.uni$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
MDR.arctic.prior.gamma.fits = apply(MDR.arctic.prior.cf.dists, 2, 
                                      function(df) fitdistr(df, "gamma")$estimate)


MDR.hypers <- MDR.arctic.prior.gamma.fits
# save(MDR.hypers, file = "R-scripts/R2jags-objects/MDRhypers.bri.Rsave")


##########
###### 2D. Fit MDR thermal responses with data-informed priors (Arctic): Briere ----
##########

load("R-scripts/R2jags-objects/MDRhypers.bri.Rsave")
MDR.arctic.prior.gamma.fits <- MDR.hypers

##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)

##### Set data
data <- data.MDR.arctic
hypers <- MDR.arctic.prior.gamma.fits * 0.1

##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, 
                 N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS -----
MDR.arctic.bri.inf <- jags(data = jag.data,
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
# save(MDR.arctic.bri.inf, file = "R-scripts/R2jags-objects/MDR.arctic.bri.inf.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/MDR.arctic.bri.inf.Rdata")


## Diagnostics ----
##### Examine output
MDR.arctic.bri.inf$BUGSoutput$summary[1:5,]
mcmcplot(MDR.arctic.bri.inf)

# Extract the DIC for future model comparisons
MDR.arctic.bri.inf$BUGSoutput$DIC

## Plot data + fit ----
df.MDR.arctic.bri.inf <- data.frame(MDR.arctic.bri.inf$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.MDR.arctic.bri.inf)

##### Plot
# plot.MDR.arctic.bri.inf <- df.MDR.arctic.bri.inf %>%
#   ggplot(aes(x = temp)) +
#   geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "pink", alpha = 0.5) +
#   geom_line(aes(y = mean), color = "red", linewidth = 1) +
#   geom_point(data = data, aes(x = temp, y = trait), size = 2) +
#   # Customize the axes and labels
#   #scale_x_continuous(limits = c(0, 41)) +
#   #scale_y_continuous(limits = c(-0.005, 0.19)) +
#   labs(
#     x = expression(paste("Temperature (", degree, "C)")),
#     y = "Development rate (days-1)"
#   ) +
#   theme_bw()
# 
# plot.MDR.arctic.bri.inf

plot.MDR.arctic.bri.inf <- df.MDR.arctic.bri.inf %>%
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "#4363d8", alpha = 0.5) +
  geom_line(aes(y = mean), color = "blue", linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait, color = species), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) +
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Development rate (days-1)"
  ) +
  scale_colour_manual(name = "Species", labels = "Ae. nigripes", values = "black") +
  theme_bw()

plot.MDR.arctic.bri.inf

# ggsave("figures/MDR.arctic.bri.inf.png", plot.MDR.arctic.bri.inf, 
#        width = 10.3, height = 5.6)


##########
###### 2E. Plot all TPCs in the same graph (for comparison): Briere ----
##########

# Add an identifying column in each model output dataframe
df.MDR.arctic.bri.uni <- df.MDR.arctic.bri.uni %>% 
  mutate(type = "Arctic uniform")

df.MDR.arctic.bri.inf <- df.MDR.arctic.bri.inf %>% 
  mutate(type = "Arctic informative")

# Combine the three dataframes
df.all <- rbind(df.MDR.arctic.bri.uni, df.MDR.arctic.bri.inf)

# Plot
plot.all <- df.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = mean, color = type), linewidth = 1) +
  geom_point(data = data.MDR.arctic, aes(x = temp, y = trait), size = 2) +
  #geom_point(data = data.MDR.nonarctic, aes(x = temp, y = trait), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Development rate (days-1)"
  ) +
  # Customize the colours
  ## ribbon
  scale_fill_manual(values = c("Arctic uniform" = "#4363d8",
                             "Arctic informative" = "pink")) +
  ## line
  scale_color_manual(values = c("Arctic uniform" = "blue",
                               "Arctic informative" = "red")) +
  theme_bw()

plot.all

# ggsave("figures/MDR.arctic.bri.all.png", plot.all, width = 10.3, height = 5.6)

MDR.arctic.bri.uni$BUGSoutput$DIC
MDR.arctic.bri.inf$BUGSoutput$DIC


##########
###### 3A. Fit MDR thermal responses with uniform priors (Arctic): Quadratic ----
##########

##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)

##### Set data
data <- data.MDR.arctic

## Set priors
prior <- data.frame(q = c(0, 1),
                    T0 = c(0, 20),
                    Tm = c(20, 45)
)

##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, 
                 N.Temp.xs = N.Temp.xs, prior = prior)

##### Run JAGS -----
MDR.arctic.quad.uni <- jags(data = jag.data,
                            inits = inits,
                            parameters.to.save = parameters,
                            model.file = "R-scripts/quad_T.txt",
                            n.thin = nt,
                            n.chains = nc,
                            n.burnin = nb,
                            n.iter = ni,
                            DIC = T,
                            working.directory = getwd()
)

## Save the model as Rdata 
# save(MDR.arctic.quad.uni, file = "R-scripts/R2jags-objects/MDR.arctic.quad.uni.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/MDR.arctic.quad.uni.Rdata")


## Diagnostics ----
##### Examine output
MDR.arctic.quad.uni$BUGSoutput$summary[1:5,]
mcmcplot(MDR.arctic.quad.uni)

# Extract the DIC for future model comparisons
MDR.arctic.quad.uni$BUGSoutput$DIC

## Plot data + fit ----
df.MDR.arctic.quad.uni <- data.frame(MDR.arctic.quad.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.MDR.arctic.quad.uni)

##### Plot
plot.MDR.arctic.quad.uni <- df.MDR.arctic.quad.uni %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "#4363d8", alpha = 0.5) +
  geom_line(aes(y = mean), color = "blue", linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait), position = "jitter", size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Development rate (days-1)"
  ) +
  theme_bw()

plot.MDR.arctic.quad.uni

# ggsave("figures/MDR.arctic.quad.uni.png", plot.MDR.arctic.quad.uni, 
#        width = 10.3, height = 5.6)

##########
###### 3B. Fit MDR thermal responses for priors (non-Arctic): Quadratic ----
##########

##### Temp sequence for derived quantity calculations
# For priors - fewer temps for derived calculations makes it go faster
Temp.xs <- seq(0, 45, 0.5)
N.Temp.xs <-length(Temp.xs)


##### Set data
data <- data.MDR.nonarctic

## Set priors
prior <- data.frame(q = c(0, 1),
                    T0 = c(0, 20),
                    Tm = c(20, 45)
)

##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, 
                 N.Temp.xs = N.Temp.xs, prior = prior)

##### Run JAGS -----

MDR.nonarctic.quad.uni <- jags(data = jag.data,
                             inits = inits,
                             parameters.to.save = parameters,
                             model.file = "R-scripts/quad_T.txt",
                             n.thin = nt,
                             n.chains = nc,
                             n.burnin = nb,
                             n.iter = ni,
                             DIC = T,
                             working.directory = getwd()
)

## Save the model as Rdata 
# save(MDR.nonarctic.quad.uni, file = "R-scripts/R2jags-objects/MDR.nonarctic.quad.uni.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/MDR.nonarctic.quad.uni.Rdata")


## Diagnostics ----
##### Examine output
MDR.nonarctic.quad.uni$BUGSoutput$summary[1:5,]
mcmcplot(MDR.nonarctic.quad.uni)

# Extract the DIC for future model comparisons
MDR.nonarctic.quad.uni$BUGSoutput$DIC

## Plot data + fit ----
df.MDR.nonarctic.quad.uni <- data.frame(MDR.nonarctic.quad.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.MDR.nonarctic.quad.uni)

##### Plot
plot.df.MDR.nonarctic.quad.uni <- df.MDR.nonarctic.quad.uni %>% 
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

plot.df.MDR.nonarctic.quad.uni

# ggsave("figures/MDR.nonarctic.quad.uni.png", plot.df.MDR.nonarctic.quad.uni, 
#        width = 10.3, height = 5.6)


##########
###### 3C. Fit gamma distributions to MDR prior thermal responses: Quadratic ----
##########

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
MDR.arctic.prior.cf.dists <- data.frame(q = as.vector(MDR.nonarctic.quad.uni$BUGSoutput$sims.list$cf.q),
                                          T0 = as.vector(MDR.nonarctic.quad.uni$BUGSoutput$sims.list$cf.T0),
                                          Tm = as.vector(MDR.nonarctic.quad.uni$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
MDR.arctic.prior.gamma.fits = apply(MDR.arctic.prior.cf.dists, 2, 
                                      function(df) fitdistr(df, "gamma")$estimate)


MDR.hypers <- MDR.arctic.prior.gamma.fits
save(MDR.hypers, file = "R-scripts/R2jags-objects/MDRhypers.quad.Rsave")


##########
###### 3D. Fit MDR thermal responses with data-informed priors (Arctic): Quadratic ----
##########

load("R-scripts/R2jags-objects/MDRhypers.quad.Rsave")
MDR.arctic.prior.gamma.fits <- MDR.hypers

##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)

##### Set data
data <- data.MDR.arctic
hypers <- MDR.arctic.prior.gamma.fits * 0.1

##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, 
                 N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS -----
MDR.arctic.quad.inf <- jags(data = jag.data,
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
# save(MDR.arctic.quad.inf, file = "R-scripts/R2jags-objects/MDR.arctic.quad.inf.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/MDR.arctic.quad.inf.Rdata")


## Diagnostics ----
##### Examine output
MDR.arctic.quad.inf$BUGSoutput$summary[1:5,]
mcmcplot(MDR.arctic.quad.inf)

# Extract the DIC for future model comparisons
MDR.arctic.quad.inf$BUGSoutput$DIC

## Plot data + fit ----
df.MDR.arctic.quad.inf <- data.frame(MDR.arctic.quad.inf$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.MDR.arctic.quad.inf)

##### Plot
plot.MDR.arctic.quad.inf <- df.MDR.arctic.quad.inf %>% 
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

plot.MDR.arctic.quad.inf

# ggsave("figures/MDR.arctic.quad.inf.png", plot.MDR.arctic.quad.inf, 
#        width = 10.3, height = 5.6)


##########
###### 3E. Plot all three TPCs in the same graph (for comparison) ----
##########

# Add an identifying column in each model output dataframe
df.MDR.arctic.quad.uni <- df.MDR.arctic.quad.uni %>% 
  mutate(type = "Arctic uniform")

df.MDR.arctic.quad.inf <- df.MDR.arctic.quad.inf %>% 
  mutate(type = "Arctic informative")

# Combine the three dataframes
df.all <- rbind(df.MDR.arctic.quad.uni, df.MDR.arctic.quad.inf)

##### Plot
plot.all <- df.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = mean, color = type), linewidth = 1) +
  geom_point(data = data.MDR.arctic, aes(x = temp, y = trait), size = 2) +
  #geom_point(data = data.MDR.nonarctic, aes(x = temp, y = trait), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Development rate (days-1)"
  ) +
  # Customize the colours
  ## ribbon
  scale_fill_manual(values = c("Arctic uniform" = "#4363d8", 
                               "Arctic informative" = "pink")) +
  ## line
  scale_color_manual(values = c("Arctic uniform" = "blue", 
                                "Arctic informative" = "red")) +
  theme_bw()

plot.all

# ggsave("figures/MDR.all.quad.png", plot.all,
#        width = 10.3, height = 5.6)


##### Plot all arctic TPCs for comparison ----
# Add an identifying column in each model output dataframe
df.MDR.arctic.bri.uni <- df.MDR.arctic.bri.uni %>% 
  mutate(type = "Briere (uni)")

df.MDR.arctic.bri.inf <- df.MDR.arctic.bri.inf %>% 
  mutate(type = "Briere (inf)")

df.MDR.arctic.quad.uni <- df.MDR.arctic.quad.uni %>% 
  mutate(type = "Quadratic (uni)")

df.MDR.arctic.quad.inf <- df.MDR.arctic.quad.inf %>% 
  mutate(type = "Quadratic (inf)")

# Combine the three dataframes
df.all <- rbind(df.MDR.arctic.bri.uni, df.MDR.arctic.bri.inf, 
                df.MDR.arctic.quad.uni, df.MDR.arctic.quad.inf)

##### Plot
plot.all <- df.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = mean, color = type), linewidth = 1) +
  geom_point(data = data.MDR.arctic, aes(x = temp, y = trait), size = 2) +
  #geom_point(data = data.MDR.nonarctic, aes(x = temp, y = trait), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Development rate (days-1)"
  ) +
  # Customize the colours
  scale_fill_jco() +
  scale_color_jco() +
  theme_bw()

plot.all

# ggsave("figures/MDR.all.arctic.png", plot.all,
#        width = 10.3, height = 5.6)


#### DIC ----
MDR.arctic.bri.uni$BUGSoutput$DIC
MDR.arctic.bri.inf$BUGSoutput$DIC # This is the best fitting TPC
MDR.arctic.quad.uni$BUGSoutput$DIC
MDR.arctic.quad.inf$BUGSoutput$DIC 

##########
###### 4. Process and save model output for plotting ----
##########

## Analyze TPC model
MDR.TPC.analysis <- extractTPC(MDR.arctic.bri.inf, "MDR", Temp.xs)
MDR.predictions.summary <- MDR.TPC.analysis[[1]]
MDR.params.summary <- MDR.TPC.analysis[[2]]
MDR.params.fullposts <- MDR.TPC.analysis[[3]]

write_csv(MDR.predictions.summary, "data-processed/MDR.predictions.summary.csv")
write_csv(MDR.params.summary, "data-processed/MDR.params.summary.csv")
write_csv(MDR.params.fullposts, "data-processed/MDR.params.fullposts.csv")
