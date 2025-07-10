## Lilian Chan, University of Guelph
## Arctic vector-borne disease transmission suitability model
##
## Purpose: use Bayesian inference (JAGS) to fit TPCs for larval-to-adult
## survival (pLA) for Arctic mosquito species using data from Aedes vexans (Brust 1967)
## and from 3 non-Arctic mosquito species (Ae. nigromaculis, Ae. sollicitans, Ae. triseriatus)
##     1) with uniform priors; and 
##     2) with data-informed priors from the non-Arctic species.
##
## 
## Table of content:
##    0. Set-up workspace
##
##    1. MCMC settings for all models
##
##    2. Fitting TPC (Briere)
##        A. Fit pLA thermal responses with uniform priors (Arctic species)
##        B. Fit pLA thermal responses for priors (non-Arctic species)
##            i. All non-Arctic species
##           ii. only Ae. nigromaculis
##          iii. only Ae. sollicitans
##           iv. only Ae. triseriatus
##
##        C. Fit gamma distributions to pLA prior thermal responses
##        D. Fit pLA thermal responses with data-informed priors (Arctic)
##        E. Plot all three TPCs in the same graph (for comparison)
##
##    3. Fitting TPC (Quadratic)
##        A. Fit pLA thermal responses with uniform priors (Arctic)
##        B. Fit pLA thermal responses for priors (non-Arctic species)
##            i. All non-Arctic species
##           ii. only Ae. nigromaculis
##          iii. only Ae. sollicitans
##           iv. only Ae. triseriatus
##
##        C. Fit gamma distributions to pLA prior thermal responses
##        D. Fit pLA thermal responses with data-informed priors (Arctic)


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
data <- read_csv("data/data-processed/TraitData_pLA.csv")
unique(data$species)

     
# Subset data
## Arctic species
data.pLA.arctic <- subset(data, species  == "vexans")

## Non-rctic species
data.pLA.nonarctic <- subset(data, species != "vexans")

## Ae. nigromaculis
data.pLA.nigromaculis <- subset(data, species == "nigromaculis")

## Ae. sollicitans
data.pLA.sollicitans <- subset(data, species == "sollicitans")

## Ae. triseriatus
data.pLA.triseriatus <- subset(data, species == "triseriatus")

## Plot raw data
plot.data.pLA <- data %>% 
  mutate(type = c(rep("Arctic", 6), rep("non-Arctic", 31))) %>% 
  ggplot(aes(x = temp, y = trait)) +
  geom_point(aes(colour = species, shape = citation)) +
  geom_line(data = ~filter(.x, species != "triseriatus"),
            aes(colour = species)) +
  labs(y = "Larval survival", x = "Temperature ºC") +
  scale_colour_discrete(name = "Species", labels = c("Ae. nigromaculis",
                                                     "Ae. sollicitans",
                                                     "Ae. triseriatus",
                                                     "Ae. vexans"
  )) +
  facet_grid(rows = vars(type)) +
  theme_bw()

plot.data.pLA

# ggsave("figures/raw_data/plot.data.pLA.png", plot.data.pLA, width = 9.83, height = 6.17)



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
###### 2A. Fit pLA thermal responses with uniform priors (Arctic): Briere ----
##########

##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)


##### Set data
data <- data.pLA.arctic

##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
pLA.arctic.bri.uni <- jags(
  data = jag.data,
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
# save(pLA.arctic.bri.uni, file = "R-scripts/R2jags-objects/pLA.arctic.bri.uni.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/pLA.arctic.bri.uni.Rdata")


## Diagnostics ----
##### Examine output
pLA.arctic.bri.uni$BUGSoutput$summary[1:5,]
mcmcplot(pLA.arctic.bri.uni)

# Extract the DIC for future model comparisons
pLA.arctic.bri.uni$BUGSoutput$DIC

## Plot data + fit ----
df.pLA.arctic.bri.uni <- data.frame(pLA.arctic.bri.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.pLA.arctic.bri.uni)

##### Plot
plot.pLA.arctic.bri.uni <- df.pLA.arctic.bri.uni %>%
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.),
              fill = "#4363d8",
              alpha = 0.5) +
  geom_line(aes(y = mean), color = "blue", linewidth = 1) +
  geom_point(data = data,
             aes(x = temp, y = trait),
             size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) +
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(x = expression(paste("Temperature (", degree, "C)")), y = "Larval survival (pLA)") +
  # Customize legend
  # scale_color_discrete(name = "Species",
  #                      labels = c("V. eleguneniensis", "S. tundra")) +
  theme_bw()

plot.pLA.arctic.bri.uni

# ggsave("figures/pLA.arctic.bri.uni.png", plot.pLA.arctic.bri.uni,
#        width = 10.3, height = 5.6)


##########
###### 2B i. Fit pLA thermal responses for priors (all non-arctic species): Briere ----
##########

##### Temp sequence for derived quantity calculations
# For priors - fewer temps for derived calculations makes it go faster
Temp.xs <- seq(0, 45, 0.5)
N.Temp.xs <-length(Temp.xs)


##### Set data
data <- data.pLA.nonarctic


##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
pLA.nonarctic.bri.uni <- jags(
  data = jag.data,
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

## Add random effects ----

## Create a unique id for each species-study combination
data <- data %>% 
  group_by(species, host.species, citation) %>% 
  mutate(unique_id = cur_group_id())

##### inits Function
inits <- function(){list(
  cf.q = 0.01,
  cf.Tm = 35,
  cf.T0 = 5,
  cf.sigma = rlnorm(1),
  sigma_q = rlnorm(0.1))}

##### Parameters to Estimate
parameters <- c("cf.q", "cf.T0", "cf.Tm", "cf.sigma", "sigma_q", "z.trait.mu.pred")

##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp
unique.id <- as.integer(as.factor(data$unique_id))
Nids <- max(unique.id)

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, Nids = Nids, unique.id = unique.id)

##### Run JAGS
pLA.nonarctic.bri.uni <- jags(
  data = jag.data,
  inits = inits,
  parameters.to.save = parameters,
  model.file = "R-scripts/briere_T_randeff.txt",
  n.thin = nt,
  n.chains = nc,
  n.burnin = nb,
  n.iter = ni,
  DIC = T,
  working.directory = getwd()
)

## Random effects END ----

## Save the model as Rdata 
# save(pLA.nonarctic.bri.uni, file = "R-scripts/R2jags-objects/pLA.nonarctic.bri.uni.Rdata")

# Read the .Rdata
load("R-scripts/R2jags-objects/pLA.nonarctic.bri.uni.Rdata")


## Diagnostics ----
##### Examine output
pLA.nonarctic.bri.uni$BUGSoutput$summary[1:5,]
mcmcplot(pLA.nonarctic.bri.uni)

# Extract the DIC for future model comparisons
pLA.nonarctic.bri.uni$BUGSoutput$DIC

## Plot data + fit ----
df.pLA.nonarctic.bri.uni <- data.frame(pLA.nonarctic.bri.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.pLA.nonarctic.bri.uni)

##### Plot
plot.pLA.nonarctic.bri.uni <- df.pLA.nonarctic.bri.uni %>%
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.),
              fill = "#4363d8",
              alpha = 0.5) +
  geom_line(aes(y = mean), color = "blue", linewidth = 1) +
  geom_point(data = data,
             aes(x = temp, y = trait, colour = species),
             size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) +
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(x = expression(paste("Temperature (", degree, "C)")), y = "Larval survival") +
  # Customize legend
  scale_color_discrete(name = "Species",
                       labels = c("Ae. sollicitans", "Ae. triseriatus", "Ae. nigromaculis")) +
  theme_bw()

plot.pLA.nonarctic.bri.uni

# ggsave("figures/pLA.nonarctic.bri.uni.png", plot.pLA.nonarctic.bri.uni,
#        width = 10.3, height = 5.6)


##########
###### 2B ii. Fit pLA thermal responses for priors in Ae. nigromaculis: Briere ----
##########

##### Temp sequence for derived quantity calculations
# For priors - fewer temps for derived calculations makes it go faster
Temp.xs <- seq(0, 45, 0.5)
N.Temp.xs <-length(Temp.xs)

##### Set data
data <- data.pLA.nigromaculis

##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS

pLA.nigromaculis.bri.uni <- jags(
  data = jag.data,
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
# save(pLA.nigromaculis.bri.uni, file = "R-scripts/R2jags-objects/pLA.nigromaculis.bri.uni.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/pLA.nigromaculis.bri.uni.Rdata")


## Diagnostics ----
##### Examine output
pLA.nigromaculis.bri.uni$BUGSoutput$summary[1:5,]
mcmcplot(pLA.nigromaculis.bri.uni)

# Extract the DIC for future model comparisons
pLA.nigromaculis.bri.uni$BUGSoutput$DIC

## Plot data + fit ----
df.pLA.nigromaculis.bri.uni <- data.frame(pLA.nigromaculis.bri.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.pLA.nigromaculis.bri.uni)

##### Plot
plot.df.pLA.nigromaculis.bri.uni <- df.pLA.nigromaculis.bri.uni %>% 
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
    y = "Larval survival (pLA)"
  ) +
  theme_bw()

plot.df.pLA.nigromaculis.bri.uni

# ggsave("figures/pLA.nigromaculis.bri.uni.png", plot.df.pLA.nigromaculis.bri.uni, 
#        width = 10.3, height = 5.6)



##########
###### 2B iii. Fit pLA thermal responses for priors in Ae. sollicitans: Briere ----
##########

##### Temp sequence for derived quantity calculations
# For priors - fewer temps for derived calculations makes it go faster
Temp.xs <- seq(0, 45, 0.5)
N.Temp.xs <-length(Temp.xs)

##### Set data
data <- data.pLA.sollicitans

##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS

pLA.sollicitans.bri.uni <- jags(
  data = jag.data,
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
# save(pLA.sollicitans.bri.uni, file = "R-scripts/R2jags-objects/pLA.sollicitans.bri.uni.Rdata")

# Read the .Rdata
load("R-scripts/R2jags-objects/pLA.sollicitans.bri.uni.Rdata")


## Diagnostics ----
##### Examine output
pLA.sollicitans.bri.uni$BUGSoutput$summary[1:5,]
mcmcplot(pLA.sollicitans.bri.uni)

# Extract the DIC for future model comparisons
pLA.sollicitans.bri.uni$BUGSoutput$DIC

## Plot data + fit ----
df.pLA.sollicitans.bri.uni <- data.frame(pLA.sollicitans.bri.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.pLA.sollicitans.bri.uni)

##### Plot
plot.df.pLA.sollicitans.bri.uni <- df.pLA.sollicitans.bri.uni %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "grey", alpha = 0.5) +
  geom_line(aes(y = mean), color = "#868686FF", linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait), size = 2 
             , position = "jitter"
  ) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 1)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Larval survival (pLA)"
  ) +
  theme_bw()

plot.df.pLA.sollicitans.bri.uni

# ggsave("figures/pLA.sollicitans.bri.uni.png", plot.df.pLA.sollicitans.bri.uni, 
#        width = 10.3, height = 5.6)



##########
###### 2B iv. Fit pLA thermal responses for priors (Ae. triseriatus): Briere ----
##########

##### Temp sequence for derived quantity calculations
# For priors - fewer temps for derived calculations makes it go faster
Temp.xs <- seq(0, 45, 0.5)
N.Temp.xs <-length(Temp.xs)

##### Set data
data <- data.pLA.triseriatus

##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS

pLA.triseriatus.bri.uni <- jags(
  data = jag.data,
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
# save(pLA.triseriatus.bri.uni, file = "R-scripts/R2jags-objects/pLA.triseriatus.bri.uni.Rdata")

# Read the .Rdata
load("R-scripts/R2jags-objects/pLA.triseriatus.bri.uni.Rdata")


## Diagnostics ----
##### Examine output
pLA.triseriatus.bri.uni$BUGSoutput$summary[1:5,]
mcmcplot(pLA.triseriatus.bri.uni)

# Extract the DIC for future model comparisons
pLA.triseriatus.bri.uni$BUGSoutput$DIC

## Plot data + fit ----
df.pLA.triseriatus.bri.uni <- data.frame(pLA.triseriatus.bri.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.pLA.triseriatus.bri.uni)

##### Plot
plot.df.pLA.triseriatus.bri.uni <- df.pLA.triseriatus.bri.uni %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "grey", alpha = 0.5) +
  geom_line(aes(y = mean), color = "#868686FF", linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait), size = 2 
             , position = "jitter"
  ) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 1)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Larval survival (pLA)"
  ) +
  theme_bw()

plot.df.pLA.triseriatus.bri.uni

# ggsave("figures/pLA.triseriatus.bri.uni.png", plot.df.pLA.triseriatus.bri.uni, 
#        width = 10.3, height = 5.6)



##########
###### 2C. Fit gamma distributions to pLA prior thermal responses: Briere ----
##########

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
pLA.arctic.prior.cf.dists <- data.frame(q = as.vector(pLA.nonarctic.bri.uni$BUGSoutput$sims.list$cf.q),
                                          T0 = as.vector(pLA.nonarctic.bri.uni$BUGSoutput$sims.list$cf.T0),
                                          Tm = as.vector(pLA.nonarctic.bri.uni$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
pLA.arctic.prior.gamma.fits = apply(pLA.arctic.prior.cf.dists, 2, 
                                      function(df) fitdistr(df, "gamma")$estimate)


pLA.hypers <- pLA.arctic.prior.gamma.fits
# save(pLA.hypers, file = "R-scripts/R2jags-objects/pLAhypers.bri.Rsave")


##########
###### 2D. Fit pLA thermal responses with data-informed priors (Arctic): Briere ----
##########

load("R-scripts/R2jags-objects/pLAhypers.bri.Rsave")
pLA.arctic.prior.gamma.fits <- pLA.hypers

##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)

##### Set data
data <- data.pLA.arctic
hypers <- pLA.arctic.prior.gamma.fits * 0.1

##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, 
                 N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS -----
pLA.arctic.bri.inf <- jags(data = jag.data,
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
# save(pLA.arctic.bri.inf, file = "R-scripts/R2jags-objects/pLA.arctic.bri.inf.Rdata")

# Read the .Rdata
load("R-scripts/R2jags-objects/pLA.arctic.bri.inf.Rdata")


## Diagnostics ----
##### Examine output
pLA.arctic.bri.inf$BUGSoutput$summary[1:5,]
mcmcplot(pLA.arctic.bri.inf)

# Extract the DIC for future model comparisons
pLA.arctic.bri.inf$BUGSoutput$DIC

## Plot data + fit ----
df.pLA.arctic.bri.inf <- data.frame(pLA.arctic.bri.inf$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.pLA.arctic.bri.inf)

##### Plot
plot.pLA.arctic.bri.inf <- df.pLA.arctic.bri.inf %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "pink", alpha = 0.5) +
  geom_line(aes(y = mean), color = "red", linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Larval survival (pLA)"
  ) +
  theme_bw()

plot.pLA.arctic.bri.inf

# ggsave("figures/pLA.arctic.bri.inf.png", plot.pLA.arctic.bri.inf, 
#        width = 10.3, height = 5.6)


##########
###### 2E. Plot all three TPCs in the same graph (for comparison): Briere ----
##########

# Add an identifying column in each model output dataframe
df.pLA.arctic.bri.uni <- df.pLA.arctic.bri.uni %>% 
  mutate(type = "Arctic uniform")

df.pLA.nonarctic.bri.uni <- df.pLA.nonarctic.bri.uni %>% 
  mutate(type = "non-Arctic uniform")

df.pLA.arctic.bri.inf <- df.pLA.arctic.bri.inf %>% 
  mutate(type = "Arctic informative")

# Combine the three dataframes
df.all <- rbind(df.pLA.arctic.bri.uni, df.pLA.nonarctic.bri.uni, df.pLA.arctic.bri.inf)

# Plot
plot.all <- df.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = mean, color = type), linewidth = 1) +
  geom_point(data = data.pLA.arctic, aes(x = temp, y = trait), size = 2) +
  #geom_point(data = data.pLA.nonarctic, aes(x = temp, y = trait), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Larval survival (pLA)"
  ) +
  # Customize the colours
  ## ribbon
  scale_fill_manual(values = c("Arctic uniform" = "#4363d8", 
                               "non-Arctic uniform" = "grey",
                               "Arctic informative" = "pink")) +
  ## line
  scale_color_manual(values = c("Arctic uniform" = "blue", 
                                "non-Arctic uniform" = "#868686FF",
                                "Arctic informative" = "red")) +
  theme_bw()

plot.all

# ggsave("figures/pLA.all.bri.png", plot.all, 
#        width = 10.3, height = 5.6)



##########
###### 3A. Fit pLA thermal responses with uniform priors (Arctic): Quadratic ----
##########

##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)

##### Set data
data <- data.pLA.arctic

##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# ##### Run JAGS -----
pLA.arctic.quad.uni <- jags(data = jag.data,
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
# save(pLA.arctic.quad.uni, file = "R-scripts/R2jags-objects/pLA.arctic.quad.uni.Rdata")

# Read the .Rdata
load("R-scripts/R2jags-objects/pLA.arctic.quad.uni.Rdata")


## Diagnostics ----
##### Examine output
pLA.arctic.quad.uni$BUGSoutput$summary[1:5,]
mcmcplot(pLA.arctic.quad.uni)

# Extract the DIC for future model comparisons
pLA.arctic.quad.uni$BUGSoutput$DIC

## Plot data + fit ----
df.pLA.arctic.quad.uni <- data.frame(pLA.arctic.quad.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.pLA.arctic.quad.uni)

##### Plot
plot.pLA.arctic.quad.uni <- df.pLA.arctic.quad.uni %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "#4363d8", alpha = 0.5) +
  geom_line(aes(y = mean), color = "blue", linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Larval survival (pLA)"
  ) +
  theme_bw()

plot.pLA.arctic.quad.uni

# ggsave("figures/pLA.arctic.quad.uni.png", plot.pLA.arctic.quad.uni, 
#        width = 10.3, height = 5.6)


##########
###### 3B i. Fit pLA thermal responses for priors (all non-arctic species): Quadratic ----
##########

##### Temp sequence for derived quantity calculations
# For priors - fewer temps for derived calculations makes it go faster
Temp.xs <- seq(0, 45, 0.5)
N.Temp.xs <-length(Temp.xs)


##### Set data
data <- data.pLA.nonarctic


##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
pLA.nonarctic.quad.uni <- jags(
  data = jag.data,
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
# save(pLA.nonarctic.quad.uni, file = "R-scripts/R2jags-objects/pLA.nonarctic.quad.uni.Rdata")

# Read the .Rdata
load("R-scripts/R2jags-objects/pLA.nonarctic.quad.uni.Rdata")


## Diagnostics ----
##### Examine output
pLA.nonarctic.quad.uni$BUGSoutput$summary[1:5,]
mcmcplot(pLA.nonarctic.quad.uni)

# Extract the DIC for future model comparisons
pLA.nonarctic.quad.uni$BUGSoutput$DIC

## Plot data + fit ----
df.pLA.nonarctic.quad.uni <- data.frame(pLA.nonarctic.quad.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.pLA.nonarctic.quad.uni)

##### Plot
plot.df.pLA.nonarctic.quad.uni <- df.pLA.nonarctic.quad.uni %>% 
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
    y = "Larval survival (pLA)"
  ) +
  theme_bw()

plot.df.pLA.nonarctic.quad.uni

# ggsave("figures/pLA.nonarctic.quad.uni.png", plot.df.pLA.nonarctic.quad.uni, 
#        width = 10.3, height = 5.6)


##########
###### 3C. Fit gamma distributions to pLA prior thermal responses: Quadratic ----
##########

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
pLA.arctic.prior.cf.dists <- data.frame(q = as.vector(pLA.nonarctic.quad.uni$BUGSoutput$sims.list$cf.q),
                                          T0 = as.vector(pLA.nonarctic.quad.uni$BUGSoutput$sims.list$cf.T0),
                                          Tm = as.vector(pLA.nonarctic.quad.uni$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
pLA.arctic.prior.gamma.fits = apply(pLA.arctic.prior.cf.dists, 2, 
                                      function(df) fitdistr(df, "gamma")$estimate)


pLA.hypers <- pLA.arctic.prior.gamma.fits
save(pLA.hypers, file = "R-scripts/R2jags-objects/pLAhypers.quad.Rsave")


##########
###### 3D. Fit pLA thermal responses with data-informed priors (Arctic): Quadratic ----
##########

load("R-scripts/R2jags-objects/pLAhypers.quad.Rsave")
pLA.arctic.prior.gamma.fits <- pLA.hypers

##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)

##### Set data
data <- data.pLA.arctic
hypers <- pLA.arctic.prior.gamma.fits * 0.1

##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, 
                 N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS -----
pLA.arctic.quad.inf <- jags(data = jag.data,
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
# save(pLA.arctic.quad.inf, file = "R-scripts/R2jags-objects/pLA.arctic.quad.inf.Rdata")

# Read the .Rdata
load("R-scripts/R2jags-objects/pLA.arctic.quad.inf.Rdata")


## Diagnostics ----
##### Examine output
pLA.arctic.quad.inf$BUGSoutput$summary[1:5,]
mcmcplot(pLA.arctic.quad.inf)

# Extract the DIC for future model comparisons
pLA.arctic.quad.inf$BUGSoutput$DIC

## Plot data + fit ----
df.pLA.arctic.quad.inf <- data.frame(pLA.arctic.quad.inf$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.pLA.arctic.quad.inf)

##### Plot
plot.pLA.arctic.quad.inf <- df.pLA.arctic.quad.inf %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "pink", alpha = 0.5) +
  geom_line(aes(y = mean), color = "red", linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Larval survival (pLA)"
  ) +
  theme_bw()

plot.pLA.arctic.quad.inf

# ggsave("figures/pLA.arctic.quad.inf.png", plot.pLA.arctic.quad.inf, 
#        width = 10.3, height = 5.6)


##########
###### 3E. Plot all three TPCs in the same graph (for comparison) ----
##########

# Add an identifying column in each model output dataframe
df.pLA.arctic.quad.uni <- df.pLA.arctic.quad.uni %>% 
  mutate(type = "Arctic uniform")

df.pLA.nonarctic.quad.uni <- df.pLA.nonarctic.quad.uni %>% 
  mutate(type = "Ae. nonarctic uniform")

df.pLA.arctic.quad.inf <- df.pLA.arctic.quad.inf %>% 
  mutate(type = "Arctic informative")

# Combine the three dataframes
df.all <- rbind(df.pLA.arctic.quad.uni, df.pLA.nonarctic.quad.uni, df.pLA.arctic.quad.inf)

##### Plot
plot.all <- df.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = mean, color = type), linewidth = 1) +
  geom_point(data = data.pLA.arctic, aes(x = temp, y = trait), size = 2) +
  #geom_point(data = data.pLA.nonarctic, aes(x = temp, y = trait), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Larval survival (pLA)"
  ) +
  # Customize the colours
  ## ribbon
  scale_fill_manual(values = c("Arctic uniform" = "#4363d8", 
                               "Ae. nonarctic uniform" = "grey",
                               "Arctic informative" = "pink")) +
  ## line
  scale_color_manual(values = c("Arctic uniform" = "blue", 
                                "Ae. nonarctic uniform" = "#868686FF",
                                "Arctic informative" = "red")) +
  theme_bw()

plot.all

# ggsave("figures/pLA.all.quad.png", plot.all,
#        width = 10.3, height = 5.6)

##### Plot all arctic TPCs for comparison ----
# Add an identifying column in each model output dataframe
df.pLA.arctic.bri.uni <- df.pLA.arctic.bri.uni %>% 
  mutate(type = "Briere (uni)")

df.pLA.arctic.bri.inf <- df.pLA.arctic.bri.inf %>% 
  mutate(type = "Briere (inf)")

df.pLA.arctic.quad.uni <- df.pLA.arctic.quad.uni %>% 
  mutate(type = "Quadratic (uni)")

df.pLA.arctic.quad.inf <- df.pLA.arctic.quad.inf %>% 
  mutate(type = "Quadratic (inf)")

# Combine the three dataframes
df.all <- rbind(df.pLA.arctic.bri.uni, df.pLA.arctic.bri.inf, 
                df.pLA.arctic.quad.uni, df.pLA.arctic.quad.inf)

##### Plot
plot.all <- df.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = mean, color = type), linewidth = 1) +
  geom_point(data = data.pLA.arctic, aes(x = temp, y = trait), size = 2) +
  #geom_point(data = data.pLA.nonarctic, aes(x = temp, y = trait), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Larval survival (pLA)"
  ) +
  # Customize the colours
  scale_fill_jco() +
  scale_color_jco() +
  theme_bw()

plot.all

# ggsave("figures/pLA.all.arctic.png", plot.all,
#        width = 10.3, height = 5.6)