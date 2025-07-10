## Lilian Chan, University of Guelph
## Arctic vector-borne disease transmission suitability model
##
## Purpose: use Bayesian inference (JAGS) to fit TPCs for parasite development 
## rate (PDR) for Arctic nematode using data from Varestrongylus eleguneniensis 
## (Kafle et al. 2018) and from Setaria tundra (Laaksonen et al. 2009).
##     1) with uniform priors; and 
##     2) with data-informed priors from Dirofilaria immitis and wuchereria bancrofti
##
## (V. eleguneniensis is a nematode infecting caribou and muskoxen in the Canadian
## Arctic. It is transmitted by gastropod.)
## 
## Table of content:
##    0. Set-up workspace
##
##    1. MCMC settings for all models
##
##    2. Fitting TPC (Briere)
##        A. Fit PDR thermal responses with uniform priors (Arctic species)
##        B. Fit PDR thermal responses for priors (non-Arctic species)
##            i. All non-Arctic species
##           ii. only D. immitis in Ae. trivittatus
##          iii. only D. immitis in Ae. aegypti
##           iv. only W. bancrofti in Ae. polynesiensis
##
##        C. Fit gamma distributions to PDR prior thermal responses
##        D. Fit PDR thermal responses with data-informed priors (Arctic)
##        E. Plot all three TPCs in the same graph (for comparison)
##
##    3. Fitting TPC (Quadratic)
##        A. Fit PDR thermal responses with uniform priors (Arctic)
##        B. Fit PDR thermal responses for priors (non-Arctic species)
##            i. All non-Arctic species
##           ii. only D. immitis in Ae. trivittatus
##          iii. only D. immitis in Ae. aegypti
##           iv. only W. bancrofti in Ae. polynesiensis
##        C. Fit gamma distributions to PDR prior thermal responses
##        D. Fit PDR thermal responses with data-informed priors (Arctic)


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
data <- read_csv("data/data-processed/TraitData_PDR.csv")
unique(data$species)

## Convert development time (1/PDR) to development rate (PDR)
data <- data %>% 
  mutate(trait = 1/trait) %>% 
  mutate(trait_name = "PDR") 
     
# Subset data
## Arctic species
data.PDR.arctic <- subset(data, species %in% c("eleguneniensis", "tundra"))

## Non-rctic species
data.PDR.nonarctic <- subset(data, !(species %in% c("eleguneniensis", "tundra")))

## D. immitis in Ae. trivittatus
data.PDR.immitis.trivittatus <- subset(data, species == "immitis" & 
                                         host.species == "trivittatus")
## D. immitis in Ae. aegypti
data.PDR.immitis.aegypti <- subset(data, species == "immitis" & 
                                     host.species == "aegypti")

## W. bancrofti in Ae. polynesiensis
data.PDR.bancrofti <- subset(data, species == "bancrofti")

# Plot the data
plot.data.PDR <- data %>% 
  mutate(type = c(rep("Arctic", 7), rep("non-Arctic", 18))) %>% 
  ggplot(aes(x = temp, y = trait)) +
  geom_point(aes(colour = citation)) +
  geom_line(aes(colour = citation)) +
  labs(y = "Parasite development rate (1/days)", x = "Temperature ºC") +
  scale_colour_discrete(name = "Species", labels = c("D. immitis (in Ae. Trivittatus)",
                                                     "V. eleguneniensis",
                                                     "S. tundra",
                                                     "W. bancrofti (in Ae. polynesiensis)",
                                                     "D. immitis (in Ae. aegypti)"
  )) +
  facet_grid(rows = vars(type)) +
  theme_bw()

plot.data.PDR

# ggsave("figures/raw_data/plot.data.PDR.png", plot.data.PDR, , width = 9.83, height = 6.17)



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
ni <- 25000 # number of iterations in each chain
nb <- 5000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 3 # number of chains


##########
###### 2A. Fit PDR thermal responses with uniform priors (Arctic): Briere ----
##########

##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)


##### Set data
data <- data.PDR.arctic

##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
PDR.arctic.bri.uni <- jags(
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
# save(PDR.arctic.bri.uni, file = "R-scripts/R2jags-objects/PDR.arctic.bri.uni.Rdata")

# Read the .Rdata
load("R-scripts/R2jags-objects/PDR.arctic.bri.uni.Rdata")


## Diagnostics ----
##### Examine output
PDR.arctic.bri.uni$BUGSoutput$summary[1:5,]
mcmcplot(PDR.arctic.bri.uni)

# Extract the DIC for future model comparisons
PDR.arctic.bri.uni$BUGSoutput$DIC

## Plot data + fit ----
df.PDR.arctic.bri.uni <- data.frame(PDR.arctic.bri.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.PDR.arctic.bri.uni)

##### Plot
plot.PDR.arctic.bri.uni <- df.PDR.arctic.bri.uni %>%
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
  labs(x = expression(paste("Temperature (", degree, "C)")), y = "Development rate (days-1)") +
  # Customize legend
  scale_color_discrete(name = "Species",
                       labels = c("V. eleguneniensis", "S. tundra")) +
  theme_bw()

plot.PDR.arctic.bri.uni

# ggsave("figures/PDR.arctic.bri.uni.png", plot.PDR.arctic.bri.uni,
#        width = 10.3, height = 5.6)


##########
###### 2B i. Fit PDR thermal responses for priors (all non-arctic species): Briere ----
##########

##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)


##### Set data
data <- data.PDR.nonarctic

## Create a unique id for each species-study combination
data <- data %>% 
  group_by(species, host.species, citation) %>% 
  mutate(unique_id = cur_group_id())

##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
PDR.nonarctic.bri.uni <- jags(
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
# save(PDR.nonarctic.bri.uni, file = "R-scripts/R2jags-objects/PDR.nonarctic.bri.uni.Rdata")

# Read the .Rdata
load("R-scripts/R2jags-objects/PDR.nonarctic.bri.uni.Rdata")


## Diagnostics ----
##### Examine output
PDR.nonarctic.bri.uni$BUGSoutput$summary[1:5,]
mcmcplot(PDR.nonarctic.bri.uni)

# Extract the DIC for future model comparisons
PDR.nonarctic.bri.uni$BUGSoutput$DIC

## Plot data + fit ----
df.PDR.nonarctic.bri.uni <- data.frame(PDR.nonarctic.bri.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.PDR.nonarctic.bri.uni)

##### Plot
plot.PDR.nonarctic.bri.uni <- df.PDR.nonarctic.bri.uni %>%
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.),
              fill = "#4363d8",
              alpha = 0.5) +
  geom_line(aes(y = mean), color = "blue", linewidth = 1) +
  geom_point(data = data,
             aes(x = temp, y = trait, colour = species, shape = host.species),
             size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) +
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(x = expression(paste("Temperature (", degree, "C)")), y = "Development rate (days-1)") +
  # Customize legend
  scale_color_discrete(name = "Parasite Species",
                       labels = c("W. bancrofti", "D. immitis")) +
  scale_shape_discrete(name = "Host Species",
                       labels = c("Ae. aegypti", "Ae. polynesiensis", "Ae. trivittatus")) +
  theme_bw()

plot.PDR.nonarctic.bri.uni

# ggsave("figures/PDR.nonarctic.bri.uni.png", plot.PDR.nonarctic.bri.uni,
#        width = 10.3, height = 5.6)


##########
###### 2B ii. Fit PDR thermal responses for priors D. immitis in Ae. trivittatus: Briere ----
##########

##### Temp sequence for derived quantity calculations
# For priors - fewer temps for derived calculations makes it go faster
Temp.xs <- seq(0, 45, 0.5)
N.Temp.xs <-length(Temp.xs)

##### Set data
data <- data.PDR.immitis.trivittatus

##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS

PDR.immitis.trivittatus.bri.uni <- jags(
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
#save(PDR.immitis.trivittatus.bri.uni, file = "R-scripts/R2jags-objects/PDR.immitis.trivittatus.bri.uni.Rdata")

# Read the .Rdata
load("R-scripts/R2jags-objects/PDR.immitis.trivittatus.bri.uni.Rdata")


## Diagnostics ----
##### Examine output
PDR.immitis.trivittatus.bri.uni$BUGSoutput$summary[1:5,]
mcmcplot(PDR.immitis.trivittatus.bri.uni)

# Extract the DIC for future model comparisons
PDR.immitis.trivittatus.bri.uni$BUGSoutput$DIC

## Plot data + fit ----
df.PDR.immitis.trivittatus.bri.uni <- data.frame(PDR.immitis.trivittatus.bri.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.PDR.immitis.trivittatus.bri.uni)

##### Plot
plot.df.PDR.immitis.trivittatus.bri.uni <- df.PDR.immitis.trivittatus.bri.uni %>% 
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

plot.df.PDR.immitis.trivittatus.bri.uni

# ggsave("figures/PDR.immitis.trivittatus.bri.uni.png", plot.df.PDR.immitis.trivittatus.bri.uni, 
#        width = 10.3, height = 5.6)



##########
###### 2B iii. Fit PDR thermal responses for priors (D. immitis in Ae. aegypti): Briere ----
##########

##### Temp sequence for derived quantity calculations
# For priors - fewer temps for derived calculations makes it go faster
Temp.xs <- seq(0, 45, 0.5)
N.Temp.xs <-length(Temp.xs)

##### Set data
data <- data.PDR.immitis.aegypti

##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS

PDR.immitis.aegypti.bri.uni <- jags(
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
#save(PDR.immitis.aegypti.bri.uni, file = "R-scripts/R2jags-objects/PDR.immitis.aegypti.bri.uni.Rdata")

# Read the .Rdata
load("R-scripts/R2jags-objects/PDR.immitis.aegypti.bri.uni.Rdata")


## Diagnostics ----
##### Examine output
PDR.immitis.aegypti.bri.uni$BUGSoutput$summary[1:5,]
mcmcplot(PDR.immitis.aegypti.bri.uni)

# Extract the DIC for future model comparisons
PDR.immitis.aegypti.bri.uni$BUGSoutput$DIC

## Plot data + fit ----
df.PDR.immitis.aegypti.bri.uni <- data.frame(PDR.immitis.aegypti.bri.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.PDR.immitis.aegypti.bri.uni)

##### Plot
plot.df.PDR.immitis.aegypti.bri.uni <- df.PDR.immitis.aegypti.bri.uni %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "grey", alpha = 0.5) +
  geom_line(aes(y = mean), color = "#868686FF", linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait), size = 2 
             , position = "jitter"
  ) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  scale_y_continuous(limits = c(-0.005, 1)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Development rate (days-1)"
  ) +
  theme_bw()

plot.df.PDR.immitis.aegypti.bri.uni

# ggsave("figures/PDR.immitis.aegypti.bri.uni.png", plot.df.PDR.immitis.aegypti.bri.uni, 
#        width = 10.3, height = 5.6)



##########
###### 2B iv. Fit PDR thermal responses for priors (W. bancrofti in Ae. polynesiensis): Briere ----
##########

##### Temp sequence for derived quantity calculations
# For priors - fewer temps for derived calculations makes it go faster
Temp.xs <- seq(0, 45, 0.5)
N.Temp.xs <-length(Temp.xs)

##### Set data
data <- data.PDR.bancrofti

##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS

PDR.bancrofti.bri.uni <- jags(
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
#save(PDR.bancrofti.bri.uni, file = "R-scripts/R2jags-objects/PDR.bancrofti.bri.uni.Rdata")

# Read the .Rdata
load("R-scripts/R2jags-objects/PDR.bancrofti.bri.uni.Rdata")


## Diagnostics ----
##### Examine output
PDR.bancrofti.bri.uni$BUGSoutput$summary[1:5,]
mcmcplot(PDR.bancrofti.bri.uni)

# Extract the DIC for future model comparisons
PDR.bancrofti.bri.uni$BUGSoutput$DIC

## Plot data + fit ----
df.PDR.bancrofti.bri.uni <- data.frame(PDR.bancrofti.bri.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.PDR.bancrofti.bri.uni)

##### Plot
plot.df.PDR.bancrofti.bri.uni <- df.PDR.bancrofti.bri.uni %>% 
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
    y = "Development rate (days-1)"
  ) +
  theme_bw()

plot.df.PDR.bancrofti.bri.uni

# ggsave("figures/PDR.bancrofti.bri.uni.png", plot.df.PDR.bancrofti.bri.uni, 
#        width = 10.3, height = 5.6)



##########
###### 2C. Fit gamma distributions to PDR prior thermal responses: Briere ----
##########

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
PDR.arctic.prior.cf.dists <- data.frame(q = as.vector(PDR.nonarctic.bri.uni$BUGSoutput$sims.list$cf.q),
                                          T0 = as.vector(PDR.nonarctic.bri.uni$BUGSoutput$sims.list$cf.T0),
                                          Tm = as.vector(PDR.nonarctic.bri.uni$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
PDR.arctic.prior.gamma.fits = apply(PDR.arctic.prior.cf.dists, 2, 
                                      function(df) fitdistr(df, "gamma")$estimate)


PDR.hypers <- PDR.arctic.prior.gamma.fits
#save(PDR.hypers, file = "R-scripts/R2jags-objects/PDRhypers.bri.Rsave")


##########
###### 2D. Fit PDR thermal responses with data-informed priors (Arctic): Briere ----
##########

load("R-scripts/R2jags-objects/PDRhypers.bri.Rsave")
PDR.arctic.prior.gamma.fits <- PDR.hypers

##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)

##### Set data
data <- data.PDR.arctic
hypers <- PDR.arctic.prior.gamma.fits * 0.1

##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, 
                 N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS -----
PDR.arctic.bri.inf <- jags(data = jag.data,
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
# save(PDR.arctic.bri.inf, file = "R-scripts/R2jags-objects/PDR.arctic.bri.inf.Rdata")

# Read the .Rdata
load("R-scripts/R2jags-objects/PDR.arctic.bri.inf.Rdata")


## Diagnostics ----
##### Examine output
PDR.arctic.bri.inf$BUGSoutput$summary[1:5,]
mcmcplot(PDR.arctic.bri.inf)

# Extract the DIC for future model comparisons
PDR.arctic.bri.inf$BUGSoutput$DIC

## Plot data + fit ----
df.PDR.arctic.bri.inf <- data.frame(PDR.arctic.bri.inf$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.PDR.arctic.bri.inf)

##### Plot
plot.PDR.arctic.bri.inf <- df.PDR.arctic.bri.inf %>% 
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

plot.PDR.arctic.bri.inf

# ggsave("figures/PDR.arctic.bri.inf.png", plot.PDR.arctic.bri.inf, 
#        width = 10.3, height = 5.6)


##########
###### 2E. Plot all three TPCs in the same graph (for comparison): Briere ----
##########

# Add an identifying column in each model output dataframe
df.PDR.arctic.bri.uni <- df.PDR.arctic.bri.uni %>% 
  mutate(type = "Arctic uniform")

df.PDR.nonarctic.bri.uni <- df.PDR.nonarctic.bri.uni %>% 
  mutate(type = "non-Arctic uniform")

df.PDR.arctic.bri.inf <- df.PDR.arctic.bri.inf %>% 
  mutate(type = "Arctic informative")

# Combine the three dataframes
df.all <- rbind(df.PDR.arctic.bri.uni, df.PDR.nonarctic.bri.uni, df.PDR.arctic.bri.inf)

# Plot
plot.all <- df.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = mean, color = type), linewidth = 1) +
  geom_point(data = data.PDR.arctic, aes(x = temp, y = trait), size = 2) +
  #geom_point(data = data.PDR.nonarctic, aes(x = temp, y = trait), size = 2) +
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
                               "non-Arctic uniform" = "grey",
                               "Arctic informative" = "pink")) +
  ## line
  scale_color_manual(values = c("Arctic uniform" = "blue", 
                                "non-Arctic uniform" = "#868686FF",
                                "Arctic informative" = "red")) +
  theme_bw()

plot.all

#ggsave("figures/PDR.all.bri.png", plot.all, 
#        width = 10.3, height = 5.6)



##########
###### 3A. Fit PDR thermal responses with uniform priors (Arctic): Quadratic ----
##########

##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)

##### Set data
data <- data.PDR.arctic

##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# ##### Run JAGS -----
# PDR.arctic.quad.uni <- jags(data = jag.data,
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
#save(PDR.arctic.quad.uni, file = "R-scripts/R2jags-objects/PDR.arctic.quad.uni.Rdata")

# Read the .Rdata
load("R-scripts/R2jags-objects/PDR.arctic.quad.uni.Rdata")


## Diagnostics ----
##### Examine output
PDR.arctic.quad.uni$BUGSoutput$summary[1:5,]
mcmcplot(PDR.arctic.quad.uni)

# Extract the DIC for future model comparisons
PDR.arctic.quad.uni$BUGSoutput$DIC

## Plot data + fit ----
df.PDR.arctic.quad.uni <- data.frame(PDR.arctic.quad.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.PDR.arctic.quad.uni)

##### Plot
plot.PDR.arctic.quad.uni <- df.PDR.arctic.quad.uni %>% 
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

plot.PDR.arctic.quad.uni

# ggsave("figures/PDR.arctic.quad.uni.png", plot.PDR.arctic.quad.uni, 
#        width = 10.3, height = 5.6)

##########
###### 3B. Fit PDR thermal responses for priors (Ae. sierrensis): Quadratic ----
##########

##### Temp sequence for derived quantity calculations
# For priors - fewer temps for derived calculations makes it go faster
Temp.xs <- seq(0, 45, 0.5)
N.Temp.xs <-length(Temp.xs)


##### Set data
data <- data.PDR.sierrensis

##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS -----

# This code took an hour to run!
# PDR.sierrensis.quad.uni <- jags(data = jag.data, 
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
# save(PDR.sierrensis.quad.uni, file = "R-scripts/R2jags-objects/PDR.sierrensis.quad.uni.Rdata")

# Read the .Rdata
load("R-scripts/R2jags-objects/PDR.sierrensis.quad.uni.Rdata")


## Diagnostics ----
##### Examine output
PDR.sierrensis.quad.uni$BUGSoutput$summary[1:5,]
mcmcplot(PDR.sierrensis.quad.uni)

# Extract the DIC for future model comparisons
PDR.sierrensis.quad.uni$BUGSoutput$DIC

## Plot data + fit ----
df.PDR.sierrensis.quad.uni <- data.frame(PDR.sierrensis.quad.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.PDR.sierrensis.quad.uni)

##### Plot
plot.df.PDR.sierrensis.quad.uni <- df.PDR.sierrensis.quad.uni %>% 
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

plot.df.PDR.sierrensis.quad.uni

# ggsave("figures/PDR.sierrensis.quad.uni.png", plot.df.PDR.sierrensis.quad.uni, 
#        width = 10.3, height = 5.6)


##########
###### 3C. Fit gamma distributions to PDR prior thermal responses: Quadratic ----
##########

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
PDR.arctic.prior.cf.dists <- data.frame(q = as.vector(PDR.sierrensis.quad.uni$BUGSoutput$sims.list$cf.q),
                                          T0 = as.vector(PDR.sierrensis.quad.uni$BUGSoutput$sims.list$cf.T0),
                                          Tm = as.vector(PDR.sierrensis.quad.uni$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
PDR.arctic.prior.gamma.fits = apply(PDR.arctic.prior.cf.dists, 2, 
                                      function(df) fitdistr(df, "gamma")$estimate)


PDR.hypers <- PDR.arctic.prior.gamma.fits
save(PDR.hypers, file = "R-scripts/R2jags-objects/PDRhypers.quad.Rsave")


##########
###### 3D. Fit PDR thermal responses with data-informed priors (Arctic): Quadratic ----
##########

load("R-scripts/R2jags-objects/PDRhypers.quad.Rsave")
PDR.arctic.prior.gamma.fits <- PDR.hypers

##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)

##### Set data
data <- data.PDR.arctic
hypers <- PDR.arctic.prior.gamma.fits * 0.1

##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, 
                 N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS -----
PDR.arctic.quad.inf <- jags(data = jag.data,
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
#save(PDR.arctic.quad.inf, file = "R-scripts/R2jags-objects/PDR.arctic.quad.inf.Rdata")

# Read the .Rdata
load("R-scripts/R2jags-objects/PDR.arctic.quad.inf.Rdata")


## Diagnostics ----
##### Examine output
PDR.arctic.quad.inf$BUGSoutput$summary[1:5,]
mcmcplot(PDR.arctic.quad.inf)

# Extract the DIC for future model comparisons
PDR.arctic.quad.inf$BUGSoutput$DIC

## Plot data + fit ----
df.PDR.arctic.quad.inf <- data.frame(PDR.arctic.quad.inf$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.PDR.arctic.quad.inf)

##### Plot
plot.PDR.arctic.quad.inf <- df.PDR.arctic.quad.inf %>% 
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

plot.PDR.arctic.quad.inf

# ggsave("figures/PDR.arctic.quad.inf.png", plot.PDR.arctic.quad.inf, 
#        width = 10.3, height = 5.6)


##########
###### 3E. Plot all three TPCs in the same graph (for comparison) ----
##########

# Add an identifying column in each model output dataframe
df.PDR.arctic.quad.uni <- df.PDR.arctic.quad.uni %>% 
  mutate(type = "Arctic uniform")

df.PDR.sierrensis.quad.uni <- df.PDR.sierrensis.quad.uni %>% 
  mutate(type = "Ae. sierrensis uniform")

df.PDR.arctic.quad.inf <- df.PDR.arctic.quad.inf %>% 
  mutate(type = "Arctic informative")

# Combine the three dataframes
df.all <- rbind(df.PDR.arctic.quad.uni, df.PDR.sierrensis.quad.uni, df.PDR.arctic.quad.inf)

##### Plot
plot.all <- df.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = mean, color = type), linewidth = 1) +
  geom_point(data = data.PDR.arctic, aes(x = temp, y = trait), size = 2) +
  #geom_point(data = data.PDR.sierrensis, aes(x = temp, y = trait), size = 2) +
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
                               "Ae. sierrensis uniform" = "grey",
                               "Arctic informative" = "pink")) +
  ## line
  scale_color_manual(values = c("Arctic uniform" = "blue", 
                                "Ae. sierrensis uniform" = "#868686FF",
                                "Arctic informative" = "red")) +
  theme_bw()

plot.all

# ggsave("figures/PDR.all.quad.png", plot.all,
#        width = 10.3, height = 5.6)

##### Plot all arctic TPCs for comparison ----
# Add an identifying column in each model output dataframe
df.PDR.arctic.bri.uni <- df.PDR.arctic.bri.uni %>% 
  mutate(type = "Briere (uni)")

df.PDR.arctic.bri.inf <- df.PDR.arctic.bri.inf %>% 
  mutate(type = "Briere (inf)")

df.PDR.arctic.quad.uni <- df.PDR.arctic.quad.uni %>% 
  mutate(type = "Quadratic (uni)")

df.PDR.arctic.quad.inf <- df.PDR.arctic.quad.inf %>% 
  mutate(type = "Quadratic (inf)")

# Combine the three dataframes
df.all <- rbind(df.PDR.arctic.bri.uni, df.PDR.arctic.bri.inf, 
                df.PDR.arctic.quad.uni, df.PDR.arctic.quad.inf)

##### Plot
plot.all <- df.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = mean, color = type), linewidth = 1) +
  geom_point(data = data.PDR.arctic, aes(x = temp, y = trait), size = 2) +
  #geom_point(data = data.PDR.sierrensis, aes(x = temp, y = trait), size = 2) +
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

# ggsave("figures/PDR.all.arctic.png", plot.all,
#        width = 10.3, height = 5.6)