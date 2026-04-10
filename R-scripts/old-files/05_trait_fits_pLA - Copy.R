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
##        C. Fit gamma distributions to pLA prior thermal responses
##        D. Fit pLA thermal responses with data-informed priors (Arctic)
##        E. Plot all three TPCs in the same graph (for comparison)
##
##    3. Fitting TPC (Quadratic)
##        A. Fit pLA thermal responses with uniform priors (Arctic)
##        B. Fit pLA thermal responses for priors (non-Arctic species)
##        C. Fit gamma distributions to pLA prior thermal responses
##        D. Fit pLA thermal responses with data-informed priors (Arctic)
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
library(RColorBrewer) # colour palette


# Load functions
source("R-scripts/00_Functions.R")


# Load data
data <- read_csv("data-processed/TraitData_pLA.csv")
unique(data$species)


# Subset data
## Arctic species
data.pLA.arctic <- subset(data, type == "Arctic")

## Non-Arctic species
data.pLA.nonarctic <- subset(data, type == "non-Arctic")


# Plot the raw data
plot.data.pLA <- data %>% 
  ggplot(aes(x = temp, y = trait)) +
  geom_point(aes(colour = species, shape = citation)) +
  labs(y = "Larval survival (%)", x = expression(paste("Temperature (", degree, "C)"))) +
  scale_colour_discrete(name = "Species", labels = c("Ae. nigromaculis",
                                                     "Ae. sollicitans",
                                                     "Ae. triseriatus",
                                                     "Ae. vexans"
  )) +
  scale_shape_discrete(name = "Citation", labels = c("Brust 1967",
                                                     "Shelton 1973",
                                                     "Teng 2000",
                                                     "Trpis 1970")) +
  facet_grid(rows = vars(type)) +
  theme_bw()

plot.data.pLA

# ggsave("figures/raw_data/plot.data.pLA.png", plot.data.pLA, , width = 9.83, height = 6.17)



##########
###### 1. MCMC settings for all models ----
##########

# Number of posterior dist elements = [(ni - nb) / nt] * nc = [(45000 - 5000) / 8] * 3 = 15000
ni <- 45000 # number of iterations in each chain
nb <- 5000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 3 # number of chains


##########
###### 2A. Fit pLA thermal responses with uniform priors (Arctic): Briere ----
##########

##### inits Function
inits<-function(){list(
  cf.q = 0.01,
  cf.Tm = 35,
  cf.T0 = 5,
  cf.sigma = rlnorm(1))}

##### Parameters to Estimate
parameters <- c("cf.q", "cf.T0", "cf.Tm","cf.sigma", "z.trait.mu.pred")


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
  model.file = "R-scripts/briereprob.txt",
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
  labs(x = expression(paste("Temperature (", degree, "C)")), y = "Larval survival (%)") +
  theme_bw()

plot.pLA.arctic.bri.uni

# ggsave("figures/pLA.arctic.bri.uni.png", plot.pLA.arctic.bri.uni,
#        width = 10.3, height = 5.6)



##########
###### 2B. Fit pLA thermal responses (with random effects) for priors (non-Arctic species): Briere ----
##########

##### Temp sequence for derived quantity calculations
# For priors - fewer temps for derived calculations makes it go faster
Temp.xs <- seq(0, 45, 0.5)
N.Temp.xs <-length(Temp.xs)


##### Set data
data <- data.pLA.nonarctic

## Create a unique id for each species-study combination
data <- data %>% 
  group_by(species, citation) %>% 
  mutate(unique_id = cur_group_id())


##### inits Function
inits <- function(){list(
  cf.q = 0.001,
  cf.Tm = 35,
  cf.T0 = 5,
  cf.sigma = rlnorm(1),
  sigma_q = 0.001,
  sigma_T0 = rlnorm(1),
  sigma_Tm = rlnorm(1))}


##### Parameters to Estimate
parameters <- c("cf.q", "cf.T0", "cf.Tm", "cf.sigma", "sigma_q", "sigma_T0", 
                "sigma_Tm", "z.trait.mu.pred.pop", "z.trait.mu.pred.id")


##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp
unique.id <- as.integer(data$unique_id)
Nids <- max(unique.id)

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, 
                 N.Temp.xs = N.Temp.xs, Nids = Nids, unique.id = unique.id,
                 prior = prior)

##### Run JAGS
pLA.nonarctic.bri.uni.raneff <- jags(
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


## Save the model as Rdata 
# save(pLA.nonarctic.bri.uni.raneff, file = "R-scripts/R2jags-objects/pLA.nonarctic.bri.uni.raneff.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/pLA.nonarctic.bri.uni.raneff.Rdata")


## Diagnostics ----
##### Examine output
pLA.nonarctic.bri.uni.raneff$BUGSoutput$summary[1:8,]
mcmcplot(pLA.nonarctic.bri.uni.raneff)

# Extract the DIC for future model comparisons
pLA.nonarctic.bri.uni.raneff$BUGSoutput$DIC


## Plot data + fit ----
df.pLA.nonarctic.bri.uni.raneff <- data.frame(pLA.nonarctic.bri.uni.raneff$BUGSoutput$summary)[-(1:8),]

## Extract the model prediction
## Overall curve
df.pLA.nonarctic.bri.uni.raneff.pop <- df.pLA.nonarctic.bri.uni.raneff %>% 
  filter(grepl("z.trait.mu.pred.pop", rownames(df.pLA.nonarctic.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)


## Unique ID 1: Ae. nigromaculis
df.pLA.nonarctic.bri.uni.1 <- df.pLA.nonarctic.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[1,*]"), rownames(df.pLA.nonarctic.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 1)

## Unique ID 2: Ae. sollicitans
df.pLA.nonarctic.bri.uni.2 <- df.pLA.nonarctic.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[2,*]"), rownames(df.pLA.nonarctic.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 2)

## Unique ID 3: Ae. triseriatus (Shelton 1973)
df.pLA.nonarctic.bri.uni.3 <- df.pLA.nonarctic.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[3,*]"), rownames(df.pLA.nonarctic.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 3)


## Unique ID 4: Ae. triseriatus (Teng 2000)
df.pLA.nonarctic.bri.uni.4 <- df.pLA.nonarctic.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[3,*]"), rownames(df.pLA.nonarctic.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 4)

## Combine the model prediciton of all three unique groups into a dataframe
df.pLA.nonarctic.bri.uni.raneff.sp <- rbind(df.pLA.nonarctic.bri.uni.1,
                                            df.pLA.nonarctic.bri.uni.2,
                                            df.pLA.nonarctic.bri.uni.3,
                                            df.pLA.nonarctic.bri.uni.4) 

## Change unique_id into factor type
df.pLA.nonarctic.bri.uni.raneff.sp$unique_id <- as.factor(df.pLA.nonarctic.bri.uni.raneff.sp$unique_id)


##### Plot
plot.pLA.nonarctic.bri.uni.raneff <- ggplot(data = df.pLA.nonarctic.bri.uni.raneff.pop, 
                                            aes(x = temp)) +
  ## Overall TPC
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.),
              fill = "grey",
              alpha = 0.5) +
  ## a separate TPC (and credible interval) for each unique group
  geom_ribbon(data = df.pLA.nonarctic.bri.uni.raneff.sp, aes(ymin = X2.5., ymax = X97.5., fill = unique_id),
              alpha = 0.5) +
  geom_line(aes(y = mean), color = "black", linewidth = 1) +
  geom_line(data = df.pLA.nonarctic.bri.uni.raneff.sp, aes(y = mean, color = unique_id)) +
  geom_point(data = data,
             aes(x = temp, y = trait, colour = as.factor(unique_id)),
             size = 2) +
  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), y = "Larval survival (%)") +
  # Customize legend
  scale_colour_discrete(name = element_blank(),
                        labels = c("Ae. nigromaculis",
                                   "Ae. sollicitans",
                                   "Ae. triseriatus (Shelton 1973)",
                                   "Ae. triseriatus (Teng 2000)")) +
  theme_bw()


plot.pLA.nonarctic.bri.uni.raneff

# ggsave("figures/pLA.nonarctic.bri.uni.raneff.png", plot.pLA.nonarctic.bri.uni.raneff,
#        width = 10.3, height = 5.6)



##########
###### 2C. Fit gamma distributions to pLA prior thermal responses: Briere ----
##########

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
pLA.arctic.prior.cf.dists <- data.frame(q = as.vector(pLA.nonarctic.bri.uni.raneff$BUGSoutput$sims.list$cf.q),
                                        T0 = as.vector(pLA.nonarctic.bri.uni.raneff$BUGSoutput$sims.list$cf.T0),
                                        Tm = as.vector(pLA.nonarctic.bri.uni.raneff$BUGSoutput$sims.list$cf.Tm))

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


##### Set data
data <- data.pLA.arctic
hypers <- pLA.arctic.prior.gamma.fits * 0.1


##### inits Function
inits<-function(){list(
  cf.q = 0.01,
  cf.Tm = 35,
  cf.T0 = 5,
  cf.sigma = rlnorm(1))}

##### Parameters to Estimate
parameters <- c("cf.q", "cf.T0", "cf.Tm","cf.sigma", "z.trait.mu.pred")


##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)


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
                           model.file = "R-scripts/briereprob_inf.txt",
                           n.thin = nt,
                           n.chains = nc,
                           n.burnin = nb,
                           n.iter = ni,
                           DIC = T,
                           working.directory = getwd()
)

## Save the model as Rdata 
save(pLA.arctic.bri.inf, file = "R-scripts/R2jags-objects/pLA.arctic.bri.inf.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/pLA.arctic.bri.inf.Rdata")


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
    y = "Larval survival (%)"
  ) +
  theme_bw()

plot.pLA.arctic.bri.inf

# ggsave("figures/pLA.arctic.bri.inf.png", plot.pLA.arctic.bri.inf, 
#        width = 10.3, height = 5.6)




##########
###### 2E. Plot all TPCs for Arctic species in the same graph (for comparison): Briere ----
##########

# Add an identifying column in each model output dataframe
df.pLA.arctic.bri.uni <- df.pLA.arctic.bri.uni %>% 
  mutate(type = "Briere uniform")

df.pLA.arctic.bri.inf <- df.pLA.arctic.bri.inf %>% 
  mutate(type = "Briere informative")


# Combine the three dataframes
df.all <- rbind(df.pLA.arctic.bri.uni, df.pLA.arctic.bri.inf)

df.all$type <- factor(df.all$type, levels = c( "Briere uniform", "Briere informative"))


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
    y = "Larval survival (%)"
  ) +
  # Customize the colours
  ## ribbon
  scale_fill_manual(values = c("Briere uniform" = "grey",
                               "Briere informative" = "#4363d8")) +
  
  ## line
  scale_color_manual(values = c("Briere uniform" = "#868686FF",
                                "Briere informative" = "blue")) +
  theme_bw()

plot.all

#ggsave("figures/pLA.arctic.bri.all.png", plot.all, width = 10.3, height = 5.6)

pLA.arctic.bri.uni$BUGSoutput$DIC
pLA.arctic.bri.inf$BUGSoutput$DIC



##########
###### 3A. Fit pLA thermal responses with uniform priors (Arctic): Quadratic ----
##########

##### inits Function
inits<-function(){list(
  cf.q = 0.01,
  cf.Tm = 35,
  cf.T0 = 5,
  cf.sigma = rlnorm(1))}

##### Parameters to Estimate
parameters <- c("cf.q", "cf.T0", "cf.Tm","cf.sigma", "z.trait.mu.pred")


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
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, 
                 N.Temp.xs = N.Temp.xs)

# ##### Run JAGS -----
# pLA.arctic.quad.uni <- jags(data = jag.data,
#                               inits = inits,
#                               parameters.to.save = parameters,
#                               model.file = "R-scripts/quadprob.txt",
#                               n.chains = nc,
#                               n.burnin = nb,
#                               n.iter = ni,
#                               DIC = T,
#                               working.directory = getwd()
# )

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
    y = "Larval survival (%)"
  ) +
  theme_bw()

plot.pLA.arctic.quad.uni

# ggsave("figures/pLA.arctic.quad.uni.png", plot.pLA.arctic.quad.uni, 
#        width = 10.3, height = 5.6)


##########
###### 3B. Fit pLA thermal responses (with random effects) for priors (non-Arctic species): Quadratic ----
##########

##### Temp sequence for derived quantity calculations
# For priors - fewer temps for derived calculations makes it go faster
Temp.xs <- seq(0, 45, 0.5)
N.Temp.xs <-length(Temp.xs)


##### Set data
data <- data.pLA.nonarctic

## Create a unique id for each species-study combination
data <- data %>% 
  group_by(species, citation) %>% 
  mutate(unique_id = cur_group_id())


##### inits Function
inits <- function(){list(
  cf.q = 0.001,
  cf.Tm = 35,
  cf.T0 = 5,
  cf.sigma = rlnorm(1),
  sigma_q = 0.0001,
  sigma_T0 = rlnorm(1),
  sigma_Tm = rlnorm(1))}


##### Parameters to Estimate
parameters <- c("cf.q", "cf.T0", "cf.Tm", "cf.sigma", "sigma_q", "sigma_T0", 
                "sigma_Tm", "z.trait.mu.pred.pop", "z.trait.mu.pred.id")


##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp
unique.id <- as.integer(data$unique_id)
Nids <- max(unique.id)

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, Nids = Nids, unique.id = unique.id)

##### Run JAGS
# pLA.nonarctic.quad.uni.raneff <- jags(
#   data = jag.data,
#   inits = inits,
#   parameters.to.save = parameters,
#   model.file = "R-scripts/quadprob_randeff.txt",
#   n.thin = nt,
#   n.chains = nc,
#   n.burnin = nb,
#   n.iter = ni,
#   DIC = T,
#   working.directory = getwd()
# )


## Save the model as Rdata 
# save(pLA.nonarctic.quad.uni.raneff, file = "R-scripts/R2jags-objects/pLA.nonarctic.quad.uni.raneff.Rdata")

# Read the .Rdata
load("R-scripts/R2jags-objects/pLA.nonarctic.quad.uni.raneff.Rdata")


## Diagnostics ----
##### Examine output
pLA.nonarctic.quad.uni.raneff$BUGSoutput$summary[1:8,]
mcmcplot(pLA.nonarctic.quad.uni.raneff)

# Extract the DIC for future model comparisons
pLA.nonarctic.quad.uni.raneff$BUGSoutput$DIC


## Plot data + fit ----
df.pLA.nonarctic.quad.uni.raneff <- data.frame(pLA.nonarctic.quad.uni.raneff$BUGSoutput$summary)[-(1:8),]

## Extract the model prediction
## Overall curve
df.pLA.nonarctic.quad.uni.raneff.pop <- df.pLA.nonarctic.quad.uni.raneff %>% 
  filter(grepl("z.trait.mu.pred.pop", rownames(df.pLA.nonarctic.quad.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)


## Unique ID 1: Ae. nigromaculis
df.pLA.nonarctic.quad.uni.1 <- df.pLA.nonarctic.quad.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[1,*]"), rownames(df.pLA.nonarctic.quad.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 1)

## Unique ID 2: Ae. sollicitans
df.pLA.nonarctic.quad.uni.2 <- df.pLA.nonarctic.quad.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[2,*]"), rownames(df.pLA.nonarctic.quad.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 2)

## Unique ID 3: Ae. triseriatus (Shelton 1973)
df.pLA.nonarctic.quad.uni.3 <- df.pLA.nonarctic.quad.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[3,*]"), rownames(df.pLA.nonarctic.quad.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 3)

## Unique ID 4: Ae. triseriatus (Teng 2000)
df.pLA.nonarctic.quad.uni.4 <- df.pLA.nonarctic.quad.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[4,*]"), rownames(df.pLA.nonarctic.quad.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 4)

## Combine the model prediciton of all three unique groups into a dataframe
df.pLA.nonarctic.quad.uni.raneff.sp <- rbind(df.pLA.nonarctic.quad.uni.1,
                                             df.pLA.nonarctic.quad.uni.2,
                                             df.pLA.nonarctic.quad.uni.3,
                                             df.pLA.nonarctic.quad.uni.4
                                             ) 

## Change unique_id into factor type
df.pLA.nonarctic.quad.uni.raneff.sp$unique_id <- as.factor(df.pLA.nonarctic.quad.uni.raneff.sp$unique_id)


##### Plot
plot.pLA.nonarctic.quad.uni.raneff <- ggplot() +
  ## Overall TPC
  geom_ribbon(data = df.pLA.nonarctic.quad.uni.raneff.pop, 
              aes(x = temp, ymin = X2.5., ymax = X97.5.),
              fill = "grey",
              alpha = 0.5) +
  
  ## a separate TPC (and credible interval) for each unique group
  # geom_ribbon(data = df.pLA.nonarctic.quad.uni.raneff.sp, 
  #             aes(x = temp, ymin = X2.5., ymax = X97.5., fill = unique_id), alpha = 0.5) +
  geom_point(data = data,
             aes(x = temp, y = trait, colour = as.factor(unique_id)),
             size = 2) +
  geom_line(data = df.pLA.nonarctic.quad.uni.raneff.sp, 
            aes(x = temp, y = mean, color = unique_id)) +
  geom_line(data = df.pLA.nonarctic.quad.uni.raneff.pop,
            aes(x = temp, y = mean), color = "black", linewidth = 1.5) +
  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), y = "Larval survival (%)") +
  # Customize legend
  scale_colour_discrete(name = element_blank(),
                        labels = c("Ae. nigromaculis",
                                   "Ae. sollicitans",
                                   "Ae. triseriatus (Shelton 1973)",
                                   "Ae. triseriatus (Teng 2000)")) +
  theme_bw()


plot.pLA.nonarctic.quad.uni.raneff

# ggsave("figures/pLA.nonarctic.quad.uni.raneff.png", plot.pLA.nonarctic.quad.uni.raneff,
#        width = 10.3, height = 5.6)



##########
###### 3C. Fit gamma distributions to pLA prior thermal responses: Quadratic ----
##########

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
pLA.arctic.prior.cf.dists <- data.frame(q = as.vector(pLA.nonarctic.quad.uni.raneff$BUGSoutput$sims.list$cf.q),
                                        T0 = as.vector(pLA.nonarctic.quad.uni.raneff$BUGSoutput$sims.list$cf.T0),
                                        Tm = as.vector(pLA.nonarctic.quad.uni.raneff$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
pLA.arctic.prior.gamma.fits = apply(pLA.arctic.prior.cf.dists, 2, 
                                    function(df) fitdistr(df, "gamma")$estimate)


pLA.hypers <- pLA.arctic.prior.gamma.fits
# save(pLA.hypers, file = "R-scripts/R2jags-objects/pLAhypers.quad.Rsave")

q <-rgamma(1000000, shape = 4.311846, rate = 707.790282)
T0 <-rgamma(1000000, shape = 60.505941, rate = 5.343882)
Tm <-rgamma(1000000, shape = 119.62490, rate = 3.32988)

df <- data.frame(T0 = T0)

gamma.dist <- ggplot() +
  geom_histogram(data = pLA.arctic.prior.cf.dists, aes(x = T0, y = ..density..), 
                 bins = 50, fill = "grey", color = "black", alpha = 0.6) +
  stat_function(fun = dgamma, args = list(shape = 60.505941, rate = 5.343882), 
                color = "black", linewidth = 1.2) +
  labs(x = "Tmin",
       y = "Density") +
  theme_bw()

gamma.dist

# ggsave("figures/pLA.T0.gamma.dist.png", gamma.dist,
#        width = 10.3, height = 5.6)


mean(T0)
mean(Tm)
plot(density(q))
plot(density(T0))
plot(density(Tm))

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
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "#4363d8", alpha = 0.5) +
  geom_line(aes(y = mean), color = "blue", linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait, colour = species), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Larval survival (%)"
  ) +
  scale_colour_manual(name = "Species", labels = "Ae. vexans", values = "black") +
  theme_bw()

plot.pLA.arctic.quad.inf

# ggsave("figures/pLA.arctic.quad.inf.png", plot.pLA.arctic.quad.inf, 
#        width = 10.3, height = 5.6)




##########
###### 3E. Plot all three TPCs in the same graph (for comparison) ----
##########

# Add an identifying column in each model output dataframe
df.pLA.arctic.quad.uni <- df.pLA.arctic.quad.uni %>% 
  mutate(type = "Quadratic uniform")


df.pLA.arctic.quad.inf <- df.pLA.arctic.quad.inf %>% 
  mutate(type = "Quadratic informative")


# Combine the three dataframes
df.all <- rbind(df.pLA.arctic.quad.uni, df.pLA.arctic.quad.inf)

##### Plot
plot.all <- df.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = mean, color = type), linewidth = 1) +
  geom_point(data = data.pLA.arctic, aes(x = temp, y = trait), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Larval survival (%)"
  ) +
  # Customize the colours
  ## ribbon
  scale_fill_manual(values = c("Quadratic uniform" = "grey", 
                               "Quadratic informative" = "#4363d8")) +
  ## line
  scale_color_manual(values = c("Quadratic uniform" = "#868686FF", 
                                "Quadratic informative" = "blue")) +
  theme_bw()

plot.all

# ggsave("figures/pLA.arctic.quad.all.png", plot.all, width = 10.3, height = 5.6)


##### Plot all best fitting TPCs for comparison ----


# Combine the three dataframes
df.all <- rbind(df.pLA.arctic.bri.uni, 
                df.pLA.arctic.bri.inf, 
                df.pLA.arctic.quad.uni,
                df.pLA.arctic.quad.inf
                )



##### Plot
plot.all <- df.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = mean, color = type), linewidth = 1) +
  geom_point(data = data.pLA.arctic, aes(x = temp, y = trait), size = 2) +
  #geom_point(data = data.pLA.sierrensis, aes(x = temp, y = trait), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Larval survival (%)"
  ) +
  # Customize the colours
  # scale_fill_jco() +
  # scale_color_jco() +
  # scale_fill_brewer(palette = "Accent") +
  # scale_color_brewer(palette = "Accent") +
  theme_bw()

plot.all

# ggsave("figures/pLA.arctic.all.png", plot.all, width = 10.3, height = 5.6)

#### DIC ----
pLA.arctic.bri.uni$BUGSoutput$DIC
pLA.arctic.bri.inf$BUGSoutput$DIC
pLA.arctic.quad.uni$BUGSoutput$DIC
pLA.arctic.quad.inf$BUGSoutput$DIC # This is the best fitting TPC


##########
###### 4. Process and save model output for plotting ----
##########

## Analyze TPC model
pLA.TPC.analysis <- extractTPC(pLA.arctic.quad.inf, "pLA", Temp.xs)
pLA.predictions.summary <- pLA.TPC.analysis[[1]]
pLA.params.summary <- pLA.TPC.analysis[[2]]
pLA.params.fullposts <- pLA.TPC.analysis[[3]]

write_csv(pLA.predictions.summary, "data-processed/pLA.predictions.summary.csv")
write_csv(pLA.params.summary, "data-processed/pLA.params.summary.csv")
write_csv(pLA.params.fullposts, "data-processed/pLA.params.fullposts.csv")
