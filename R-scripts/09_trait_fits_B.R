## Lilian Chan, University of Guelph
## Arctic vector-borne disease transmission suitability model
##
## Purpose: use Bayesian inference (JAGS) to fit TPCs for lifetime egg 
## production (B). Since  I only have Eggs per female per day (EFD) data for these 
## species, I will fit the TPC for EFD and adult lifespan to obtain 
## lifetime egg production. Then, I will combine the TPC fits and data from Arctic 
## species (Aedes hexodontus (Barlow 1955), Aedes cinereus, Aedes 
## communis, Aedes impiger, and Aedes punctor (Sommerman 1969)) and fit another 
## TPC.
##
##
## Table of content:
##    0. Set-up workspace
##
##    1. MCMC settings for all models
##
##    2. Fit EFD thermal responses (with random effects)
##        A. Briere
##        B. Quadratic
##        C. Compare the two TPC fits
##
##    3. Fit lf thermal responses (with random effects) 
##        A. Briere
##        B. Quadratic
##        C. Compare the two TPC fits
##
##    4. Calculate B and fit TPC
##       A. Calculate B for non-Arctic species
##       B. Combine data from non-Arctic and Arctic species
##       C. Fit thermal responses to B (all data): Briere
##       D. Fit thermal responses to B (all data): Quadratic



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

##### Load functions
source("R-scripts/00_Functions.R")

# Load data
data <- read_csv("data-processed/TraitData_B.csv")
unique(data$species)


# Subset data
## Arctic species
data.B.arctic <- subset(data, !(species %in% c("aegypti")))

## Non-Arctic species
data.B.nonarctic <- subset(data, species %in% "aegypti")


# Plot the raw data
plot.data.B <- data %>% 
  mutate(type = c(rep("Arctic",6), rep("non-Arctic: EFD", 30), rep("non-Arctic: lf", 907))) %>% 
  ggplot() +
  geom_point(aes(x = temp, y = trait, colour = species)) +
  labs(y = "eggs or adult lifespan", x = "Temperature ºC") +
  scale_colour_discrete(name = "species", labels = c("Ae. aegypti",
                                                     "Ae. cinereus",
                                                     "Ae. communis",
                                                     "Ae. hexodontus",
                                                     "Ae. impiger",
                                                     "Ae. punctor"
  )) +
  facet_grid(rows = vars(type), scales = "free_y") +
  theme_bw()


plot.data.B

# ggsave("figures/raw_data/plot.data.fecundity.png", plot.data.B, , width = 9.83, height = 6.17)

plot.data.B.nonarctic <- data.B.nonarctic %>% 
  mutate(type = c(rep("non-Arctic: EFD", 30), rep("non-Arctic: lf", 907))) %>% 
  ggplot() +
  geom_point(aes(x = temp, y = trait, colour = citation)) +
  labs(y = "eggs per days or adult lifespan", x = "Temperature ºC") +
  # scale_colour_discrete(name = "species", labels = c("Ae. aegypti"
  # )) +
  facet_grid(rows = vars(type), scales = "free_y") +
  scale_colour_discrete(name = "citation", labels = c()) +
  theme_bw()


plot.data.B.nonarctic

# ggsave("figures/raw_data/plot.data.B.nonarctic.png", plot.data.B.nonarctic, width = 9.83, height = 6.17)


##########
###### 1. MCMC settings for all models ----
##########

# Number of posterior dist elements = [(ni - nb) / nt] * nc = [(45000 - 5000) / 8] * 3 = 15000
ni <- 45000 # number of iterations in each chain
nb <- 5000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 3 # number of chains



##########
###### 2A. Fit EFD thermal responses (with random effects): Briere ----
##########

##### Temp sequence for derived quantity calculations
# Fewer temps for derived calculations makes it go faster
Temp.xs <- seq(0, 45, 0.5)
N.Temp.xs <-length(Temp.xs)

##### Set data
data <- data.B.nonarctic %>% 
  filter(trait_name == "EFD")


## Create a unique id for each species-study combination
data <- data %>% 
  group_by(species, citation) %>% 
  mutate(unique_id = cur_group_id())

## Set priors
prior <- data.frame(q = c(0, 0.1),
                    T0 = c(0, 20),
                    Tm = c(30, 45),
                    sigma_q = c(0, 0.0001),
                    sigma_T0 = c(0, 10),
                    sigma_Tm = c(0, 10)
)

##### inits Function
inits <- function(){list(
  cf.q = 0.01,
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
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, 
                 N.Temp.xs = N.Temp.xs, Nids = Nids, unique.id = unique.id,
                 prior = prior)

##### Run JAGS
B.EFD.nonarctic.bri.uni.raneff <- jags(
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
# save(B.EFD.nonarctic.bri.uni.raneff, file = "R-scripts/R2jags-objects/B.EFD.nonarctic.bri.uni.raneff.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/B.EFD.nonarctic.bri.uni.raneff.Rdata")


## Diagnostics ----
##### Examine output
B.EFD.nonarctic.bri.uni.raneff$BUGSoutput$summary[1:8,]
mcmcplot(B.EFD.nonarctic.bri.uni.raneff)

# Extract the DIC for future model comparisons
B.EFD.nonarctic.bri.uni.raneff$BUGSoutput$DIC


## Plot data + fit ----
df.B.EFD.nonarctic.bri.uni.raneff <- data.frame(B.EFD.nonarctic.bri.uni.raneff$BUGSoutput$summary)[-(1:8),]

## Extract the model prediction
## Overall curve
df.B.EFD.nonarctic.bri.uni.raneff.pop <- df.B.EFD.nonarctic.bri.uni.raneff %>% 
  filter(grepl("z.trait.mu.pred.pop", rownames(df.B.EFD.nonarctic.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)


## Unique ID 1: Ae. aegypti (Beserra 2009)
df.B.EFD.nonarctic.bri.uni.raneff.1 <- df.B.EFD.nonarctic.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[1,*]"), rownames(df.B.EFD.nonarctic.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 1)

## Unique ID 2: Ae. aegypti (Yang et al 2008)
df.B.EFD.nonarctic.bri.uni.raneff.2 <- df.B.EFD.nonarctic.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[2,*]"), rownames(df.B.EFD.nonarctic.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 2)


## Combine the model prediciton of all three unique groups into a dataframe
df.B.EFD.nonarctic.bri.uni.raneff.sp <- rbind(df.B.EFD.nonarctic.bri.uni.raneff.1,
                                            df.B.EFD.nonarctic.bri.uni.raneff.2) 

## Change unique_id into factor type
df.B.EFD.nonarctic.bri.uni.raneff.sp$unique_id <- as.factor(df.B.EFD.nonarctic.bri.uni.raneff.sp$unique_id)


##### Plot
plot.B.EFD.nonarctic.bri.uni.raneff <- ggplot(data = df.B.EFD.nonarctic.bri.uni.raneff.pop, 
                                        aes(x = temp)) +
  ## Overall TPC
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.),
              fill = "grey",
              alpha = 0.5) +
  ## a separate TPC (and credible interval) for each unique group
  # geom_ribbon(data = df.B.EFD.nonarctic.bri.uni.raneff.sp, aes(ymin = X2.5., ymax = X97.5., fill = unique_id),
  #             alpha = 0.5) +
  geom_line(data = df.B.EFD.nonarctic.bri.uni.raneff.sp, aes(y = mean, color = unique_id)) +
  geom_line(aes(y = mean), color = "black", linewidth = 1.5) +
  geom_point(data = data,
             aes(x = temp, y = trait, colour = as.factor(unique_id)),
             size = 2) +
  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), y = "Eggs per female per day") +
  # Customize legend
  scale_colour_discrete(name = element_blank(),
                        labels = c("Ae. aegypti (Beserra 2009)",
                                   "Ae. aegypti (Yang et al 2008)")) +
  theme_bw()


plot.B.EFD.nonarctic.bri.uni.raneff

# ggsave("figures/B.EFD.nonarctic.bri.uni.raneff.png", plot.B.EFD.nonarctic.bri.uni.raneff,
#        width = 10.3, height = 5.6)


##########
###### 2B. Fit EFD thermal responses (with random effects): Quadratic ----
##########

##### Temp sequence for derived quantity calculations
# Fewer temps for derived calculations makes it go faster
Temp.xs <- seq(0, 45, 0.5)
N.Temp.xs <-length(Temp.xs)

##### Set data
data <- data.B.nonarctic %>% 
  filter(trait_name == "EFD")


## Create a unique id for each species-study combination
data <- data %>% 
  group_by(species, citation) %>% 
  mutate(unique_id = cur_group_id())

## Set priors
prior <- data.frame(q = c(0, 0.1),
                    T0 = c(0, 20),
                    Tm = c(30, 45),
                    sigma_q = c(0, 0.0002),
                    sigma_T0 = c(0, 10),
                    sigma_Tm = c(0, 10)
)

##### inits Function
inits <- function(){list(
  cf.q = 0.1,
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
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, 
                 N.Temp.xs = N.Temp.xs, Nids = Nids, unique.id = unique.id,
                 prior = prior)

##### Run JAGS
B.EFD.nonarctic.quad.uni.raneff <- jags(
  data = jag.data,
  inits = inits,
  parameters.to.save = parameters,
  model.file = "R-scripts/quad_T_randeff.txt",
  n.thin = nt,
  n.chains = nc,
  n.burnin = nb,
  n.iter = ni,
  DIC = T,
  working.directory = getwd()
)


## Save the model as Rdata 
# save(B.EFD.nonarctic.quad.uni.raneff, file = "R-scripts/R2jags-objects/B.EFD.nonarctic.quad.uni.raneff.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/B.EFD.nonarctic.quad.uni.raneff.Rdata")


## Diagnostics ----
##### Examine output
B.EFD.nonarctic.quad.uni.raneff$BUGSoutput$summary[1:8,]
mcmcplot(B.EFD.nonarctic.quad.uni.raneff)

# Extract the DIC for future model comparisons
B.EFD.nonarctic.quad.uni.raneff$BUGSoutput$DIC


## Plot data + fit ----
df.B.EFD.nonarctic.quad.uni.raneff <- data.frame(B.EFD.nonarctic.quad.uni.raneff$BUGSoutput$summary)[-(1:8),]

## Extract the model prediction
## Overall curve
df.B.EFD.nonarctic.quad.uni.raneff.pop <- df.B.EFD.nonarctic.quad.uni.raneff %>% 
  filter(grepl("z.trait.mu.pred.pop", rownames(df.B.EFD.nonarctic.quad.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)


## Unique ID 1: Ae. aegypti (Beserra 2009)
df.B.EFD.nonarctic.quad.uni.raneff.1 <- df.B.EFD.nonarctic.quad.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[1,*]"), rownames(df.B.EFD.nonarctic.quad.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 1)

## Unique ID 2: Ae. aegypti (Yang et al 2008)
df.B.EFD.nonarctic.quad.uni.raneff.2 <- df.B.EFD.nonarctic.quad.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[2,*]"), rownames(df.B.EFD.nonarctic.quad.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 2)


## Combine the model prediciton of all three unique groups into a dataframe
df.B.EFD.nonarctic.quad.uni.raneff.sp <- rbind(df.B.EFD.nonarctic.quad.uni.raneff.1,
                                              df.B.EFD.nonarctic.quad.uni.raneff.2) 

## Change unique_id into factor type
df.B.EFD.nonarctic.quad.uni.raneff.sp$unique_id <- as.factor(df.B.EFD.nonarctic.quad.uni.raneff.sp$unique_id)


##### Plot
plot.B.EFD.nonarctic.quad.uni.raneff <- ggplot(data = df.B.EFD.nonarctic.quad.uni.raneff.pop, 
                                              aes(x = temp)) +
  ## Overall TPC
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.),
              fill = "grey",
              alpha = 0.5) +
  ## a separate TPC (and credible interval) for each unique group
  # geom_ribbon(data = df.B.EFD.nonarctic.quad.uni.raneff.sp, aes(ymin = X2.5., ymax = X97.5., fill = unique_id),
  #             alpha = 0.5) +
  geom_line(data = df.B.EFD.nonarctic.quad.uni.raneff.sp, aes(y = mean, color = unique_id)) +
  geom_line(aes(y = mean), color = "black", linewidth = 1.5) +
  geom_point(data = data,
             aes(x = temp, y = trait, colour = as.factor(unique_id)),
             size = 2) +
  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), y = "Eggs per female per day") +
  # Customize legend
  scale_colour_discrete(name = element_blank(),
                        labels = c("Ae. aegypti (Beserra 2009)",
                                   "Ae. aegypti (Yang et al 2008)")) +
  theme_bw()


plot.B.EFD.nonarctic.quad.uni.raneff

# ggsave("figures/B.EFD.nonarctic.quad.uni.raneff.png", plot.B.EFD.nonarctic.quad.uni.raneff,
#        width = 10.3, height = 5.6)


##########
###### 2C. Fit EFD thermal responses (NO random effects): Briere ----
##########

##### Temp sequence for derived quantity calculations
# Fewer temps for derived calculations makes it go faster
Temp.xs <- seq(0, 45, 0.5)
N.Temp.xs <-length(Temp.xs)


#### Set data
data <- data.B.nonarctic %>% 
  filter(trait_name == "EFD")


## Create a unique id for each species-study combination
data <- data %>% 
  group_by(species, citation) %>% 
  mutate(unique_id = cur_group_id())


## Set priors
prior <- data.frame(q = c(0, 0.01),
                    T0 = c(0, 20),
                    Tm = c(30, 45)
)


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
                 N.Temp.xs = N.Temp.xs, prior = prior)

##### Run JAGS -----
B.EFD.nonarctic.bri.uni <- jags(data = jag.data,
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
# save(B.EFD.nonarctic.bri.uni, file = "R-scripts/R2jags-objects/B.EFD.nonarctic.bri.uni.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/B.EFD.nonarctic.bri.uni.Rdata")

## Diagnostics ----
##### Examine output
B.EFD.nonarctic.bri.uni$BUGSoutput$summary[1:5,]
mcmcplot(B.EFD.nonarctic.bri.uni)

# Extract the DIC for future model comparisons
B.EFD.nonarctic.bri.uni$BUGSoutput$DIC

## Plot data + fit ----
df.B.EFD.nonarctic.bri.uni <- data.frame(B.EFD.nonarctic.bri.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.B.EFD.nonarctic.bri.uni)

##### Plot
plot.B.EFD.nonarctic.bri.uni <- df.B.EFD.nonarctic.bri.uni %>% 
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

plot.B.EFD.nonarctic.bri.uni

# ggsave("figures/B.EFD.nonarctic.bri.uni.png", plot.B.EFD.nonarctic.bri.uni, 
#        width = 10.3, height = 5.6)



##########
###### 2D. Fit EFD thermal responses (NO random effects): Quadratic ----
##########

##### Temp sequence for derived quantity calculations
# Fewer temps for derived calculations makes it go faster
Temp.xs <- seq(0, 45, 0.5)
N.Temp.xs <-length(Temp.xs)


#### Set data
data <- data.B.nonarctic %>% 
  filter(trait_name == "EFD")


## Create a unique id for each species-study combination
data <- data %>% 
  group_by(species, citation) %>% 
  mutate(unique_id = cur_group_id())


## Set priors
prior <- data.frame(q = c(0, 1),
                    T0 = c(0, 20),
                    Tm = c(30, 45)
)

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
                 N.Temp.xs = N.Temp.xs, prior = prior)

##### Run JAGS -----
B.EFD.nonarctic.quad.uni <- jags(data = jag.data,
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
# save(B.EFD.nonarctic.quad.uni, file = "R-scripts/R2jags-objects/B.EFD.nonarctic.quad.uni.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/B.EFD.nonarctic.quad.uni.Rdata")

## Diagnostics ----
##### Examine output
B.EFD.nonarctic.quad.uni$BUGSoutput$summary[1:5,]
mcmcplot(B.EFD.nonarctic.quad.uni)

# Extract the DIC for future model comparisons
B.EFD.nonarctic.quad.uni$BUGSoutput$DIC

## Plot data + fit ----
df.B.EFD.nonarctic.quad.uni <- data.frame(B.EFD.nonarctic.quad.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.B.EFD.nonarctic.quad.uni)

##### Plot
plot.B.EFD.nonarctic.quad.uni <- df.B.EFD.nonarctic.quad.uni %>% 
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

plot.B.EFD.nonarctic.quad.uni

# ggsave("figures/B.EFD.nonarctic.quad.uni.png", plot.B.EFD.nonarctic.quad.uni, 
#        width = 10.3, height = 5.6)



##########
###### 2E. Compare the TPC fits ----
##########

## DIC
B.EFD.nonarctic.bri.uni$BUGSoutput$DIC
B.EFD.nonarctic.bri.uni.raneff$BUGSoutput$DIC
B.EFD.nonarctic.quad.uni$BUGSoutput$DIC
B.EFD.nonarctic.quad.uni.raneff$BUGSoutput$DIC


df.B.EFD.nonarctic.bri.uni <- df.B.EFD.nonarctic.bri.uni %>% 
  mutate(type = "Briere")

df.B.EFD.nonarctic.bri.uni.raneff.pop <- df.B.EFD.nonarctic.bri.uni.raneff.pop %>% 
  mutate(type = "Briere w/ random effects")

df.B.EFD.nonarctic.quad.uni <- df.B.EFD.nonarctic.quad.uni %>% 
  mutate(type = "Quadratic")

df.B.EFD.nonarctic.quad.uni.raneff.pop <- df.B.EFD.nonarctic.quad.uni.raneff.pop %>% 
  mutate(type = "Quadratic w/ random effects")


# Combine the three dataframes
df.all <- rbind(df.B.EFD.nonarctic.bri.uni,
                df.B.EFD.nonarctic.bri.uni.raneff.pop, 
                df.B.EFD.nonarctic.quad.uni,
                df.B.EFD.nonarctic.quad.uni.raneff.pop
                )

# Plot
plot.all <- df.all %>% 
  ggplot(aes(x = temp)) +
  # geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = mean, color = type), linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait, colour = as.factor(unique_id)), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Eggs per female per day"
  ) +
  # Customize the colours
  ## ribbon
  # scale_fill_manual(name = element_blank(),
  #                   values = c("Briere" = "#4363d8",
  #                              "Quadratic" = "grey")) +

  ## line
  # scale_color_manual(name = element_blank(),
  #                    values = c("Briere" = "blue",
  #                               "Quadratic" = "#868686FF",
  #                               "1" = "#F8766D",
  #                               "2" = "#00BFC4"),
  #                    label = c("Beserra 2009", "Yang et al 2008", "Briere", "Quadratic")) +
  theme_bw()

plot.all

# ggsave("figures/B.EFD.nonarctic.all.png", plot.all, width = 10.3, height = 5.6)





##########
###### 3A. Fit lf thermal responses (with random effects): Briere ----
##########

##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.5)
N.Temp.xs <-length(Temp.xs)

##### Set data
data <- data.B.nonarctic %>% 
  filter(trait_name == "lf")


## Create a unique id for each species-study combination
data <- data %>% 
  group_by(species, citation) %>% 
  mutate(unique_id = cur_group_id())


## Set priors
prior <- data.frame(q = c(0, 0.01),
                    T0 = c(0, 20),
                    Tm = c(30, 45),
                    sigma_q = c(0, 0.0001),
                    sigma_T0 = c(0, 10),
                    sigma_Tm = c(0, 10)
)

##### inits Function
inits <- function(){list(
  cf.q = 0.01,
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
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, 
                 N.Temp.xs = N.Temp.xs, Nids = Nids, unique.id = unique.id,
                 prior = prior)

##### Run JAGS
# B.lf.nonarctic.bri.uni.raneff <- jags(
#   data = jag.data,
#   inits = inits,
#   parameters.to.save = parameters,
#   model.file = "R-scripts/briere_T_randeff_B.txt",
#   n.thin = nt,
#   n.chains = nc,
#   n.burnin = nb,
#   n.iter = ni,
#   DIC = T,
#   working.directory = getwd()
# )


## Save the model as Rdata 
# save(B.lf.nonarctic.bri.uni.raneff, file = "R-scripts/R2jags-objects/B.lf.nonarctic.bri.uni.raneff.Rdata")

# Read the .Rdata
load("R-scripts/R2jags-objects/B.lf.nonarctic.bri.uni.raneff.Rdata")


## Diagnostics ----
##### Examine output
B.lf.nonarctic.bri.uni.raneff$BUGSoutput$summary[1:8,]
mcmcplot(lf.nonarctic.bri.uni.raneff)

# Extract the DIC for future model comparisons
B.lf.nonarctic.bri.uni.raneff$BUGSoutput$DIC


## Plot data + fit ----
df.B.lf.nonarctic.bri.uni.raneff <- data.frame(B.lf.nonarctic.bri.uni.raneff$BUGSoutput$summary)[-(1:8),]

## Extract the model prediction
## Overall curve
df.B.lf.nonarctic.bri.uni.raneff.pop <- df.B.lf.nonarctic.bri.uni.raneff %>% 
  filter(grepl("z.trait.mu.pred.pop", rownames(df.B.lf.nonarctic.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)


## Unique ID 1: Ae. aegypti (Beserra 2009)
df.B.lf.nonarctic.bri.uni.raneff.1 <- df.B.lf.nonarctic.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[1,*]"), rownames(df.B.lf.nonarctic.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 1)

## Unique ID 2: Ae. aegypti (Goindin et al. 2015)
df.B.lf.nonarctic.bri.uni.raneff.2 <- df.B.lf.nonarctic.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[2,*]"), rownames(df.B.lf.nonarctic.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 2)

## Unique ID 3: Ae. aegypti (Huxley et al. 2021)
df.B.lf.nonarctic.bri.uni.raneff.3 <- df.B.lf.nonarctic.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[3,*]"), rownames(df.B.lf.nonarctic.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 3)


## Unique ID 4: Ae. aegypti (Huxley et al. 2022)
df.B.lf.nonarctic.bri.uni.raneff.4 <- df.B.lf.nonarctic.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[4,*]"), rownames(df.B.lf.nonarctic.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 4)


## Unique ID 5: Ae. aegypti (Marinho et al. 2016)
df.B.lf.nonarctic.bri.uni.raneff.5 <- df.B.lf.nonarctic.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[5,*]"), rownames(df.B.lf.nonarctic.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 5)


## Unique ID 6: Ae. aegypti (Rocha-Santos et al. 2021)
df.B.lf.nonarctic.bri.uni.raneff.6 <- df.B.lf.nonarctic.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[6,*]"), rownames(df.B.lf.nonarctic.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 6)


## Unique ID 7: Ae. aegypti (Yang et al. 2009)
df.B.lf.nonarctic.bri.uni.raneff.7 <- df.B.lf.nonarctic.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[7,*]"), rownames(df.B.lf.nonarctic.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 7)


## Combine the model prediciton of all three unique groups into a dataframe
df.B.lf.nonarctic.bri.uni.raneff.sp <- rbind(df.B.lf.nonarctic.bri.uni.raneff.1,
                                            df.B.lf.nonarctic.bri.uni.raneff.2,
                                            df.B.lf.nonarctic.bri.uni.raneff.3,
                                            df.B.lf.nonarctic.bri.uni.raneff.4,
                                            df.B.lf.nonarctic.bri.uni.raneff.5,
                                            df.B.lf.nonarctic.bri.uni.raneff.6,
                                            df.B.lf.nonarctic.bri.uni.raneff.7) 

## Change unique_id into factor type
df.B.lf.nonarctic.bri.uni.raneff.sp$unique_id <- as.factor(df.B.lf.nonarctic.bri.uni.raneff.sp$unique_id)


##### Plot
plot.B.lf.nonarctic.bri.uni.raneff <- ggplot(data = df.B.lf.nonarctic.bri.uni.raneff.pop, 
                                            aes(x = temp)) +
  ## Overall TPC
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.),
              fill = "grey",
              alpha = 0.5) +
  ## a separate TPC (and credible interval) for each unique group
  geom_point(data = data,
             aes(x = temp, y = trait, colour = as.factor(unique_id)),
             size = 2) +
  geom_line(data = df.B.lf.nonarctic.bri.uni.raneff.sp, aes(y = mean, color = unique_id)) +
  geom_line(aes(y = mean), color = "black", linewidth = 1.5) +
  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), y = "lifespan (days)") +
  # Customize legend
  scale_colour_discrete(name = element_blank(),
                        labels = c("Ae. aegypti (Beserra 2009)",
                                   "Ae. aegypti (Goindin et al. 2015)",
                                   "Ae. aegypti (Huxley et al. 2021)",
                                   "Ae. aegypti (Huxley et al. 2022)",
                                   "Ae. aegypti (Marinho et al. 2016)",
                                   "Ae. aegypti (Rocha-Santos et al. 2021)",
                                   "Ae. aegypti (Yang et al. 2009)")) +
  theme_bw()


plot.B.lf.nonarctic.bri.uni.raneff

# ggsave("figures/B.lf.nonarctic.bri.uni.raneff.png", plot.B.lf.nonarctic.bri.uni.raneff,
#        width = 10.3, height = 5.6)




##########
###### 3B. Fit lf thermal responses (with random effects): Quadratic ----
##########

##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.5)
N.Temp.xs <-length(Temp.xs)

##### Set data
data <- data.B.nonarctic %>% 
  filter(trait_name == "lf")


## Create a unique id for each species-study combination
data <- data %>% 
  group_by(species, citation) %>% 
  mutate(unique_id = cur_group_id())


## Set priors
prior <- data.frame(q = c(0, 0.1),
                    T0 = c(0, 20),
                    Tm = c(30, 45),
                    sigma_q = c(0, 0.0002),
                    sigma_T0 = c(0, 10),
                    sigma_Tm = c(0, 10)
)

##### inits Function
inits <- function(){list(
  cf.q = 0.1,
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
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, 
                 N.Temp.xs = N.Temp.xs, Nids = Nids, unique.id = unique.id,
                 prior = prior)

##### Run JAGS
B.lf.nonarctic.quad.uni.raneff <- jags(
  data = jag.data,
  inits = inits,
  parameters.to.save = parameters,
  model.file = "R-scripts/quad_T_randeff.txt",
  n.thin = nt,
  n.chains = nc,
  n.burnin = nb,
  n.iter = ni,
  DIC = T,
  working.directory = getwd()
)


## Save the model as Rdata 
# save(B.lf.nonarctic.quad.uni.raneff, file = "R-scripts/R2jags-objects/B.lf.nonarctic.quad.uni.raneff.Rdata")

# Read the .Rdata
load("R-scripts/R2jags-objects/B.lf.nonarctic.quad.uni.raneff.Rdata")


## Diagnostics ----
##### Examine output
B.lf.nonarctic.quad.uni.raneff$BUGSoutput$summary[1:8,]
mcmcplot(lf.nonarctic.quad.uni.raneff)

# Extract the DIC for future model comparisons
B.lf.nonarctic.quad.uni.raneff$BUGSoutput$DIC


## Plot data + fit ----
df.B.lf.nonarctic.quad.uni.raneff <- data.frame(B.lf.nonarctic.quad.uni.raneff$BUGSoutput$summary)[-(1:8),]

## Extract the model prediction
## Overall curve
df.B.lf.nonarctic.quad.uni.raneff.pop <- df.B.lf.nonarctic.quad.uni.raneff %>% 
  filter(grepl("z.trait.mu.pred.pop", rownames(df.B.lf.nonarctic.quad.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)


## Unique ID 1: Ae. aegypti (Beserra 2009)
df.B.lf.nonarctic.quad.uni.raneff.1 <- df.B.lf.nonarctic.quad.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[1,*]"), rownames(df.B.lf.nonarctic.quad.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 1)

## Unique ID 2: Ae. aegypti (Goindin et al. 2015)
df.B.lf.nonarctic.quad.uni.raneff.2 <- df.B.lf.nonarctic.quad.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[2,*]"), rownames(df.B.lf.nonarctic.quad.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 2)

## Unique ID 3: Ae. aegypti (Huxley et al. 2021)
df.B.lf.nonarctic.quad.uni.raneff.3 <- df.B.lf.nonarctic.quad.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[3,*]"), rownames(df.B.lf.nonarctic.quad.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 3)


## Unique ID 4: Ae. aegypti (Huxley et al. 2022)
df.B.lf.nonarctic.quad.uni.raneff.4 <- df.B.lf.nonarctic.quad.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[4,*]"), rownames(df.B.lf.nonarctic.quad.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 4)


## Unique ID 5: Ae. aegypti (Marinho et al. 2016)
df.B.lf.nonarctic.quad.uni.raneff.5 <- df.B.lf.nonarctic.quad.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[5,*]"), rownames(df.B.lf.nonarctic.quad.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 5)


## Unique ID 6: Ae. aegypti (Rocha-Santos et al. 2021)
df.B.lf.nonarctic.quad.uni.raneff.6 <- df.B.lf.nonarctic.quad.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[6,*]"), rownames(df.B.lf.nonarctic.quad.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 6)


## Unique ID 7: Ae. aegypti (Yang et al. 2009)
df.B.lf.nonarctic.quad.uni.raneff.7 <- df.B.lf.nonarctic.quad.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[7,*]"), rownames(df.B.lf.nonarctic.quad.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 7)


## Combine the model prediciton of all three unique groups into a dataframe
df.B.lf.nonarctic.quad.uni.raneff.sp <- rbind(df.B.lf.nonarctic.quad.uni.raneff.1,
                                             df.B.lf.nonarctic.quad.uni.raneff.2,
                                             df.B.lf.nonarctic.quad.uni.raneff.3,
                                             df.B.lf.nonarctic.quad.uni.raneff.4,
                                             df.B.lf.nonarctic.quad.uni.raneff.5,
                                             df.B.lf.nonarctic.quad.uni.raneff.6,
                                             df.B.lf.nonarctic.quad.uni.raneff.7) 

## Change unique_id into factor type
df.B.lf.nonarctic.quad.uni.raneff.sp$unique_id <- as.factor(df.B.lf.nonarctic.quad.uni.raneff.sp$unique_id)


##### Plot
plot.B.lf.nonarctic.quad.uni.raneff <- ggplot(data = df.B.lf.nonarctic.quad.uni.raneff.pop, 
                                             aes(x = temp)) +
  ## Overall TPC
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.),
              fill = "grey",
              alpha = 0.5) +
  ## a separate TPC (and credible interval) for each unique group
  # geom_ribbon(data = df.B.lf.nonarctic.quad.uni.raneff.sp, aes(ymin = X2.5., ymax = X97.5., fill = unique_id),
  #             alpha = 0.5) +
  geom_line(data = df.B.lf.nonarctic.quad.uni.raneff.sp, aes(y = mean, color = unique_id)) +
  geom_line(aes(y = mean), color = "black", linewidth = 1.5) +
  geom_point(data = data,
             aes(x = temp, y = trait, colour = as.factor(unique_id)),
             size = 2) +
  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), y = "lifespan (days)") +
  # Customize legend
  scale_colour_discrete(name = element_blank(),
                        labels = c("Ae. aegypti (Beserra 2009)",
                                   "Ae. aegypti (Goindin et al. 2015)",
                                   "Ae. aegypti (Huxley et al. 2021)",
                                   "Ae. aegypti (Huxley et al. 2022)",
                                   "Ae. aegypti (Marinho et al. 2016)",
                                   "Ae. aegypti (Rocha-Santos et al. 2021)",
                                   "Ae. aegypti (Yang et al. 2009)")) +
  theme_bw()


plot.B.lf.nonarctic.quad.uni.raneff

# ggsave("figures/B.lf.nonarctic.quad.uni.raneff.png", plot.B.lf.nonarctic.quad.uni.raneff,
#        width = 10.3, height = 5.6)




##########
###### 3C. Compare TPC fits ----
##########

## DIC
B.lf.nonarctic.bri.uni.raneff$BUGSoutput$DIC
B.lf.nonarctic.quad.uni.raneff$BUGSoutput$DIC


df.B.lf.nonarctic.bri.uni.raneff.pop <- df.B.lf.nonarctic.bri.uni.raneff.pop %>% 
  mutate(type = "Briere w/ random effects")

df.B.lf.nonarctic.quad.uni.raneff.pop <- df.B.lf.nonarctic.quad.uni.raneff.pop %>% 
  mutate(type = "Quadratic w/ random effects")


# Combine the three dataframes
df.all <- rbind(df.B.lf.nonarctic.bri.uni.raneff.pop,
                df.B.lf.nonarctic.quad.uni.raneff.pop
)

# Plot
plot.all <- df.all %>% 
  ggplot(aes(x = temp)) +
  # geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_point(data = data, aes(x = temp, y = trait, colour = as.factor(unique_id)), size = 2) +
  geom_line(aes(y = mean, color = type), linewidth = 1) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "lifespan (days)"
  ) +
  # Customize the colours
  ## ribbon
  # scale_fill_manual(name = element_blank(),
  #                   values = c("Briere" = "#4363d8",
  #                              "Quadratic" = "grey")) +
  
  ## line
  # scale_color_manual(name = element_blank(),
  #                    values = c("Briere" = "blue",
  #                               "Quadratic" = "#868686FF",
  #                               "1" = "#F8766D",
  #                               "2" = "#00BFC4"),
  #                    label = c("Beserra 2009", "Yang et al 2008", "Briere", "Quadratic")) +
  theme_bw()

plot.all

# ggsave("figures/B.lf.nonarctic.all.png", plot.all, width = 10.3, height = 5.6)


##########
###### 4A. Calculate B for non-Arctic species ----
##########

## Load the models
load("R-scripts/R2jags-objects/B.EFD.nonarctic.bri.uni.Rdata")
load("R-scripts/R2jags-objects/B.lf.nonarctic.quad.uni.raneff.Rdata")

## Pull out the derived/predicted values:
B.EFD.pred <- B.EFD.nonarctic.bri.uni$BUGSoutput$sims.list$z.trait.mu.pred
B.lf.pred <- B.lf.nonarctic.quad.uni.raneff$BUGSoutput$sims.list$z.trait.mu.pred.pop


## Calculate B for Ae. aegypti
B.nonarctic.calc <- B.EFD.pred * B.lf.pred

##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.5)
N.Temp.xs <-length(Temp.xs)


## Get the mean, median, CI, and upper and lower quartile
df.B.nonarctic <- calcPostQuants(B.nonarctic.calc, Temp.xs)


## plot 
##### Plot
plot.B.nonarctic <- df.B.nonarctic %>%
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = lowerCI, ymax = upperCI),
              fill = "#4363d8",
              alpha = 0.5) +
  geom_line(aes(y = mean), color = "blue", linewidth = 1) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) +
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(x = expression(paste("Temperature (", degree, "C)")), y = "Lifetime egg production") +
  theme_bw()

plot.B.nonarctic

# ggsave("figures/B.nonarctic.png", plot.B.nonarctic,
#        width = 10.3, height = 5.6)


##########
###### 4B. Combine data from non-Arctic and Arctic species ----
##########

## Get the mean values for B
data <- df.B.nonarctic[, c("temp", "mean")]

colnames(data) <- c("temp", "trait")
data$species <- "aegypti"

## Combined with the arctic data
data <- bind_rows(data, data.B.arctic[, c("temp", "trait", "species")])


B.alldata <- ggplot(data = data) +
  geom_point(aes(x = temp, y = trait, colour = species)) +
  labs(y = "Lifetime eggs", x = "Temperature ºC") +
  scale_colour_discrete(name = "species", labels = c("Ae. aegypti",
                                                     "Ae. cinereus",
                                                     "Ae. communis",
                                                     "Ae. hexodontus",
                                                     "Ae. impiger",
                                                     "Ae. punctor"
  )) +
  theme_bw()

B.alldata

# ggsave("figures/raw_data/B.alldata.png", B.alldata, , width = 9.83, height = 6.17)


##########
###### 4C. Fit thermal responses to B (all data): Briere ----
##########

##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)


## Set priors
prior <- data.frame(q = c(0, 1),
                    T0 = c(0, 20),
                    Tm = c(30, 45)
)

##### inits Function
inits<-function(){list(
  cf.q = 0.1,
  cf.Tm = 35,
  cf.T0 = 5,
  cf.sigma = rlnorm(1))}


##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp


##### Parameters to Estimate
parameters <- c("cf.q", "cf.T0", "cf.Tm","cf.sigma", "z.trait.mu.pred")


##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, 
                 N.Temp.xs = N.Temp.xs, prior = prior)

##### Run JAGS
B.alldata.bri.uni <- jags(
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
# save(B.alldata.bri.uni, file = "R-scripts/R2jags-objects/B.alldata.bri.uni.Rdata")

# Read the .Rdata
load("R-scripts/R2jags-objects/B.alldata.bri.uni.Rdata")


## Diagnostics ----
##### Examine output
B.alldata.bri.uni$BUGSoutput$summary[1:5,]
mcmcplot(B.alldata.bri.uni)

# Extract the DIC for future model comparisons
B.alldata.bri.uni$BUGSoutput$DIC

## Plot data + fit ----
df.B.alldata.bri.uni <- data.frame(B.alldata.bri.uni$BUGSoutput$summary)[-(1:5),] %>%
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.B.alldata.bri.uni)

##### Plot
plot.B.alldata.bri.uni <- df.B.alldata.bri.uni%>%
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
  labs(x = expression(paste("Temperature (", degree, "C)")), y = "Lifetime egg production") +
  scale_colour_discrete(name = "species", labels = c("Ae. aegypti",
                                                     "Ae. cinereus",
                                                     "Ae. communis",
                                                     "Ae. hexodontus",
                                                     "Ae. impiger",
                                                     "Ae. punctor"
  )) +
  theme_bw()

plot.B.alldata.bri.uni

# ggsave("figures/B.alldata.bri.uni.png", plot.B.alldata.bri.uni,
#        width = 10.3, height = 5.6)



##########
###### 4D. Fit thermal responses to B (all data): Quadratic ----
##########

##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)


## Set priors
prior <- data.frame(q = c(0, 1),
                    T0 = c(0, 20),
                    Tm = c(30, 45)
)

##### inits Function
inits<-function(){list(
  cf.q = 0.1,
  cf.Tm = 35,
  cf.T0 = 5,
  cf.sigma = rlnorm(1))}


##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp


##### Parameters to Estimate
parameters <- c("cf.q", "cf.T0", "cf.Tm","cf.sigma", "z.trait.mu.pred")


##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, 
                 N.Temp.xs = N.Temp.xs, prior = prior)

##### Run JAGS
B.alldata.quad.uni <- jags(
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
# save(B.alldata.quad.uni, file = "R-scripts/R2jags-objects/B.alldata.quad.uni.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/B.alldata.quad.uni.Rdata")


## Diagnostics ----
##### Examine output
B.alldata.quad.uni$BUGSoutput$summary[1:5,]
mcmcplot(B.alldata.quad.uni)

# Extract the DIC for future model comparisons
B.alldata.quad.uni$BUGSoutput$DIC

## Plot data + fit ----
df.B.alldata.quad.uni <- data.frame(B.alldata.quad.uni$BUGSoutput$summary)[-(1:5),] %>%
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.B.alldata.quad.uni)

##### Plot
plot.B.alldata.quad.uni <- df.B.alldata.quad.uni%>%
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
  labs(x = expression(paste("Temperature (", degree, "C)")), y = "Lifetime egg production") +
  scale_colour_discrete(name = "species", labels = c("Ae. aegypti",
                                                     "Ae. cinereus",
                                                     "Ae. communis",
                                                     "Ae. hexodontus",
                                                     "Ae. impiger",
                                                     "Ae. punctor"
  )) +
  theme_bw()

plot.B.alldata.quad.uni

# ggsave("figures/B.alldata.quad.uni.png", plot.B.alldata.quad.uni,
#        width = 10.3, height = 5.6)

##########
###### 4E. Compare TPC fits ----
##########

## DIC
B.alldata.bri.uni$BUGSoutput$DIC
B.alldata.quad.uni$BUGSoutput$DIC


df.B.alldata.bri.uni <- df.B.alldata.bri.uni %>% 
  mutate(type = "Briere")

df.B.alldata.quad.uni <- df.B.alldata.quad.uni %>% 
  mutate(type = "Quadratic")


# Combine the three dataframes
df.all <- rbind(df.B.alldata.bri.uni,
                df.B.alldata.quad.uni
)

# Plot
plot.all <- df.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_point(data = data, aes(x = temp, y = trait), size = 2) +
  geom_line(aes(y = mean, color = type), linewidth = 1) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "lifespan (days)"
  ) +
  # Customize the colours
  ## ribbon
  scale_fill_manual(name = element_blank(),
                    values = c("Briere" = "#4363d8",
                               "Quadratic" = "grey")) +

  ## line
  scale_color_manual(name = element_blank(),
                     values = c("Briere" = "blue",
                                "Quadratic" = "#868686FF")) +
  theme_bw()

plot.all

# ggsave("figures/B.alldata.all.png", plot.all, width = 10.3, height = 5.6)
