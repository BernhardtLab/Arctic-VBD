## Lilian Chan, University of Guelph
## Arctic vector-borne disease transmission suitability model
##
## Purpose: Fit thermal performance curves (TPCs) for adult mosquito lifespan
## (lf) using Bayesian inference (JAGS).
## 
## Table of content:
##    0. Set-up workspace
##
##    1. MCMC settings for all models
##
##    2. Fitting TPC (Briere)
##        A. Fit Arctic TPC using uniform priors
##        B. Fit non-Arctic TPC to generate priors
##        C. Fit gamma distributions to non-Arctic TPC parameters
##        D. Fit Arctic TPC using data-informed priors
##
##    3. Fitting TPC (Quadratic)
##        A. Fit Arctic TPC using uniform priors
##        B. Fit non-Arctic TPC to generate priors
##        C. Fit gamma distributions to non-Arctic TPC parameters
##        D. Fit Arctic TPC using data-informed priors
##
##    4. Compare model fit between Briere and Quadratic models
##    5. Process and save model output for visualization
##
##
## Inputs:
## data-processed/TraitData_lf.csv - 
##     Synthesized published trait data for lf
##
## Outputs: 
## R-scripts/R2jags-objects/best-fitting-mods/lf.arctic.mod.Rdata - 
##     Best-fitting TPC models for Arctic species
##
## R-scripts/R2jags-objects/best-fitting-mods/lf.nonarctic.mod.Rdata -
##     Best-fitting TPC models for non-Arctic species
##
## data-processed/lf/lf.arctic.predictions.summary.csv -
##     Posterior summary of TPC predictions for Arctic species across temperatures
##
## data-processed/lf/lf.arctic.params.summary.csv -
##     Summary statistics of TPC parameters (Arctic TPC)
##
## data-processed/lf/lf.arctic.params.fullposts.csv -
##     Full posterior distributions for TPC parameters (Arctic TPC)
##
## data-processed/lf/lf.nonarctic.predictions.summary.csv -
##     Posterior summary of TPC predictions for non-Arctic species
##
## data-processed/lf/lf.nonarctic.params.summary.csv -
##     Summary statistics of TPC parameters (non-Arctic TPC)
##
## data-processed/lf/lf.nonarctic.params.fullposts.csv -
##     Full posterior distributions for TPC parameters (non-Arctic TPC)


# 0. Set-up workspace ----------------------------------------------------------

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
data.all <- read_csv("data-processed/TraitData_lf.csv")
unique(data.all$species)

# Subset data
## Arctic species
data.lf.arctic <- subset(data.all, type == "Arctic")

## Non-Arctic species
data.lf.nonarctic <- subset(data.all, type == "non-Arctic")


## Plot raw data
plot.data.lf <- data.all %>% 
  ggplot(aes(x = temp, y = trait, colour = species)) +
  geom_point(aes(colour = species)) +
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Time (days)"
       ) +
  scale_colour_discrete(name = "Species", labels = c("Ae. albopictus",
                                                     "Ae. cinereus", 
                                                     "Ae. communis",
                                                     "Ae. impiger",
                                                     "Ae. punctor",
                                                     "Ae. vexans"
  )) +
  facet_grid(rows = vars(type)) +
  theme_bw()

plot.data.lf



# 1. MCMC Settings for all models ----------------------------------------------

# Number of posterior dist elements = [(ni - nb) / nt] * nc = [(450000 - 50000) / 100] * 3 = 12000
ni <- 450000 # number of iterations in each chain
nb <- 50000 # number of 'burn in' iterations to discard
nt <- 100 # thinning rate - jags saves every nt iterations in each chain
nc <- 3 # number of chains


# 2. Fitting TPC (Briere) ------------------------------------------------------

## 2A. Fit Arctic TPC using uniform priors -------------------------------------

##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)

##### Set data
data <- data.lf.arctic

# Since this dataset has contains data from multiple species or multiple studies
# of the same species, we incorporated random effects on each thermal response
# parameter (q, T0, Tm) to addressed non-independence among observations 

## Create a unique id for each species-study combination
data <- data %>% 
  group_by(species, citation) %>% 
  mutate(unique_id = cur_group_id())


## Set priors
prior <- data.frame(q = c(0, 1),
                    T0 = c(0, 20),
                    Tm = c(21, 45),
                    sigma_q = c(0, 0.01),
                    sigma_T0 = c(0, 10),
                    sigma_Tm = c(0, 10)
)

##### inits Function
inits <- function(){list(
  cf.q = 0.1,
  cf.Tm = 35,
  cf.T0 = 5,
  cf.sigma = rlnorm(1),
  sigma_q = 0.001,
  sigma_T0 = rlnorm(1),
  sigma_Tm = rlnorm(1))}


##### Parameters to Estimate
parameters <- c("cf.q", "cf.T0", "cf.Tm", "cf.sigma", 
                "q", "T0", "Tm", "sigma_q", "sigma_T0", 
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
set.seed(123) # for reproducibility
lf.arctic.bri.uni <- jags(data = jag.data,
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
save(lf.arctic.bri.uni, file = "R-scripts/R2jags-objects/all-mods/lf.arctic.bri.uni.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/all-mods/lf.arctic.bri.uni.Rdata")


## Diagnostics
##### Examine output
lf.arctic.bri.uni$BUGSoutput$summary[c("cf.T0", "cf.Tm", "cf.q", "cf.sigma", "deviance", "sigma_T0", "sigma_Tm", "sigma_q"),]
mcmcplot(lf.arctic.bri.uni, parms = c("cf.T0", "cf.Tm", "cf.q", "cf.sigma", "deviance", "sigma_T0", "sigma_Tm", "sigma_q"))

# Extract the DIC for future model comparisons
lf.arctic.bri.uni$BUGSoutput$DIC


## Plot data + fit
df.lf.arctic.bri.uni <- data.frame(lf.arctic.bri.uni$BUGSoutput$summary)

## Extract the model prediction
## Overall curve
df.lf.arctic.bri.uni.pop <- df.lf.arctic.bri.uni %>% 
  filter(grepl("z.trait.mu.pred.pop", rownames(df.lf.arctic.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)


## Unique ID 1: Ae. cinereus
df.lf.arctic.bri.uni.1 <- df.lf.arctic.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[1,*]"), rownames(df.lf.arctic.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 1)

## Unique ID 2: Ae. communis
df.lf.arctic.bri.uni.2 <- df.lf.arctic.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[2,*]"), rownames(df.lf.arctic.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 2)

## Unique ID 3: Ae. impiger
df.lf.arctic.bri.uni.3 <- df.lf.arctic.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[3,*]"), rownames(df.lf.arctic.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 3)

## Unique ID 4: Ae. punctor
df.lf.arctic.bri.uni.4 <- df.lf.arctic.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[4,*]"), rownames(df.lf.arctic.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 4)

## Unique ID 5: Ae. vexans
df.lf.arctic.bri.uni.5 <- df.lf.arctic.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[5,*]"), rownames(df.lf.arctic.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 5)



## Combine the model prediciton of all three unique groups into a dataframe
df.lf.arctic.bri.uni.sp <- rbind(df.lf.arctic.bri.uni.1,
                                 df.lf.arctic.bri.uni.2,
                                 df.lf.arctic.bri.uni.3,
                                 df.lf.arctic.bri.uni.4,
                                 df.lf.arctic.bri.uni.5
                                 ) 

## Change unique_id into factor type
df.lf.arctic.bri.uni.sp$unique_id <- as.factor(df.lf.arctic.bri.uni.sp$unique_id)


##### Plot
plot.lf.arctic.bri.uni <- ggplot() +
  ## data
  geom_point(data = data,
             aes(x = temp, y = trait, colour = as.factor(unique_id)),
             size = 2) +
  
  ## a separate TPC for each unique group
  geom_line(data = df.lf.arctic.bri.uni.sp, 
            aes(x = temp, y = X50., color = unique_id)) +
  
  ## Overall TPC
  geom_ribbon(data = df.lf.arctic.bri.uni.pop,
              aes(x = temp, ymin = X2.5., ymax = X97.5.),
              fill = "grey",
              alpha = 0.5) +
  geom_line(data = df.lf.arctic.bri.uni.pop,
            aes(x = temp, y = X50.), color = "black", linewidth = 1) +
  
  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Time (days)",
       title = "A) lf, arctic sp., briere, uniform priors"
       ) +
  
  # Customize legend
  scale_colour_discrete(name = "",
                        labels = c("Ae. cinereus",
                                   "Ae. communis",
                                   "Ae. impiger",
                                   "Ae. punctor",
                                   "Ae. vexans")) +
  theme_bw()


plot.lf.arctic.bri.uni

ggsave("figures/lf.arctic.bri.uni.png", plot.lf.arctic.bri.uni,
       width = 10.3, height = 5.6)


## 2B. Fit non-Arctic TPC to generate priors -----------------------------------

##### Temp sequence for derived quantity calculations
# For priors - fewer temps for derived calculations makes it go faster
Temp.xs <- seq(0, 45, 0.5)
N.Temp.xs <-length(Temp.xs)


##### Set data
data <- data.lf.nonarctic

# Since this dataset has contains data from multiple species or multiple studies
# of the same species, we incorporated random effects on each thermal response
# parameter (q, T0, Tm) to addressed non-independence among observations 

## Create a unique id for each species-study combination
data <- data %>% 
  group_by(species, citation) %>% 
  mutate(unique_id = cur_group_id())


## Set priors
prior <- data.frame(q = c(0, 1),
                    T0 = c(0, 20),
                    Tm = c(21, 45),
                    sigma_q = c(0, 0.01),
                    sigma_T0 = c(0, 10),
                    sigma_Tm = c(0, 10)
)


##### inits Function
inits <- function(){list(
  cf.q = 0.01,
  cf.Tm = 35,
  cf.T0 = 5,
  cf.sigma = rlnorm(1))}


##### Parameters to Estimate
parameters <- c("cf.q", "cf.T0", "cf.Tm", "cf.sigma", 
                "q", "T0", "Tm", "sigma_q", "sigma_T0", 
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
set.seed(123) # for reproducibility
lf.nonarctic.bri.uni <- jags(data = jag.data,
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
save(lf.nonarctic.bri.uni, file = "R-scripts/R2jags-objects/all-mods/lf.nonarctic.bri.uni.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/all-mods/lf.nonarctic.bri.uni.Rdata")


## Diagnostics
##### Examine output
lf.nonarctic.bri.uni$BUGSoutput$summary[c("cf.T0", "cf.Tm", "cf.q", "cf.sigma", "deviance", "sigma_T0", "sigma_Tm", "sigma_q"),]
mcmcplot(lf.nonarctic.bri.uni, parms = c("cf.T0", "cf.Tm", "cf.q", "cf.sigma", "deviance", "sigma_T0", "sigma_Tm", "sigma_q"))

# Extract the DIC for future model comparisons
lf.nonarctic.bri.uni$BUGSoutput$DIC


## Plot data + fit
df.lf.nonarctic.bri.uni <- data.frame(lf.nonarctic.bri.uni$BUGSoutput$summary)

## Extract the model prediction
## Overall curve
df.lf.nonarctic.bri.uni.pop <- df.lf.nonarctic.bri.uni %>% 
  filter(grepl("z.trait.mu.pred.pop", rownames(df.lf.nonarctic.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)


## unique ID 1: Ae. albopictus (Alto 2001)
df.lf.nonarctic.bri.uni.1 <- df.lf.nonarctic.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[1,*]"), rownames(df.lf.nonarctic.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 1)

## unique ID 2: Ae. albopictus (Calado 2002)
df.lf.nonarctic.bri.uni.2 <- df.lf.nonarctic.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[2,*]"), rownames(df.lf.nonarctic.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 2)

## unique ID 3: Ae.albpictus (Ezeakacha 2015)
df.lf.nonarctic.bri.uni.3 <- df.lf.nonarctic.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[3,*]"), rownames(df.lf.nonarctic.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 3)

## unique ID 4: Ae. albopictus (Marini 2020)
df.lf.nonarctic.bri.uni.4 <- df.lf.nonarctic.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[4,*]"), rownames(df.lf.nonarctic.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 4)


## unique ID 5: Ae. albopictus (Tsuda 1994)
df.lf.nonarctic.bri.uni.5 <- df.lf.nonarctic.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[5,*]"), rownames(df.lf.nonarctic.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 5)



## Combine the model prediciton of all three unique groups into a dataframe
df.lf.nonarctic.bri.uni.sp <- rbind(df.lf.nonarctic.bri.uni.1,
                                    df.lf.nonarctic.bri.uni.2,
                                    df.lf.nonarctic.bri.uni.3,
                                    df.lf.nonarctic.bri.uni.4,
                                    df.lf.nonarctic.bri.uni.5
                                    ) 

## Change unique_id into factor type
df.lf.nonarctic.bri.uni.sp$unique_id <- as.factor(df.lf.nonarctic.bri.uni.sp$unique_id)



##### Plot
plot.lf.nonarctic.bri.uni <- ggplot() +
  ## data
  geom_point(data = data,
             aes(x = temp, y = trait, colour = as.factor(unique_id)),
             size = 2) +
  
  ## a separate TPC for each unique group
  geom_line(data = df.lf.nonarctic.bri.uni.sp, 
            aes(x = temp, y = X50., color = unique_id)) +
  
  ## Overall TPC
  geom_ribbon(data = df.lf.nonarctic.bri.uni.pop,
              aes(x = temp, ymin = X2.5., ymax = X97.5.),
              fill = "grey",
              alpha = 0.5) +
  geom_line(data = df.lf.nonarctic.bri.uni.pop,
            aes(x = temp, y = X50.), color = "black", linewidth = 1) +
  
  # Customize the axes and labels
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Time (days)"
  ) +
  
  # Customize legend
  scale_colour_discrete(name = "",
                        labels = c("Ae. albopictus 1",
                                   "Ae. albopictus 2",
                                   "Ae. albopictus 3",
                                   "Ae. albopictus 4",
                                   "Ae. albopictus 5")) +
  theme_bw()


plot.lf.nonarctic.bri.uni

ggsave("figures/lf.nonarctic.bri.uni.png", plot.lf.nonarctic.bri.uni,
       width = 10.3, height = 5.6)


## 2C. Fit gamma distributions to non-Arctic TPC parameters --------------------

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
lf.bri.prior.cf.dists <- data.frame(q = as.vector(lf.nonarctic.bri.uni$BUGSoutput$sims.list$cf.q),
                                    T0 = as.vector(lf.nonarctic.bri.uni$BUGSoutput$sims.list$cf.T0),
                                    Tm = as.vector(lf.nonarctic.bri.uni$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
lf.bri.prior.gamma.fits = apply(lf.bri.prior.cf.dists, 2, 
                                function(df) fitdistr(df, "gamma")$estimate)


save(lf.bri.prior.gamma.fits, file = "R-scripts/R2jags-objects/priors/lf.bri.priors.Rsave")


## 2D. Fit Arctic TPC using data-informed priors -------------------------------

# load("R-scripts/R2jags-objects/priors/lf.bri.priors.Rsave")


##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)


##### Set data
data <- data.lf.arctic
# Since this dataset has contains data from multiple species or multiple studies
# of the same species, we incorporated random effects on each thermal response
# parameter (q, T0, Tm) to addressed non-independence among observations 

## Create a unique id for each species-study combination
data <- data %>% 
  group_by(species, citation) %>% 
  mutate(unique_id = cur_group_id())

hypers <- lf.bri.prior.gamma.fits * 0.2
hypers[,2] <- lf.bri.prior.gamma.fits[,2]

# Gamma prior density for Tm

df.lf.Tm.prior <- tibble(
  temp = seq(0, 45, 0.1),
  density = dgamma(
    temp,
    shape = hypers["shape", "Tm"],
    rate  = hypers["rate", "Tm"]
  )
)

TPC.apex <- max(df.lf.nonarctic.bri.uni.pop$X50., na.rm = TRUE)

df.lf.Tm.prior <- df.lf.Tm.prior %>%
  mutate(
    density.scaled = density / max(density, na.rm = TRUE) * TPC.apex * 0.5
  )

df.lf.T0.prior <- tibble(
  temp = seq(0, 45, 0.1),
  density = dgamma(
    temp,
    shape = hypers["shape", "T0"],
    rate  = hypers["rate", "T0"]
  )
)

df.lf.T0.prior <- df.lf.T0.prior %>%
  mutate(
    max_density = max(density[is.finite(density)]), 
    density.scaled = density / max_density  * TPC.apex * 0.5
  ) %>% 
  dplyr::select(-max_density)


##### inits Function
inits<-function(){list(
  cf.q = 0.01,
  cf.Tm = 35,
  cf.T0 = 5,
  cf.sigma = rlnorm(1))}


##### Parameters to Estimate
parameters <- c("cf.q", "cf.T0", "cf.Tm", "cf.sigma", 
                "q", "T0", "Tm", "sigma_q", "sigma_T0", 
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
                 hypers = hypers)

##### Run JAGS
set.seed(123) # for reproducibility
lf.arctic.bri.inf <- jags(data = jag.data,
                          inits = inits,
                          parameters.to.save = parameters,
                          model.file = "R-scripts/briere_inf_raneff.txt",
                          n.thin = nt,
                          n.chains = nc,
                          n.burnin = nb,
                          n.iter = ni,
                          DIC = T,
                          working.directory = getwd()
                          )

## Save the model as Rdata 
save(lf.arctic.bri.inf, file = "R-scripts/R2jags-objects/all-mods/lf.arctic.bri.inf.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/all-mods/lf.arctic.bri.inf.Rdata")


## Diagnostics
##### Examine output
lf.nonarctic.bri.uni$BUGSoutput$summary[c("cf.T0", "cf.Tm", "cf.q", "cf.sigma", "deviance", "sigma_T0", "sigma_Tm", "sigma_q"),]
mcmcplot(lf.arctic.bri.inf, parms = c("cf.T0", "cf.Tm", "cf.q", "cf.sigma", "deviance", "sigma_T0", "sigma_Tm", "sigma_q"))

# Extract the DIC for future model comparisons
lf.arctic.bri.inf$BUGSoutput$DIC

## Plot data + fit
df.lf.arctic.bri.inf <- data.frame(lf.arctic.bri.inf$BUGSoutput$summary)

## Extract the model prediction
## Overall curve
df.lf.arctic.bri.inf.pop <- df.lf.arctic.bri.inf %>% 
  filter(grepl("z.trait.mu.pred.pop", rownames(df.lf.arctic.bri.inf))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)


## unique ID 1: Ae. cinereus
df.lf.arctic.bri.inf.1 <- df.lf.arctic.bri.inf %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[1,*]"), rownames(df.lf.arctic.bri.inf))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 1)

## unique ID 2: Ae. communis
df.lf.arctic.bri.inf.2 <- df.lf.arctic.bri.inf %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[2,*]"), rownames(df.lf.arctic.bri.inf))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 2)

## unique ID 3: Ae. impiger
df.lf.arctic.bri.inf.3 <- df.lf.arctic.bri.inf %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[3,*]"), rownames(df.lf.arctic.bri.inf))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 3)

## unique ID 4: Ae. punctor
df.lf.arctic.bri.inf.4 <- df.lf.arctic.bri.inf %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[4,*]"), rownames(df.lf.arctic.bri.inf))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 4)

## unique ID 5: Ae. vexans
df.lf.arctic.bri.inf.5 <- df.lf.arctic.bri.inf %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[5,*]"), rownames(df.lf.arctic.bri.inf))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 5)



## Combine the model prediciton of all three unique groups into a dataframe
df.lf.arctic.bri.inf.sp <- rbind(df.lf.arctic.bri.inf.1,
                                 df.lf.arctic.bri.inf.2,
                                 df.lf.arctic.bri.inf.3,
                                 df.lf.arctic.bri.inf.4,
                                 df.lf.arctic.bri.inf.5
                                 ) 

## Change unique_id into factor type
df.lf.arctic.bri.inf.sp$unique_id <- as.factor(df.lf.arctic.bri.inf.sp$unique_id)


##### Plot
plot.lf.arctic.bri.inf <- ggplot() +
  ## data
  geom_point(data = data,
             aes(x = temp, y = trait, colour = as.factor(unique_id)),
             size = 2) +
  
  # a separate TPC for each unique group
  geom_line(data = df.lf.arctic.bri.inf.sp,
            aes(x = temp, y = X50., color = unique_id)) +

  ## Overall TPC
  geom_ribbon(data = df.lf.arctic.bri.inf.pop,
              aes(x = temp, ymin = X2.5., ymax = X97.5.),
              fill = "grey",
              alpha = 0.5) +
  geom_line(data = df.lf.arctic.bri.inf.pop,
            aes(x = temp, y = X50.), color = "black", linewidth = 1) +

  geom_line(
    data = df.lf.Tm.prior,
    aes(x = temp, y = density.scaled),
    colour = "firebrick3",
    linetype = "dashed",
    linewidth = 1
  ) +
  
  geom_line(
    data = df.lf.T0.prior,
    aes(x = temp, y = density.scaled),
    colour = "skyblue",
    linetype = "dashed",
    linewidth = 1
  ) +
  
  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Time (days)",
       title = "C) lf, arctic sp. briere w/ gamma post"
       ) +
  
  # Customize legend
  scale_colour_discrete(name = "",
                        labels = c("Ae. cinereus",
                                   "Ae. communis",
                                   "Ae. impiger",
                                   "Ae. punctor",
                                   "Ae. vexans")) +
  theme_bw()


plot.lf.arctic.bri.inf

ggsave("figures/lf.arctic.bri.inf0.2.png", plot.lf.arctic.bri.inf,
       width = 10.3, height = 5.6)



# 3. Fitting TPC (quadratic) ---------------------------------------------------

## 3A. Fit Arctic TPC using uniform priors -------------------------------------

##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)

##### Set data
data <- data.lf.arctic

# Since this dataset has contains data from multiple species or multiple studies
# of the same species, we incorporated random effects on each thermal response
# parameter (q, T0, Tm) to addressed non-independence among observations 

## Create a unique id for each species-study combination
data <- data %>% 
  group_by(species, citation) %>% 
  mutate(unique_id = cur_group_id())


## Set priors
prior <- data.frame(q = c(0, 1),
                    T0 = c(0, 20),
                    Tm = c(21, 45),
                    sigma_q = c(0, 0.01),
                    sigma_T0 = c(0, 10),
                    sigma_Tm = c(0, 10)
)

##### inits Function
inits <- function(){list(
  cf.q = 0.1,
  cf.Tm = 35,
  cf.T0 = 5,
  cf.sigma = rlnorm(1),
  sigma_q = 0.001,
  sigma_T0 = rlnorm(1),
  sigma_Tm = rlnorm(1))}


##### Parameters to Estimate
parameters <- c("cf.q", "cf.T0", "cf.Tm", "cf.sigma", 
                "q", "T0", "Tm", "sigma_q", "sigma_T0", 
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
set.seed(123) # for reproducibility
lf.arctic.quad.uni <- jags(data = jag.data,
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
save(lf.arctic.quad.uni, file = "R-scripts/R2jags-objects/all-mods/lf.arctic.quad.uni.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/all-mods/lf.arctic.quad.uni.Rdata")


## Diagnostics
##### Examine output
lf.arctic.quad.uni$BUGSoutput$summary[c("cf.T0", "cf.Tm", "cf.q", "cf.sigma", "deviance", "sigma_T0", "sigma_Tm", "sigma_q"),]
mcmcplot(lf.arctic.quad.uni, parms = c("cf.T0", "cf.Tm", "cf.q", "cf.sigma", "deviance", "sigma_T0", "sigma_Tm", "sigma_q"))

# Extract the DIC for future model comparisons
lf.arctic.quad.uni$BUGSoutput$DIC


## Plot data + fit
df.lf.arctic.quad.uni <- data.frame(lf.arctic.quad.uni$BUGSoutput$summary)

## Extract the model prediction
## Overall curve
df.lf.arctic.quad.uni.pop <- df.lf.arctic.quad.uni %>% 
  filter(grepl("z.trait.mu.pred.pop", rownames(df.lf.arctic.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)


## Unique ID 1: Ae. cinereus
df.lf.arctic.quad.uni.1 <- df.lf.arctic.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[1,*]"), rownames(df.lf.arctic.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 1)

## Unique ID 2: Ae. communis
df.lf.arctic.quad.uni.2 <- df.lf.arctic.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[2,*]"), rownames(df.lf.arctic.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 2)

## Unique ID 3: Ae. impiger
df.lf.arctic.quad.uni.3 <- df.lf.arctic.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[3,*]"), rownames(df.lf.arctic.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 3)

## Unique ID 4: Ae. punctor
df.lf.arctic.quad.uni.4 <- df.lf.arctic.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[4,*]"), rownames(df.lf.arctic.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 4)

## Unique ID 5: Ae. vexans
df.lf.arctic.quad.uni.5 <- df.lf.arctic.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[5,*]"), rownames(df.lf.arctic.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 5)



## Combine the model prediciton of all three unique groups into a dataframe
df.lf.arctic.quad.uni.sp <- rbind(df.lf.arctic.quad.uni.1,
                                  df.lf.arctic.quad.uni.2,
                                  df.lf.arctic.quad.uni.3,
                                  df.lf.arctic.quad.uni.4,
                                  df.lf.arctic.quad.uni.5
                                  ) 

## Change unique_id into factor type
df.lf.arctic.quad.uni.sp$unique_id <- as.factor(df.lf.arctic.quad.uni.sp$unique_id)


##### Plot
plot.lf.arctic.quad.uni <- ggplot() +
  ## data
  geom_point(data = data,
             aes(x = temp, y = trait, colour = as.factor(unique_id)),
             size = 2) +
  
  ## a separate TPC for each unique group
  geom_line(data = df.lf.arctic.quad.uni.sp, 
            aes(x = temp, y = X50., color = unique_id)) +
  
  ## Overall TPC
  geom_ribbon(data = df.lf.arctic.quad.uni.pop,
              aes(x = temp, ymin = X2.5., ymax = X97.5.),
              fill = "grey",
              alpha = 0.5) +
  geom_line(data = df.lf.arctic.quad.uni.pop,
            aes(x = temp, y = X50.), color = "black", linewidth = 1) +
  
  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Time (days)",
       title = "B) lf, arctic sp., quadratic, uniform priors"
       ) +
  
  # Customize legend
  scale_colour_discrete(name = "",
                        labels = c("Ae. cinereus",
                                   "Ae. communis",
                                   "Ae. impiger",
                                   "Ae. punctor",
                                   "Ae. vexans")) +
  theme_bw()


plot.lf.arctic.quad.uni

ggsave("figures/lf.arctic.quad.uni.png", plot.lf.arctic.quad.uni,
       width = 10.3, height = 5.6)


## 3B. Fit non-Arctic TPC to generate priors -----------------------------------

##### Temp sequence for derived quantity calculations
# For priors - fewer temps for derived calculations makes it go faster
Temp.xs <- seq(0, 45, 0.5)
N.Temp.xs <-length(Temp.xs)


##### Set data
data <- data.lf.nonarctic

# Since this dataset has contains data from multiple species or multiple studies
# of the same species, we incorporated random effects on each thermal response
# parameter (q, T0, Tm) to addressed non-independence among observations 

## Create a unique id for each species-study combination
data <- data %>% 
  group_by(species, citation) %>% 
  mutate(unique_id = cur_group_id())


## Set priors
prior <- data.frame(q = c(0, 2),
                    T0 = c(0, 20),
                    Tm = c(21, 45),
                    sigma_q = c(0, 0.1),
                    sigma_T0 = c(0, 10),
                    sigma_Tm = c(0, 10)
)


##### inits Function
inits <- function(){list(
  cf.q = 0.01,
  cf.Tm = 35,
  cf.T0 = 5,
  cf.sigma = rlnorm(1))}


##### Parameters to Estimate
parameters <- c("cf.q", "cf.T0", "cf.Tm", "cf.sigma", 
                "q", "T0", "Tm", "sigma_q", "sigma_T0", 
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
set.seed(123) # for reproducibility
lf.nonarctic.quad.uni <- jags(data = jag.data,
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
save(lf.nonarctic.quad.uni, file = "R-scripts/R2jags-objects/all-mods/lf.nonarctic.quad.uni.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/all-mods/lf.nonarctic.quad.uni.Rdata")


## Diagnostics
##### Examine output
lf.nonarctic.quad.uni$BUGSoutput$summary[c("cf.T0", "cf.Tm", "cf.q", "cf.sigma", "deviance", "sigma_T0", "sigma_Tm", "sigma_q"),]
mcmcplot(lf.nonarctic.quad.uni, parms = c("cf.T0", "cf.Tm", "cf.q", "cf.sigma", "deviance", "sigma_T0", "sigma_Tm", "sigma_q"))

# Extract the DIC for future model comparisons
lf.nonarctic.quad.uni$BUGSoutput$DIC


## Plot data + fit
df.lf.nonarctic.quad.uni <- data.frame(lf.nonarctic.quad.uni$BUGSoutput$summary)

## Extract the model prediction
## Overall curve
df.lf.nonarctic.quad.uni.pop <- df.lf.nonarctic.quad.uni %>% 
  filter(grepl("z.trait.mu.pred.pop", rownames(df.lf.nonarctic.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)


## unique ID 1: Ae. albopictus (Alto 2001)
df.lf.nonarctic.quad.uni.1 <- df.lf.nonarctic.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[1,*]"), rownames(df.lf.nonarctic.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 1)

## unique ID 2: Ae. albopictus (Calado 2002)
df.lf.nonarctic.quad.uni.2 <- df.lf.nonarctic.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[2,*]"), rownames(df.lf.nonarctic.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 2)

## unique ID 3: Ae.albpictus (Ezeakacha 2015)
df.lf.nonarctic.quad.uni.3 <- df.lf.nonarctic.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[3,*]"), rownames(df.lf.nonarctic.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 3)

## unique ID 4: Ae. albopictus (Marini 2020)
df.lf.nonarctic.quad.uni.4 <- df.lf.nonarctic.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[4,*]"), rownames(df.lf.nonarctic.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 4)


## unique ID 5: Ae. albopictus (Tsuda 1994)
df.lf.nonarctic.quad.uni.5 <- df.lf.nonarctic.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[5,*]"), rownames(df.lf.nonarctic.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 5)



## Combine the model prediciton of all three unique groups into a dataframe
df.lf.nonarctic.quad.uni.sp <- rbind(df.lf.nonarctic.quad.uni.1,
                                    df.lf.nonarctic.quad.uni.2,
                                    df.lf.nonarctic.quad.uni.3,
                                    df.lf.nonarctic.quad.uni.4,
                                    df.lf.nonarctic.quad.uni.5
) 

## Change unique_id into factor type
df.lf.nonarctic.quad.uni.sp$unique_id <- as.factor(df.lf.nonarctic.quad.uni.sp$unique_id)



##### Plot
plot.lf.nonarctic.quad.uni <- ggplot() +
  ## data
  geom_point(data = data,
             aes(x = temp, y = trait, colour = as.factor(unique_id)),
             size = 2) +
  
  ## a separate TPC for each unique group
  geom_line(data = df.lf.nonarctic.quad.uni.sp, 
            aes(x = temp, y = X50., color = unique_id)) +
  
  ## Overall TPC
  geom_ribbon(data = df.lf.nonarctic.quad.uni.pop,
              aes(x = temp, ymin = X2.5., ymax = X97.5.),
              fill = "grey",
              alpha = 0.5) +
  geom_line(data = df.lf.nonarctic.quad.uni.pop,
            aes(x = temp, y = X50.), color = "black", linewidth = 1) +
  
  # Customize the axes and labels
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Time (days)"
  ) +
  
  # Customize legend
  scale_colour_discrete(name = "",
                        labels = c("Ae. albopictus 1",
                                   "Ae. albopictus 2",
                                   "Ae. albopictus 3",
                                   "Ae. albopictus 4",
                                   "Ae. albopictus 5")) +
  theme_bw()


plot.lf.nonarctic.quad.uni

ggsave("figures/lf.nonarctic.quad.uni.png", plot.lf.nonarctic.quad.uni,
       width = 10.3, height = 5.6)


## 3C. Fit gamma distributions to non-Arctic TPC parameters --------------------

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
lf.quad.prior.cf.dists <- data.frame(q = as.vector(lf.nonarctic.quad.uni$BUGSoutput$sims.list$cf.q),
                                    T0 = as.vector(lf.nonarctic.quad.uni$BUGSoutput$sims.list$cf.T0),
                                    Tm = as.vector(lf.nonarctic.quad.uni$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
lf.quad.prior.gamma.fits = apply(lf.quad.prior.cf.dists, 2, 
                                function(df) fitdistr(df, "gamma")$estimate)


save(lf.quad.prior.gamma.fits, file = "R-scripts/R2jags-objects/priors/lf.quad.priors.Rsave")


## 3D. Fit Arctic TPC using data-informed priors -------------------------------

# load("R-scripts/R2jags-objects/priors/lf.quad.priors.Rsave")


##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)


##### Set data
data <- data.lf.arctic
# Since this dataset has contains data from multiple species or multiple studies
# of the same species, we incorporated random effects on each thermal response
# parameter (q, T0, Tm) to addressed non-independence among observations 

## Create a unique id for each species-study combination
data <- data %>% 
  group_by(species, citation) %>% 
  mutate(unique_id = cur_group_id())

hypers <- lf.quad.prior.gamma.fits * 0.5
hypers[,2] <- lf.quad.prior.gamma.fits[,2]


# Gamma prior density for Tm

df.lf.Tm.prior <- tibble(
  temp = seq(0, 45, 0.1),
  density = dgamma(
    temp,
    shape = hypers["shape", "Tm"],
    rate  = hypers["rate", "Tm"]
  )
)

TPC.apex <- max(df.lf.nonarctic.quad.uni.pop$X50., na.rm = TRUE)

df.lf.Tm.prior <- df.lf.Tm.prior %>%
  mutate(
    density.scaled = density / max(density, na.rm = TRUE) * TPC.apex * 0.5
  )

df.lf.T0.prior <- tibble(
  temp = seq(0, 45, 0.1),
  density = dgamma(
    temp,
    shape = hypers["shape", "T0"],
    rate  = hypers["rate", "T0"]
  )
)

df.lf.T0.prior <- df.lf.T0.prior %>%
  mutate(
    max_density = max(density[is.finite(density)]), 
    density.scaled = density / max_density  * TPC.apex * 0.5
  ) %>% 
  dplyr::select(-max_density)
  


##### inits Function
inits<-function(){list(
  cf.q = 0.01,
  cf.Tm = 35,
  cf.T0 = 5,
  cf.sigma = rlnorm(1))}


##### Parameters to Estimate
parameters <- c("cf.q", "cf.T0", "cf.Tm", "cf.sigma", 
                "q", "T0", "Tm", "sigma_q", "sigma_T0", 
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
                 hypers = hypers)

##### Run JAGS
set.seed(123) # for reproducibility
lf.arctic.quad.inf <- jags(data = jag.data,
                          inits = inits,
                          parameters.to.save = parameters,
                          model.file = "R-scripts/quad_inf_raneff.txt",
                          n.thin = nt,
                          n.chains = nc,
                          n.burnin = nb,
                          n.iter = ni,
                          DIC = T,
                          working.directory = getwd()
)

## Save the model as Rdata 
save(lf.arctic.quad.inf, file = "R-scripts/R2jags-objects/all-mods/lf.arctic.quad.inf.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/all-mods/lf.arctic.quad.inf.Rdata")


## Diagnostics
##### Examine output
lf.nonarctic.quad.uni$BUGSoutput$summary[c("cf.T0", "cf.Tm", "cf.q", "cf.sigma", "deviance", "sigma_T0", "sigma_Tm", "sigma_q"),]
mcmcplot(lf.arctic.quad.inf, parms = c("cf.T0", "cf.Tm", "cf.q", "cf.sigma", "deviance", "sigma_T0", "sigma_Tm", "sigma_q"))

# Extract the DIC for future model comparisons
lf.arctic.quad.inf$BUGSoutput$DIC

## Plot data + fit
df.lf.arctic.quad.inf <- data.frame(lf.arctic.quad.inf$BUGSoutput$summary)

## Extract the model prediction
## Overall curve
df.lf.arctic.quad.inf.pop <- df.lf.arctic.quad.inf %>% 
  filter(grepl("z.trait.mu.pred.pop", rownames(df.lf.arctic.quad.inf))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)


## unique ID 1: Ae. cinereus
df.lf.arctic.quad.inf.1 <- df.lf.arctic.quad.inf %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[1,*]"), rownames(df.lf.arctic.quad.inf))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 1)

## unique ID 2: Ae. communis
df.lf.arctic.quad.inf.2 <- df.lf.arctic.quad.inf %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[2,*]"), rownames(df.lf.arctic.quad.inf))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 2)

## unique ID 3: Ae. impiger
df.lf.arctic.quad.inf.3 <- df.lf.arctic.quad.inf %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[3,*]"), rownames(df.lf.arctic.quad.inf))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 3)

## unique ID 4: Ae. punctor
df.lf.arctic.quad.inf.4 <- df.lf.arctic.quad.inf %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[4,*]"), rownames(df.lf.arctic.quad.inf))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 4)

## unique ID 5: Ae. vexans
df.lf.arctic.quad.inf.5 <- df.lf.arctic.quad.inf %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[5,*]"), rownames(df.lf.arctic.quad.inf))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 5)



## Combine the model prediciton of all three unique groups into a dataframe
df.lf.arctic.quad.inf.sp <- rbind(df.lf.arctic.quad.inf.1,
                                  df.lf.arctic.quad.inf.2,
                                  df.lf.arctic.quad.inf.3,
                                  df.lf.arctic.quad.inf.4,
                                  df.lf.arctic.quad.inf.5
                                  ) 

## Change unique_id into factor type
df.lf.arctic.quad.inf.sp$unique_id <- as.factor(df.lf.arctic.quad.inf.sp$unique_id)


##### Plot
plot.lf.arctic.quad.inf <- ggplot() +
  ## data
  geom_point(data = data,
             aes(x = temp, y = trait, colour = as.factor(unique_id)),
             size = 2) +
  
  # a separate TPC for each unique group
  geom_line(data = df.lf.arctic.quad.inf.sp,
            aes(x = temp, y = X50., color = unique_id)) +

  ## Overall TPC
  geom_ribbon(data = df.lf.arctic.quad.inf.pop,
              aes(x = temp, ymin = X2.5., ymax = X97.5.),
              fill = "grey",
              alpha = 0.5) +
  geom_line(data = df.lf.arctic.quad.inf.pop,
            aes(x = temp, y = X50.), color = "black", linewidth = 1) +
  
  geom_line(
    data = df.lf.Tm.prior,
    aes(x = temp, y = density.scaled),
    colour = "firebrick3",
    linetype = "dashed",
    linewidth = 1
  ) +
  
  geom_line(
    data = df.lf.T0.prior,
    aes(x = temp, y = density.scaled),
    colour = "skyblue",
    linetype = "dashed",
    linewidth = 1
  ) +
  
  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Time (days)",
       title = "D) lf, arctic sp. quadratic w/ gamma post") +
  
  # Customize legend
  scale_colour_discrete(name = "",
                        labels = c("Ae. cinereus",
                                   "Ae. communis",
                                   "Ae. impiger",
                                   "Ae. punctor",
                                   "Ae. vexans")) +
  theme_bw()


plot.lf.arctic.quad.inf

ggsave("figures/lf.arctic.quad.inf.png", plot.lf.arctic.quad.inf,
       width = 10.3, height = 5.6)



# 4. Compare model fit between Briere and Quadratic models ---------------------
plot.lf.arctic <- plot_grid(plot.lf.arctic.bri.uni, plot.lf.arctic.quad.uni,
                             plot.lf.arctic.bri.inf, plot.lf.arctic.quad.inf, ncol = 2)

plot.lf.arctic

ggsave("figures/lf.arctic.png", plot.lf.arctic,
       width = 10.3, height = 5.6)

##### Find best fitting model #####
# Add an identifying column in each model output dataframe
df.lf.arctic.bri.inf.pop <- df.lf.arctic.bri.inf.pop %>% 
  mutate(type = "briere")

df.lf.arctic.quad.inf.pop <- df.lf.arctic.quad.inf.pop %>% 
  mutate(type = "quadratic")

# Combine the two dataframes
df.all <- bind_rows(df.lf.arctic.bri.inf.pop, df.lf.arctic.quad.inf.pop)

##### Plot
plot.all <- df.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = X50., color = type), linewidth = 1) +
  geom_point(data = data.lf.arctic, aes(x = temp, y = trait), size = 2) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Time (days)"
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

ggsave("figures/lf.bri.quad.png", plot.all, width = 10.3, height = 5.6)


## DIC
lf.arctic.bri.uni$BUGSoutput$DIC 
lf.arctic.bri.inf$BUGSoutput$DIC 
lf.arctic.quad.uni$BUGSoutput$DIC
lf.arctic.quad.inf$BUGSoutput$DIC # This is the best fitting TPC



# Save best-fitting TPC in a separate folder
lf.arctic.mod <- lf.arctic.quad.inf
lf.nonarctic.mod <- lf.nonarctic.quad.uni

## Save the model as Rdata 
save(lf.arctic.mod, file = "R-scripts/R2jags-objects/best-fitting-mods/lf.arctic.mod.Rdata")
save(lf.nonarctic.mod, file = "R-scripts/R2jags-objects/best-fitting-mods/lf.nonarctic.mod.Rdata")


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
lf.TPC.analysis <- extractTPC_raneff(lf.arctic.mod, "lf", Temp.xs)
lf.arctic.predictions.summary <- lf.TPC.analysis[[1]]
lf.arctic.params.summary <- lf.TPC.analysis[[2]]
lf.arctic.params.fullposts <- lf.TPC.analysis[[3]]

write_csv(lf.arctic.predictions.summary, "data-processed/lf/lf.arctic.predictions.summary.csv")
write_csv(lf.arctic.params.summary, "data-processed/lf/lf.arctic.params.summary.csv")
write_csv(lf.arctic.params.fullposts, "data-processed/lf/lf.arctic.params.fullposts.csv")

##### non-Arctic #####
Temp.xs <- seq(0, 45, 0.5)
lf.TPC.analysis <- extractTPC_raneff(lf.nonarctic.mod, "lf", Temp.xs)
lf.nonarctic.predictions.summary <- lf.TPC.analysis[[1]]
lf.nonarctic.params.summary <- lf.TPC.analysis[[2]]
lf.nonarctic.params.fullposts <- lf.TPC.analysis[[3]]

write_csv(lf.nonarctic.predictions.summary, "data-processed/lf/lf.nonarctic.predictions.summary.csv")
write_csv(lf.nonarctic.params.summary, "data-processed/lf/lf.nonarctic.params.summary.csv")
write_csv(lf.nonarctic.params.fullposts, "data-processed/lf/lf.nonarctic.params.fullposts.csv")


##### Briere model #####
Temp.xs <- seq(0, 45, 0.1)
lf.TPC.analysis <- extractTPC_raneff(lf.arctic.bri.inf, "lf", Temp.xs)
lf.arctic.predictions.summary <- lf.TPC.analysis[[1]]
lf.arctic.params.summary <- lf.TPC.analysis[[2]]
lf.arctic.params.fullposts <- lf.TPC.analysis[[3]]

write_csv(lf.arctic.predictions.summary, "data-processed/supplemental-analysis/briere-only/lf.arctic.predictions.summary.csv")
write_csv(lf.arctic.params.summary, "data-processed/supplemental-analysis/briere-only/lf.arctic.params.summary.csv")
write_csv(lf.arctic.params.fullposts, "data-processed/supplemental-analysis/briere-only/lf.arctic.params.fullposts.csv")

Temp.xs <- seq(0, 45, 0.5)
lf.TPC.analysis <- extractTPC_raneff(lf.nonarctic.bri.uni, "lf", Temp.xs)
lf.nonarctic.predictions.summary <- lf.TPC.analysis[[1]]
lf.nonarctic.params.summary <- lf.TPC.analysis[[2]]
lf.nonarctic.params.fullposts <- lf.TPC.analysis[[3]]

write_csv(lf.nonarctic.predictions.summary, "data-processed/supplemental-analysis/briere-only/lf.nonarctic.predictions.summary.csv")
write_csv(lf.nonarctic.params.summary, "data-processed/supplemental-analysis/briere-only/lf.nonarctic.params.summary.csv")
write_csv(lf.nonarctic.params.fullposts, "data-processed/supplemental-analysis/briere-only/lf.nonarctic.params.fullposts.csv")
