## Lilian Chan, University of Guelph
## Arctic vector-borne disease transmission suitability model
##
## Purpose: Fit thermal performance curves (TPCs) for pathogen development rate 
## (PDR) using Bayesian inference (JAGS). Arctic species models are fit 
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
## data-processed/TraitData_PDR.csv -
##     Synthesized published trait data for PDR
##
## Outputs: 
## R-scripts/R2jags-objects/best-fitting-mods/PDR.arctic.mod.Rdata - 
##     Best-fitting TPC models for Arctic species
##
## R-scripts/R2jags-objects/best-fitting-mods/PDR.nonarctic.mod.Rdata -
##     Best-fitting TPC models for non-Arctic species
##
## data-processed/PDR/PDR.arctic.predictions.summary.csv -
##     Posterior summary of TPC predictions for Arctic species across temperatures
##
## data-processed/PDR/PDR.arctic.params.summary.csv -
##     Summary statistics of TPC parameters (Arctic TPC)
##
## data-processed/PDR/PDR.arctic.params.fullposts.csv -
##     Full posterior distributions for TPC parameters (Arctic TPC)
##
## data-processed/PDR/PDR.nonarctic.predictions.summary.csv -
##     Posterior summary of TPC predictions for non-Arctic species
##
## data-processed/PDR/PDR.nonarctic.params.summary.csv -
##     Summary statistics of TPC parameters (non-Arctic TPC)
##
## data-processed/PDR/PDR.nonarctic.params.fullposts.csv - 
##     Full posterior distributions for TPC parameters (non-Arctic TPC)


# 0. Set-up workspace -----------------------------------------------------

library(tidyverse)
library(readxl)
library(janitor)
library(R2jags)
library(mcmcplots) # Diagnostic plots for fits
library(MASS)
library(ggsci)
library(cowplot)

# Load functions
source("R-scripts/00_Functions.R")


# Load data
data.all <- read_csv("data-processed/TraitData_PDR.csv")
unique(data.all$species)

     
# Subset data
## Arctic species
data.PDR.arctic <- subset(data.all, type == "Arctic")

## Non-Arctic species
data.PDR.nonarctic <- subset(data.all, type == "non-Arctic")


# Plot the raw data
plot.data.PDR <- data.all %>% 
  ggplot(aes(x = temp, y = trait)) +
  geom_point(aes(colour = host_species)) +
  labs(y = "Parasite development rate (1/days)", x = "Temperature ºC") +
  scale_colour_discrete(name = "Species", labels = c("D. immitis (in Ae. canadensis)",
                                                     "V. eleguneniensis",
                                                     "W. bancrofti (in Ae. polynesiensis)",
                                                     "D. immitis (in Ae. triseriatus)",
                                                     "D. immitis (in Ae. trivittatus)",
                                                     "D. immitis (in Ae. vexans)",
                                                     "S. tundra"
  )) +
  facet_grid(rows = vars(type)) +
  theme_bw()

plot.data.PDR



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
data <- data.PDR.nonarctic

# Since this dataset has contains data from multiple species or multiple studies
# of the same species, we incorporated random effects on each thermal response
# parameter (q, T0, Tm) to addressed non-independence among observations 

## Create a unique id for each species-study combination
data <- data %>% 
  group_by(paras_species, host_species, citation) %>% 
  mutate(unique_id = cur_group_id())


##### Set priors
prior <- data.frame(q = c(0, 1),
                    T0 = c(0, 20),
                    Tm = c(30, 45),
                    sigma_q = c(0, 0.001),
                    sigma_T0 = c(0, 10),
                    sigma_Tm = c(0, 10)
                    )


##### inits Function
inits <- function(){list(
  cf.q = 0.01,
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
PDR.nonarctic.bri.uni <- jags(data = jag.data,
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
save(PDR.nonarctic.bri.uni, file = "R-scripts/R2jags-objects/all-mods/PDR.nonarctic.bri.uni.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/all-mods/PDR.nonarctic.bri.uni.Rdata")


## Diagnostics
##### Examine output
PDR.nonarctic.bri.uni$BUGSoutput$summary[1:8,]
mcmcplot(PDR.nonarctic.bri.uni)

# Extract the DIC for future model comparisons
PDR.nonarctic.bri.uni$BUGSoutput$DIC

## Plot data + fit
df.PDR.nonarctic.bri.uni <- data.frame(PDR.nonarctic.bri.uni$BUGSoutput$summary)[-(1:8),]

## Extract the model prediction
## Overall curve
df.PDR.nonarctic.bri.uni.pop <- df.PDR.nonarctic.bri.uni %>% 
  filter(grepl("z.trait.mu.pred.pop", rownames(df.PDR.nonarctic.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)


## unique ID 1: Wuchereria bancrofti
df.PDR.nonarctic.bri.uni.1 <- df.PDR.nonarctic.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[1,*]"), rownames(df.PDR.nonarctic.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 1)

## unique ID 2: Dirofilaria immitis (in Aedes canadensis)
df.PDR.nonarctic.bri.uni.2 <- df.PDR.nonarctic.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[2,*]"), rownames(df.PDR.nonarctic.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 2)

## unique ID 3: Dirofilaria immitis (in Aedes triseriatus)
df.PDR.nonarctic.bri.uni.3 <- df.PDR.nonarctic.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[3,*]"), rownames(df.PDR.nonarctic.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 3)

## unique ID 4: Dirofilaria immitis (in Aedes trivittatus)
df.PDR.nonarctic.bri.uni.4 <- df.PDR.nonarctic.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[4,*]"), rownames(df.PDR.nonarctic.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 4)


## unique ID 5: Dirofilaria immitis (in Aedes vexans 1)
df.PDR.nonarctic.bri.uni.5 <- df.PDR.nonarctic.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[5,*]"), rownames(df.PDR.nonarctic.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 5)


## unique ID 6: Dirofilaria immitis (in Aedes vexans 2)
df.PDR.nonarctic.bri.uni.6 <- df.PDR.nonarctic.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[6,*]"), rownames(df.PDR.nonarctic.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 6)


## Combine the model prediciton of all three unique groups into a dataframe
df.PDR.nonarctic.bri.uni.sp <- rbind(df.PDR.nonarctic.bri.uni.1,
                                     df.PDR.nonarctic.bri.uni.2,
                                     df.PDR.nonarctic.bri.uni.3,
                                     df.PDR.nonarctic.bri.uni.4,
                                     df.PDR.nonarctic.bri.uni.5,
                                     df.PDR.nonarctic.bri.uni.6
                                     ) 

## Change unique_id into factor type
df.PDR.nonarctic.bri.uni.sp$unique_id <- as.factor(df.PDR.nonarctic.bri.uni.sp$unique_id)


head(df.PDR.nonarctic.bri.uni)


##### Plot
plot.PDR.nonarctic.bri.uni <- ggplot() +
  ## data
  geom_point(data = data,
             aes(x = temp, y = trait, colour = as.factor(unique_id)),
             size = 2) +
  
  ## a separate TPC for each unique group
  geom_line(data = df.PDR.nonarctic.bri.uni.sp, 
            aes(x = temp, y = X50., color = unique_id)) +
  
  ## Overall TPC
  geom_ribbon(data = df.PDR.nonarctic.bri.uni.pop,
              aes(x = temp, ymin = X2.5., ymax = X97.5.),
              fill = "grey",
              alpha = 0.5) +
  geom_line(data = df.PDR.nonarctic.bri.uni.pop,
            aes(x = temp, y = X50.), color = "black", linewidth = 1) +
  
  # Customize the axes and labels
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Development rate (days-1)"
  ) +
  
  # Customize legend
  scale_colour_discrete(name = "",
                        labels = c("Wuchereria bancrofti",
                                   "Dirofilaria immitis (in Ae. canadensis)",
                                   "Dirofilaria immitis (in Ae. triseriatus)",
                                   "Dirofilaria immitis (in Ae. trivittatus)",
                                   "Dirofilaria immitis (in Ae. vexans 1)",
                                   "Dirofilaria immitis (in Ae. vexans 2)")) +
  theme_bw()


plot.PDR.nonarctic.bri.uni

ggsave("figures/PDR.nonarctic.bri.uni.png", plot.PDR.nonarctic.bri.uni,
       width = 10.3, height = 5.6)



## 2B. Fit gamma distributions to non-Arctic TPC parameters --------------------

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
PDR.bri.prior.cf.dists <- data.frame(q = as.vector(PDR.nonarctic.bri.uni$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(PDR.nonarctic.bri.uni$BUGSoutput$sims.list$cf.T0),
                                    Tm = as.vector(PDR.nonarctic.bri.uni$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
PDR.bri.prior.gamma.fits <- apply(PDR.bri.prior.cf.dists, 2, 
                                  function(df) fitdistr(df, "gamma")$estimate)

save(PDR.bri.prior.gamma.fits, file = "R-scripts/R2jags-objects/priors/PDR.bri.priors.Rsave")



## 2C. Fit Arctic TPC using data-informed priors -------------------------------

# load("R-scripts/R2jags-objects/priors/PDR.bri.priors.Rsave")


##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)


##### Set data
data <- data.PDR.arctic
hypers <- PDR.bri.prior.gamma.fits * 0.1
hypers[,3] <- PDR.bri.prior.gamma.fits[,3]

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
ggsave("figures/PDR.bri.priors.png", priors, width = 10.3, height = 5.6)

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
save(PDR.arctic.bri.inf, file = "R-scripts/R2jags-objects/all-mods/PDR.arctic.bri.inf.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/all-mods/PDR.arctic.bri.inf.Rdata")


## Diagnostics
##### Examine output
PDR.arctic.bri.inf$BUGSoutput$summary[1:5,]
mcmcplot(PDR.arctic.bri.inf)

# Extract the DIC for future model comparisons
PDR.arctic.bri.inf$BUGSoutput$DIC

## Plot data + fit
df.PDR.arctic.bri.inf <- data.frame(PDR.arctic.bri.inf$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)

head(df.PDR.arctic.bri.inf)

##### Plot
plot.PDR.arctic.bri.inf <- df.PDR.arctic.bri.inf %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "#4363d8", alpha = 0.5) +
  geom_line(aes(y = X50.), color = "blue", linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait), size = 2) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Parasite development rate (days-1)"
  ) +
  theme_bw()

plot.PDR.arctic.bri.inf

ggsave("figures/PDR.arctic.bri.inf.png", plot.PDR.arctic.bri.inf,
       width = 10.3, height = 5.6)




# 3. Fitting TPC (quadratic) ---------------------------------------------------

## 3A. Fit non-Arctic TPC to generate priors -----------------------------------

##### Temp sequence for derived quantity calculations
# For priors - fewer temps for derived calculations makes it go faster
Temp.xs <- seq(0, 45, 0.5)
N.Temp.xs <-length(Temp.xs)

##### Set data
data <- data.PDR.nonarctic

# Since this dataset has contains data from multiple species or multiple studies
# of the same species, we incorporated random effects on each thermal response
# parameter (q, T0, Tm) to addressed non-independence among observations 

## Create a unique id for each species-study combination
data <- data %>% 
  group_by(paras_species, host_species, citation) %>% 
  mutate(unique_id = cur_group_id())


##### Set priors
prior <- data.frame(q = c(0, 1),
                    T0 = c(0, 20),
                    Tm = c(30, 45),
                    sigma_q = c(0, 0.001),
                    sigma_T0 = c(0, 10),
                    sigma_Tm = c(0, 10)
)


##### inits Function
inits <- function(){list(
  cf.q = 0.01,
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
PDR.nonarctic.quad.uni <- jags(data = jag.data,
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
save(PDR.nonarctic.quad.uni, file = "R-scripts/R2jags-objects/all-mods/PDR.nonarctic.quad.uni.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/all-mods/PDR.nonarctic.quad.uni.Rdata")


## Diagnostics
##### Examine output
PDR.nonarctic.quad.uni$BUGSoutput$summary[1:8,]
mcmcplot(PDR.nonarctic.quad.uni)

# Extract the DIC for future model comparisons
PDR.nonarctic.quad.uni$BUGSoutput$DIC

## Plot data + fit
df.PDR.nonarctic.quad.uni <- data.frame(PDR.nonarctic.quad.uni$BUGSoutput$summary)[-(1:8),]

## Extract the model prediction
## Overall curve
df.PDR.nonarctic.quad.uni.pop <- df.PDR.nonarctic.quad.uni %>% 
  filter(grepl("z.trait.mu.pred.pop", rownames(df.PDR.nonarctic.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)


## unique ID 1: Wuchereria bancrofti
df.PDR.nonarctic.quad.uni.1 <- df.PDR.nonarctic.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[1,*]"), rownames(df.PDR.nonarctic.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 1)

## unique ID 2: Dirofilaria immitis (in Aedes canadensis)
df.PDR.nonarctic.quad.uni.2 <- df.PDR.nonarctic.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[2,*]"), rownames(df.PDR.nonarctic.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 2)

## unique ID 3: Dirofilaria immitis (in Aedes triseriatus)
df.PDR.nonarctic.quad.uni.3 <- df.PDR.nonarctic.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[3,*]"), rownames(df.PDR.nonarctic.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 3)

## unique ID 4: Dirofilaria immitis (in Aedes trivittatus)
df.PDR.nonarctic.quad.uni.4 <- df.PDR.nonarctic.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[4,*]"), rownames(df.PDR.nonarctic.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 4)


## unique ID 5: Dirofilaria immitis (in Aedes vexans 1)
df.PDR.nonarctic.quad.uni.5 <- df.PDR.nonarctic.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[5,*]"), rownames(df.PDR.nonarctic.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 5)


## unique ID 6: Dirofilaria immitis (in Aedes vexans 2)
df.PDR.nonarctic.quad.uni.6 <- df.PDR.nonarctic.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[6,*]"), rownames(df.PDR.nonarctic.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 6)


## Combine the model prediciton of all three unique groups into a dataframe
df.PDR.nonarctic.quad.uni.sp <- rbind(df.PDR.nonarctic.quad.uni.1,
                                     df.PDR.nonarctic.quad.uni.2,
                                     df.PDR.nonarctic.quad.uni.3,
                                     df.PDR.nonarctic.quad.uni.4,
                                     df.PDR.nonarctic.quad.uni.5,
                                     df.PDR.nonarctic.quad.uni.6
) 

## Change unique_id into factor type
df.PDR.nonarctic.quad.uni.sp$unique_id <- as.factor(df.PDR.nonarctic.quad.uni.sp$unique_id)


head(df.PDR.nonarctic.quad.uni)


##### Plot
plot.PDR.nonarctic.quad.uni <- ggplot() +
  ## data
  geom_point(data = data,
             aes(x = temp, y = trait, colour = as.factor(unique_id)),
             size = 2) +
  
  ## a separate TPC for each unique group
  geom_line(data = df.PDR.nonarctic.quad.uni.sp, 
            aes(x = temp, y = X50., color = unique_id)) +
  
  ## Overall TPC
  geom_ribbon(data = df.PDR.nonarctic.quad.uni.pop,
              aes(x = temp, ymin = X2.5., ymax = X97.5.),
              fill = "grey",
              alpha = 0.5) +
  geom_line(data = df.PDR.nonarctic.quad.uni.pop,
            aes(x = temp, y = X50.), color = "black", linewidth = 1) +
  
  # Customize the axes and labels
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Development rate (days-1)"
  ) +
  
  # Customize legend
  scale_colour_discrete(name = "",
                        labels = c("Wuchereria bancrofti",
                                   "Dirofilaria immitis (in Aedes canadensis)",
                                   "Dirofilaria immitis (in Aedes triseriatus)",
                                   "Dirofilaria immitis (in Aedes trivittatus)",
                                   "Dirofilaria immitis (in Aedes vexans 1)",
                                   "Dirofilaria immitis (in Aedes vexans 2)")) +
  theme_bw()


plot.PDR.nonarctic.quad.uni

ggsave("figures/PDR.nonarctic.quad.uni.png", plot.PDR.nonarctic.quad.uni,
       width = 10.3, height = 5.6)


## 3B. Fit gamma distributions to non-Arctic TPC parameters --------------------

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
PDR.quad.prior.cf.dists <- data.frame(q = as.vector(PDR.nonarctic.quad.uni$BUGSoutput$sims.list$cf.q),
                                      T0 = as.vector(PDR.nonarctic.quad.uni$BUGSoutput$sims.list$cf.T0),
                                      Tm = as.vector(PDR.nonarctic.quad.uni$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
PDR.quad.prior.gamma.fits = apply(PDR.quad.prior.cf.dists, 2, 
                                      function(df) fitdistr(df, "gamma")$estimate)


save(PDR.quad.prior.gamma.fits, file = "R-scripts/R2jags-objects/priors/PDR.quad.priors.Rsave")


## 3C. Fit Arctic TPC using data-informed priors -------------------------------

# load("R-scripts/R2jags-objects/priors/PDR.quad.priors.Rsave")

##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)


##### Set data
data <- data.PDR.arctic
hypers <- PDR.quad.prior.gamma.fits * 0.1
hypers[,3] <- PDR.quad.prior.gamma.fits[,3]

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
ggsave("figures/PDR.quad.priors.png", priors, width = 10.3, height = 5.6)

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
save(PDR.arctic.quad.inf, file = "R-scripts/R2jags-objects/all-mods/PDR.arctic.quad.inf.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/all-mods/PDR.arctic.quad.inf.Rdata")


## Diagnostics
##### Examine output
PDR.arctic.quad.inf$BUGSoutput$summary[1:5,]
mcmcplot(PDR.arctic.quad.inf)

# Extract the DIC for future model comparisons
PDR.arctic.quad.inf$BUGSoutput$DIC

## Plot data + fit
df.PDR.arctic.quad.inf <- data.frame(PDR.arctic.quad.inf$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)

head(df.PDR.arctic.quad.inf)


##### Plot
plot.PDR.arctic.quad.inf <- df.PDR.arctic.quad.inf %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "#4363d8", alpha = 0.5) +
  geom_line(aes(y = X50.), color = "blue", linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait), size = 2) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Development rate (days-1)"
  ) +
  theme_bw()

plot.PDR.arctic.quad.inf

ggsave("figures/PDR.arctic.quad.inf.png", plot.PDR.arctic.quad.inf,
       width = 10.3, height = 5.6)



# 4. Compare model fit between Briere and Quadratic models ---------------------

##### Find best fitting model #####
# Add an identifying column in each model output dataframe
df.PDR.arctic.bri.inf <- df.PDR.arctic.bri.inf %>% 
  mutate(type = "briere")

df.PDR.arctic.quad.inf <- df.PDR.arctic.quad.inf %>% 
  mutate(type = "quadratic")

# Combine the three dataframes
df.all <- bind_rows(df.PDR.arctic.bri.inf, df.PDR.arctic.quad.inf)


##### Plot
plot.all <- df.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = X50., color = type), linewidth = 1) +
  geom_point(data = data.PDR.arctic, aes(x = temp, y = trait), size = 2) +

  # Customize the axes and labels
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Development rate (days-1)"
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

ggsave("figures/PDR.bri.quad.png", plot.all, width = 10.3, height = 5.6)


## DIC
PDR.arctic.bri.inf$BUGSoutput$DIC # This is the best fitting TPC
PDR.arctic.quad.inf$BUGSoutput$DIC 


##### Plot Arctic vs. non-Arctic TPCs for the best fitting TPC #####
df.PDR.nonarctic.bri.uni.pop <- df.PDR.nonarctic.bri.uni.pop %>% 
  mutate(type = "non-Arctic")

df.PDR.arctic.bri.inf <- df.PDR.arctic.bri.inf %>% 
  mutate(type = "Arctic")

df.arctic.nonarctic <- bind_rows(df.PDR.nonarctic.bri.uni.pop, df.PDR.arctic.bri.inf)

plot.arctic.nonarctic <- df.arctic.nonarctic %>% 
  ggplot(aes(x = temp)) +
  geom_point(data = data.all, aes(x = temp, y = trait, colour = type), size = 2) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = X50., color = type), linewidth = 1) +
  
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Development rate (days-1)"
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

ggsave("figures/PDR.arctic.nonarctic.png", plot.arctic.nonarctic, width = 10.3, height = 5.6)


# Save best-fitting TPC in a separate folder
PDR.arctic.mod <- PDR.arctic.bri.inf
PDR.nonarctic.mod <- PDR.nonarctic.bri.uni

## Save the model as Rdata 
save(PDR.arctic.mod, file = "R-scripts/R2jags-objects/best-fitting-mods/PDR.arctic.mod.Rdata")
save(PDR.nonarctic.mod, file = "R-scripts/R2jags-objects/best-fitting-mods/PDR.nonarctic.mod.Rdata")



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
PDR.TPC.analysis <- extractTPC(PDR.arctic.bri.inf, "PDR", Temp.xs)
PDR.arctic.predictions.summary <- PDR.TPC.analysis[[1]]
PDR.arctic.params.summary <- PDR.TPC.analysis[[2]]
PDR.arctic.params.fullposts <- PDR.TPC.analysis[[3]]

write_csv(PDR.arctic.predictions.summary, "data-processed/PDR/PDR.arctic.predictions.summary.csv")
write_csv(PDR.arctic.params.summary, "data-processed/PDR/PDR.arctic.params.summary.csv")
write_csv(PDR.arctic.params.fullposts, "data-processed/PDR/PDR.arctic.params.fullposts.csv")

##### non-Arctic #####
Temp.xs <- seq(0, 45, 0.5)
PDR.TPC.analysis <- extractTPC_raneff(PDR.nonarctic.bri.uni, "PDR", Temp.xs)
PDR.nonarctic.predictions.summary <- PDR.TPC.analysis[[1]]
PDR.nonarctic.params.summary <- PDR.TPC.analysis[[2]]
PDR.nonarctic.params.fullposts <- PDR.TPC.analysis[[3]]

write_csv(PDR.nonarctic.predictions.summary, "data-processed/PDR/PDR.nonarctic.predictions.summary.csv")
write_csv(PDR.nonarctic.params.summary, "data-processed/PDR/PDR.nonarctic.params.summary.csv")
write_csv(PDR.nonarctic.params.fullposts, "data-processed/PDR/PDR.nonarctic.params.fullposts.csv")


