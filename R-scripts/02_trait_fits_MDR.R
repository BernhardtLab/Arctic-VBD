## Lilian Chan, University of Guelph
## Arctic vector-borne disease transmission suitability model
##
## Purpose: Fit thermal performance curves (TPCs) for mosquito development rate 
## (MDR) using Bayesian inference (JAGS). Arctic species models are fit 
## using data-informed priors derived from non-Arctic species.
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
## data-processed/TraitData_MDR.csv - Synthesized published trait data for MDR
##
## Outputs: 
## data-processed/MDR/MDR.arctic.predictions.summary.csv - Posterior summary of
## TPC predictions for Arctic species across temperatures
##
## data-processed/MDR/MDR.arctic.params.summary.csv - Summary statistics of TPC 
## parameters (Arctic TPC)
##
## data-processed/MDR/MDR.arctic.params.fullposts.csv - Full posterior 
## distributions for q, Tmin, Tmax, and Tbreadth (Arctic TPC)
##
## data-processed/MDR/MDR.nonarctic.predictions.summary.csv - Posterior summary 
## of TPC predictions for non-Arctic species
##
## data-processed/MDR/MDR.nonarctic.params.summary.csv -  Summary statistics of 
## TPC parameters (non-Arctic TPC)
##
## data-processed/MDR/MDR.nonarctic.params.fullposts.csv - Full posterior 
## distributions for q, Tmin, Tmax, and Tbreadth (non-Arctic TPC)


# 0. Set-up workspace ----------------------------------------------------------

library(tidyverse)
library(readxl)
library(janitor)
library(R2jags)
library(mcmcplots) # Diagnostic plots for fits
library(MASS)

# Load functions
source("R-scripts/00_Functions.R")

# Load data
data.all <- read_csv("data-processed/TraitData_MDR.csv")
unique(data.all$species)


# Subset data
## Arctic species
data.MDR.arctic <- subset(data.all, type == "Arctic")

## Non-Arctic species
data.MDR.nonarctic <- subset(data.all, type == "non-Arctic")


## Plot raw data
plot.data.MDR <- data.all %>% 
  ggplot(aes(x = temp, y = trait, colour = species)) +
  geom_point() +
  labs(y = "Mosquito development rate (days)", x = "Temperature ºC") +
  facet_grid(rows = vars(type), scales = "free_y") +
  theme_bw()

plot.data.MDR



# 1. MCMC Settings for all models ----------------------------------------------

# Number of posterior dist elements = [(ni - nb) / nt] * nc = [(45000 - 5000) / 8] * 3 = 15000
ni <- 45000 # number of iterations in each chain
nb <- 5000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 3 # number of chains

set.seed(123) # for reproducibility


# 2. Fitting TPC (Briere) ------------------------------------------------------

## 2A. Fit non-Arctic TPC to generate priors -----------------------------------

##### Temp sequence for derived quantity calculations
# For priors - fewer temps for derived calculations makes it go faster
Temp.xs <- seq(0, 45, 0.5)
N.Temp.xs <-length(Temp.xs)


##### Set data
data <- data.MDR.nonarctic

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
                    Tm = c(25, 45),
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
MDR.nonarctic.bri.uni <- jags(data = jag.data,
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
save(MDR.nonarctic.bri.uni, file = "R-scripts/R2jags-objects/all-mods/MDR.nonarctic.bri.uni.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/all-mods/MDR.nonarctic.bri.uni.Rdata")


## Diagnostics
##### Examine output
MDR.nonarctic.bri.uni$BUGSoutput$summary[1:8,]
# mcmcplot(MDR.nonarctic.bri.uni)

# Extract the DIC for future model comparisons
MDR.nonarctic.bri.uni$BUGSoutput$DIC

## Plot data + fit
df.MDR.nonarctic.bri.uni <- data.frame(MDR.nonarctic.bri.uni$BUGSoutput$summary)[-(1:8),]

## Extract the model prediction
## Overall curve
df.MDR.nonarctic.bri.uni.pop <- df.MDR.nonarctic.bri.uni %>% 
  filter(grepl("z.trait.mu.pred.pop", rownames(df.MDR.nonarctic.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)


## unique ID 1: Ae. albopictus (Alto 2001)
df.MDR.nonarctic.bri.uni.1 <- df.MDR.nonarctic.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[1,*]"), rownames(df.MDR.nonarctic.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 1)

## unique ID 2: Ae. albopictus (Briegel 2001)
df.MDR.nonarctic.bri.uni.2 <- df.MDR.nonarctic.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[2,*]"), rownames(df.MDR.nonarctic.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 2)

## unique ID 3: Ae. albopictus (Calado 2002)
df.MDR.nonarctic.bri.uni.3 <- df.MDR.nonarctic.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[3,*]"), rownames(df.MDR.nonarctic.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 3)

## unique ID 4: Ae. albopictus (Delatte 2009)
df.MDR.nonarctic.bri.uni.4 <- df.MDR.nonarctic.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[4,*]"), rownames(df.MDR.nonarctic.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 4)


## unique ID 5: Ae. albopictus (Ezeakacha 2015)
df.MDR.nonarctic.bri.uni.5 <- df.MDR.nonarctic.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[5,*]"), rownames(df.MDR.nonarctic.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 5)


## unique ID 6: Ae. albopictus (Muturi 2011)
df.MDR.nonarctic.bri.uni.6 <- df.MDR.nonarctic.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[6,*]"), rownames(df.MDR.nonarctic.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 6)


## unique ID 7: Ae. albopictus (Westbrook thesis 2010)
df.MDR.nonarctic.bri.uni.7 <- df.MDR.nonarctic.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[7,*]"), rownames(df.MDR.nonarctic.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 7)


## unique ID 8: Ae. albopictus (Westbrook 2010)
df.MDR.nonarctic.bri.uni.8 <- df.MDR.nonarctic.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[8,*]"), rownames(df.MDR.nonarctic.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 8)


## unique ID 9: Ae. albopictus (Witwatana 2006)
df.MDR.nonarctic.bri.uni.9 <- df.MDR.nonarctic.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[9,*]"), rownames(df.MDR.nonarctic.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 9)


## unique ID 10: Ae. albopictus (Yee 2016)
df.MDR.nonarctic.bri.uni.10 <- df.MDR.nonarctic.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[10,*]"), rownames(df.MDR.nonarctic.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 10)


## unique ID 11: Ae. nigromaculis
df.MDR.nonarctic.bri.uni.11 <- df.MDR.nonarctic.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[11,*]"), rownames(df.MDR.nonarctic.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 11)


## unique ID 12: Ae. sollicitans
df.MDR.nonarctic.bri.uni.12 <- df.MDR.nonarctic.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[12,*]"), rownames(df.MDR.nonarctic.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 12)


## unique ID 13: Ae. triseriatus (Jalil 1972)
df.MDR.nonarctic.bri.uni.13 <- df.MDR.nonarctic.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[13,*]"), rownames(df.MDR.nonarctic.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 13)


## unique ID 14: Ae. triseriatus (Shelton 1973)
df.MDR.nonarctic.bri.uni.14 <- df.MDR.nonarctic.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[14,*]"), rownames(df.MDR.nonarctic.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 14)



## Combine the model prediciton of all three unique groups into a dataframe
df.MDR.nonarctic.bri.uni.sp <- rbind(df.MDR.nonarctic.bri.uni.1,
                                     df.MDR.nonarctic.bri.uni.2,
                                     df.MDR.nonarctic.bri.uni.3,
                                     df.MDR.nonarctic.bri.uni.4,
                                     df.MDR.nonarctic.bri.uni.5,
                                     df.MDR.nonarctic.bri.uni.6,
                                     df.MDR.nonarctic.bri.uni.7,
                                     df.MDR.nonarctic.bri.uni.8,
                                     df.MDR.nonarctic.bri.uni.9,
                                     df.MDR.nonarctic.bri.uni.10,
                                     df.MDR.nonarctic.bri.uni.11,
                                     df.MDR.nonarctic.bri.uni.12,
                                     df.MDR.nonarctic.bri.uni.13,
                                     df.MDR.nonarctic.bri.uni.14
) 

## Change unique_id into factor type
df.MDR.nonarctic.bri.uni.sp$unique_id <- as.factor(df.MDR.nonarctic.bri.uni.sp$unique_id)


head(df.MDR.nonarctic.bri.uni)

##### Plot
plot.MDR.nonarctic.bri.uni <- ggplot() +
  ## data
  geom_point(data = data,
             aes(x = temp, y = trait, colour = as.factor(unique_id)),
             size = 2) +
  
  ## a separate TPC for each unique group
  geom_line(data = df.MDR.nonarctic.bri.uni.sp, 
            aes(x = temp, y = X50., color = unique_id)) +
  
  ## Overall TPC
  geom_ribbon(data = df.MDR.nonarctic.bri.uni.pop,
              aes(x = temp, ymin = X2.5., ymax = X97.5.),
              fill = "grey",
              alpha = 0.5) +
  geom_line(data = df.MDR.nonarctic.bri.uni.pop,
            aes(x = temp, y = X50.), color = "black", linewidth = 1) +
  
  # Customize the axes and labels
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Development rate (days-1)"
    ) +
  
  # Customize legend
  scale_colour_discrete(name = "",
                        labels = c("Ae. albopictus 1",
                                   "Ae. albopictus 2",
                                   "Ae. albopictus 3",
                                   "Ae. albopictus 4",
                                   "Ae. albopictus 5",
                                   "Ae. albopictus 6",
                                   "Ae. albopictus 7",
                                   "Ae. albopictus 8",
                                   "Ae. albopictus 9",
                                   "Ae. albopictus 10",
                                   "Ae. nigromaculis",
                                   "Ae. sollicitans",
                                   "Ae. triseriatus 1",
                                   "Ae. triseriatus 2")) +
  theme_bw()


plot.MDR.nonarctic.bri.uni

ggsave("figures/MDR.nonarctic.bri.uni.png", plot.MDR.nonarctic.bri.uni,
       width = 10.3, height = 5.6)




## 2B. Fit gamma distributions to non-Arctic TPC parameters --------------------


# Get the posterior dists for 3 main parameters (not sigma) into a data frame
MDR.bri.prior.cf.dists <- data.frame(q = as.vector(MDR.nonarctic.bri.uni$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(MDR.nonarctic.bri.uni$BUGSoutput$sims.list$cf.T0),
                                     Tm = as.vector(MDR.nonarctic.bri.uni$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
MDR.bri.prior.gamma.fits <- apply(MDR.bri.prior.cf.dists, 2, 
                                 function(df) fitdistr(df, "gamma")$estimate)


save(MDR.bri.prior.gamma.fits, file = "R-scripts/R2jags-objects/priors/MDR.bri.priors.Rsave")



## 2C. Fit Arctic TPC using data-informed priors -------------------------------

# load("R-scripts/R2jags-objects/priors/MDR.bri.priors.Rsave")

##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)

##### Set data
data <- data.MDR.arctic

# Although this dataset contains data from multiple species or multiple studies
# of the same species, there are less than 5 species-study combination. 
# Thus we will not incorporate random effects 

hypers <- MDR.bri.prior.gamma.fits

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
                 N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
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
save(MDR.arctic.bri.inf, file = "R-scripts/R2jags-objects/all-mods/MDR.arctic.bri.inf.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/all-mods/MDR.arctic.bri.inf.Rdata")


## Diagnostics
##### Examine output
MDR.arctic.bri.inf$BUGSoutput$summary[1:5,]
# mcmcplot(MDR.arctic.bri.inf)

# Extract the DIC for future model comparisons
MDR.arctic.bri.inf$BUGSoutput$DIC

## Plot data + fit
df.MDR.arctic.bri.inf <- data.frame(MDR.arctic.bri.inf$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, X50., sd, X2.5., X50., X97.5.)

head(df.MDR.arctic.bri.inf)

##### Plot
plot.MDR.arctic.bri.inf <- df.MDR.arctic.bri.inf %>%
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "#4363d8", alpha = 0.5) +
  geom_line(aes(y = X50.), color = "blue", linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait), size = 2) +

  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Development rate (days-1)"
    ) +
  theme_bw()

plot.MDR.arctic.bri.inf

ggsave("figures/MDR.arctic.bri.inf.png", plot.MDR.arctic.bri.inf,
       width = 10.3, height = 5.6)




# 3. Fitting TPC (quadratic) ---------------------------------------------------

## 3A. Fit non-Arctic TPC to generate priors -----------------------------------


##### Temp sequence for derived quantity calculations
# For priors - fewer temps for derived calculations makes it go faster
Temp.xs <- seq(0, 45, 0.5)
N.Temp.xs <-length(Temp.xs)


##### Set data
data <- data.MDR.nonarctic

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
                    Tm = c(25, 45),
                    sigma_q = c(0, 0.01),
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
MDR.nonarctic.quad.uni <- jags(data = jag.data,
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
save(MDR.nonarctic.quad.uni, file = "R-scripts/R2jags-objects/all-mods/MDR.nonarctic.quad.uni.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/all-mods/MDR.nonarctic.quad.uni.Rdata")


## Diagnostics
##### Examine output
MDR.nonarctic.quad.uni$BUGSoutput$summary[1:8,]
# mcmcplot(MDR.nonarctic.quad.uni)

# Extract the DIC for future model comparisons
MDR.nonarctic.quad.uni$BUGSoutput$DIC

## Plot data + fit
df.MDR.nonarctic.quad.uni <- data.frame(MDR.nonarctic.quad.uni$BUGSoutput$summary)[-(1:8),]

## Extract the model prediction
## Overall curve
df.MDR.nonarctic.quad.uni.pop <- df.MDR.nonarctic.quad.uni %>% 
  filter(grepl("z.trait.mu.pred.pop", rownames(df.MDR.nonarctic.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)


## unique ID 1: Ae. albopictus (Alto 2001)
df.MDR.nonarctic.quad.uni.1 <- df.MDR.nonarctic.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[1,*]"), rownames(df.MDR.nonarctic.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 1)

## unique ID 2: Ae. albopictus (quadegel 2001)
df.MDR.nonarctic.quad.uni.2 <- df.MDR.nonarctic.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[2,*]"), rownames(df.MDR.nonarctic.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 2)

## unique ID 3: Ae. albopictus (Calado 2002)
df.MDR.nonarctic.quad.uni.3 <- df.MDR.nonarctic.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[3,*]"), rownames(df.MDR.nonarctic.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 3)

## unique ID 4: Ae. albopictus (Delatte 2009)
df.MDR.nonarctic.quad.uni.4 <- df.MDR.nonarctic.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[4,*]"), rownames(df.MDR.nonarctic.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 4)


## unique ID 5: Ae. albopictus (Ezeakacha 2015)
df.MDR.nonarctic.quad.uni.5 <- df.MDR.nonarctic.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[5,*]"), rownames(df.MDR.nonarctic.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 5)


## unique ID 6: Ae. albopictus (Muturi 2011)
df.MDR.nonarctic.quad.uni.6 <- df.MDR.nonarctic.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[6,*]"), rownames(df.MDR.nonarctic.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 6)


## unique ID 7: Ae. albopictus (Westbrook thesis 2010)
df.MDR.nonarctic.quad.uni.7 <- df.MDR.nonarctic.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[7,*]"), rownames(df.MDR.nonarctic.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 7)


## unique ID 8: Ae. albopictus (Westbrook 2010)
df.MDR.nonarctic.quad.uni.8 <- df.MDR.nonarctic.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[8,*]"), rownames(df.MDR.nonarctic.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 8)


## unique ID 9: Ae. albopictus (Witwatana 2006)
df.MDR.nonarctic.quad.uni.9 <- df.MDR.nonarctic.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[9,*]"), rownames(df.MDR.nonarctic.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 9)


## unique ID 10: Ae. albopictus (Yee 2016)
df.MDR.nonarctic.quad.uni.10 <- df.MDR.nonarctic.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[10,*]"), rownames(df.MDR.nonarctic.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 10)


## unique ID 11: Ae. nigromaculis
df.MDR.nonarctic.quad.uni.11 <- df.MDR.nonarctic.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[11,*]"), rownames(df.MDR.nonarctic.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 11)


## unique ID 12: Ae. sollicitans
df.MDR.nonarctic.quad.uni.12 <- df.MDR.nonarctic.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[12,*]"), rownames(df.MDR.nonarctic.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 12)


## unique ID 13: Ae. triseriatus (Jalil 1972)
df.MDR.nonarctic.quad.uni.13 <- df.MDR.nonarctic.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[13,*]"), rownames(df.MDR.nonarctic.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 13)


## unique ID 14: Ae. triseriatus (Shelton 1973)
df.MDR.nonarctic.quad.uni.14 <- df.MDR.nonarctic.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[14,*]"), rownames(df.MDR.nonarctic.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 14)



## Combine the model prediciton of all three unique groups into a dataframe
df.MDR.nonarctic.quad.uni.sp <- rbind(df.MDR.nonarctic.quad.uni.1,
                                      df.MDR.nonarctic.quad.uni.2,
                                      df.MDR.nonarctic.quad.uni.3,
                                      df.MDR.nonarctic.quad.uni.4,
                                      df.MDR.nonarctic.quad.uni.5,
                                      df.MDR.nonarctic.quad.uni.6,
                                      df.MDR.nonarctic.quad.uni.7,
                                      df.MDR.nonarctic.quad.uni.8,
                                      df.MDR.nonarctic.quad.uni.9,
                                      df.MDR.nonarctic.quad.uni.10,
                                      df.MDR.nonarctic.quad.uni.11,
                                      df.MDR.nonarctic.quad.uni.12,
                                      df.MDR.nonarctic.quad.uni.13,
                                      df.MDR.nonarctic.quad.uni.14
                                      ) 

## Change unique_id into factor type
df.MDR.nonarctic.quad.uni.sp$unique_id <- as.factor(df.MDR.nonarctic.quad.uni.sp$unique_id)



##### Plot
plot.MDR.nonarctic.quad.uni <- ggplot() +
  ## data
  geom_point(data = data,
             aes(x = temp, y = trait, colour = as.factor(unique_id)),
             size = 2) +
  
  ## a separate TPC for each unique group
  geom_line(data = df.MDR.nonarctic.quad.uni.sp, 
            aes(x = temp, y = X50., color = unique_id)) +
  
  ## Overall TPC
  geom_ribbon(data = df.MDR.nonarctic.quad.uni.pop,
              aes(x = temp, ymin = X2.5., ymax = X97.5.),
              fill = "grey",
              alpha = 0.5) +
  geom_line(data = df.MDR.nonarctic.quad.uni.pop,
            aes(x = temp, y = X50.), color = "black", linewidth = 1) +
  
  # Customize the axes and labels
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Development rate (days-1)"
  ) +
  
  # Customize legend
  scale_colour_discrete(name = "",
                        labels = c("Ae. albopictus 1",
                                   "Ae. albopictus 2",
                                   "Ae. albopictus 3",
                                   "Ae. albopictus 4",
                                   "Ae. albopictus 5",
                                   "Ae. albopictus 6",
                                   "Ae. albopictus 7",
                                   "Ae. albopictus 8",
                                   "Ae. albopictus 9",
                                   "Ae. albopictus 10",
                                   "Ae. nigromaculis",
                                   "Ae. sollicitans",
                                   "Ae. triseriatus 1",
                                   "Ae. triseriatus 2")) +
  theme_bw()


plot.MDR.nonarctic.quad.uni

ggsave("figures/MDR.nonarctic.quad.uni.png", plot.MDR.nonarctic.quad.uni,
       width = 10.3, height = 5.6)




## 3B. Fit gamma distributions to non-Arctic TPC parameters --------------------

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
MDR.quad.prior.cf.dists <- data.frame(q = as.vector(MDR.nonarctic.quad.uni$BUGSoutput$sims.list$cf.q),
                                      T0 = as.vector(MDR.nonarctic.quad.uni$BUGSoutput$sims.list$cf.T0),
                                      Tm = as.vector(MDR.nonarctic.quad.uni$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
MDR.quad.prior.gamma.fits <- apply(MDR.quad.prior.cf.dists, 2, 
                                    function(df) fitdistr(df, "gamma")$estimate)


save(MDR.quad.prior.gamma.fits, file = "R-scripts/R2jags-objects/priors/MDR.quad.priors.Rsave")



## 3C. Fit Arctic TPC using data-informed priors -------------------------------

# load("R-scripts/R2jags-objects/priors/MDR.quad.priors.Rsave")

##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)

##### Set data
data <- data.MDR.arctic
hypers <- MDR.quad.prior.gamma.fits

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
                 N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
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

# Save the model as Rdata
save(MDR.arctic.quad.inf, file = "R-scripts/R2jags-objects/all-mods/MDR.arctic.quad.inf.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/all-mods/MDR.arctic.quad.inf.Rdata")


## Diagnostics
##### Examine output
MDR.arctic.quad.inf$BUGSoutput$summary[1:5,]
# mcmcplot(MDR.arctic.quad.inf)

# Extract the DIC for future model comparisons
MDR.arctic.quad.inf$BUGSoutput$DIC

## Plot data + fit
df.MDR.arctic.quad.inf <- data.frame(MDR.arctic.quad.inf$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, X50., sd, X2.5., X97.5.)

head(df.MDR.arctic.quad.inf)

##### Plot
plot.MDR.arctic.quad.inf <- df.MDR.arctic.quad.inf %>%
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "#4363d8", alpha = 0.5) +
  geom_line(aes(y = X50.), color = "blue", linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait), size = 2) +

  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Development rate (days-1)"
  ) +
  theme_bw()

plot.MDR.arctic.quad.inf

ggsave("figures/MDR.arctic.quad.inf.png", plot.MDR.arctic.quad.inf,
       width = 10.3, height = 5.6)




# 4. Compare model fit between Briere and Quadratic models ---------------------

##### Find best fitting model #####
# Add an identifying column in each model output dataframe
df.MDR.arctic.bri.inf <- df.MDR.arctic.bri.inf %>% 
  mutate(type = "briere")

df.MDR.arctic.quad.inf <- df.MDR.arctic.quad.inf %>% 
  mutate(type = "quadratic")

# Combine the two dataframes
df.all <- bind_rows(df.MDR.arctic.bri.inf, df.MDR.arctic.quad.inf)

##### Plot
plot.all <- df.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = mean, color = type), linewidth = 1) +
  geom_point(data = data.MDR.arctic, aes(x = temp, y = trait), size = 2) +
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

ggsave("figures/MDR.bri.quad.png", plot.all, width = 10.3, height = 5.6)


## DIC
MDR.arctic.bri.inf$BUGSoutput$DIC  # This is the best-fitting TPC
MDR.arctic.quad.inf$BUGSoutput$DIC 


##### Plot Arctic vs. non-Arctic TPCs for the best fitting TPC #####
df.MDR.nonarctic.bri.uni.pop <- df.MDR.nonarctic.bri.uni.pop %>% 
  mutate(type = "non-Arctic")

df.MDR.arctic.bri.inf <- df.MDR.arctic.bri.inf %>% 
  mutate(type = "Arctic")

df.arctic.nonarctic <- bind_rows(df.MDR.nonarctic.bri.uni.pop, df.MDR.arctic.bri.inf)

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

ggsave("figures/MDR.arctic.nonarctic.png", plot.arctic.nonarctic, width = 10.3, height = 5.6)


# Save the Briere TPCs in a separate folder
MDR.arctic.mod <- MDR.arctic.bri.inf
MDR.nonarctic.mod <- MDR.nonarctic.bri.uni

## Save the model as Rdata 
save(MDR.arctic.mod, file = "R-scripts/R2jags-objects/best-fitting-mods/MDR.arctic.mod.Rdata")
save(MDR.nonarctic.mod, file = "R-scripts/R2jags-objects/best-fitting-mods/MDR.nonarctic.mod.Rdata")


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
MDR.TPC.analysis <- extractTPC(MDR.arctic.bri.inf, "MDR", Temp.xs)
MDR.arctic.predictions.summary <- MDR.TPC.analysis[[1]]
MDR.arctic.params.summary <- MDR.TPC.analysis[[2]]
MDR.arctic.params.fullposts <- MDR.TPC.analysis[[3]]

write_csv(MDR.arctic.predictions.summary, "data-processed/MDR/MDR.arctic.predictions.summary.csv")
write_csv(MDR.arctic.params.summary, "data-processed/MDR/MDR.arctic.params.summary.csv")
write_csv(MDR.arctic.params.fullposts, "data-processed/MDR/MDR.arctic.params.fullposts.csv")

##### non-Arctic #####
Temp.xs <- seq(0, 45, 0.5)
MDR.TPC.analysis <- extractTPC_raneff(MDR.nonarctic.bri.uni, "MDR", Temp.xs)
MDR.nonarctic.predictions.summary <- MDR.TPC.analysis[[1]]
MDR.nonarctic.params.summary <- MDR.TPC.analysis[[2]]
MDR.nonarctic.params.fullposts <- MDR.TPC.analysis[[3]]

write_csv(MDR.nonarctic.predictions.summary, "data-processed/MDR/MDR.nonarctic.predictions.summary.csv")
write_csv(MDR.nonarctic.params.summary, "data-processed/MDR/MDR.nonarctic.params.summary.csv")
write_csv(MDR.nonarctic.params.fullposts, "data-processed/MDR/MDR.nonarctic.params.fullposts.csv")

