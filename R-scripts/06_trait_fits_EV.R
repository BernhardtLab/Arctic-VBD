## Lilian Chan, University of Guelph
## Arctic vector-borne disease transmission suitability model
##
## Purpose: Fit thermal performance curves (TPCs) for egg viability (EV) using 
## Bayesian inference (JAGS). Arctic species models are fit using data-informed 
## priors derived from non-Arctic species.
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
## data-processed/TraitData_EV.csv - 
##     Synthesized published trait data for EV
##
## Outputs: 
## R-scripts/R2jags-objects/best-fitting-mods/EV.arctic.mod.Rdata - 
##     Best-fitting TPC models for Arctic species
##
## R-scripts/R2jags-objects/best-fitting-mods/EV.nonarctic.mod.Rdata -
##     Best-fitting TPC models for non-Arctic species
##
## data-processed/EV/EV.arctic.predictions.summary.csv -
##     Posterior summary of TPC predictions for Arctic species across temperatures
##
## data-processed/EV/EV.arctic.params.summary.csv -
##     Summary statistics of TPC parameters (Arctic TPC)
##
## data-processed/EV/EV.arctic.params.fullposts.csv -
##     Full posterior distributions for TPC parameters (Arctic TPC)
##
## data-processed/EV/EV.nonarctic.predictions.summary.csv -
##     Posterior summary of TPC predictions for non-Arctic species
##
## data-processed/EV/EV.nonarctic.params.summary.csv -
##     Summary statistics of TPC parameters (non-Arctic TPC)
##
## data-processed/EV/EV.nonarctic.params.fullposts.csv -
##     Full posterior distributions for TPC parameters (non-Arctic TPC)


# 0. Set-up workspace -----------------------------------------------------

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
data.all <- read_csv("data-processed/TraitData_EV.csv")
unique(data.all$species)


# Subset data
## Arctic species
data.EV.arctic <- subset(data.all, type == "Arctic")

## Non-Arctic species
data.EV.nonarctic <- subset(data.all, type == "non-Arctic")


## Plot raw data
plot.data.EV <- data.all %>% 
  ggplot() +
  geom_point(aes(x = temp, y = trait, colour = species
  )) +
  labs(y = "Proportion hatching", x = "Temperature ºC") +
  scale_colour_discrete(name = "species", labels = c("Ae. albopictus",
                                                     "Ae. dorsalis",
                                                     "Ae. nigromaculis",
                                                     "Ae. triseriatus",
                                                     "Ae. vexans"
  )) +
  facet_grid(rows = vars(type)) +
  theme_bw()

plot.data.EV



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
data <- data.EV.arctic

##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp


prior <- data.frame(q = c(0, 1),
                    T0 = c(0, 20),
                    Tm = c(25, 45)
)

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
                 N.Temp.xs = N.Temp.xs, prior = prior)

##### Run JAGS
set.seed(123) # for reproducibility
EV.arctic.bri.uni <- jags(data = jag.data,
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
save(EV.arctic.bri.uni, file = "R-scripts/R2jags-objects/all-mods/EV.arctic.bri.uni.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/all-mods/EV.arctic.bri.uni.Rdata")


## Diagnostics
##### Examine output
EV.arctic.bri.uni$BUGSoutput$summary[c("cf.T0", "cf.Tm", "cf.q", "cf.sigma", "deviance"),]
mcmcplot(EV.arctic.bri.uni, parms = c("cf.T0", "cf.Tm", "cf.q", "cf.sigma", "deviance"))

# Extract the DIC for future model comparisons
EV.arctic.bri.uni$BUGSoutput$DIC

## Plot data + fit
df.EV.arctic.bri.uni <- data.frame(EV.arctic.bri.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)

head(df.EV.arctic.bri.uni)

##### Plot
plot.EV.arctic.bri.uni <- df.EV.arctic.bri.uni %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "#4363d8", alpha = 0.5) +
  geom_line(aes(y = X50.), color = "blue", linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Proportion hatching",
    title = "A) EV, arctic sp., briere, uniform priors"
  ) +
  theme_bw()

plot.EV.arctic.bri.uni

ggsave("figures/EV.arctic.bri.uni.png", plot.EV.arctic.bri.uni,
       width = 10.3, height = 5.6)



## 2B. Fit non-Arctic TPC to generate priors -----------------------------------


##### Temp sequence for derived quantity calculations
# For priors - fewer temps for derived calculations makes it go faster
Temp.xs <- seq(0, 45, 0.5)
N.Temp.xs <-length(Temp.xs)


##### Set data
data <- data.EV.nonarctic

# Since this dataset has contains data from multiple species or multiple studies
# of the same species, we incorporated random effects on each thermal response
# parameter (q, T0, Tm) to addressed non-independence among observations 

## Create a unique id for each species-study combination
data <- data %>% 
  group_by(species, citation) %>% 
  mutate(unique_id = cur_group_id())

                    
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
EV.nonarctic.bri.uni <- jags(data = jag.data,
                                    inits = inits,
                                    parameters.to.save = parameters,
                                    model.file = "R-scripts/briereprob_randeff.txt",
                                    n.thin = nt,
                                    n.chains = nc,
                                    n.burnin = nb,
                                    n.iter = ni,
                                    DIC = T,
                                    working.directory = getwd()
)


## Save the model as Rdata 
save(EV.nonarctic.bri.uni, file = "R-scripts/R2jags-objects/all-mods/EV.nonarctic.bri.uni.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/all-mods/EV.nonarctic.bri.uni.Rdata")


## Diagnostics
##### Examine output
EV.nonarctic.bri.uni$BUGSoutput$summary[c("cf.T0", "cf.Tm", "cf.q", "cf.sigma", "deviance", "sigma_T0", "sigma_Tm", "sigma_q"),]
mcmcplot(EV.nonarctic.bri.uni, parms = c("cf.T0", "cf.Tm", "cf.q", "cf.sigma", "deviance", "sigma_T0", "sigma_Tm", "sigma_q"))

# Extract the DIC for future model comparisons
EV.nonarctic.bri.uni$BUGSoutput$DIC


## Plot data + fit
df.EV.nonarctic.bri.uni <- data.frame(EV.nonarctic.bri.uni$BUGSoutput$summary)[-(1:35),]

## Extract the model prediction
## Overall curve
df.EV.nonarctic.bri.uni.pop <- df.EV.nonarctic.bri.uni %>% 
  filter(grepl("z.trait.mu.pred.pop", rownames(df.EV.nonarctic.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)


## Unique ID 1: Ae. albopictus (Blagrove et al. 2013)
df.EV.nonarctic.bri.uni.1 <- df.EV.nonarctic.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[1,*]"), rownames(df.EV.nonarctic.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 1)

## Unique ID 2: Ae. albopictus (Calado et al. 2002)
df.EV.nonarctic.bri.uni.2 <- df.EV.nonarctic.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[2,*]"), rownames(df.EV.nonarctic.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 2)

## Unique ID 3: Ae. albopictus (Delatte et al 2009)
df.EV.nonarctic.bri.uni.3 <- df.EV.nonarctic.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[3,*]"), rownames(df.EV.nonarctic.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 3)

## Unique ID 4: Ae. albopictus (Li et al 2021)
df.EV.nonarctic.bri.uni.4 <- df.EV.nonarctic.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[4,*]"), rownames(df.EV.nonarctic.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 4)

## Unique ID 5: Ae. albopictus (Monteiro et al 2007)
df.EV.nonarctic.bri.uni.5 <- df.EV.nonarctic.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[5,*]"), rownames(df.EV.nonarctic.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 5)

## Unique ID 6: Ae. albopictus (Zhang et al 2015)
df.EV.nonarctic.bri.uni.6 <- df.EV.nonarctic.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[6,*]"), rownames(df.EV.nonarctic.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 6)

## Unique ID 7: Ae. dorsalis
df.EV.nonarctic.bri.uni.7 <- df.EV.nonarctic.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[7,*]"), rownames(df.EV.nonarctic.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 7)

## Unique ID 8: Ae. nigromaculis
df.EV.nonarctic.bri.uni.8 <- df.EV.nonarctic.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[8,*]"), rownames(df.EV.nonarctic.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 8)

## Unique ID 9: Ae. triseriatus
df.EV.nonarctic.bri.uni.9 <- df.EV.nonarctic.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[9,*]"), rownames(df.EV.nonarctic.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 9)




## Combine the model prediciton of all three unique groups into a dataframe
df.EV.nonarctic.bri.uni.sp <- rbind(df.EV.nonarctic.bri.uni.1,
                                    df.EV.nonarctic.bri.uni.2,
                                    df.EV.nonarctic.bri.uni.3,
                                    df.EV.nonarctic.bri.uni.4,
                                    df.EV.nonarctic.bri.uni.5,
                                    df.EV.nonarctic.bri.uni.6,
                                    df.EV.nonarctic.bri.uni.7,
                                    df.EV.nonarctic.bri.uni.8,
                                    df.EV.nonarctic.bri.uni.9
                                    ) 


## Change unique_id into factor type
df.EV.nonarctic.bri.uni.sp$unique_id <- as.factor(df.EV.nonarctic.bri.uni.sp$unique_id)


##### Plot
plot.EV.nonarctic.bri.uni <- ggplot() +
  ## data
  geom_point(data = data,
             aes(x = temp, y = trait, colour = as.factor(unique_id)),
             size = 2) +
  
  ## a separate TPC for each unique group
  geom_line(data = df.EV.nonarctic.bri.uni.sp, aes(x = temp, y = X50., 
                                                   color = unique_id)) +
  
  ## Overall TPC
  geom_ribbon(data = df.EV.nonarctic.bri.uni.pop,
              aes(x = temp, ymin = X2.5., ymax = X97.5.),
              fill = "grey",
              alpha = 0.5) +
  geom_line(data = df.EV.nonarctic.bri.uni.pop,
            aes(x = temp, y = X50.), color = "black", linewidth = 1) +
  
  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), y = "Proportion hatching") +
  
  # Customize legend
  scale_colour_discrete(name = "",
                        labels = c("Ae. albopictus (Blagrove et al. 2013)",
                                   "Ae. albopictus (Calado et al. 2002)",
                                   "Ae. albopictus (Delatte et al 2009)",
                                   "Ae. albopictus (Li et al 2021)",
                                   "Ae. albopictus (Monteiro et al 2007)",
                                   "Ae. albopictus (Zhang et al 2015)",
                                   "Ae. dorsalis",
                                   "Ae. nigromaculis",
                                   "Ae. triseriatus")) +
  theme_bw()


plot.EV.nonarctic.bri.uni

ggsave("figures/EV.nonarctic.bri.uni.png", plot.EV.nonarctic.bri.uni,
       width = 10.3, height = 5.6)



## 2C. Fit gamma distributions to non-Arctic TPC parameters --------------------


# Get the posterior dists for 3 main parameters (not sigma) into a data frame
EV.bri.prior.cf.dists <- data.frame(q = as.vector(EV.nonarctic.bri.uni$BUGSoutput$sims.list$cf.q),
                                        T0 = as.vector(EV.nonarctic.bri.uni$BUGSoutput$sims.list$cf.T0),
                                        Tm = as.vector(EV.nonarctic.bri.uni$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
EV.bri.prior.gamma.fits = apply(EV.bri.prior.cf.dists, 2, 
                                function(df) fitdistr(df, "gamma")$estimate)

save(EV.bri.prior.gamma.fits, file = "R-scripts/R2jags-objects/priors/EV.bri.priors.Rsave")



## 2D. Fit Arctic TPC using data-informed priors -------------------------------

# load("R-scripts/R2jags-objects/priors/EV.bri.priors.Rsave")

##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)

##### Set data
data <- data.EV.arctic
hypers <- EV.bri.prior.gamma.fits * 0.1

# Gamma prior density for Tm

df.EV.Tm.prior <- tibble(
  temp = seq(0, 45, 0.1),
  density = dgamma(
    temp,
    shape = hypers["shape", "Tm"],
    rate  = hypers["rate", "Tm"]
  )
)

TPC.apex <- max(df.EV.nonarctic.bri.uni.pop$X50., na.rm = TRUE)

df.EV.Tm.prior <- df.EV.Tm.prior %>%
  mutate(
    density.scaled = density / max(density, na.rm = TRUE) * TPC.apex * 0.5
  )

df.EV.T0.prior <- tibble(
  temp = seq(0, 45, 0.1),
  density = dgamma(
    temp,
    shape = hypers["shape", "T0"],
    rate  = hypers["rate", "T0"]
  )
)

df.EV.T0.prior <- df.EV.T0.prior %>%
  mutate(
    density.scaled = density / max(density, na.rm = TRUE) * TPC.apex * 0.5
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
                 N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
set.seed(123) # for reproducibility
EV.arctic.bri.inf <- jags(data = jag.data,
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
save(EV.arctic.bri.inf, file = "R-scripts/R2jags-objects/all-mods/EV.arctic.bri.inf.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/all-mods/EV.arctic.bri.inf.Rdata")


## Diagnostics
##### Examine output
EV.arctic.bri.inf$BUGSoutput$summary[c("cf.T0", "cf.Tm", "cf.q", "cf.sigma", "deviance"),]
mcmcplot(EV.arctic.bri.inf, parms = c("cf.T0", "cf.Tm", "cf.q", "cf.sigma", "deviance"))

# Extract the DIC for future model comparisons
EV.arctic.bri.inf$BUGSoutput$DIC

## Plot data + fit
df.EV.arctic.bri.inf <- data.frame(EV.arctic.bri.inf$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)

head(df.EV.arctic.bri.inf)

##### Plot
plot.EV.arctic.bri.inf <- df.EV.arctic.bri.inf %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "#4363d8", alpha = 0.5) +
  geom_line(aes(y = X50.), color = "blue", linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait), size = 2) +
  
  geom_line(
    data = df.EV.Tm.prior,
    aes(x = temp, y = density.scaled),
    colour = "firebrick3",
    linetype = "dashed",
    linewidth = 1
  ) +
  
  geom_line(
    data = df.EV.T0.prior,
    aes(x = temp, y = density),
    colour = "skyblue",
    linetype = "dashed",
    linewidth = 1
  ) +
  
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Proportion hatching"
  ) +
  theme_bw()

plot.EV.arctic.bri.inf

ggsave("figures/EV.arctic.bri.inf.png", plot.EV.arctic.bri.inf, 
       width = 10.3, height = 5.6)




# 3. Fitting TPC (quadratic) ---------------------------------------------------

## 3A. Fit Arctic TPC using uniform priors -------------------------------------

##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)


##### Set data
data <- data.EV.arctic

##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp


prior <- data.frame(q = c(0, 1),
                    T0 = c(0, 20),
                    Tm = c(25, 45)
                    )

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
                 N.Temp.xs = N.Temp.xs, prior = prior)

##### Run JAGS
set.seed(123) # for reproducibility
EV.arctic.quad.uni <- jags(data = jag.data,
                            inits = inits,
                            parameters.to.save = parameters,
                            model.file = "R-scripts/quadprob.txt",
                            n.thin = nt,
                            n.chains = nc,
                            n.burnin = nb,
                            n.iter = ni,
                            DIC = T,
                            working.directory = getwd()
                            )

## Save the model as Rdata 
save(EV.arctic.quad.uni, file = "R-scripts/R2jags-objects/all-mods/EV.arctic.quad.uni.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/all-mods/EV.arctic.quad.uni.Rdata")


## Diagnostics
##### Examine output
EV.arctic.quad.uni$BUGSoutput$summary[c("cf.T0", "cf.Tm", "cf.q", "cf.sigma", "deviance"),]
mcmcplot(EV.arctic.quad.uni, parms = c("cf.T0", "cf.Tm", "cf.q", "cf.sigma", "deviance"))

# Extract the DIC for future model comparisons
EV.arctic.quad.uni$BUGSoutput$DIC

## Plot data + fit
df.EV.arctic.quad.uni <- data.frame(EV.arctic.quad.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)

head(df.EV.arctic.quad.uni)

##### Plot
plot.EV.arctic.quad.uni <- df.EV.arctic.quad.uni %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "#4363d8", alpha = 0.5) +
  geom_line(aes(y = X50.), color = "blue", linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Proportion hatching",
    title = "B) EV, arctic sp., quadratic, uniform priors"
  ) +
  theme_bw()

plot.EV.arctic.quad.uni

ggsave("figures/EV.arctic.quad.uni.png", plot.EV.arctic.quad.uni,
       width = 10.3, height = 5.6)


## 3B. Fit non-Arctic TPC to generate priors -----------------------------------


##### Temp sequence for derived quantity calculations
# For priors - fewer temps for derived calculations makes it go faster
Temp.xs <- seq(0, 45, 0.5)
N.Temp.xs <-length(Temp.xs)


##### Set data
data <- data.EV.nonarctic

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
EV.nonarctic.quad.uni <- jags(data = jag.data,
                              inits = inits,
                              parameters.to.save = parameters,
                              model.file = "R-scripts/quadprob_randeff.txt",
                              n.thin = nt,
                              n.chains = nc,
                              n.burnin = nb,
                              n.iter = ni,
                              DIC = T,
                              working.directory = getwd()
                              )


## Save the model as Rdata 
save(EV.nonarctic.quad.uni, file = "R-scripts/R2jags-objects/all-mods/EV.nonarctic.quad.uni.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/all-mods/EV.nonarctic.quad.uni.Rdata")


## Diagnostics
##### Examine output
EV.nonarctic.quad.uni$BUGSoutput$summary[c("cf.T0", "cf.Tm", "cf.q", "cf.sigma", "deviance", "sigma_T0", "sigma_Tm", "sigma_q"),]
mcmcplot(EV.nonarctic.quad.uni, parms = c("cf.T0", "cf.Tm", "cf.q", "cf.sigma", "deviance", "sigma_T0", "sigma_Tm", "sigma_q"))


# Extract the DIC for future model comparisons
EV.nonarctic.quad.uni$BUGSoutput$DIC


## Plot data + fit
df.EV.nonarctic.quad.uni <- data.frame(EV.nonarctic.quad.uni$BUGSoutput$summary)

## Extract the model prediction
## Overall curve
df.EV.nonarctic.quad.uni.pop <- df.EV.nonarctic.quad.uni %>% 
  filter(grepl("z.trait.mu.pred.pop", rownames(df.EV.nonarctic.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)


## Unique ID 1: Ae. albopictus (Blagrove et al. 2013)
df.EV.nonarctic.quad.uni.1 <- df.EV.nonarctic.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[1,*]"), rownames(df.EV.nonarctic.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 1)

## Unique ID 2: Ae. albopictus (Calado et al. 2002)
df.EV.nonarctic.quad.uni.2 <- df.EV.nonarctic.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[2,*]"), rownames(df.EV.nonarctic.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 2)

## Unique ID 3: Ae. albopictus (Delatte et al 2009)
df.EV.nonarctic.quad.uni.3 <- df.EV.nonarctic.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[3,*]"), rownames(df.EV.nonarctic.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 3)

## Unique ID 4: Ae. albopictus (Li et al 2021)
df.EV.nonarctic.quad.uni.4 <- df.EV.nonarctic.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[4,*]"), rownames(df.EV.nonarctic.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 4)

## Unique ID 5: Ae. albopictus (Monteiro et al 2007)
df.EV.nonarctic.quad.uni.5 <- df.EV.nonarctic.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[5,*]"), rownames(df.EV.nonarctic.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 5)

## Unique ID 6: Ae. albopictus (Zhang et al 2015)
df.EV.nonarctic.quad.uni.6 <- df.EV.nonarctic.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[6,*]"), rownames(df.EV.nonarctic.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 6)

## Unique ID 7: Ae. dorsalis
df.EV.nonarctic.quad.uni.7 <- df.EV.nonarctic.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[7,*]"), rownames(df.EV.nonarctic.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 7)

## Unique ID 8: Ae. nigromaculis
df.EV.nonarctic.quad.uni.8 <- df.EV.nonarctic.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[8,*]"), rownames(df.EV.nonarctic.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 8)

## Unique ID 9: Ae. triseriatus
df.EV.nonarctic.quad.uni.9 <- df.EV.nonarctic.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[9,*]"), rownames(df.EV.nonarctic.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 9)



## Combine the model prediciton of all three unique groups into a dataframe
df.EV.nonarctic.quad.uni.sp <- rbind(df.EV.nonarctic.quad.uni.1,
                                     df.EV.nonarctic.quad.uni.2,
                                     df.EV.nonarctic.quad.uni.3,
                                     df.EV.nonarctic.quad.uni.4,
                                     df.EV.nonarctic.quad.uni.5,
                                     df.EV.nonarctic.quad.uni.6,
                                     df.EV.nonarctic.quad.uni.7,
                                     df.EV.nonarctic.quad.uni.8) 

## Change unique_id into factor type
df.EV.nonarctic.quad.uni.sp$unique_id <- as.factor(df.EV.nonarctic.quad.uni.sp$unique_id)


##### Plot
plot.EV.nonarctic.quad.uni <- ggplot() +
  ## data
  geom_point(data = data,
             aes(x = temp, y = trait, colour = as.factor(unique_id)),
             size = 2) +
  
  ## a separate TPC for each unique group
  geom_line(data = df.EV.nonarctic.quad.uni.sp, 
            aes(x = temp, y = X50., color = unique_id)) +
  
  ## Overall TPC
  geom_ribbon(data = df.EV.nonarctic.quad.uni.pop,
              aes(x = temp, ymin = X2.5., ymax = X97.5.),
              fill = "grey",
              alpha = 0.5) +
  geom_line(data = df.EV.nonarctic.quad.uni.pop,
            aes(x = temp, y = X50.), color = "black", linewidth = 1) +
  

  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Proportion hatching") +
  # Customize legend
  scale_colour_discrete(name = element_blank(),
                        labels = c("Ae. albopictus (Blagrove et al. 2013)",
                                   "Ae. albopictus (Calado et al. 2002)",
                                   "Ae. albopictus (Delatte et al 2009)",
                                   "Ae. albopictus (Li et al 2021)",
                                   "Ae. albopictus (Monteiro et al 2007)",
                                   "Ae. albopictus (Zhang et al 2015)",
                                   "Ae. dorsalis",
                                   "Ae. nigromaculis",
                                   "Ae. triseriatus")) +
  theme_bw()


plot.EV.nonarctic.quad.uni

ggsave("figures/EV.nonarctic.quad.uni.png", plot.EV.nonarctic.quad.uni,
       width = 10.3, height = 5.6)



## 3C. Fit gamma distributions to non-Arctic TPC parameters --------------------

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
EV.quad.prior.cf.dists <- data.frame(q = as.vector(EV.nonarctic.quad.uni$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(EV.nonarctic.quad.uni$BUGSoutput$sims.list$cf.T0),
                                     Tm = as.vector(EV.nonarctic.quad.uni$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
EV.quad.prior.gamma.fits = apply(EV.quad.prior.cf.dists, 2, 
                                 function(df) fitdistr(df, "gamma")$estimate)


save(EV.quad.prior.gamma.fits, file = "R-scripts/R2jags-objects/priors/EV.quad.priors.Rsave")



## 3D. Fit Arctic TPC using data-informed priors -------------------------------

# load("R-scripts/R2jags-objects/EVhypers.quad.Rsave")

##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)

##### Set data
data <- data.EV.arctic
hypers <- EV.quad.prior.gamma.fits * 0.1

# Gamma prior density for Tm

df.EV.Tm.prior <- tibble(
  temp = seq(0, 45, 0.1),
  density = dgamma(
    temp,
    shape = hypers["shape", "Tm"],
    rate  = hypers["rate", "Tm"]
  )
)

TPC.apex <- max(df.EV.nonarctic.quad.uni.pop$X50., na.rm = TRUE)

df.EV.Tm.prior <- df.EV.Tm.prior %>%
  mutate(
    density.scaled = density / max(density, na.rm = TRUE) * TPC.apex * 0.5
  )

df.EV.T0.prior <- tibble(
  temp = seq(0, 45, 0.1),
  density = dgamma(
    temp,
    shape = hypers["shape", "T0"],
    rate  = hypers["rate", "T0"]
  )
)

df.EV.T0.prior <- df.EV.T0.prior %>%
  mutate(
    density.scaled = density / max(density, na.rm = TRUE) * TPC.apex * 0.5
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
                 N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
set.seed(123) # for reproducibility
EV.arctic.quad.inf <- jags(data = jag.data,
                            inits = inits,
                            parameters.to.save = parameters,
                            model.file = "R-scripts/quadprob_inf.txt",
                            n.thin = nt,
                            n.chains = nc,
                            n.burnin = nb,
                            n.iter = ni,
                            DIC = T,
                            working.directory = getwd()
                           )

## Save the model as Rdata 
save(EV.arctic.quad.inf, file = "R-scripts/R2jags-objects/all-mods/EV.arctic.quad.inf.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/all-mods/EV.arctic.quad.inf.Rdata")


## Diagnostics
##### Examine output
EV.arctic.quad.inf$BUGSoutput$summary[c("cf.T0", "cf.Tm", "cf.q", "cf.sigma", "deviance"),]
mcmcplot(EV.arctic.quad.inf, parms = c("cf.T0", "cf.Tm", "cf.q", "cf.sigma", "deviance"))

# Extract the DIC for future model comparisons
EV.arctic.quad.inf$BUGSoutput$DIC

## Plot data + fit
df.EV.arctic.quad.inf <- data.frame(EV.arctic.quad.inf$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)

head(df.EV.arctic.quad.inf)

##### Plot
plot.EV.arctic.quad.inf <- df.EV.arctic.quad.inf %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "#4363d8", alpha = 0.5) +
  geom_line(aes(y = X50.), color = "blue", linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait), size = 2) +

  geom_line(
    data = df.EV.Tm.prior,
    aes(x = temp, y = density.scaled),
    colour = "firebrick3",
    linetype = "dashed",
    linewidth = 1
  ) +
  
  geom_line(
    data = df.EV.T0.prior,
    aes(x = temp, y = density),
    colour = "skyblue",
    linetype = "dashed",
    linewidth = 1
  ) +
  
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Proportion hatching"
  ) +
  theme_bw()

plot.EV.arctic.quad.inf

ggsave("figures/EV.arctic.quad.inf.png", plot.EV.arctic.quad.inf,
       width = 10.3, height = 5.6)


plot.EV.arctic <- plot_grid(plot.EV.arctic.bri.uni, plot.EV.arctic.quad.uni,
                             plot.EV.arctic.bri.inf, plot.EV.arctic.quad.inf, ncol = 2)

plot.EV.arctic

ggsave("figures/EV.arctic.png", plot.EV.arctic,
       width = 10.3, height = 5.6)

# 4. Compare model fit between Briere and Quadratic models ---------------------

##### Find best fitting model #####
# Add an identifying column in each model output dataframe
df.EV.arctic.bri.inf <- df.EV.arctic.bri.inf %>% 
  mutate(type = "briere")

df.EV.arctic.quad.inf <- df.EV.arctic.quad.inf %>% 
  mutate(type = "quadratic")

# Combine the two dataframes
df.all <- bind_rows(df.EV.arctic.bri.inf, df.EV.arctic.quad.inf)

##### Plot
plot.all <- df.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = X50., color = type), linewidth = 1) +
  geom_point(data = data.EV.arctic, aes(x = temp, y = trait), size = 2) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Proportion hatching"
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

ggsave("figures/EV.bri.quad.png", plot.all, width = 10.3, height = 5.6)


## DIC
EV.arctic.bri.uni$BUGSoutput$DIC 
EV.arctic.bri.inf$BUGSoutput$DIC 
EV.arctic.quad.uni$BUGSoutput$DIC
EV.arctic.quad.inf$BUGSoutput$DIC # This is the best fitting TPC


##### Plot Arctic vs. non-Arctic TPCs for the best fitting TPC #####
df.EV.nonarctic.quad.uni.pop <- df.EV.nonarctic.quad.uni.pop %>% 
  mutate(type = "non-Arctic")

df.EV.arctic.quad.inf <- df.EV.arctic.quad.inf %>% 
  mutate(type = "Arctic")

df.arctic.nonarctic <- bind_rows(df.EV.nonarctic.quad.uni.pop, df.EV.arctic.quad.inf)

plot.arctic.nonarctic <- df.arctic.nonarctic %>% 
  ggplot(aes(x = temp)) +
  geom_point(data = data.all, aes(x = temp, y = trait, colour = type), size = 2) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = X50., color = type), linewidth = 1) +
  
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Proportion hatching"
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

ggsave("figures/EV.arctic.nonarctic.png", plot.arctic.nonarctic, width = 10.3, height = 5.6)


# Save best-fitting TPC in a separate folder
EV.arctic.mod <- EV.arctic.quad.inf
EV.nonarctic.mod <- EV.nonarctic.quad.uni

## Save the model as Rdata 
save(EV.arctic.mod, file = "R-scripts/R2jags-objects/best-fitting-mods/EV.arctic.mod.Rdata")
save(EV.nonarctic.mod, file = "R-scripts/R2jags-objects/best-fitting-mods/EV.nonarctic.mod.Rdata")



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
EV.TPC.analysis <- extractTPC(EV.arctic.quad.inf, "EV", Temp.xs)
EV.arctic.predictions.summary <- EV.TPC.analysis[[1]]
EV.arctic.params.summary <- EV.TPC.analysis[[2]]
EV.arctic.params.fullposts <- EV.TPC.analysis[[3]]

write_csv(EV.arctic.predictions.summary, "data-processed/EV/EV.arctic.predictions.summary.csv")
write_csv(EV.arctic.params.summary, "data-processed/EV/EV.arctic.params.summary.csv")
write_csv(EV.arctic.params.fullposts, "data-processed/EV/EV.arctic.params.fullposts.csv")

##### non-Arctic #####
Temp.xs <- seq(0, 45, 0.5)
EV.TPC.analysis <- extractTPC_raneff(EV.nonarctic.quad.uni, "EV", Temp.xs)
EV.nonarctic.predictions.summary <- EV.TPC.analysis[[1]]
EV.nonarctic.params.summary <- EV.TPC.analysis[[2]]
EV.nonarctic.params.fullposts <- EV.TPC.analysis[[3]]

write_csv(EV.nonarctic.predictions.summary, "data-processed/EV/EV.nonarctic.predictions.summary.csv")
write_csv(EV.nonarctic.params.summary, "data-processed/EV/EV.nonarctic.params.summary.csv")
write_csv(EV.nonarctic.params.fullposts, "data-processed/EV/EV.nonarctic.params.fullposts.csv")



##### Briere model #####
Temp.xs <- seq(0, 45, 0.1)
EV.TPC.analysis <- extractTPC(EV.arctic.bri.inf, "EV", Temp.xs) 
EV.arctic.predictions.summary <- EV.TPC.analysis[[1]]
EV.arctic.params.summary <- EV.TPC.analysis[[2]]
EV.arctic.params.fullposts <- EV.TPC.analysis[[3]]

write_csv(EV.arctic.predictions.summary, "data-processed/supplemental-analysis/briere-only/EV.arctic.predictions.summary.csv")
write_csv(EV.arctic.params.summary, "data-processed/supplemental-analysis/briere-only/EV.arctic.params.summary.csv")
write_csv(EV.arctic.params.fullposts, "data-processed/supplemental-analysis/briere-only/EV.arctic.params.fullposts.csv")

Temp.xs <- seq(0, 45, 0.5)
EV.TPC.analysis <- extractTPC_raneff(EV.nonarctic.bri.uni, "EV", Temp.xs)
EV.nonarctic.predictions.summary <- EV.TPC.analysis[[1]]
EV.nonarctic.params.summary <- EV.TPC.analysis[[2]]
EV.nonarctic.params.fullposts <- EV.TPC.analysis[[3]]

write_csv(EV.nonarctic.predictions.summary, "data-processed/supplemental-analysis/briere-only/EV.nonarctic.predictions.summary.csv")
write_csv(EV.nonarctic.params.summary, "data-processed/supplemental-analysis/briere-only/EV.nonarctic.params.summary.csv")
write_csv(EV.nonarctic.params.fullposts, "data-processed/supplemental-analysis/briere-only/EV.nonarctic.params.fullposts.csv")
