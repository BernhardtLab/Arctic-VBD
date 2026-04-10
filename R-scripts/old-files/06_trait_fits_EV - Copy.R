## Lilian Chan, University of Guelph
## Arctic vector-borne disease transmission suitability model
##
## Purpose: use Bayesian inference (JAGS) to fit TPCs for Egg Viability (EV) for 
## Arctic mosquito species using data from Aedes vexans (McHaffey 1972)
## and from 3 non-Arctic mosquito species (Ae. dorsalis, Ae. nigromaculis, 
## Ae. triseriatus)
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
##        A. Fit EV thermal responses with uniform priors (Arctic species)
##        B. Fit EV thermal responses for priors (non-Arctic species)
##        C. Fit gamma distributions to EV prior thermal responses
##        D. Fit EV thermal responses with data-informed priors (Arctic)
##        E. Plot all three TPCs in the same graph (for comparison)
##
##    3. Fitting TPC (Quadratic)
##        A. Fit EV thermal responses with uniform priors (Arctic)
##        B. Fit EV thermal responses for priors (non-Arctic species)
##        C. Fit gamma distributions to EV prior thermal responses
##        D. Fit EV thermal responses with data-informed priors (Arctic)
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
data <- read_csv("data-processed/TraitData_EV.csv")
unique(data$species)


# Subset data
## Arctic species
data.EV.arctic <- subset(data, type == "Arctic")

## Non-Arctic species
data.EV.nonarctic <- subset(data, type == "non-Arctic")


## Plot raw data
plot.data.EV <- data %>% 
  ggplot() +
  geom_point(aes(x = temp, y = trait, colour = species
  )) +
  labs(y = "Egg viability (%)", x = "Temperature ºC") +
  scale_colour_discrete(name = "species", labels = c("Ae. albopictus",
                                                     "Ae. dorsalis",
                                                     "Ae. nigromaculis",
                                                     "Ae. triseriatus",
                                                     "Ae. vexans"
  )) +
  facet_grid(rows = vars(type)) +
  theme_bw()

plot.data.EV

# ggsave("figures/raw_data/plot.data.EV.png", plot.data.EV, , width = 9.83, height = 6.17)



##########
###### 1. MCMC settings for all models ----
##########

# Number of posterior dist elements = [(ni - nb) / nt] * nc = [(45000 - 5000) / 8] * 3 = 15000
ni <- 45000 # number of iterations in each chain
nb <- 5000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 3 # number of chains


##########
###### 2A. Fit EV thermal responses with uniform priors (Arctic): Briere ----
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
data <- data.EV.arctic


prior <- data.frame(q = c(0, 1),
                    T0 = c(0, 24),
                    Tm = c(25, 50)
)


##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, 
                 N.Temp.xs = N.Temp.xs, prior = prior)

##### Run JAGS
EV.arctic.bri.uni <- jags(
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
save(EV.arctic.bri.uni, file = "R-scripts/R2jags-objects/EV.arctic.bri.uni.Rdata")

# Read the .Rdata
load("R-scripts/R2jags-objects/EV.arctic.bri.uni.Rdata")


## Diagnostics ----
##### Examine output
EV.arctic.bri.uni$BUGSoutput$summary[1:5,]
mcmcplot(EV.arctic.bri.uni)

# Extract the DIC for future model comparisons
EV.arctic.bri.uni$BUGSoutput$DIC

## Plot data + fit ----
df.EV.arctic.bri.uni <- data.frame(EV.arctic.bri.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.EV.arctic.bri.uni)

##### Plot
plot.EV.arctic.bri.uni <- df.EV.arctic.bri.uni %>%
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
  labs(x = expression(paste("Temperature (", degree, "C)")), y = "Egg viability (%)") +
  theme_bw()

plot.EV.arctic.bri.uni

# ggsave("figures/EV.arctic.bri.uni.png", plot.EV.arctic.bri.uni,
#        width = 10.3, height = 5.6)



##########
###### 2B. Fit EV thermal responses (with random effects) for priors (non-Arctic species): Briere ----
##########

##### Temp sequence for derived quantity calculations
# For priors - fewer temps for derived calculations makes it go faster
Temp.xs <- seq(0, 45, 0.5)
N.Temp.xs <-length(Temp.xs)


##### Set data
data <- data.EV.nonarctic

## Create a unique id for each species-study combination
data <- data %>% 
  group_by(species, citation) %>% 
  mutate(unique_id = cur_group_id())

                    
prior <- data.frame(q = c(0, 1),
                    T0 = c(0, 20),
                    Tm = c(20, 45),
                    sigma_q = c(0, 0.001),
                    sigma_T0 = c(0, 10),
                    sigma_Tm = c(0, 10)
)

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
EV.nonarctic.bri.uni.raneff <- jags(
  data = jag.data,
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
# save(EV.nonarctic.bri.uni.raneff, file = "R-scripts/R2jags-objects/EV.nonarctic.bri.uni.raneff.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/EV.nonarctic.bri.uni.raneff.Rdata")


## Diagnostics ----
##### Examine output
EV.nonarctic.bri.uni.raneff$BUGSoutput$summary[1:8,]
mcmcplot(EV.nonarctic.bri.uni.raneff)

# Extract the DIC for future model comparisons
EV.nonarctic.bri.uni.raneff$BUGSoutput$DIC


## Plot data + fit ----
df.EV.nonarctic.bri.uni.raneff <- data.frame(EV.nonarctic.bri.uni.raneff$BUGSoutput$summary)[-(1:8),]

## Extract the model prediction
## Overall curve
df.EV.nonarctic.bri.uni.raneff.pop <- df.EV.nonarctic.bri.uni.raneff %>% 
  filter(grepl("z.trait.mu.pred.pop", rownames(df.EV.nonarctic.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)


## Unique ID 1: Ae. albopictus (Blagrove et al. 2013)
df.EV.nonarctic.bri.uni.1 <- df.EV.nonarctic.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[1,*]"), rownames(df.EV.nonarctic.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 1)

## Unique ID 2: Ae. albopictus (Calado et al. 2002)
df.EV.nonarctic.bri.uni.2 <- df.EV.nonarctic.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[2,*]"), rownames(df.EV.nonarctic.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 2)

## Unique ID 3: Ae. albopictus (Delatte et al 2009)
df.EV.nonarctic.bri.uni.3 <- df.EV.nonarctic.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[3,*]"), rownames(df.EV.nonarctic.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 3)

## Unique ID 4: Ae. albopictus (Li et al 2021)
df.EV.nonarctic.bri.uni.4 <- df.EV.nonarctic.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[4,*]"), rownames(df.EV.nonarctic.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 4)

## Unique ID 5: Ae. albopictus (Monteiro et al 2007)
df.EV.nonarctic.bri.uni.5 <- df.EV.nonarctic.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[5,*]"), rownames(df.EV.nonarctic.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 5)

## Unique ID 6: Ae. albopictus (Zhang et al 2015)
df.EV.nonarctic.bri.uni.6 <- df.EV.nonarctic.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[6,*]"), rownames(df.EV.nonarctic.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 6)

## Unique ID 7: Ae. nigromaculis
df.EV.nonarctic.bri.uni.7 <- df.EV.nonarctic.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[7,*]"), rownames(df.EV.nonarctic.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 7)

## Unique ID 8: Ae. triseriatus
df.EV.nonarctic.bri.uni.8 <- df.EV.nonarctic.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[8,*]"), rownames(df.EV.nonarctic.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 8)





## Combine the model prediciton of all three unique groups into a dataframe
df.EV.nonarctic.bri.uni.raneff.sp <- rbind(df.EV.nonarctic.bri.uni.1,
                                           df.EV.nonarctic.bri.uni.2,
                                           df.EV.nonarctic.bri.uni.3,
                                           df.EV.nonarctic.bri.uni.4,
                                           df.EV.nonarctic.bri.uni.5,
                                           df.EV.nonarctic.bri.uni.6,
                                           df.EV.nonarctic.bri.uni.7,
                                           df.EV.nonarctic.bri.uni.8) 

## Change unique_id into factor type
df.EV.nonarctic.bri.uni.raneff.sp$unique_id <- as.factor(df.EV.nonarctic.bri.uni.raneff.sp$unique_id)


##### Plot
plot.EV.nonarctic.bri.uni.raneff <- ggplot(data = df.EV.nonarctic.bri.uni.raneff.pop, 
                                            aes(x = temp)) +
  ## Overall TPC
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.),
              fill = "grey",
              alpha = 0.5) +
  ## a separate TPC (and credible interval) for each unique group
  # geom_ribbon(data = df.EV.nonarctic.bri.uni.raneff.sp, aes(ymin = X2.5., ymax = X97.5., fill = unique_id),
  #             alpha = 0.5) +
  geom_line(aes(y = mean), color = "black", linewidth = 1) +
  geom_line(data = df.EV.nonarctic.bri.uni.raneff.sp, aes(y = mean, color = unique_id)) +
  geom_point(data = data,
             aes(x = temp, y = trait, colour = as.factor(unique_id)),
             size = 2) +
  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), y = "Egg viability (%)") +
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


plot.EV.nonarctic.bri.uni.raneff

# ggsave("figures/EV.nonarctic.bri.uni.raneff.png", plot.EV.nonarctic.bri.uni.raneff,
#        width = 10.3, height = 5.6)



##########
###### 2C. Fit gamma distributions to EV prior thermal responses: Briere ----
##########

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
EV.arctic.prior.cf.dists <- data.frame(q = as.vector(EV.nonarctic.bri.uni.raneff$BUGSoutput$sims.list$cf.q),
                                        T0 = as.vector(EV.nonarctic.bri.uni.raneff$BUGSoutput$sims.list$cf.T0),
                                        Tm = as.vector(EV.nonarctic.bri.uni.raneff$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
EV.arctic.prior.gamma.fits = apply(EV.arctic.prior.cf.dists, 2, 
                                    function(df) fitdistr(df, "gamma")$estimate)


EV.hypers <- EV.arctic.prior.gamma.fits
# save(EV.hypers, file = "R-scripts/R2jags-objects/EVhypers.bri.Rsave")



##########
###### 2D. Fit EV thermal responses with data-informed priors (Arctic): Briere ----
##########

load("R-scripts/R2jags-objects/EVhypers.bri.Rsave")
EV.arctic.prior.gamma.fits <- EV.hypers


##### Set data
data <- data.EV.arctic
hypers <- EV.arctic.prior.gamma.fits * 0.1


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
# save(EV.arctic.bri.inf, file = "R-scripts/R2jags-objects/EV.arctic.bri.inf.Rdata")

# Read the .Rdata
load("R-scripts/R2jags-objects/EV.arctic.bri.inf.Rdata")


## Diagnostics ----
##### Examine output
EV.arctic.bri.inf$BUGSoutput$summary[1:5,]
mcmcplot(EV.arctic.bri.inf)

# Extract the DIC for future model comparisons
EV.arctic.bri.inf$BUGSoutput$DIC

## Plot data + fit ----
df.EV.arctic.bri.inf <- data.frame(EV.arctic.bri.inf$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.EV.arctic.bri.inf)

##### Plot
plot.EV.arctic.bri.inf <- df.EV.arctic.bri.inf %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "#4363d8", alpha = 0.5) +
  geom_line(aes(y = mean), color = "blue", linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait, colour = species), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Egg viability (%)"
  ) +
  scale_color_discrete(name = element_blank(),
                     labels = c("Ae. vexans")) +
  theme_bw()

plot.EV.arctic.bri.inf

# ggsave("figures/EV.arctic.bri.inf.png", plot.EV.arctic.bri.inf, 
#        width = 10.3, height = 5.6)





##########
###### 2E. Plot all TPCs for Arctic species in the same graph (for comparison): Briere ----
##########

# Add an identifying column in each model output dataframe
df.EV.arctic.bri.uni <- df.EV.arctic.bri.uni %>% 
  mutate(type = "Briere uniform")

df.EV.arctic.bri.inf <- df.EV.arctic.bri.inf %>% 
  mutate(type = "Briere informative")


# Combine the three dataframes
df.all <- rbind(df.EV.arctic.bri.uni, df.EV.arctic.bri.inf)

df.all$type <- factor(df.all$type, levels = c( "Briere uniform", "Briere informative"))


# Plot
plot.all <- df.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = mean, color = type), linewidth = 1) +
  geom_point(data = data.EV.arctic, aes(x = temp, y = trait), size = 2) +
  #geom_point(data = data.EV.nonarctic, aes(x = temp, y = trait), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Egg viability (%)"
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

#ggsave("figures/EV.arctic.bri.all.png", plot.all, width = 10.3, height = 5.6)

EV.arctic.bri.uni$BUGSoutput$DIC
EV.arctic.bri.inf$BUGSoutput$DIC



##########
###### 3A. Fit EV thermal responses with uniform priors (Arctic): Quadratic ----
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
data <- data.EV.arctic

prior <- data.frame(q = c(0, 1),
                    T0 = c(0, 24),
                    Tm = c(26, 50)
)

##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, 
                 N.Temp.xs = N.Temp.xs, prior = prior)

# ##### Run JAGS -----
EV.arctic.quad.uni <- jags(data = jag.data,
                           inits = inits,
                           parameters.to.save = parameters,
                           model.file = "R-scripts/quadprob.txt",
                           n.chains = nc,
                           n.burnin = nb,
                           n.iter = ni,
                           DIC = T,
                           working.directory = getwd()
)

## Save the model as Rdata 
# save(EV.arctic.quad.uni, file = "R-scripts/R2jags-objects/EV.arctic.quad.uni.Rdata")

# Read the .Rdata
load("R-scripts/R2jags-objects/EV.arctic.quad.uni.Rdata")


## Diagnostics ----
##### Examine output
EV.arctic.quad.uni$BUGSoutput$summary[1:5,]
mcmcplot(EV.arctic.quad.uni)

# Extract the DIC for future model comparisons
EV.arctic.quad.uni$BUGSoutput$DIC

## Plot data + fit ----
df.EV.arctic.quad.uni <- data.frame(EV.arctic.quad.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.EV.arctic.quad.uni)

##### Plot
plot.EV.arctic.quad.uni <- df.EV.arctic.quad.uni %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "#4363d8", alpha = 0.5) +
  geom_line(aes(y = mean), color = "blue", linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait, colour = species), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Egg viability (%)"
  ) +
  scale_color_discrete(name = element_blank(),
                       labels = c("Ae. vexans")) +
  theme_bw()

plot.EV.arctic.quad.uni

# ggsave("figures/EV.arctic.quad.uni.png", plot.EV.arctic.quad.uni, 
#        width = 10.3, height = 5.6)


##########
###### 3B. Fit EV thermal responses for priors (non-Arctic species): Quadratic ----
##########

##### Temp sequence for derived quantity calculations
# For priors - fewer temps for derived calculations makes it go faster
Temp.xs <- seq(0, 45, 0.5)
N.Temp.xs <-length(Temp.xs)


##### Set data
data <- data.EV.nonarctic

## Create a unique id for each species-study combination
data <- data %>% 
  group_by(species, citation) %>% 
  mutate(unique_id = cur_group_id())

## Set priors
prior <- data.frame(q = c(0, 1),
                    T0 = c(0, 20),
                    Tm = c(20, 45),
                    sigma_q = c(0, 0.1),
                    sigma_T0 = c(0, 10),
                    sigma_Tm = c(0, 10)
)


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
EV.nonarctic.quad.uni.raneff <- jags(
  data = jag.data,
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
save(EV.nonarctic.quad.uni.raneff, file = "R-scripts/R2jags-objects/EV.nonarctic.quad.uni.raneff.Rdata")

# Read the .Rdata
load("R-scripts/R2jags-objects/EV.nonarctic.quad.uni.raneff.Rdata")


## Diagnostics ----
##### Examine output
EV.nonarctic.quad.uni.raneff$BUGSoutput$summary[1:8,]
mcmcplot(EV.nonarctic.quad.uni.raneff)

# Extract the DIC for future model comparisons
EV.nonarctic.quad.uni.raneff$BUGSoutput$DIC


## Plot data + fit ----
df.EV.nonarctic.quad.uni.raneff <- data.frame(EV.nonarctic.quad.uni.raneff$BUGSoutput$summary)[-(1:8),]

## Extract the model prediction
## Overall curve
df.EV.nonarctic.quad.uni.raneff.pop <- df.EV.nonarctic.quad.uni.raneff %>% 
  filter(grepl("z.trait.mu.pred.pop", rownames(df.EV.nonarctic.quad.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)


## Unique ID 1: Ae. albopictus (Blagrove et al. 2013)
df.EV.nonarctic.quad.uni.1 <- df.EV.nonarctic.quad.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[1,*]"), rownames(df.EV.nonarctic.quad.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 1)

## Unique ID 2: Ae. albopictus (Calado et al. 2002)
df.EV.nonarctic.quad.uni.2 <- df.EV.nonarctic.quad.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[2,*]"), rownames(df.EV.nonarctic.quad.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 2)

## Unique ID 3: Ae. albopictus (Delatte et al 2009)
df.EV.nonarctic.quad.uni.3 <- df.EV.nonarctic.quad.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[3,*]"), rownames(df.EV.nonarctic.quad.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 3)

## Unique ID 4: Ae. albopictus (Li et al 2021)
df.EV.nonarctic.quad.uni.4 <- df.EV.nonarctic.quad.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[4,*]"), rownames(df.EV.nonarctic.quad.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 4)

## Unique ID 5: Ae. albopictus (Monteiro et al 2007)
df.EV.nonarctic.quad.uni.5 <- df.EV.nonarctic.quad.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[5,*]"), rownames(df.EV.nonarctic.quad.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 5)

## Unique ID 6: Ae. albopictus (Zhang et al 2015)
df.EV.nonarctic.quad.uni.6 <- df.EV.nonarctic.quad.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[6,*]"), rownames(df.EV.nonarctic.quad.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 6)

## Unique ID 7: Ae. nigromaculis
df.EV.nonarctic.quad.uni.7 <- df.EV.nonarctic.quad.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[7,*]"), rownames(df.EV.nonarctic.quad.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 7)

## Unique ID 8: Ae. triseriatus
df.EV.nonarctic.quad.uni.8 <- df.EV.nonarctic.quad.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[8,*]"), rownames(df.EV.nonarctic.quad.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 8)





## Combine the model prediciton of all three unique groups into a dataframe
df.EV.nonarctic.quad.uni.raneff.sp <- rbind(df.EV.nonarctic.quad.uni.1,
                                           df.EV.nonarctic.quad.uni.2,
                                           df.EV.nonarctic.quad.uni.3,
                                           df.EV.nonarctic.quad.uni.4,
                                           df.EV.nonarctic.quad.uni.5,
                                           df.EV.nonarctic.quad.uni.6,
                                           df.EV.nonarctic.quad.uni.7,
                                           df.EV.nonarctic.quad.uni.8) 

## Change unique_id into factor type
df.EV.nonarctic.quad.uni.raneff.sp$unique_id <- as.factor(df.EV.nonarctic.quad.uni.raneff.sp$unique_id)


##### Plot
plot.EV.nonarctic.quad.uni.raneff <- ggplot(data = df.EV.nonarctic.quad.uni.raneff.pop, 
                                           aes(x = temp)) +
  ## Overall TPC
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.),
              fill = "grey",
              alpha = 0.5) +
  ## a separate TPC (and credible interval) for each unique group
  # geom_ribbon(data = df.EV.nonarctic.quad.uni.raneff.sp, aes(ymin = X2.5., ymax = X97.5., fill = unique_id),
  #             alpha = 0.5) +
  geom_point(data = data,
             aes(x = temp, y = trait, colour = as.factor(unique_id)),
             size = 2) +
  geom_line(data = df.EV.nonarctic.quad.uni.raneff.sp, aes(y = mean, color = unique_id)) +
  geom_line(aes(y = mean), color = "black", linewidth = 1.5) +

  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), y = "Egg viability (%)") +
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


plot.EV.nonarctic.quad.uni.raneff

# ggsave("figures/EV.nonarctic.quad.uni.raneff.png", plot.EV.nonarctic.quad.uni.raneff,
#        width = 10.3, height = 5.6)


##########
###### 3C. Fit gamma distributions to EV prior thermal responses: Quadratic ----
##########

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
EV.arctic.prior.cf.dists <- data.frame(q = as.vector(EV.nonarctic.quad.uni.raneff$BUGSoutput$sims.list$cf.q),
                                       T0 = as.vector(EV.nonarctic.quad.uni.raneff$BUGSoutput$sims.list$cf.T0),
                                       Tm = as.vector(EV.nonarctic.quad.uni.raneff$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
EV.arctic.prior.gamma.fits = apply(EV.arctic.prior.cf.dists, 2, 
                                    function(df) fitdistr(df, "gamma")$estimate)


EV.hypers <- EV.arctic.prior.gamma.fits
# save(EV.hypers, file = "R-scripts/R2jags-objects/EVhypers.quad.Rsave")


##########
###### 3D. Fit EV thermal responses with data-informed priors (Arctic): Quadratic ----
##########

load("R-scripts/R2jags-objects/EVhypers.quad.Rsave")
EV.arctic.prior.gamma.fits <- EV.hypers

##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)

##### Set data
data <- data.EV.arctic
hypers <- EV.arctic.prior.gamma.fits * 0.1


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
# save(EV.arctic.quad.inf, file = "R-scripts/R2jags-objects/EV.arctic.quad.inf.Rdata")

# Read the .Rdata
load("R-scripts/R2jags-objects/EV.arctic.quad.inf.Rdata")


## Diagnostics ----
##### Examine output
EV.arctic.quad.inf$BUGSoutput$summary[1:5,]
mcmcplot(EV.arctic.quad.inf)

# Extract the DIC for future model comparisons
EV.arctic.quad.inf$BUGSoutput$DIC

## Plot data + fit ----
df.EV.arctic.quad.inf <- data.frame(EV.arctic.quad.inf$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.EV.arctic.quad.inf)

##### Plot
plot.EV.arctic.quad.inf <- df.EV.arctic.quad.inf %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "#4363d8", alpha = 0.5) +
  geom_line(aes(y = mean), color = "blue", linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait, colour = species), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Egg viability (%)"
  ) +
  # Customize legend
  scale_colour_discrete(name = element_blank(),
                        labels = c("Ae. dorsalis",
                                   "Ae. vexanss")) +
  theme_bw()

plot.EV.arctic.quad.inf

# ggsave("figures/EV.arctic.quad.inf.png", plot.EV.arctic.quad.inf,
#        width = 10.3, height = 5.6)




##########
###### 3E. Plot all three TPCs in the same graph (for comparison) ----
##########

# Add an identifying column in each model output dataframe
df.EV.arctic.quad.uni <- df.EV.arctic.quad.uni %>% 
  mutate(type = "Quadratic uniform")


df.EV.arctic.quad.inf <- df.EV.arctic.quad.inf %>% 
  mutate(type = "Quadratic informative")


# Combine the three dataframes
df.all <- rbind(df.EV.arctic.quad.uni, df.EV.arctic.quad.inf)

##### Plot
plot.all <- df.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = mean, color = type), linewidth = 1) +
  geom_point(data = data.EV.arctic, aes(x = temp, y = trait), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Egg viability (%)"
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

# ggsave("figures/EV.arctic.quad.all.png", plot.all, width = 10.3, height = 5.6)


##### Plot all best fitting TPCs for comparison ----

# Combine the three dataframes
df.all <- rbind(#df.EV.arctic.bri.uni, 
                #df.EV.arctic.bri.inf, 
                df.EV.arctic.quad.uni,
                df.EV.arctic.quad.inf)



##### Plot
plot.all <- df.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = mean, color = type), linewidth = 1) +
  geom_point(data = data.EV.arctic, aes(x = temp, y = trait), size = 2) +
  #geom_point(data = data.EV.sierrensis, aes(x = temp, y = trait), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Egg viability (%)"
  ) +
  # Customize the colours
  # scale_fill_jco() +
  # scale_color_jco() +
  # scale_fill_brewer(palette = "Accent") +
  # scale_color_brewer(palette = "Accent") +
  theme_bw()

plot.all

# ggsave("figures/EV.arctic.all.png", plot.all, width = 10.3, height = 5.6)


#### DIC ----
EV.arctic.bri.uni$BUGSoutput$DIC
EV.arctic.bri.inf$BUGSoutput$DIC
EV.arctic.quad.uni$BUGSoutput$DIC
EV.arctic.quad.inf$BUGSoutput$DIC # This is the best fitting TPC


##########
###### 4. Process and save model output for plotting ----
##########

## Analyze TPC model
EV.TPC.analysis <- extractTPC(EV.arctic.quad.inf, "EV", Temp.xs)
EV.predictions.summary <- EV.TPC.analysis[[1]]
EV.params.summary <- EV.TPC.analysis[[2]]
EV.params.fullposts <- EV.TPC.analysis[[3]]

write_csv(EV.predictions.summary, "data-processed/EV.predictions.summary.csv")
write_csv(EV.params.summary, "data-processed/EV.params.summary.csv")
write_csv(EV.params.fullposts, "data-processed/EV.params.fullposts.csv")
