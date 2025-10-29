## Lilian Chan, University of Guelph
## Arctic vector-borne disease transmission suitability model
##
## Purpose: use Bayesian inference (JAGS) to fit TPCs for lifetime egg 
## production (B) for Aedes hexodontus (Barlow 1955), Aedes cinereus, Aedes 
## communis, Aedes impiger, and Aedes punctor (Sommerman 1969)
##     1) with uniform priors; and 
##     2) with data-informed priors from Aedes aegypti data
##        Since I only have Eggs per female per day (EFD) data for these 
##        species, I will fit the TPC for EFD and adult lifespan to obtain 
##        lifetime egg production
##
## Table of content:
##    0. Set-up workspace
##
##    1. MCMC settings for all models
##
##    2. Fitting TPC (Briere)
##        A. Fit EFD and lf thermal responses (with random effects) for priors (non-Arctic species)
##        B. Calculate B for non-Arctic species
##        C. Fit gamma distributions to B prior thermal responses
##        D. Fit B thermal responses with data-informed priors (Arctic species)
##        E. Plot all TPCs for Arctic species in the same graph (for comparison)
##
##    3. Fitting TPC (Quadratic)
##        A. Fit EFD and lf thermal responses (with random effects) for priors (non-Arctic species)
##        B. Calculate B for non-Arctic species
##        C. Fit gamma distributions to B prior thermal responses
##        D. Fit B thermal responses with data-informed priors (Arctic species)
##        E. Plot all TPCs for Arctic species in the same graph (for comparison)


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


sink("R-scripts/briere_T_randeff_B.txt")
cat("
    model{

    ## Priors
    cf.q ~ dunif(prior[1,1], prior[2,1])
    cf.T0 ~ dunif(prior[1,2], prior[2,2])
    cf.Tm ~ dunif(prior[1,3], prior[2,3])
    cf.sigma ~ dunif(0, 1000)
    cf.tau <- 1 / (cf.sigma * cf.sigma)
    
    ## Random effect priors
    sigma_q ~ dunif(prior[1,4], prior[2,4])
    tau_q <- 1 / (sigma_q * sigma_q)
    
    sigma_T0 ~ dunif(prior[1,5], prior[2,5])
    tau_T0 <- 1 / (sigma_T0 * sigma_T0)
    
    sigma_Tm ~ dunif(prior[1,6], prior[2,6])
    tau_Tm <- 1 / (sigma_Tm * sigma_Tm)
    
    ## Random effects for each species-study combination (unique_id)
     
    for (j in 1:Nids) {
    q[j] ~ dnorm(0, tau_q)
    T0[j] ~ dnorm(0, tau_T0)
    Tm[j] ~ dnorm(0, tau_Tm)
    }
		
    ## Likelihood
    for(i in 1:N.obs){
    trait.mu[i] <- (cf.q + q[unique.id[i]]) * temp[i] * (temp[i] - (cf.T0 + T0[unique.id[i]])) * sqrt(((cf.Tm + Tm[unique.id[i]]) - temp[i]) * ((cf.Tm + Tm[unique.id[i]]) > temp[i])) * ((cf.T0 + T0[unique.id[i]]) < temp[i])
    trait[i] ~ dnorm(trait.mu[i], cf.tau)T(0,)
    }

    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred.pop[i] <- cf.q * Temp.xs[i] * (Temp.xs[i] - cf.T0) * sqrt((cf.Tm - Temp.xs[i]) * (cf.Tm > Temp.xs[i])) * (cf.T0 < Temp.xs[i])}
    
    for (j in 1:Nids) {
      for(i in 1:N.Temp.xs){
        z.trait.mu.pred.id[j,i] <- (cf.q + q[j]) * Temp.xs[i] * (Temp.xs[i] - (cf.T0 + T0[j])) * sqrt(((cf.Tm + Tm[j]) - Temp.xs[i]) * ((cf.Tm + Tm[j]) > Temp.xs[i])) * ((cf.T0 + T0[j]) < Temp.xs[i])
      }
    }
    
    } # close model
    ",fill=T)
sink()


## Set priors
prior <- data.frame(q = c(0, 0.1),
                    T0 = c(0, 20),
                    Tm = c(30, 45),
                    sigma_q = c(0, 0.0001),
                    sigma_T0 = c(0, 10),
                    sigma_Tm = c(0, 10)
)

##########
###### 2A. Fit EFD and lf thermal responses (with random effects) for priors (non-Arctic species): Briere ----
##########

##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)

## EFD ----
##### Set data
data <- data.B.nonarctic %>% 
  filter(trait_name == "EFD")


## Create a unique id for each species-study combination
data <- data %>% 
  group_by(species, citation) %>% 
  mutate(unique_id = cur_group_id())


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
  model.file = "R-scripts/briere_T_randeff_B.txt",
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


## lf ----

# ##### Set data
# data <- data.B.nonarctic %>% 
#   filter(trait_name == "lf")
# 
# ##### Organize data for JAGS
# trait <- data$trait
# N.obs <- length(trait)
# temp <- data$temp
# 
# ##### define data for JAGS in a list object
# jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)
# 
# ##### Run JAGS
# lf.nonarctic.bri.uni <- jags(
#   data = jag.data,
#   inits = inits,
#   parameters.to.save = parameters,
#   model.file = "R-scripts/briere_T_B.txt",
#   n.thin = nt,
#   n.chains = nc,
#   n.burnin = nb,
#   n.iter = ni,
#   DIC = T,
#   working.directory = getwd()
# )
# 
# ## Save the model as Rdata 
# # save(lf.nonarctic.bri.uni, file = "R-scripts/R2jags-objects/lf.nonarctic.bri.uni.Rdata")
# 
# # Read the .Rdata
# # load("R-scripts/R2jags-objects/lf.nonarctic.bri.uni.Rdata")
# 
# 
# ## Diagnostics ----
# ##### Examine output
# lf.nonarctic.bri.uni$BUGSoutput$summary[1:5,]
# mcmcplot(lf.nonarctic.bri.uni)
# 
# # Extract the DIC for future model comparisons
# lf.nonarctic.bri.uni$BUGSoutput$DIC
# 
# ## Plot data + fit ----
# df.lf.nonarctic.bri.uni <- data.frame(lf.nonarctic.bri.uni$BUGSoutput$summary)[-(1:5),] %>% 
#   mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
#   dplyr::select(temp, mean, sd, X2.5., X97.5.)
# 
# head(df.lf.nonarctic.bri.uni)
# 
# ##### Plot
# plot.lf.nonarctic.bri.uni <- df.lf.nonarctic.bri.uni %>%
#   ggplot(aes(x = temp)) +
#   geom_ribbon(aes(ymin = X2.5., ymax = X97.5.),
#               fill = "#4363d8",
#               alpha = 0.5) +
#   geom_line(aes(y = mean), color = "blue", linewidth = 1) +
#   geom_point(data = data,
#              aes(x = temp, y = trait, colour = species),
#              size = 2) +
#   # Customize the axes and labels
#   #scale_x_continuous(limits = c(0, 41)) +
#   #scale_y_continuous(limits = c(-0.005, 0.19)) +
#   labs(x = expression(paste("Temperature (", degree, "C)")), y = "Adult mosquito lifespan (days)") +
#   # Customize legend
#   # scale_color_discrete(name = "Species",
#   #                      labels = c("Ae. cinereus",
#   #                                 "Ae. communis",
#   #                                 "Ae. impiger",
#   #                                 "Ae. punctor",
#   #                                 "Ae. vexans")) +
#   theme_bw()
# 
# plot.lf.nonarctic.bri.uni
# 
# # ggsave("figures/B.lf.nonarctic.bri.uni.png", plot.lf.nonarctic.bri.uni,
# #        width = 10.3, height = 5.6)


## Random effects START ----

## Set priors
prior <- data.frame(q = c(0, 0.01),
                    T0 = c(0, 20),
                    Tm = c(30, 45),
                    sigma_q = c(0, 0.0001),
                    sigma_T0 = c(0, 10),
                    sigma_Tm = c(0, 10)
)

##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)

##### Set data
data <- data.B.nonarctic %>% 
  filter(trait_name == "lf")


## Create a unique id for each species-study combination
data <- data %>% 
  group_by(species, citation) %>% 
  mutate(unique_id = cur_group_id())


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
# load("R-scripts/R2jags-objects/B.lf.nonarctic.bri.uni.raneff.Rdata")


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
  # geom_ribbon(data = df.B.lf.nonarctic.bri.uni.raneff.sp, aes(ymin = X2.5., ymax = X97.5., fill = unique_id),
  #             alpha = 0.5) +
  geom_line(data = df.B.lf.nonarctic.bri.uni.raneff.sp, aes(y = mean, color = unique_id)) +
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


plot.B.lf.nonarctic.bri.uni.raneff

# ggsave("figures/B.lf.nonarctic.bri.uni.raneff.png", plot.B.lf.nonarctic.bri.uni.raneff,
#        width = 10.3, height = 5.6)


##########
###### 2B. Calculate B for non-Arctic species ----
##########

## Pull out the derived/predicted values:
B.EFD.nonarctic.bri.uni.raneff.pred <- B.EFD.nonarctic.bri.uni.raneff$BUGSoutput$sims.list$z.trait.mu.pred.pop
B.lf.nonarctic.bri.uni.raneff.pred <- B.lf.nonarctic.bri.uni.raneff$BUGSoutput$sims.list$z.trait.mu.pred.pop



## Specify function to calculate mean & quantiles
calcPostQuants = function(input, grad.xs) {
  
  # Get length of gradient
  N.grad.xs <- length(grad.xs)
  
  # Create output dataframe
  output.df <- data.frame("mean" = numeric(N.Temp.xs), "median" = numeric(N.Temp.xs), 
                          "lowerCI" = numeric(N.Temp.xs), "upperCI" = numeric(N.Temp.xs), 
                          "lowerQuartile" = numeric(N.Temp.xs), "upperQuartile" = numeric(N.Temp.xs), temp = grad.xs)
  
  # Calculate mean & quantiles
  for(i in 1:N.grad.xs){
    output.df$mean[i] <- mean(input[ ,i])
    output.df$median[i] <- quantile(input[ ,i], 0.5, na.rm = TRUE)
    output.df$lowerCI[i] <- quantile(input[ ,i], 0.025, na.rm = TRUE)
    output.df$upperCI[i] <- quantile(input[ ,i], 0.975, na.rm = TRUE)
    output.df$lowerQuartile[i] <- quantile(input[ ,i], 0.25, na.rm = TRUE)
    output.df$upperQuartile[i] <- quantile(input[ ,i], 0.75, na.rm = TRUE)
  }
  
  output.df # return output
  
}


## Calculate B for Ae. aegypti
B.nonarctic.bri.uni.calc <- B.EFD.nonarctic.bri.uni.raneff.pred * B.lf.nonarctic.bri.uni.raneff.pred

## Get the mean, median, CI, and upper and lower quartile
df.B.nonarctic.bri.uni <- calcPostQuants(B.nonarctic.bri.uni.calc, Temp.xs)


## plot 
##### Plot
plot.B.nonarctic.bri.uni <- df.B.nonarctic.bri.uni %>%
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

plot.B.nonarctic.bri.uni

# ggsave("figures/B.nonarctic.bri.uni.png", plot.B.nonarctic.bri.uni,
#        width = 10.3, height = 5.6)


##########
###### 2C. Fit gamma distributions to B prior thermal responses: Briere ----
##########
## First fit a Briere to the mean
data <- df.B.nonarctic.bri.uni[, c("temp", "mean")]

colnames(data) <- c("temp", "trait")

## Include the arctic data
data <-rbind(data, data.B.arctic[, c("temp", "trait")])

ggplot(data = data) +
  geom_point(aes(x = temp, y = trait)) +
  theme_bw()

##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp


##### inits Function
inits<-function(){list(
  cf.q = 0.01,
  cf.Tm = 35,
  cf.T0 = 5,
  cf.sigma = rlnorm(1))}

##### Parameters to Estimate
parameters <- c("cf.q", "cf.T0", "cf.Tm","cf.sigma", "z.trait.mu.pred")


##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
B.nonarctic.bri.uni <- jags(
  data = jag.data,
  inits = inits,
  parameters.to.save = parameters,
  model.file = "R-scripts/briere_T_B.txt",
  n.thin = nt,
  n.chains = nc,
  n.burnin = nb,
  n.iter = ni,
  DIC = T,
  working.directory = getwd()
)


## Save the model as Rdata
# save(B.nonarctic.bri.uni, file = "R-scripts/R2jags-objects/B.nonarctic.bri.uni.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/B.nonarctic.bri.uni.Rdata")


## Diagnostics ----
##### Examine output
B.nonarctic.bri.uni$BUGSoutput$summary[1:5,]
mcmcplot(B.nonarctic.bri.uni)

# Extract the DIC for future model comparisons
B.nonarctic.bri.uni$BUGSoutput$DIC

## Plot data + fit ----
df.B.nonarctic.bri.uni <- data.frame(B.nonarctic.bri.uni$BUGSoutput$summary)[-(1:5),] %>%
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.B.nonarctic.bri.uni)

##### Plot
plot.B.nonarctic.bri.uni <- df.B.nonarctic.bri.uni%>%
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
  labs(x = expression(paste("Temperature (", degree, "C)")), y = "Lifetime egg production") +
  theme_bw()

plot.B.nonarctic.bri.uni

# ggsave("figures/B.nonarctic.bri.uni.png", plot.B.nonarctic.bri.uni,
#        width = 10.3, height = 5.6)




# Get the posterior dists for 3 main parameters (not sigma) into a data frame
B.arctic.prior.cf.dists <- data.frame(q = as.vector(B.nonarctic.bri.uni$BUGSoutput$sims.list$cf.q),
                                       T0 = as.vector(B.nonarctic.bri.uni$BUGSoutput$sims.list$cf.T0),
                                       Tm = as.vector(B.nonarctic.bri.uni$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
B.arctic.prior.gamma.fits = apply(B.arctic.prior.cf.dists, 2, 
                                    function(df) fitdistr(df, "gamma")$estimate)


B.hypers <- B.arctic.prior.gamma.fits
# save(B.hypers, file = "R-scripts/R2jags-objects/Bhypers.bri.Rsave")



##########
###### 2D. Fit B thermal responses with data-informed priors (Arctic species): Briere ----
##########

load("R-scripts/R2jags-objects/Bhypers.bri.Rsave")
B.arctic.prior.gamma.fits <- B.hypers


##### Set data
data <- data.B.arctic
hypers <- B.arctic.prior.gamma.fits * 0.1

##### No random effect ----
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
B.arctic.bri.inf <- jags(data = jag.data,
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
# save(B.arctic.bri.inf, file = "R-scripts/R2jags-objects/B.arctic.bri.inf.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/B.arctic.bri.inf.Rdata")


## Diagnostics ----
##### Examine output
B.arctic.bri.inf$BUGSoutput$summary[1:5,]
mcmcplot(B.arctic.bri.inf)

# Extract the DIC for future model comparisons
B.arctic.bri.inf$BUGSoutput$DIC

## Plot data + fit ----
df.B.arctic.bri.inf <- data.frame(B.arctic.bri.inf$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.B.arctic.bri.inf)

##### Plot
plot.B.arctic.bri.inf <- df.B.arctic.bri.inf %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "#4363d8", alpha = 0.5) +
  geom_line(aes(y = mean), color = "blue", linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait, color = species), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Mosquito adult lifespan (days)"
  ) +
  # Customize legend
  scale_color_discrete(name = "Species",
                       labels = c("Ae. cinereus",
                                  "Ae. communis",
                                  "Ae. impiger",
                                  "Ae. punctor",
                                  "Ae. vexans")) +
  theme_bw()

plot.B.arctic.bri.inf

# ggsave("figures/B.arctic.bri.inf.png", plot.B.arctic.bri.inf,
#        width = 10.3, height = 5.6)



##########
###### 2E. Plot all TPCs for Arctic species in the same graph (for comparison): Briere ----
##########

# Add an identifying column in each model output dataframe
df.B.arctic.bri.uni <- df.B.arctic.bri.uni %>% 
  mutate(type = "Briere uniform")

df.B.arctic.bri.inf <- df.B.arctic.bri.inf %>% 
  mutate(type = "Briere informative")


# Combine the three dataframes
df.all <- rbind(df.B.arctic.bri.uni, df.B.arctic.bri.inf)

df.all$type <- factor(df.all$type, levels = c( "Briere uniform", "Briere informative"))


# Plot
plot.all <- df.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = mean, color = type), linewidth = 1) +
  geom_point(data = data.B.arctic, aes(x = temp, y = trait), size = 2) +
  #geom_point(data = data.B.nonarctic, aes(x = temp, y = trait), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Mosquito adult lifespan (days)"
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

# ggsave("figures/B.arctic.bri.all.png", plot.all, width = 10.3, height = 5.6)

B.arctic.bri.uni$BUGSoutput$DIC
B.arctic.bri.inf$BUGSoutput$DIC



##########
###### 3A. Fit EFD and lf thermal responses (with random effects) for priors (non-Arctic species): Quadratic ----
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

## EFD ----
##### Set data
data <- data.B.nonarctic %>% 
  filter(trait_name == "EFD")



##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
EFD.nonarctic.quad.uni <- jags(
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
# save(EFD.nonarctic.quad.uni, file = "R-scripts/R2jags-objects/EFD.nonarctic.quad.uni.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/EFD.nonarctic.quad.uni.Rdata")


## Diagnostics ----
##### Examine output
EFD.nonarctic.quad.uni$BUGSoutput$summary[1:5,]
mcmcplot(EFD.nonarctic.quad.uni)

# Extract the DIC for future model comparisons
EFD.nonarctic.quad.uni$BUGSoutput$DIC

## Plot data + fit ----
df.EFD.nonarctic.quad.uni <- data.frame(EFD.nonarctic.quad.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.EFD.nonarctic.quad.uni)

##### Plot
plot.EFD.nonarctic.quad.uni <- df.EFD.nonarctic.quad.uni %>%
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
  labs(x = expression(paste("Temperature (", degree, "C)")), y = "Eggs per female per day") +
  # Customize legend
  # scale_color_discrete(name = "Species",
  #                      labels = c("Ae. cinereus",
  #                                 "Ae. communis",
  #                                 "Ae. impiger",
  #                                 "Ae. punctor",
  #                                 "Ae. vexans")) +
  theme_bw()

plot.EFD.nonarctic.quad.uni

# ggsave("figures/B.EFD.nonarctic.quad.uni.png", plot.EFD.nonarctic.quad.uni,
#        width = 10.3, height = 5.6)


## lf ----

##### Set data
data <- data.B.nonarctic %>% 
  filter(trait_name == "lf")

##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
lf.nonarctic.quad.uni <- jags(
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
# save(lf.nonarctic.quad.uni, file = "R-scripts/R2jags-objects/lf.nonarctic.quad.uni.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/lf.nonarctic.quad.uni.Rdata")


## Diagnostics ----
##### Examine output
lf.nonarctic.quad.uni$BUGSoutput$summary[1:5,]
mcmcplot(lf.nonarctic.quad.uni)

# Extract the DIC for future model comparisons
lf.nonarctic.quad.uni$BUGSoutput$DIC

## Plot data + fit ----
df.lf.nonarctic.quad.uni <- data.frame(lf.nonarctic.quad.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.lf.nonarctic.quad.uni)

##### Plot
plot.lf.nonarctic.quad.uni <- df.lf.nonarctic.quad.uni %>%
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
  labs(x = expression(paste("Temperature (", degree, "C)")), y = "Adult mosquito lifespan (days)") +
  # Customize legend
  # scale_color_discrete(name = "Species",
  #                      labels = c("Ae. cinereus",
  #                                 "Ae. communis",
  #                                 "Ae. impiger",
  #                                 "Ae. punctor",
  #                                 "Ae. vexans")) +
  theme_bw()

plot.lf.nonarctic.quad.uni

# ggsave("figures/B.lf.nonarctic.quad.uni.png", plot.lf.nonarctic.quad.uni,
#        width = 10.3, height = 5.6)

##########
###### 3B. Calculate B for non-Arctic species ----
##########

## Pull out the derived/predicted values:
EFD.nonarctic.quad.uni.pred <- EFD.nonarctic.quad.uni$BUGSoutput$sims.list$z.trait.mu.pred
lf.nonarctic.quad.uni.pred <- lf.nonarctic.quad.uni$BUGSoutput$sims.list$z.trait.mu.pred



## Specify function to calculate mean & quantiles
calcPostQuants = function(input, grad.xs) {
  
  # Get length of gradient
  N.grad.xs <- length(grad.xs)
  
  # Create output dataframe
  output.df <- data.frame("mean" = numeric(N.Temp.xs), "median" = numeric(N.Temp.xs), 
                          "lowerCI" = numeric(N.Temp.xs), "upperCI" = numeric(N.Temp.xs), 
                          "lowerQuartile" = numeric(N.Temp.xs), "upperQuartile" = numeric(N.Temp.xs), temp = grad.xs)
  
  # Calculate mean & quantiles
  for(i in 1:N.grad.xs){
    output.df$mean[i] <- mean(input[ ,i])
    output.df$median[i] <- quantile(input[ ,i], 0.5, na.rm = TRUE)
    output.df$lowerCI[i] <- quantile(input[ ,i], 0.025, na.rm = TRUE)
    output.df$upperCI[i] <- quantile(input[ ,i], 0.975, na.rm = TRUE)
    output.df$lowerQuartile[i] <- quantile(input[ ,i], 0.25, na.rm = TRUE)
    output.df$upperQuartile[i] <- quantile(input[ ,i], 0.75, na.rm = TRUE)
  }
  
  output.df # return output
  
}


## Calculate B for Ae. aegypti
B.nonarctic.quad.uni.calc <- EFD.nonarctic.quad.uni.pred * lf.nonarctic.quad.uni.pred

## Get the mean, median, CI, and upper and lower quartile
df.B.nonarctic.quad.uni <- calcPostQuants(B.nonarctic.quad.uni.calc, Temp.xs)


## plot 
##### Plot
plot.B.nonarctic.quad.uni <- df.B.nonarctic.quad.uni %>%
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

plot.B.nonarctic.quad.uni

# ggsave("figures/B.nonarctic.quad.uni.png", plot.B.nonarctic.quad.uni,
#        width = 10.3, height = 5.6)


##########
###### 3C. Fit gamma distributions to B prior thermal responses: Quadratic ----
##########

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
B.arctic.prior.cf.dists <- data.frame(q = as.vector(B.nonarctic.quad.uni$BUGSoutput$sims.list$cf.q),
                                        T0 = as.vector(B.nonarctic.quad.uni$BUGSoutput$sims.list$cf.T0),
                                        Tm = as.vector(B.nonarctic.quad.uni$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
B.arctic.prior.gamma.fits = apply(B.arctic.prior.cf.dists, 2, 
                                    function(df) fitdistr(df, "gamma")$estimate)


B.hypers <- B.arctic.prior.gamma.fits
save(B.hypers, file = "R-scripts/R2jags-objects/Bhypers.quad.Rsave")



##########
###### 3D. Fit B thermal responses with data-informed priors (Arctic species): Quadratic ----
##########

load("R-scripts/R2jags-objects/Bhypers.quad.Rsave")
B.arctic.prior.gamma.fits <- B.hypers

##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)

##### Set data
data <- data.B.arctic
hypers <- B.arctic.prior.gamma.fits * 0.1


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
B.arctic.quad.inf <- jags(data = jag.data,
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
# save(B.arctic.quad.inf, file = "R-scripts/R2jags-objects/B.arctic.quad.inf.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/B.arctic.quad.inf.Rdata")


## Diagnostics ----
##### Examine output
B.arctic.quad.inf$BUGSoutput$summary[1:5,]
mcmcplot(B.arctic.quad.inf)

# Extract the DIC for future model comparisons
B.arctic.quad.inf$BUGSoutput$DIC

## Plot data + fit ----
df.B.arctic.quad.inf <- data.frame(B.arctic.quad.inf$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.B.arctic.quad.inf)

##### Plot
plot.B.arctic.quad.inf <- df.B.arctic.quad.inf %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "pink", alpha = 0.5) +
  geom_line(aes(y = mean), color = "red", linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Mosquito adult lifespan (days)"
  ) +
  theme_bw()

plot.B.arctic.quad.inf

# ggsave("figures/B.arctic.quad.inf.png", plot.B.arctic.quad.inf, 
#        width = 10.3, height = 5.6)


##        A. Fit EFD and lf thermal responses (with random effects) for priors (non-Arctic species)
##        B. Calculate B for non-Arctic species
##        C. Fit gamma distributions to B prior thermal responses
##        D. 
##        E. Plot all TPCs in the same graph (for comparison)

##########
###### 3E. Plot all TPCs in the same graph (for comparison) ----
##########

# Add an identifying column in each model output dataframe
df.B.arctic.quad.uni <- df.B.arctic.quad.uni %>% 
  mutate(type = "Quadratic uniform")


df.B.arctic.quad.inf <- df.B.arctic.quad.inf %>% 
  mutate(type = "Quadratic informative")


# Combine the three dataframes
df.all <- rbind(df.B.arctic.quad.uni, df.B.arctic.quad.inf)

##### Plot
plot.all <- df.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = mean, color = type), linewidth = 1) +
  geom_point(data = data.B.arctic, aes(x = temp, y = trait), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Mosquito adult lifespan (days)"
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

# ggsave("figures/B.arctic.quad.all.png", plot.all, width = 10.3, height = 5.6)


##### Plot all best fitting TPCs for comparison ----

#### DIC ----
B.arctic.bri.uni$BUGSoutput$DIC
B.arctic.bri.inf$BUGSoutput$DIC
B.arctic.quad.uni$BUGSoutput$DIC
B.arctic.quad.inf$BUGSoutput$DIC

# Combine the three dataframes
df.all <- rbind(df.B.arctic.bri.uni, 
                df.B.arctic.bri.inf, 
                df.B.arctic.quad.uni,
                df.B.arctic.quad.inf)



##### Plot
plot.all <- df.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = mean, color = type), linewidth = 1) +
  geom_point(data = data.B.arctic, aes(x = temp, y = trait), size = 2) +
  #geom_point(data = data.B.sierrensis, aes(x = temp, y = trait), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Mosquito adult lifespan (days)"
  ) +
  # Customize the colours
  scale_fill_jco() +
  scale_color_jco() +
  # scale_fill_brewer(palette = "Accent") +
  # scale_color_brewer(palette = "Accent") +
  theme_bw()

plot.all

# ggsave("figures/B.arctic.arctic.all.png", plot.all, width = 10.3, height = 5.6)



