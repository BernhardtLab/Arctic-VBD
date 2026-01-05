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
##    2. Fit EFGC thermal responses (with random effects)
##        A. Briere
##        B. Quadratic
##        C. Compare the two TPC fits
##
##    3. Process and save model output for plotting


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
data <- read_csv("data-processed/TraitData_EFGC.csv")
unique(data$species)


# Subset data
data.EFGC.nonarctic <- data

## Ae. albopictus
data.EFGC.aalb <- subset(data, species == c("albopictus"))

## Cx. pipiens
data.EFGC.cpip <- subset(data, species == c("pipiens"))



## Plot raw data
plot.data.EFGC <- data %>% 
  ggplot() +
  geom_point(aes(x = temp, y = trait, colour = species, shape = citation)) +
  # xlim(c(10,35)) +
  labs(y = "Eggs per female per gonotrophic cycles", x = "Temperature ºC") +
  scale_colour_discrete(name = "Species", labels = c("Ae albopictus",
                                                     "Cx. pipiens"
  )) +
  scale_shape_discrete(name = "Citation", labels = c("Delatte 2009",
                                                     "Ezeakacha 2015",
                                                     "Ju-lin 2017",
                                                     "Yee 2016"
  )) +
  theme_bw()



plot.data.EFGC


# ggsave("figures/raw_data/plot.data.EFGC.png", plot.data.EFGC, width = 9.83, height = 6.17)


##########
###### 1. MCMC settings for all models ----
##########

# Number of posterior dist elements = [(ni - nb) / nt] * nc = [(45000 - 5000) / 8] * 3 = 15000
ni <- 45000 # number of iterations in each chain
nb <- 5000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 3 # number of chains




##########
###### 2A. Fit EFGC thermal responses (Cx. Pipiens): Briere ----
##########

##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)


#### Set data
data <- data.EFGC.cpip


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
EFGC.cpip.bri.uni <- jags(data = jag.data,
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
# save(EFGC.cpip.bri.uni, file = "R-scripts/R2jags-objects/EFGC.cpip.bri.uni.Rdata")

# Read the .Rdata
load("R-scripts/R2jags-objects/EFGC.cpip.bri.uni.Rdata")

## Diagnostics ----
##### Examine output
EFGC.cpip.bri.uni$BUGSoutput$summary[1:5,]
mcmcplot(EFGC.cpip.bri.uni)

# Extract the DIC for future model comparisons
EFGC.cpip.bri.uni$BUGSoutput$DIC

## Plot data + fit ----
df.EFGC.cpip.bri.uni <- data.frame(EFGC.cpip.bri.uni$BUGSoutput$summary)[-(1:5),]; %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.EFGC.cpip.bri.uni)

##### Plot
plot.EFGC.cpip.bri.uni <- df.EFGC.cpip.bri.uni %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "#4363d8", alpha = 0.5) +
  geom_line(aes(y = mean), color = "blue", linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Eggs per female per gonotrophic cycle"
  ) +
  theme_bw()


plot.EFGC.cpip.bri.uni

# ggsave("figures/EFGC.cpip.bri.uni.png", plot.EFGC.cpip.bri.uni, 
#        width = 10.3, height = 5.6)


##########
###### 2B. Fit EFGC thermal responses (Cx. pipiens): Quadratic ----
##########

##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)


#### Set data
data <- data.EFGC.cpip


## Set priors
prior <- data.frame(q = c(0, 2),
                    T0 = c(0, 15),
                    Tm = c(30, 45)
)


##### inits Function
inits<-function(){list(
  cf.q = 0.5,
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
EFGC.cpip.quad.uni <- jags(data = jag.data,
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
# save(EFGC.cpip.quad.uni, file = "R-scripts/R2jags-objects/EFGC.cpip.quad.uni.Rdata")

# Read the .Rdata
load("R-scripts/R2jags-objects/EFGC.cpip.quad.uni.Rdata")

## Diagnostics ----
##### Examine output
EFGC.cpip.quad.uni$BUGSoutput$summary[1:5,]
mcmcplot(EFGC.cpip.quad.uni)

# Extract the DIC for future model comparisons
EFGC.cpip.quad.uni$BUGSoutput$DIC

## Plot data + fit ----
df.EFGC.cpip.quad.uni <- data.frame(EFGC.cpip.quad.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.EFGC.cpip.quad.uni)

##### Plot
plot.EFGC.cpip.quad.uni <- df.EFGC.cpip.quad.uni %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "#4363d8", alpha = 0.5) +
  geom_line(aes(y = mean), color = "blue", linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait, colour = as.factor(unique_id)), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Eggs per female per gonotrophic cycle"
  ) +
  # Customize legend
  scale_colour_discrete(name = element_blank(),
                        labels = c("Ae. albopictus (Delatte 2009)",
                                   "Ae. albopictus (Ezeakacha 2015)",
                                   "Ae. albopictus (Yee 2016)",
                                   "Cx. pipiens")) +
  theme_bw()


plot.EFGC.cpip.quad.uni

# ggsave("figures/EFGC.cpip.quad.uni.png", plot.EFGC.cpip.quad.uni, 
#        width = 10.3, height = 5.6)



##########
###### 2C. Compare the TPC fits ----
##########

## DIC
EFGC.cpip.bri.uni$BUGSoutput$DIC
EFGC.cpip.quad.uni$BUGSoutput$DIC

df.EFGC.cpip.bri.uni <- df.EFGC.cpip.bri.uni %>% 
  mutate(type = "Briere")

df.EFGC.cpip.quad.uni <- df.EFGC.cpip.quad.uni %>% 
  mutate(type = "Quadratic")


# Combine the three dataframes
df.all <- rbind(df.EFGC.cpip.bri.uni,
                df.EFGC.cpip.quad.uni
                )

# Plot
plot.all <- df.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = mean, color = type), linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Eggs per female per gonotrophic cycle"
  ) +
  theme_bw()

plot.all

# ggsave("figures/EFGC.cpip.all.png", plot.all, width = 10.3, height = 5.6)


##########
###### 3A. Fit EFGC thermal responses (Ae. albopictus): Briere ----
##########

##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)

##### Set data
data <- data.EFGC.aalb


## Create a unique id for each species-study combination
data <- data %>% 
  group_by(species, citation) %>% 
  mutate(unique_id = cur_group_id())

## Set priors
prior <- data.frame(q = c(0, 1),
                    T0 = c(0, 20),
                    Tm = c(30, 45),
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
EFGC.aalb.bri.uni.raneff <- jags(
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
# save(EFGC.aalb.bri.uni.raneff, file = "R-scripts/R2jags-objects/EFGC.aalb.bri.uni.raneff.Rdata")

# Read the .Rdata
load("R-scripts/R2jags-objects/EFGC.aalb.bri.uni.raneff.Rdata")


## Diagnostics ----
##### Examine output
EFGC.aalb.bri.uni.raneff$BUGSoutput$summary[1:8,]
mcmcplot(EFGC.aalb.bri.uni.raneff)

# Extract the DIC for future model comparisons
EFGC.aalb.bri.uni.raneff$BUGSoutput$DIC


## Plot data + fit ----
df.EFGC.aalb.bri.uni.raneff <- data.frame(EFGC.aalb.bri.uni.raneff$BUGSoutput$summary)[-(1:8),]

## Extract the model prediction
## Overall curve
df.EFGC.aalb.bri.uni.raneff.pop <- df.EFGC.aalb.bri.uni.raneff %>% 
  filter(grepl("z.trait.mu.pred.pop", rownames(df.EFGC.aalb.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)


## Unique ID 1: Ae. albopictus (Delatte 2009)
df.EFGC.aalb.bri.uni.raneff.1 <- df.EFGC.aalb.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[1,*]"), rownames(df.EFGC.aalb.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 1)

## Unique ID 2: Ae. albopictus (Ezeakacha 2015)
df.EFGC.aalb.bri.uni.raneff.2 <- df.EFGC.aalb.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[2,*]"), rownames(df.EFGC.aalb.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 2)

## Unique ID 3: Ae. albopictus (Yee 2016)
df.EFGC.aalb.bri.uni.raneff.3 <- df.EFGC.aalb.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[3,*]"), rownames(df.EFGC.aalb.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 3)



## Combine the model prediciton of all three unique groups into a dataframe
df.EFGC.aalb.bri.uni.raneff.sp <- rbind(df.EFGC.aalb.bri.uni.raneff.1,
                                        df.EFGC.aalb.bri.uni.raneff.2,
                                        df.EFGC.aalb.bri.uni.raneff.3) 

## Change unique_id into factor type
df.EFGC.aalb.bri.uni.raneff.sp$unique_id <- as.factor(df.EFGC.aalb.bri.uni.raneff.sp$unique_id)


##### Plot
plot.EFGC.aalb.bri.uni.raneff <- ggplot(data = df.EFGC.aalb.bri.uni.raneff.pop, 
                                              aes(x = temp)) +
  ## Overall TPC
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.),
              fill = "grey",
              alpha = 0.5) +
  ## a separate TPC (and credible interval) for each unique group
  # geom_ribbon(data = df.EFGC.aalb.bri.uni.raneff.sp, aes(ymin = X2.5., ymax = X97.5., fill = unique_id),
  #             alpha = 0.5) +
  geom_line(data = df.EFGC.aalb.bri.uni.raneff.sp, aes(y = mean, color = unique_id)) +
  geom_line(aes(y = mean), color = "black", linewidth = 1.5) +
  geom_point(data = data,
             aes(x = temp, y = trait, colour = as.factor(unique_id)),
             size = 2) +
  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), y = "Eggs per female per gonotrophic cycle") +
  # Customize legend
  scale_colour_discrete(name = element_blank(),
                        labels = c("Ae. albopictus (Delatte 2009)",
                                   "Ae. albopictus (Ezeakacha 2015)",
                                   "Ae. albopictus (Yee 2016)")) +
  theme_bw()


plot.EFGC.aalb.bri.uni.raneff

# ggsave("figures/EFGC.aalb.bri.uni.raneff.png", plot.EFGC.aalb.bri.uni.raneff,
#        width = 10.3, height = 5.6)


##########
###### 3B. Fit EFGC thermal responses: Quadratic ----
##########

##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)

##### Set data
data <- data.EFGC.aalb


## Create a unique id for each species-study combination
data <- data %>% 
  group_by(species, citation) %>% 
  mutate(unique_id = cur_group_id())

## Set priors
prior <- data.frame(q = c(0, 1),
                    T0 = c(0, 20),
                    Tm = c(30, 45),
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
EFGC.aalb.quad.uni.raneff <- jags(
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
# save(EFGC.aalb.quad.uni.raneff, file = "R-scripts/R2jags-objects/EFGC.aalb.quad.uni.raneff.Rdata")

# Read the .Rdata
load("R-scripts/R2jags-objects/EFGC.aalb.quad.uni.raneff.Rdata")


## Diagnostics ----
##### Examine output
EFGC.aalb.quad.uni.raneff$BUGSoutput$summary[1:8,]
mcmcplot(EFGC.aalb.quad.uni.raneff)

# Extract the DIC for future model comparisons
EFGC.aalb.quad.uni.raneff$BUGSoutput$DIC


## Plot data + fit ----
df.EFGC.aalb.quad.uni.raneff <- data.frame(EFGC.aalb.quad.uni.raneff$BUGSoutput$summary)[-(1:8),]

## Extract the model prediction
## Overall curve
df.EFGC.aalb.quad.uni.raneff.pop <- df.EFGC.aalb.quad.uni.raneff %>% 
  filter(grepl("z.trait.mu.pred.pop", rownames(df.EFGC.aalb.quad.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)


## Unique ID 1: Ae. albopictus (Delatte 2009)
df.EFGC.aalb.quad.uni.raneff.1 <- df.EFGC.aalb.quad.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[1,*]"), rownames(df.EFGC.aalb.quad.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 1)

## Unique ID 2: Ae. albopictus (Ezeakacha 2015)
df.EFGC.aalb.quad.uni.raneff.2 <- df.EFGC.aalb.quad.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[2,*]"), rownames(df.EFGC.aalb.quad.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 2)

## Unique ID 3: Ae. albopictus (Yee 2016)
df.EFGC.aalb.quad.uni.raneff.3 <- df.EFGC.aalb.quad.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[3,*]"), rownames(df.EFGC.aalb.quad.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 3)



## Combine the model prediciton of all three unique groups into a dataframe
df.EFGC.aalb.quad.uni.raneff.sp <- rbind(df.EFGC.aalb.quad.uni.raneff.1,
                                         df.EFGC.aalb.quad.uni.raneff.2,
                                         df.EFGC.aalb.quad.uni.raneff.3) 

## Change unique_id into factor type
df.EFGC.aalb.quad.uni.raneff.sp$unique_id <- as.factor(df.EFGC.aalb.quad.uni.raneff.sp$unique_id)


##### Plot
plot.EFGC.aalb.quad.uni.raneff <- ggplot(data = df.EFGC.aalb.quad.uni.raneff.pop, 
                                             aes(x = temp)) +
  ## Overall TPC
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.),
              fill = "grey",
              alpha = 0.5) +
  ## a separate TPC (and credible interval) for each unique group
  # geom_ribbon(data = df.EFGC.aalb.quad.uni.raneff.sp, aes(ymin = X2.5., ymax = X97.5., fill = unique_id),
  #             alpha = 0.5) +
  geom_line(data = df.EFGC.aalb.quad.uni.raneff.sp, aes(y = mean, color = unique_id)) +
  geom_line(aes(y = mean), color = "black", linewidth = 1.5) +
  geom_point(data = data,
             aes(x = temp, y = trait, colour = as.factor(unique_id)),
             size = 2) +
  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), y = "Eggs per female per gonotrophic cycle") +
  # Customize legend
  scale_colour_discrete(name = element_blank(),
                        labels = c("Ae. albopictus (Delatte 2009)",
                                   "Ae. albopictus (Ezeakacha 2015)",
                                   "Ae. albopictus (Yee 2016)",
                                   "Cx. pipiens")) +
  theme_bw()


plot.EFGC.aalb.quad.uni.raneff

# ggsave("figures/EFGC.aalb.quad.uni.raneff.png", plot.EFGC.aalb.quad.uni.raneff,
#        width = 10.3, height = 5.6)


##########
###### 3C. Compare the TPC fits ----
##########

## DIC
EFGC.aalb.bri.uni.raneff$BUGSoutput$DIC
EFGC.aalb.quad.uni.raneff$BUGSoutput$DIC

df.EFGC.aalb.bri.uni.raneff.pop <- df.EFGC.aalb.bri.uni.raneff.pop %>% 
  mutate(type = "Briere")

df.EFGC.aalb.quad.uni.raneff.pop <- df.EFGC.aalb.quad.uni.raneff.pop %>% 
  mutate(type = "Quadratic")


# Combine the three dataframes
df.all <- rbind(df.EFGC.aalb.bri.uni.raneff.pop,
                df.EFGC.aalb.quad.uni.raneff.pop
)

# Plot
plot.all <- df.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = mean, color = type), linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait, shape = citation), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Eggs per female per gonotrophic cycle"
  ) +
  theme_bw()

plot.all

# ggsave("figures/EFGC.aalb.all.png", plot.all, width = 10.3, height = 5.6)



##########
###### 4A. Fit EFGC thermal responses (all data): Briere ----
##########

##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)

##### Set data
data <- data.EFGC.nonarctic


## Create a unique id for each species-study combination
data <- data %>% 
  group_by(species, citation) %>% 
  mutate(unique_id = cur_group_id())

## Set priors
prior <- data.frame(q = c(0, 1),
                    T0 = c(0, 20),
                    Tm = c(30, 45),
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
EFGC.nonarctic.bri.uni.raneff <- jags(
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
# save(EFGC.nonarctic.bri.uni.raneff, file = "R-scripts/R2jags-objects/EFGC.nonarctic.bri.uni.raneff.Rdata")

# Read the .Rdata
load("R-scripts/R2jags-objects/EFGC.nonarctic.bri.uni.raneff.Rdata")


## Diagnostics ----
##### Examine output
EFGC.nonarctic.bri.uni.raneff$BUGSoutput$summary[1:8,]
mcmcplot(EFGC.nonarctic.bri.uni.raneff)

# Extract the DIC for future model comparisons
EFGC.nonarctic.bri.uni.raneff$BUGSoutput$DIC


## Plot data + fit ----
df.EFGC.nonarctic.bri.uni.raneff <- data.frame(EFGC.nonarctic.bri.uni.raneff$BUGSoutput$summary)[-(1:8),]

## Extract the model prediction
## Overall curve
df.EFGC.nonarctic.bri.uni.raneff.pop <- df.EFGC.nonarctic.bri.uni.raneff %>% 
  filter(grepl("z.trait.mu.pred.pop", rownames(df.EFGC.nonarctic.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)


## Unique ID 1: Ae. albopictus (Delatte 2009)
df.EFGC.nonarctic.bri.uni.raneff.1 <- df.EFGC.nonarctic.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[1,*]"), rownames(df.EFGC.nonarctic.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 1)

## Unique ID 2: Ae. albopictus (Ezeakacha 2015)
df.EFGC.nonarctic.bri.uni.raneff.2 <- df.EFGC.nonarctic.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[2,*]"), rownames(df.EFGC.nonarctic.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 2)

## Unique ID 3: Ae. albopictus (Yee 2016)
df.EFGC.nonarctic.bri.uni.raneff.3 <- df.EFGC.nonarctic.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[3,*]"), rownames(df.EFGC.nonarctic.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 3)

## Unique ID 4: Cx. pipiens
df.EFGC.nonarctic.bri.uni.raneff.4 <- df.EFGC.nonarctic.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[4,*]"), rownames(df.EFGC.nonarctic.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 4)


## Combine the model prediciton of all three unique groups into a dataframe
df.EFGC.nonarctic.bri.uni.raneff.sp <- rbind(df.EFGC.nonarctic.bri.uni.raneff.1,
                                             df.EFGC.nonarctic.bri.uni.raneff.2,
                                             df.EFGC.nonarctic.bri.uni.raneff.3,
                                             df.EFGC.nonarctic.bri.uni.raneff.4) 

## Change unique_id into factor type
df.EFGC.nonarctic.bri.uni.raneff.sp$unique_id <- as.factor(df.EFGC.nonarctic.bri.uni.raneff.sp$unique_id)


##### Plot
plot.EFGC.nonarctic.bri.uni.raneff <- ggplot(data = df.EFGC.nonarctic.bri.uni.raneff.pop, 
                                             aes(x = temp)) +
  ## Overall TPC
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.),
              fill = "grey",
              alpha = 0.5) +
  ## a separate TPC (and credible interval) for each unique group
  # geom_ribbon(data = df.EFGC.nonarctic.bri.uni.raneff.sp, aes(ymin = X2.5., ymax = X97.5., fill = unique_id),
  #             alpha = 0.5) +
  geom_line(data = df.EFGC.nonarctic.bri.uni.raneff.sp, aes(y = mean, color = unique_id)) +
  geom_line(aes(y = mean), color = "black", linewidth = 1.5) +
  geom_point(data = data,
             aes(x = temp, y = trait, colour = as.factor(unique_id)),
             size = 2) +
  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), y = "Eggs per female per gonotrophic cycle") +
  # Customize legend
  scale_colour_discrete(name = element_blank(),
                        labels = c("Ae. albopictus (Delatte 2009)",
                                   "Ae. albopictus (Ezeakacha 2015)",
                                   "Ae. albopictus (Yee 2016)",
                                   "Cx. pipiens")) +
  theme_bw()


plot.EFGC.nonarctic.bri.uni.raneff

# ggsave("figures/EFGC.nonarctic.bri.uni.raneff.png", plot.EFGC.nonarctic.bri.uni.raneff,
#        width = 10.3, height = 5.6)


##########
###### 4B. Fit EFGC thermal responses (all data): Quadratic ----
##########

##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)

##### Set data
data <- data.EFGC.nonarctic


## Create a unique id for each species-study combination
data <- data %>% 
  group_by(species, citation) %>% 
  mutate(unique_id = cur_group_id())

## Set priors
prior <- data.frame(q = c(0, 1),
                    T0 = c(0, 20),
                    Tm = c(30, 45),
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
EFGC.nonarctic.quad.uni.raneff <- jags(
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
# save(EFGC.nonarctic.quad.uni.raneff, file = "R-scripts/R2jags-objects/EFGC.nonarctic.quad.uni.raneff.Rdata")

# Read the .Rdata
load("R-scripts/R2jags-objects/EFGC.nonarctic.quad.uni.raneff.Rdata")


## Diagnostics ----
##### Examine output
EFGC.nonarctic.quad.uni.raneff$BUGSoutput$summary[1:8,]
mcmcplot(EFGC.nonarctic.quad.uni.raneff)

# Extract the DIC for future model comparisons
EFGC.nonarctic.quad.uni.raneff$BUGSoutput$DIC


## Plot data + fit ----
df.EFGC.nonarctic.quad.uni.raneff <- data.frame(EFGC.nonarctic.quad.uni.raneff$BUGSoutput$summary)[-(1:8),]

## Extract the model prediction
## Overall curve
df.EFGC.nonarctic.quad.uni.raneff.pop <- df.EFGC.nonarctic.quad.uni.raneff %>% 
  filter(grepl("z.trait.mu.pred.pop", rownames(df.EFGC.nonarctic.quad.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)


## Unique ID 1: Ae. albopictus (Delatte 2009)
df.EFGC.nonarctic.quad.uni.raneff.1 <- df.EFGC.nonarctic.quad.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[1,*]"), rownames(df.EFGC.nonarctic.quad.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 1)

## Unique ID 2: Ae. albopictus (Ezeakacha 2015)
df.EFGC.nonarctic.quad.uni.raneff.2 <- df.EFGC.nonarctic.quad.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[2,*]"), rownames(df.EFGC.nonarctic.quad.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 2)

## Unique ID 3: Ae. albopictus (Yee 2016)
df.EFGC.nonarctic.quad.uni.raneff.3 <- df.EFGC.nonarctic.quad.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[3,*]"), rownames(df.EFGC.nonarctic.quad.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 3)

## Unique ID 4: Cx. pipiens
df.EFGC.nonarctic.quad.uni.raneff.4 <- df.EFGC.nonarctic.quad.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[4,*]"), rownames(df.EFGC.nonarctic.quad.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 4)


## Combine the model prediciton of all three unique groups into a dataframe
df.EFGC.nonarctic.quad.uni.raneff.sp <- rbind(df.EFGC.nonarctic.quad.uni.raneff.1,
                                              df.EFGC.nonarctic.quad.uni.raneff.2,
                                              df.EFGC.nonarctic.quad.uni.raneff.3,
                                              df.EFGC.nonarctic.quad.uni.raneff.4) 

## Change unique_id into factor type
df.EFGC.nonarctic.quad.uni.raneff.sp$unique_id <- as.factor(df.EFGC.nonarctic.quad.uni.raneff.sp$unique_id)


##### Plot
plot.EFGC.nonarctic.quad.uni.raneff <- ggplot(data = df.EFGC.nonarctic.quad.uni.raneff.pop, 
                                              aes(x = temp)) +
  ## Overall TPC
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.),
              fill = "grey",
              alpha = 0.5) +
  ## a separate TPC (and credible interval) for each unique group
  # geom_ribbon(data = df.EFGC.nonarctic.quad.uni.raneff.sp, aes(ymin = X2.5., ymax = X97.5., fill = unique_id),
  #             alpha = 0.5) +
  geom_line(data = df.EFGC.nonarctic.quad.uni.raneff.sp, aes(y = mean, color = unique_id)) +
  geom_line(aes(y = mean), color = "black", linewidth = 1.5) +
  geom_point(data = data,
             aes(x = temp, y = trait, colour = as.factor(unique_id)),
             size = 2) +
  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), y = "Eggs per female per gonotrophic cycle") +
  # Customize legend
  scale_colour_discrete(name = element_blank(),
                        labels = c("Ae. albopictus (Delatte 2009)",
                                   "Ae. albopictus (Ezeakacha 2015)",
                                   "Ae. albopictus (Yee 2016)",
                                   "Cx. pipiens")) +
  theme_bw()


plot.EFGC.nonarctic.quad.uni.raneff

# ggsave("figures/EFGC.nonarctic.quad.uni.raneff.png", plot.EFGC.nonarctic.quad.uni.raneff,
#        width = 10.3, height = 5.6)


##########
###### 4C. Compare the TPC fits ----
##########

## DIC
EFGC.nonarctic.bri.uni.raneff$BUGSoutput$DIC
EFGC.nonarctic.quad.uni.raneff$BUGSoutput$DIC

df.EFGC.nonarctic.bri.uni.raneff.pop <- df.EFGC.nonarctic.bri.uni.raneff.pop %>% 
  mutate(type = "Briere")

df.EFGC.nonarctic.quad.uni.raneff.pop <- df.EFGC.nonarctic.quad.uni.raneff.pop %>% 
  mutate(type = "Quadratic")


# Combine the three dataframes
df.all <- rbind(df.EFGC.nonarctic.bri.uni.raneff.pop,
                df.EFGC.nonarctic.quad.uni.raneff.pop
)

# Plot
plot.all <- df.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = mean, color = type), linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait, colour = species), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Eggs per female per gonotrophic cycle"
  ) +
  theme_bw()

plot.all

# ggsave("figures/EFGC.nonarctic.all.png", plot.all, width = 10.3, height = 5.6)



#### DIC ----
#EFGC.aalb.bri.uni.raneff$BUGSoutput$DIC
#EFGC.aalb.quad.uni.raneff$BUGSoutput$DIC # This is the best fitting TPC
EFGC.nonarctic.bri.uni.raneff$BUGSoutput$DIC 
EFGC.nonarctic.quad.uni.raneff$BUGSoutput$DIC # This is the best fitting TPC


##########
###### 4. Process and save model output for plotting ----
##########

## Analyze TPC model
EFGC.TPC.analysis <- extractTPC_raneff(EFGC.nonarctic.quad.uni.raneff, "EFGC", Temp.xs)
EFGC.predictions.summary <- EFGC.TPC.analysis[[1]]
EFGC.params.summary <- EFGC.TPC.analysis[[2]]
EFGC.params.fullposts <- EFGC.TPC.analysis[[3]]

write_csv(EFGC.predictions.summary, "data-processed/EFGC.predictions.summary.csv")
write_csv(EFGC.params.summary, "data-processed/EFGC.params.summary.csv")
write_csv(EFGC.params.fullposts, "data-processed/EFGC.params.fullposts.csv")


