## Lilian Chan, University of Guelph
## Arctic vector-borne disease transmission suitability model
##
## Purpose: Fit thermal performance curves (TPCs) for eggs per female per 
## gonotrophic cycle (EFGC) using Bayesian inference (JAGS).
## 
## Table of content:
##    0. Set-up workspace
##    1. MCMC settings for all models
##    2. Fitting TPC (Briere)
##    3. Fitting TPC (Quadratic)
##    4. Compare model fit between Briere and Quadratic models
##    5. Process and save model output for visualization
##
##
## Inputs:
## data-processed/TraitData_EFGC.csv - 
##     Synthesized published trait data for EFGC
##
## Outputs: 
## R-scripts/R2jags-objects/best-fitting-mods/EFGC.alldata.mod.Rdata -
##     Best-fitting TPC models 
##
## data-processed/EFGC/EFGC.alldata.predictions.summary.csv -
##     Posterior summary of TPC predictions across temperatures
##
## data-processed/EFGC/EFGC.alldata.params.summary.csv -
##     Summary statistics of TPC parameters
##
## data-processed/EFGC/EFGC.alldata.params.fullposts.csv -
##     Full posterior distributions for TPC parameters


# 0. Set-up workspace ----------------------------------------------------------

library(tidyverse)
library(readxl)
library(janitor)
library(R2jags)
library(mcmcplots) # Diagnostic plots for fits
library(MASS)

##### Load functions
source("R-scripts/00_Functions.R")

# Load data
data.all <- read_csv("data-processed/TraitData_EFGC.csv")
unique(data.all$species)


## Plot raw data
plot.data.EFGC.alldata <- data.all %>% 
  ggplot(aes(x = temp, y = trait, colour = species)) +
  geom_point() +
  labs(y = "Eggs", x = "Temperature ºC") +
  facet_grid(rows = vars(type)) +
  theme_bw()


plot.data.EFGC.alldata

## Since the Arctic dataset included fewer than three temperature treatments 
## (the minimum requirement for constructing TPC), we combined the Arctic and 
## non-Arctic species data and fitted a single TPC with uniform priors.

## Put all data into the same graph
plot.data.EFGC.combine <- data.all %>%
  ggplot(aes(x = temp, y = trait)) +
  geom_point(aes(colour = species)) +
  labs(y = "Eggs", x = "Temperature ºC") +
  scale_colour_discrete(name = "Species", labels = c("Ae. albopictus",
                                                     "Ae. cinereus",
                                                     "Ae. communis",
                                                     "Ae. hexodontus",
                                                     "Ae. impiger",
                                                     "Ae. punctor"
  )) +
  theme_bw()

plot.data.EFGC.combine


# 1. MCMC Settings for all models ----------------------------------------------

# Number of posterior dist elements = [(ni - nb) / nt] * nc = [(45000 - 5000) / 8] * 3 = 15000
ni <- 45000 # number of iterations in each chain
nb <- 5000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 3 # number of chains

set.seed(123) # for reproducibility


# 2. Fitting TPC (Briere) ------------------------------------------------------

##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)

##### Set data
data <- data.all

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
  cf.q = 0.1,
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
EFGC.alldata.bri.uni <- jags(data = jag.data,
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
save(EFGC.alldata.bri.uni, file = "R-scripts/R2jags-objects/all-mods/EFGC.alldata.bri.uni.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/all-mods/EFGC.alldata.bri.uni.Rdata")


## Diagnostics
##### Examine output
EFGC.alldata.bri.uni$BUGSoutput$summary[1:8,]
# mcmcplot(EFGC.alldata.bri.uni)

# Extract the DIC for future model comparisons
EFGC.alldata.bri.uni$BUGSoutput$DIC


## Plot data + fit
df.EFGC.alldata.bri.uni <- data.frame(EFGC.alldata.bri.uni$BUGSoutput$summary)[-(1:8),]

## Extract the model prediction
## Overall curve
df.EFGC.alldata.bri.uni.pop <- df.EFGC.alldata.bri.uni %>% 
  filter(grepl("z.trait.mu.pred.pop", rownames(df.EFGC.alldata.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)


## Unique ID 1: Ae. albopictus (Delatte 2009)
df.EFGC.alldata.bri.uni.1 <- df.EFGC.alldata.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[1,*]"), rownames(df.EFGC.alldata.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 1)

## Unique ID 2: Ae. albopictus (Ezeakacha 2015)
df.EFGC.alldata.bri.uni.2 <- df.EFGC.alldata.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[2,*]"), rownames(df.EFGC.alldata.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 2)

## Unique ID 3: Ae. albopictus (Yee 2016)
df.EFGC.alldata.bri.uni.3 <- df.EFGC.alldata.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[3,*]"), rownames(df.EFGC.alldata.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 3)

## Unique ID 4: Ae. cinereus
df.EFGC.alldata.bri.uni.4 <- df.EFGC.alldata.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[4,*]"), rownames(df.EFGC.alldata.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 4)

## Unique ID 5: Ae. communis
df.EFGC.alldata.bri.uni.5 <- df.EFGC.alldata.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[5,*]"), rownames(df.EFGC.alldata.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 5)


## Unique ID 6: Ae. hexodontus
df.EFGC.alldata.bri.uni.6 <- df.EFGC.alldata.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[6,*]"), rownames(df.EFGC.alldata.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 6)


## unique ID 7: Ae. impiger
df.EFGC.alldata.bri.uni.7 <- df.EFGC.alldata.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[7,*]"), rownames(df.EFGC.alldata.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 7)

## unique ID 8: Ae. punctor
df.EFGC.alldata.bri.uni.8 <- df.EFGC.alldata.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[8,*]"), rownames(df.EFGC.alldata.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 8)



## Combine the model prediciton of all three unique groups into a dataframe
df.EFGC.alldata.bri.uni.sp <- rbind(df.EFGC.alldata.bri.uni.1,
                                    df.EFGC.alldata.bri.uni.2,
                                    df.EFGC.alldata.bri.uni.3,
                                    df.EFGC.alldata.bri.uni.4,
                                    df.EFGC.alldata.bri.uni.5,
                                    df.EFGC.alldata.bri.uni.6,
                                    df.EFGC.alldata.bri.uni.7,
                                    df.EFGC.alldata.bri.uni.8
                                    ) 

## Change unique_id into factor type
df.EFGC.alldata.bri.uni.sp$unique_id <- as.factor(df.EFGC.alldata.bri.uni.sp$unique_id)


##### Plot
plot.EFGC.alldata.bri.uni <- ggplot() +
  ## data
  geom_point(data = data,
             aes(x = temp, y = trait, colour = as.factor(unique_id)),
             size = 2) +
  
  ## a separate TPC for each unique group
  geom_line(data = df.EFGC.alldata.bri.uni.sp, 
            aes(x = temp, y = X50., color = unique_id)) +
  
  ## Overall TPC
  geom_ribbon(data = df.EFGC.alldata.bri.uni.pop,
              aes(x = temp, ymin = X2.5., ymax = X97.5.),
              fill = "grey",
              alpha = 0.5) +
  geom_line(data = df.EFGC.alldata.bri.uni.pop,
            aes(x = temp, y = X50.), color = "black", linewidth = 1) +
  
  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Eggs") +
  
  # Customize legend
  scale_colour_discrete(name = "",
                        labels = c("Ae. albopictus 1",
                                   "Ae. albopictus 2",
                                   "Ae. albopictus 3",
                                   "Ae. cinereus",
                                   "Ae. communis",
                                   "Ae. hexodontus",
                                   "Ae. impiger",
                                   "Ae. punctor"
                                   )) +
  theme_bw()


plot.EFGC.alldata.bri.uni

ggsave("figures/EFGC.alldata.bri.uni.png", plot.EFGC.alldata.bri.uni,
       width = 10.3, height = 5.6)



# 3. Fitting TPC (quadratic) ---------------------------------------------------

##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)

##### Set data
data <- data.all

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
  cf.q = 0.1,
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
EFGC.alldata.quad.uni <- jags(data = jag.data,
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
save(EFGC.alldata.quad.uni, file = "R-scripts/R2jags-objects/all-mods/EFGC.alldata.quad.uni.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/all-mods/EFGC.alldata.quad.uni.Rdata")


## Diagnostics
##### Examine output
EFGC.alldata.quad.uni$BUGSoutput$summary[1:8,]
# mcmcplot(EFGC.alldata.quad.uni)

# Extract the DIC for future model comparisons
EFGC.alldata.quad.uni$BUGSoutput$DIC


## Plot data + fit
df.EFGC.alldata.quad.uni <- data.frame(EFGC.alldata.quad.uni$BUGSoutput$summary)[-(1:8),]

## Extract the model prediction
## Overall curve
df.EFGC.alldata.quad.uni.pop <- df.EFGC.alldata.quad.uni %>% 
  filter(grepl("z.trait.mu.pred.pop", rownames(df.EFGC.alldata.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)


## Unique ID 1: Ae. albopictus (Delatte 2009)
df.EFGC.alldata.quad.uni.1 <- df.EFGC.alldata.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[1,*]"), rownames(df.EFGC.alldata.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 1)

## Unique ID 2: Ae. albopictus (Ezeakacha 2015)
df.EFGC.alldata.quad.uni.2 <- df.EFGC.alldata.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[2,*]"), rownames(df.EFGC.alldata.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 2)

## Unique ID 3: Ae. albopictus (Yee 2016)
df.EFGC.alldata.quad.uni.3 <- df.EFGC.alldata.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[3,*]"), rownames(df.EFGC.alldata.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 3)

## Unique ID 4: Ae. cinereus
df.EFGC.alldata.quad.uni.4 <- df.EFGC.alldata.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[4,*]"), rownames(df.EFGC.alldata.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 4)

## Unique ID 5: Ae. communis
df.EFGC.alldata.quad.uni.5 <- df.EFGC.alldata.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[5,*]"), rownames(df.EFGC.alldata.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 5)


## Unique ID 6: Ae. hexodontus
df.EFGC.alldata.quad.uni.6 <- df.EFGC.alldata.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[6,*]"), rownames(df.EFGC.alldata.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 6)


## unique ID 7: Ae. impiger
df.EFGC.alldata.quad.uni.7 <- df.EFGC.alldata.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[7,*]"), rownames(df.EFGC.alldata.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 7)

## unique ID 8: Ae. punctor
df.EFGC.alldata.quad.uni.8 <- df.EFGC.alldata.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[8,*]"), rownames(df.EFGC.alldata.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 8)



## Combine the model prediciton of all three unique groups into a dataframe
df.EFGC.alldata.quad.uni.sp <- rbind(df.EFGC.alldata.quad.uni.1,
                                     df.EFGC.alldata.quad.uni.2,
                                     df.EFGC.alldata.quad.uni.3,
                                     df.EFGC.alldata.quad.uni.4,
                                     df.EFGC.alldata.quad.uni.5,
                                     df.EFGC.alldata.quad.uni.6,
                                     df.EFGC.alldata.quad.uni.7,
                                     df.EFGC.alldata.quad.uni.8
                                     ) 

## Change unique_id into factor type
df.EFGC.alldata.quad.uni.sp$unique_id <- as.factor(df.EFGC.alldata.quad.uni.sp$unique_id)


##### Plot
plot.EFGC.alldata.quad.uni <- ggplot() +
  ## data
  geom_point(data = data,
             aes(x = temp, y = trait, colour = as.factor(unique_id)),
             size = 2) +
  
  ## a separate TPC for each unique group
  geom_line(data = df.EFGC.alldata.quad.uni.sp, 
            aes(x = temp, y = X50., color = unique_id)) +
  
  ## Overall TPC
  geom_ribbon(data = df.EFGC.alldata.quad.uni.pop,
              aes(x = temp, ymin = X2.5., ymax = X97.5.),
              fill = "grey",
              alpha = 0.5) +
  geom_line(data = df.EFGC.alldata.quad.uni.pop,
            aes(x = temp, y = X50.), color = "black", linewidth = 1) +
  
  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Eggs") +
  
  # Customize legend
  scale_colour_discrete(name = "",
                        labels = c("Ae. albopictus 1",
                                   "Ae. albopictus 2",
                                   "Ae. albopictus 3",
                                   "Ae. cinereus",
                                   "Ae. communis",
                                   "Ae. hexodontus",
                                   "Ae. impiger",
                                   "Ae. punctor"
                        )) +
  theme_bw()


plot.EFGC.alldata.quad.uni

ggsave("figures/EFGC.alldata.quad.uni.png", plot.EFGC.alldata.quad.uni,
       width = 10.3, height = 5.6)


# 4. Compare model fit between Briere and Quadratic models ---------------------

##### Find best fitting model #####
# Add an identifying column in each model output dataframe
df.EFGC.alldata.bri.uni.pop <- df.EFGC.alldata.bri.uni.pop %>% 
  mutate(type = "briere")

df.EFGC.alldata.quad.uni.pop <- df.EFGC.alldata.quad.uni.pop %>% 
  mutate(type = "quadratic")

# Combine the two dataframes
df.all <- bind_rows(df.EFGC.alldata.bri.uni.pop, df.EFGC.alldata.quad.uni.pop)

##### Plot
plot.all <- df.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = X50., color = type), linewidth = 1) +
  geom_point(data = data.all, aes(x = temp, y = trait), size = 2) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Eggs per female per gonotrophic cycle"
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

ggsave("figures/EFGC.bri.quad.png", plot.all, width = 10.3, height = 5.6)


#### DIC
EFGC.alldata.bri.uni$BUGSoutput$DIC
EFGC.alldata.quad.uni$BUGSoutput$DIC # This is the best fitting TPC


# Save best-fitting TPC in a separate folder
EFGC.alldata.mod <- EFGC.alldata.quad.uni

## Save the model as Rdata 
save(EFGC.alldata.mod, file = "R-scripts/R2jags-objects/best-fitting-mods/EFGC.alldata.mod.Rdata")


# 5. Process and save model output for visualization ---------------------------

## Analyze TPC model
# We will create 3 files: 
# a. predictions.summary: showing the mean, median, and 95% credible interval of
#      the predicted trait value at each temp from 0 to 45ºC at a 0.1ºC
# b. params.summary: showing the showing the mean, median, and 95% credible 
#      interval of TPC parameters, Topt, and Tbreadth
# c. params.fullposts: showing the TPC parameter of each MCMC iteration


##### non-Arctic #####
Temp.xs <- seq(0, 45, 0.1)
EFGC.TPC.analysis <- extractTPC_raneff(EFGC.alldata.quad.uni, "EFGC", Temp.xs)
EFGC.alldata.predictions.summary <- EFGC.TPC.analysis[[1]]
EFGC.alldata.params.summary <- EFGC.TPC.analysis[[2]]
EFGC.alldata.params.fullposts <- EFGC.TPC.analysis[[3]]

write_csv(EFGC.alldata.predictions.summary, "data-processed/EFGC/EFGC.alldata.predictions.summary.csv")
write_csv(EFGC.alldata.params.summary, "data-processed/EFGC/EFGC.alldata.params.summary.csv")
write_csv(EFGC.alldata.params.fullposts, "data-processed/EFGC/EFGC.alldata.params.fullposts.csv")
