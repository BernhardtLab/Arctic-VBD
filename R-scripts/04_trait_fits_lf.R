## Lilian Chan, University of Guelph
## Arctic vector-borne disease transmission suitability model
##
## Purpose: use Bayesian unierence (JAGS) to fit TPCs for mosquito adult 
## lifespan (lf) for Arctic species with data-uniormed priors generated from 
## non-Arctic species data.


## 
## Table of content:
##    0. Set-up workspace
##
##    1. MCMC settings for all models
##
##    2. Fitting TPC (Briere)
##        A. Fit non-Arctic TPC for priors
##        B. Fit gamma distributions to non-Arctic TPC parameters
##        C. Fit Arctic TPC using data-uniormed priors
##
##    3. Fitting TPC (Quadratic)
##        A. Fit non-Arctic TPC for priors
##        B. Fit gamma distributions to non-Arctic TPC parameters
##        C. Fit Arctic TPC using data-uniormed priors
##
##    4. Compare model fit between Quadratic and Briere models
##    5. Process and save model output for plotting



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
data.all <- read_csv("data-processed/TraitData_lf.csv")
unique(data.all$species)


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

## Since Arctic datasets consisted of observations from 5 species, but 4 species 
## were measured at a single temperature, and the remaining species (Ae. vexans)
## had substantially shorter lifespans than the others. We believed that the 
## dataset containes insufficient information to fit a TPCs

## We combined the Arctic and non-Arctic species data and fitted a single TPC
## with uniform priors.

## Put all data into the same graph
plot.data.lf.combine <- data.all %>%
  ggplot(aes(x = temp, y = trait)) +
  geom_point(aes(colour = species)) +
  labs(y = "Time (days)", x = "Temperature ºC") +
  scale_colour_discrete(name = "Species", labels = c("Ae. albopictus",
                                                     "Ae. cinereus", 
                                                     "Ae. communis",
                                                     "Ae. impiger",
                                                     "Ae. punctor",
                                                     "Ae. vexans"
  )) +
  theme_bw()

plot.data.lf.combine


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
                    T0 = c(0, 17),
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
lf.alldata.bri.uni <- jags(data = jag.data,
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
save(lf.alldata.bri.uni, file = "R-scripts/R2jags-objects/all-mods/lf.alldata.bri.uni.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/all-mods/lf.alldata.bri.uni.Rdata")


## Diagnostics
##### Examine output
lf.alldata.bri.uni$BUGSoutput$summary[1:8,]
mcmcplot(lf.alldata.bri.uni)

# Extract the DIC for future model comparisons
lf.alldata.bri.uni$BUGSoutput$DIC


## Plot data + fit
df.lf.alldata.bri.uni <- data.frame(lf.alldata.bri.uni$BUGSoutput$summary)[-(1:8),]

## Extract the model prediction
## Overall curve
df.lf.alldata.bri.uni.pop <- df.lf.alldata.bri.uni %>% 
  filter(grepl("z.trait.mu.pred.pop", rownames(df.lf.alldata.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)


## Unique ID 1: Ae. albopictus (Alto 2001)
df.lf.alldata.bri.uni.1 <- df.lf.alldata.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[1,*]"), rownames(df.lf.alldata.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 1)

## Unique ID 2: Ae. albopictus (Calado and Navarro-Silva 2002)
df.lf.alldata.bri.uni.2 <- df.lf.alldata.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[2,*]"), rownames(df.lf.alldata.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 2)

## Unique ID 3: Ae. albopictus (Ezeakacha 2015)
df.lf.alldata.bri.uni.3 <- df.lf.alldata.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[3,*]"), rownames(df.lf.alldata.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 3)

## Unique ID 4: Ae. albopictus (Marini et al. 2020)
df.lf.alldata.bri.uni.4 <- df.lf.alldata.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[4,*]"), rownames(df.lf.alldata.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 4)

## Unique ID 5: Ae. albopictus (Tsuda et al. 1994)
df.lf.alldata.bri.uni.5 <- df.lf.alldata.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[5,*]"), rownames(df.lf.alldata.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 5)


## Unique ID 6: Ae. cinereus
df.lf.alldata.bri.uni.6 <- df.lf.alldata.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[6,*]"), rownames(df.lf.alldata.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 6)


## unique ID 7: Ae. communis
df.lf.alldata.bri.uni.7 <- df.lf.alldata.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[7,*]"), rownames(df.lf.alldata.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 7)

## unique ID 8: Ae. impiger
df.lf.alldata.bri.uni.8 <- df.lf.alldata.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[8,*]"), rownames(df.lf.alldata.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 8)

## unique ID 9: Ae. punctor
df.lf.alldata.bri.uni.9 <- df.lf.alldata.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[9,*]"), rownames(df.lf.alldata.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 9)

## unique ID 10: Ae. vexans
df.lf.alldata.bri.uni.10 <- df.lf.alldata.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[10,*]"), rownames(df.lf.alldata.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 10)


## Combine the model prediciton of all three unique groups into a dataframe
df.lf.alldata.bri.uni.sp <- rbind(df.lf.alldata.bri.uni.1,
                                  df.lf.alldata.bri.uni.2,
                                  df.lf.alldata.bri.uni.3,
                                  df.lf.alldata.bri.uni.4,
                                  df.lf.alldata.bri.uni.5,
                                  df.lf.alldata.bri.uni.6,
                                  df.lf.alldata.bri.uni.7,
                                  df.lf.alldata.bri.uni.8,
                                  df.lf.alldata.bri.uni.9,
                                  df.lf.alldata.bri.uni.10
                                  ) 

## Change unique_id into factor type
df.lf.alldata.bri.uni.sp$unique_id <- as.factor(df.lf.alldata.bri.uni.sp$unique_id)


##### Plot
plot.lf.alldata.bri.uni <- ggplot() +
  ## data
  geom_point(data = data,
             aes(x = temp, y = trait, colour = as.factor(unique_id)),
             size = 2) +
  
  ## a separate TPC for each unique group
  geom_line(data = df.lf.alldata.bri.uni.sp, 
            aes(x = temp, y = X50., color = unique_id)) +
  
  ## Overall TPC
  geom_ribbon(data = df.lf.alldata.bri.uni.pop,
              aes(x = temp, ymin = X2.5., ymax = X97.5.),
              fill = "grey",
              alpha = 0.5) +
  geom_line(data = df.lf.alldata.bri.uni.pop,
            aes(x = temp, y = X50.), color = "black", linewidth = 1) +
  
  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Time (days)") +
  
  # Customize legend
  scale_colour_discrete(name = "",
                        labels = c("Ae. albopictus 1",
                                   "Ae. albopictus 2",
                                   "Ae. albopictus 3",
                                   "Ae. albopictus 4",
                                   "Ae. albopictus 5",
                                   "Ae. cinereus",
                                   "Ae. communis",
                                   "Ae. impiger",
                                   "Ae. punctor",
                                   "Ae. vexans")) +
  theme_bw()


plot.lf.alldata.bri.uni

ggsave("figures/lf.alldata.bri.uni.png", plot.lf.alldata.bri.uni,
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
                    T0 = c(0, 17),
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
lf.alldata.quad.uni <- jags(data = jag.data,
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
save(lf.alldata.quad.uni, file = "R-scripts/R2jags-objects/all-mods/lf.alldata.quad.uni.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/all-mods/lf.alldata.quad.uni.Rdata")


## Diagnostics
##### Examine output
lf.alldata.quad.uni$BUGSoutput$summary[1:8,]
mcmcplot(lf.alldata.quad.uni)

# Extract the DIC for future model comparisons
lf.alldata.quad.uni$BUGSoutput$DIC


## Plot data + fit
df.lf.alldata.quad.uni <- data.frame(lf.alldata.quad.uni$BUGSoutput$summary)[-(1:8),]

## Extract the model prediction
## Overall curve
df.lf.alldata.quad.uni.pop <- df.lf.alldata.quad.uni %>% 
  filter(grepl("z.trait.mu.pred.pop", rownames(df.lf.alldata.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)


## Unique ID 1: Ae. albopictus (Alto 2001)
df.lf.alldata.quad.uni.1 <- df.lf.alldata.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[1,*]"), rownames(df.lf.alldata.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 1)

## Unique ID 2: Ae. albopictus (Calado and Navarro-Silva 2002)
df.lf.alldata.quad.uni.2 <- df.lf.alldata.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[2,*]"), rownames(df.lf.alldata.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 2)

## Unique ID 3: Ae. albopictus (Ezeakacha 2015)
df.lf.alldata.quad.uni.3 <- df.lf.alldata.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[3,*]"), rownames(df.lf.alldata.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 3)

## Unique ID 4: Ae. albopictus (Marini et al. 2020)
df.lf.alldata.quad.uni.4 <- df.lf.alldata.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[4,*]"), rownames(df.lf.alldata.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 4)

## Unique ID 5: Ae. albopictus (Tsuda et al. 1994)
df.lf.alldata.quad.uni.5 <- df.lf.alldata.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[5,*]"), rownames(df.lf.alldata.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 5)


## Unique ID 6: Ae. cinereus
df.lf.alldata.quad.uni.6 <- df.lf.alldata.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[6,*]"), rownames(df.lf.alldata.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 6)


## unique ID 7: Ae. communis
df.lf.alldata.quad.uni.7 <- df.lf.alldata.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[7,*]"), rownames(df.lf.alldata.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 7)

## unique ID 8: Ae. impiger
df.lf.alldata.quad.uni.8 <- df.lf.alldata.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[8,*]"), rownames(df.lf.alldata.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 8)

## unique ID 9: Ae. punctor
df.lf.alldata.quad.uni.9 <- df.lf.alldata.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[9,*]"), rownames(df.lf.alldata.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 9)

## unique ID 10: Ae. vexans
df.lf.alldata.quad.uni.10 <- df.lf.alldata.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[10,*]"), rownames(df.lf.alldata.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 10)



## Combine the model prediciton of all three unique groups into a dataframe
df.lf.alldata.quad.uni.sp <- rbind(df.lf.alldata.quad.uni.1,
                                   df.lf.alldata.quad.uni.2,
                                   df.lf.alldata.quad.uni.3,
                                   df.lf.alldata.quad.uni.4,
                                   df.lf.alldata.quad.uni.5,
                                   df.lf.alldata.quad.uni.6,
                                   df.lf.alldata.quad.uni.7,
                                   df.lf.alldata.quad.uni.8,
                                   df.lf.alldata.quad.uni.9,
                                   df.lf.alldata.quad.uni.10
                                   ) 

## Change unique_id into factor type
df.lf.alldata.quad.uni.sp$unique_id <- as.factor(df.lf.alldata.quad.uni.sp$unique_id)


##### Plot
plot.lf.alldata.quad.uni <- ggplot() +
  ## data
  geom_point(data = data,
             aes(x = temp, y = trait, colour = as.factor(unique_id)),
             size = 2) +
  
  ## a separate TPC for each unique group
  geom_line(data = df.lf.alldata.quad.uni.sp, 
            aes(x = temp, y = X50., color = unique_id)) +
  
  ## Overall TPC
  geom_ribbon(data = df.lf.alldata.quad.uni.pop,
              aes(x = temp, ymin = X2.5., ymax = X97.5.),
              fill = "grey",
              alpha = 0.5) +
  geom_line(data = df.lf.alldata.quad.uni.pop,
            aes(x = temp, y = X50.), color = "black", linewidth = 1) +
  
  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Time (days)") +
  
  # Customize legend
  scale_colour_discrete(name = "",
                        labels = c("Ae. albopictus 1",
                                   "Ae. albopictus 2",
                                   "Ae. albopictus 3",
                                   "Ae. albopictus 4",
                                   "Ae. albopictus 5",
                                   "Ae. cinereus",
                                   "Ae. communis",
                                   "Ae. impiger",
                                   "Ae. punctor",
                                   "Ae. vexans")) +
  theme_bw()


plot.lf.alldata.quad.uni

ggsave("figures/lf.alldata.quad.uni.png", plot.lf.alldata.quad.uni,
       width = 10.3, height = 5.6)



# 4. Compare model fit between Quadratic and Briere models ---------------------

##### Find best fitting model #####
# Add an identifying column in each model output dataframe
df.lf.alldata.bri.uni.pop <- df.lf.alldata.bri.uni.pop %>% 
  mutate(type = "briere")

df.lf.alldata.quad.uni.pop <- df.lf.alldata.quad.uni.pop %>% 
  mutate(type = "quadratic")

# Combine the two dataframes
df.all <- bind_rows(df.lf.alldata.bri.uni.pop, df.lf.alldata.quad.uni.pop)

##### Plot
plot.all <- df.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = X50., color = type), linewidth = 1) +
  geom_point(data = data.all, aes(x = temp, y = trait), size = 2) +
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
lf.alldata.bri.uni$BUGSoutput$DIC  
lf.alldata.quad.uni$BUGSoutput$DIC # This is the best fitting TPC



# Save best-fitting TPC in a separate folder
lf.alldata.mod <- lf.alldata.quad.uni

## Save the model as Rdata 
save(lf.alldata.mod, file = "R-scripts/R2jags-objects/best-fitting-mods/lf.alldata.mod.Rdata")


# 5. Process and save model output for plotting -------------------------------

## Analyze TPC model
# We will create 3 files: 
# a. predictions.summary: showing the mean, median, and 95% credible interval of
#      the predicted trait value at each temp from 0 to 45ºC at a 0.1ºC
# b. params.summary: showing the showing the mean, median, and 95% credible 
#      interval of TPC parameters, Topt, and Tbreadth
# c. params.fullposts: showing the TPC parameter of each MCMC iteration



Temp.xs <- seq(0, 45, 0.1)
lf.TPC.analysis <- extractTPC_raneff(lf.alldata.quad.uni, "lf", Temp.xs)
lf.alldata.predictions.summary <- lf.TPC.analysis[[1]]
lf.alldata.params.summary <- lf.TPC.analysis[[2]]
lf.alldata.params.fullposts <- lf.TPC.analysis[[3]]

write_csv(lf.alldata.predictions.summary, "data-processed/lf/lf.alldata.predictions.summary.csv")
write_csv(lf.alldata.params.summary, "data-processed/lf/lf.alldata.params.summary.csv")
write_csv(lf.alldata.params.fullposts, "data-processed/lf/lf.alldata.params.fullposts.csv")

