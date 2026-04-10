## Lilian Chan, University of Guelph
## Arctic vector-borne disease transmission suitability model
##
## Purpose: use Bayesian inference (JAGS) to fit TPCs for biting rate (a) using 
## data from Aedes aegypti, Aedes albopictus and four Aedes spp. from Alaska 
## (Ae. cinereus, Ae. communis, Ae. impiger, Ae. punctor).
##
## 
## Table of content:
##    0. Set-up workspace
##
##    1. MCMC settings for all models
##
##    2. Fitting TPC (Briere)
##        A. Fit a thermal responses with data from all species
##        B. Fit a thermal responses (non-Arctic species)
##
##    3. Fitting TPC (Quadratic)
##        A. Fit a thermal responses with data from all species
##        B. Fit a thermal responses for priors (non-Arctic species)
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
data <- read_csv("data-processed/TraitData_a.csv")
unique(data$species)

## Convert genotrophic cycle duration (1/a) to biting rate (a)
data <- data %>% 
  mutate(trait = ifelse(trait_name == "1/a", 1/trait, trait)) %>% 
  mutate(trait_name = ifelse(trait_name == "1/a", "a", trait_name))


# Subset data
## all species
data.a <- data

## Non-Arctic species only
data.a.nonarctic <- subset(data, species %in% c("aegypti", "albopictus"))


## Plot raw data
plot.data.a <- data.a %>%
  ggplot(aes(x = temp, y = trait)) +
  geom_point(aes(colour = species
                 #, shape = citation
  )) +
  labs(y = "Biting rate (1/days)", x = "Temperature ºC") +
  scale_colour_discrete(name = "Species", labels = c("Ae. aegypti",
                                                     "Ae. albopictus",
                                                     "Ae. cinereus",
                                                     "Ae. communis",
                                                     "Ae. impiger",
                                                     "Ae. punctor"
  )) +
  facet_grid(rows = vars(type)) +
  theme_bw()

plot.data.a

# ggsave("figures/raw_data/plot.data.a.png", plot.data.a, width = 9.83, height = 6.17)

## Put all data into the same graph
plot.data.a.combine <- data.a %>%
  ggplot(aes(x = temp, y = trait)) +
  geom_point(aes(colour = species
  )) +
  labs(y = "Biting rate (1/days)", x = "Temperature ºC") +
  scale_colour_discrete(name = "Species", labels = c("Ae. aegypti",
                                                     "Ae. albopictus",
                                                     "Ae. cinereus",
                                                     "Ae. communis",
                                                     "Ae. impiger",
                                                     "Ae. punctor"
  )) +
  theme_bw()

plot.data.a.combine

# ggsave("figures/raw_data/plot.data.a.combine.png", plot.data.a.combine, , width = 9.83, height = 6.17)



##########
###### 1. MCMC settings for all models ----
##########

# Number of posterior dist elements = [(ni - nb) / nt] * nc = [(45000 - 5000) / 8] * 3 = 15000
ni <- 45000 # number of iterations in each chain
nb <- 5000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 3 # number of chains


##########
###### 2A. Fit a thermal responses with data from all species ----
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
data <- data.a

#### Set priors
prior <- data.frame(q = c(0, 1),
                    T0 = c(0, 20),
                    Tm = c(20, 45)
)


##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, 
                 N.Temp.xs = N.Temp.xs, prior = prior)

##### Run JAGS
a.alldata.bri.uni <- jags(
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
# save(a.alldata.bri.uni, file = "R-scripts/R2jags-objects/a.alldata.bri.uni.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/a.alldata.bri.uni.Rdata")


## Diagnostics ----
##### Examine output
a.alldata.bri.uni$BUGSoutput$summary[1:5,]
mcmcplot(a.alldata.bri.uni)

a.alldata.bri.uni$BUGSoutput$DIC

## Plot data + fit ----
df.a.alldata.bri.uni <- data.frame(a.alldata.bri.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.a.alldata.bri.uni)

##### Plot
plot.a.alldata.bri.uni <- df.a.alldata.bri.uni %>%
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
  labs(x = expression(paste("Temperature (", degree, "C)")), y = "Biting rate (days-1)") +
  scale_colour_discrete(name = "Species", labels = c("Ae. aegypti",
                                                     "Ae. albopictus",
                                                     "Ae. cinereus",
                                                     "Ae. communis",
                                                     "Ae. impiger",
                                                     "Ae. punctor"
  )) +
  theme_bw()

plot.a.alldata.bri.uni

# ggsave("figures/a.alldata.bri.uni.png", plot.a.alldata.bri.uni,
#        width = 10.3, height = 5.6)



##########
###### 2B. Fit a thermal responses (with random effects): Briere ----
##########

##### Temp sequence for derived quantity calculations
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)


##### Set data
data <- data.a

## Create a unique id for each species-study combination
data <- data %>% 
  group_by(species, citation) %>% 
  mutate(unique_id = cur_group_id())

##### Set priors
prior <- data.frame(q = c(0, 0.001),
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
# a.alldata.bri.uni.raneff <- jags(
#   data = jag.data,
#   inits = inits,
#   parameters.to.save = parameters,
#   model.file = "R-scripts/briere_T_randeff.txt",
#   n.thin = nt,
#   n.chains = nc,
#   n.burnin = nb,
#   n.iter = ni,
#   DIC = T,
#   working.directory = getwd()
# )


## Save the model as Rdata 
# save(a.alldata.bri.uni.raneff, file = "R-scripts/R2jags-objects/a.alldata.bri.uni.raneff.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/a.alldata.bri.uni.raneff.Rdata")


## Diagnostics ----
##### Examine output
a.alldata.bri.uni.raneff$BUGSoutput$summary[1:8,]
mcmcplot(a.alldata.bri.uni.raneff)

# Extract the DIC for future model comparisons
a.alldata.bri.uni.raneff$BUGSoutput$DIC


## Plot data + fit ----
df.a.alldata.bri.uni.raneff <- data.frame(a.alldata.bri.uni.raneff$BUGSoutput$summary)[-(1:8),]

## Extract the model prediction
## Overall curve
df.a.alldata.bri.uni.raneff.pop <- df.a.alldata.bri.uni.raneff %>% 
  filter(grepl("z.trait.mu.pred.pop", rownames(df.a.alldata.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)


## Unique ID 1: Ae. aegypti (Focks 1993)
df.a.alldata.bri.uni.1 <- df.a.alldata.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[1,*]"), rownames(df.a.alldata.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 1)

## Unique ID 2: Ae. aegypti (Focks 2006)
df.a.alldata.bri.uni.2 <- df.a.alldata.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[2,*]"), rownames(df.a.alldata.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 2)

## Unique ID 3: Ae. aegypti (Goindin 2015)
df.a.alldata.bri.uni.3 <- df.a.alldata.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[3,*]"), rownames(df.a.alldata.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 3)


## Unique ID 4: Ae. aegypti (Morin 2015)
df.a.alldata.bri.uni.4 <- df.a.alldata.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[4,*]"), rownames(df.a.alldata.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 4)


## Unique ID 5: Ae. albopictus (Delatte 2009)
df.a.alldata.bri.uni.5 <- df.a.alldata.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[5,*]"), rownames(df.a.alldata.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 5)


## Unique ID 6: Ae. albopictus (Marini 2020)
df.a.alldata.bri.uni.6 <- df.a.alldata.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[6,*]"), rownames(df.a.alldata.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 6)


## Unique ID 7: Ae. cinereus
df.a.alldata.bri.uni.7 <- df.a.alldata.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[7,*]"), rownames(df.a.alldata.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 7)


## Unique ID 8: Ae. communis
df.a.alldata.bri.uni.8 <- df.a.alldata.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[8,*]"), rownames(df.a.alldata.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 8)


## Unique ID 9: Ae. impiger
df.a.alldata.bri.uni.9 <- df.a.alldata.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[9,*]"), rownames(df.a.alldata.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 9)


## Unique ID 10: Ae. punctor
df.a.alldata.bri.uni.10 <- df.a.alldata.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[10,*]"), rownames(df.a.alldata.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 10)

## Combine the model prediciton of all three unique groups into a dataframe
df.a.alldata.bri.uni.raneff.sp <- rbind(df.a.alldata.bri.uni.1,
                                   df.a.alldata.bri.uni.2,
                                   df.a.alldata.bri.uni.3,
                                   df.a.alldata.bri.uni.4,
                                   df.a.alldata.bri.uni.5,
                                   df.a.alldata.bri.uni.6,
                                   df.a.alldata.bri.uni.7,
                                   df.a.alldata.bri.uni.8,
                                   df.a.alldata.bri.uni.9,
                                   df.a.alldata.bri.uni.10) 

## Change unique_id into factor type
df.a.alldata.bri.uni.raneff.sp$unique_id <- as.factor(df.a.alldata.bri.uni.raneff.sp$unique_id)


##### Plot
plot.a.alldata.bri.uni.raneff <- ggplot(data = df.a.alldata.bri.uni.raneff.pop, 
                                            aes(x = temp)) +
  ## Overall TPC
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.),
              fill = "grey",
              alpha = 0.5) +
  ## a separate TPC (and credible interval) for each unique group
  # geom_ribbon(data = df.a.alldata.bri.uni.raneff.sp, aes(ymin = X2.5., ymax = X97.5., fill = unique_id),
  #             alpha = 0.5) +
    geom_point(data = data,
             aes(x = temp, y = trait, colour = as.factor(unique_id)),
             size = 2) +
  # geom_line(data = df.a.alldata.bri.uni.raneff.sp, aes(y = mean, color = unique_id)) +
  geom_line(aes(y = mean), color = "black", linewidth = 1.5) +

  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), y = "Biting rate (days-1)") +
  # Customize legend
  scale_colour_discrete(name = element_blank(),
                        labels = c("Ae. aegypti (Focks 1993)",
                                   "Ae. aegypti (Focks 2006)",
                                   "Ae. aegypti (Goindin 2015)",
                                   "Ae. aegypti (Morin 2015)",
                                   "Ae. albopictus (Delatte 2009)",
                                   "Ae. albopictus (Marini 2020)",
                                   "Ae. cinereus",
                                   "Ae. communis",
                                   "Ae. impiger",
                                   "Ae. punctor")) +
  theme_bw()


plot.a.alldata.bri.uni.raneff

# ggsave("figures/a.alldata.bri.uni.raneff.png", plot.a.alldata.bri.uni.raneff,
#        width = 10.3, height = 5.6)




##########
###### 2C. Plot all TPCs in the same graph (for comparison): Briere ----
##########

# Add an identifying column in each model output dataframe
df.a.alldata.bri.uni <- df.a.alldata.bri.uni %>% 
  mutate(type = "Briere uniform")

df.a.alldata.bri.uni.raneff.pop <- df.a.alldata.bri.uni.raneff.pop %>% 
  mutate(type = "Briere w/ random effects")


# Combine the three dataframes
df.all <- rbind(df.a.alldata.bri.uni, df.a.alldata.bri.uni.raneff.pop)

df.all$type <- factor(df.all$type, levels = c( "Briere uniform", "Briere w/ random effects"))


# Plot
plot.all <- df.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = mean, color = type), linewidth = 1) +
  geom_point(data = data.a, aes(x = temp, y = trait), size = 2) +
  #geom_point(data = data.a.nonarctic, aes(x = temp, y = trait), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Biting rate (days-1)"
  ) +
  # Customize the colours
  ## ribbon
  scale_fill_manual(values = c("Briere uniform" = "grey",
                               "Briere w/ random effects" = "#4363d8")) +
  
  ## line
  scale_color_manual(values = c("Briere uniform" = "#868686FF",
                                "Briere w/ random effects" = "blue")) +
  theme_bw()

plot.all

#ggsave("figures/a.bri.all.png", plot.all, width = 10.3, height = 5.6)

a.alldata.bri.uni$BUGSoutput$DIC
a.alldata.bri.uni.raneff$BUGSoutput$DIC



##########
###### 3A. Fit a thermal responses with uniform priors (Arctic): Quadratic ----
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
data <- data.a

##### Set priors 
prior <- data.frame(q = c(0, 1),
                    T0 = c(0, 20),
                    Tm = c(20, 45)
)


##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, 
                 N.Temp.xs = N.Temp.xs, prior = prior)

# ##### Run JAGS -----
# a.alldata.quad.uni <- jags(data = jag.data,
#                               inits = inits,
#                               parameters.to.save = parameters,
#                               model.file = "R-scripts/quad_T.txt",
#                               n.chains = nc,
#                               n.burnin = nb,
#                               n.iter = ni,
#                               DIC = T,
#                               working.directory = getwd()
# )

## Save the model as Rdata 
# save(a.alldata.quad.uni, file = "R-scripts/R2jags-objects/a.alldata.quad.uni.Rdata")

# Read the .Rdata
load("R-scripts/R2jags-objects/a.alldata.quad.uni.Rdata")


## Diagnostics ----
##### Examine output
a.alldata.quad.uni$BUGSoutput$summary[1:5,]
mcmcplot(a.alldata.quad.uni)

# Extract the DIC for future model comparisons
a.alldata.quad.uni$BUGSoutput$DIC

## Plot data + fit ----
df.a.alldata.quad.uni<- data.frame(a.alldata.quad.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.a.alldata.quad.uni)

##### Plot
plot.a.alldata.quad.uni <- df.a.alldata.quad.uni %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "#4363d8", alpha = 0.5) +
  geom_line(aes(y = mean), color = "blue", linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Biting rate (days-1)"
  ) +
  theme_bw()

plot.a.alldata.quad.uni

# ggsave("figures/a.alldata.quad.uni.png", plot.a.alldata.quad.uni, 
#        width = 10.3, height = 5.6)


##########
###### 3B. Fit a thermal responses (with random effects): Quadratic ----
##########

##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)


##### Set data
data <- data.a

## Create a unique id for each species-study combination
data <- data %>% 
  group_by(species, citation) %>% 
  mutate(unique_id = cur_group_id())

##### Set prior 
prior <- data.frame(q = c(0, 1),
                      T0 = c(0, 20),
                      Tm = c(20, 45),
                      sigma_q = c(0, 0.1),
                      sigma_T0 = c(0, 10),
                      sigma_Tm = c(0, 10)
                    )


##### inits Function
inits <- function(){list(
  cf.q = 0.1,
  cf.Tm = 35,
  cf.T0 = 5,
  cf.sigma = rlnorm(1),
  sigma_q = 0.1,
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
a.alldata.quad.uni.raneff <- jags(
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
# save(a.alldata.quad.uni.raneff, file = "R-scripts/R2jags-objects/a.alldata.quad.uni.raneff.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/a.alldata.quad.uni.raneff.Rdata")


## Diagnostics ----
##### Examine output
a.alldata.quad.uni.raneff$BUGSoutput$summary[1:8,]
mcmcplot(a.alldata.quad.uni.raneff)

# Extract the DIC for future model comparisons
a.alldata.quad.uni.raneff$BUGSoutput$DIC


## Plot data + fit ----
df.a.alldata.quad.uni.raneff <- data.frame(a.alldata.quad.uni.raneff$BUGSoutput$summary)[-(1:8),]

## Extract the model prediction
## Overall curve
df.a.alldata.quad.uni.raneff.pop <- df.a.alldata.quad.uni.raneff %>% 
  filter(grepl("z.trait.mu.pred.pop", rownames(df.a.alldata.quad.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)


## Unique ID 1: Ae. aegypti (Focks 1993)
df.a.alldata.quad.uni.1 <- df.a.alldata.quad.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[1,*]"), rownames(df.a.alldata.quad.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 1)

## Unique ID 2: Ae. aegypti (Focks 2006)
df.a.alldata.quad.uni.2 <- df.a.alldata.quad.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[2,*]"), rownames(df.a.alldata.quad.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 2)

## Unique ID 3: Ae. aegypti (Goindin 2015)
df.a.alldata.quad.uni.3 <- df.a.alldata.quad.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[3,*]"), rownames(df.a.alldata.quad.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 3)


## Unique ID 4: Ae. aegypti (Morin 2015)
df.a.alldata.quad.uni.4 <- df.a.alldata.quad.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[4,*]"), rownames(df.a.alldata.quad.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 4)


## Unique ID 5: Ae. albopictus (Delatte 2009)
df.a.alldata.quad.uni.5 <- df.a.alldata.quad.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[5,*]"), rownames(df.a.alldata.quad.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 5)


## Unique ID 6: Ae. albopictus (Marini 2020)
df.a.alldata.quad.uni.6 <- df.a.alldata.quad.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[6,*]"), rownames(df.a.alldata.quad.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 6)


## Unique ID 7: Ae. cinereus
df.a.alldata.quad.uni.7 <- df.a.alldata.quad.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[7,*]"), rownames(df.a.alldata.quad.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 7)


## Unique ID 8: Ae. communis
df.a.alldata.quad.uni.8 <- df.a.alldata.quad.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[8,*]"), rownames(df.a.alldata.quad.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 8)


## Unique ID 9: Ae. impiger
df.a.alldata.quad.uni.9 <- df.a.alldata.quad.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[9,*]"), rownames(df.a.alldata.quad.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 9)


## Unique ID 10: Ae. punctor
df.a.alldata.quad.uni.10 <- df.a.alldata.quad.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[10,*]"), rownames(df.a.alldata.quad.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 10)

## Combine the model prediciton of all three unique groups into a dataframe
df.a.alldata.quad.uni.raneff.sp <- rbind(df.a.alldata.quad.uni.1,
                                        df.a.alldata.quad.uni.2,
                                        df.a.alldata.quad.uni.3,
                                        df.a.alldata.quad.uni.4,
                                        df.a.alldata.quad.uni.5,
                                        df.a.alldata.quad.uni.6,
                                        df.a.alldata.quad.uni.7,
                                        df.a.alldata.quad.uni.8,
                                        df.a.alldata.quad.uni.9,
                                        df.a.alldata.quad.uni.10) 

## Change unique_id into factor type
df.a.alldata.quad.uni.raneff.sp$unique_id <- as.factor(df.a.alldata.quad.uni.raneff.sp$unique_id)


##### Plot
plot.a.alldata.quad.uni.raneff <- ggplot(data = df.a.alldata.quad.uni.raneff.pop, 
                                        aes(x = temp)) +
  ## Overall TPC
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.),
              fill = "grey",
              alpha = 0.5) +
  ## a separate TPC (and credible interval) for each unique group
  # geom_ribbon(data = df.a.alldata.quad.uni.raneff.sp, aes(ymin = X2.5., ymax = X97.5., fill = unique_id),
  #             alpha = 0.5) +
  geom_point(data = data,
             aes(x = temp, y = trait, colour = as.factor(unique_id)),
             size = 2) +
  geom_line(data = df.a.alldata.quad.uni.raneff.sp, aes(y = mean, color = unique_id)) +
  geom_line(aes(y = mean), color = "black", linewidth = 1.5) +
  
  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), y = "Biting rate (days-1)") +
  # Customize legend
  scale_colour_discrete(name = element_blank(),
                        labels = c("Ae. aegypti (Focks 1993)",
                                   "Ae. aegypti (Focks 2006)",
                                   "Ae. aegypti (Goindin 2015)",
                                   "Ae. aegypti (Morin 2015)",
                                   "Ae. albopictus (Delatte 2009)",
                                   "Ae. albopictus (Marini 2020)",
                                   "Ae. cinereus",
                                   "Ae. communis",
                                   "Ae. impiger",
                                   "Ae. punctor")) +
  theme_bw()


plot.a.alldata.quad.uni.raneff

# ggsave("figures/a.alldata.quad.uni.raneff.png", plot.a.alldata.quad.uni.raneff,
#        width = 10.3, height = 5.6)


##########
###### 3C. Plot all TPCs in the same graph (for comparison): Briere ----
##########

# Add an identifying column in each model output dataframe
df.a.alldata.quad.uni <- df.a.alldata.quad.uni %>% 
  mutate(type = "Quadratic uniform")

df.a.alldata.quad.uni.raneff.pop <- df.a.alldata.quad.uni.raneff.pop %>% 
  mutate(type = "Quadratic w/ random effects")


# Combine the three dataframes
df.all <- rbind(df.a.alldata.bri.uni, df.a.alldata.bri.uni.raneff.pop,
                df.a.alldata.quad.uni, df.a.alldata.quad.uni.raneff.pop)

df.all$type <- factor(df.all$type, levels = c("Briere uniform", 
                                              "Briere w/ random effects",
                                              "Quadratic uniform",
                                              "Quadratics w/ random effects"))


# Plot
plot.all <- df.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = mean, color = type), linewidth = 1) +
  geom_point(data = data.a, aes(x = temp, y = trait), size = 2) +
  #geom_point(data = data.a.nonarctic, aes(x = temp, y = trait), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Biting rate (days-1)"
  ) +
  # Customize the colours
  ## ribbon
  scale_fill_manual(values = c("Briere uniform" = "grey",
                               "Briere w/ random effects" = "#4363d8")) +
  
  ## line
  scale_color_manual(values = c("Briere uniform" = "#868686FF",
                                "Briere w/ random effects" = "blue")) +
  theme_bw()

plot.all

#ggsave("figures/a.bri.all.png", plot.all, width = 10.3, height = 5.6)

a.alldata.bri.uni$BUGSoutput$DIC
a.alldata.bri.uni.raneff$BUGSoutput$DIC # This is the best fitting TPC
a.alldata.quad.uni$BUGSoutput$DIC
a.alldata.quad.uni.raneff$BUGSoutput$DIC


##########
###### 4. Process and save model output for plotting ----
##########

## Analyze TPC model
a.TPC.analysis <- extractTPC_raneff(a.alldata.bri.uni.raneff, "a", Temp.xs)
a.predictions.summary <- a.TPC.analysis[[1]]
a.params.summary <- a.TPC.analysis[[2]]
a.params.fullposts <- a.TPC.analysis[[3]]

write_csv(a.predictions.summary, "data-processed/a.predictions.summary.csv")
write_csv(a.params.summary, "data-processed/a.params.summary.csv")
write_csv(a.params.fullposts, "data-processed/a.params.fullposts.csv")

