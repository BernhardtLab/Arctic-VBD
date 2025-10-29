## Lilian Chan, University of Guelph
## Arctic vector-borne disease transmission suitability model
##
## Purpose: use Bayesian inference (JAGS) to fit TPCs for mosquito adult lifespan (lf) 
## for Aedes vexans (Costello and Brust 1971) 
##     1) with uniform priors; and 
##     2) with data-informed priors from Aedes sierrensis data (Couper et al. 2024)
##
## Table of content:
##    0. Set-up workspace
##
##    1. MCMC settings for all models
##
##    2. Fitting TPC (Briere)
##        A. Fit lf thermal responses with uniform priors (Ae. vexans)
##        B. Fit lf thermal responses for priors (Ae. sierrensis)
##        C. Fit gamma distributions to lf prior thermal responses
##        D. Fit lf thermal responses with data-informed priors (Ae. vexans)
##        E. Plot all three TPCs in the same graph (for comparison)
##
##    3. Fitting TPC (Quadratic)
##        A. Fit lf thermal responses with uniform priors (Ae. vexans)
##        B. Fit lf thermal responses for priors (Ae. sierrensis)
##        C. Fit gamma distributions to lf prior thermal responses
##        D. Fit lf thermal responses with data-informed priors (Ae. vexans)


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
data <- read_csv("data-processed/TraitData_lf.csv")
unique(data$species)


# Subset data
## Arctic species
data.lf.arctic <- subset(data, species %in% c("cinereus", "communis", "impiger",
                                              "punctor", "vexans"))

## Non-Arctic species
data.lf.nonarctic <- subset(data, species %in% c("aegypti", "albopictus", "sierrensis"))



## Plot raw data
plot.data.lf <- data %>% 
  mutate(type = c(rep("Arctic", 54), rep("non-Arctic", 1712))) %>% 
  ggplot(aes(x = temp, y = trait, colour = species)) +
  geom_point(aes(colour = species)) +
  
  ## Since the Ae. aegypti, albopictus, and sierrensis has many data, I will just plot the mean±SE
  # geom_point(data = ~filter(.x, type == "Arctic")) +
  # geom_point(data = ~filter(.x, type == "Arctic")) +
  # stat_summary(data = ~filter(.x, type == "non-Arctic"),
  #              fun = mean, geom = "point") +
  # stat_summary(data = ~filter(.x, type == "non-Arctic"),
  #              fun.data = "mean_se", geom = "errorbar") +
  
  labs(y = "Mosquito adult lifespan (days)", x = "Temperature ºC") +
  scale_colour_discrete(name = "Species", labels = c("Ae. aegypti",
                                                     "Ae. albopictus",
                                                     "Ae. cinereus", 
                                                     "Ae. communis",
                                                     "Ae. impiger",
                                                     "Ae. punctor",
                                                     "Ae. sierrensis",
                                                     "Ae. vexans"
  )) +
  facet_grid(rows = vars(type)) +
  theme_bw()

plot.data.lf

# ggsave("figures/raw_data/plot.data.lf.png", plot.data.lf, , width = 9.83, height = 6.17)



##########
###### 1. MCMC settings for all models ----
##########

# Number of posterior dist elements = [(ni - nb) / nt] * nc = [(45000 - 5000) / 8] * 3 = 15000
ni <- 45000 # number of iterations in each chain
nb <- 5000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 3 # number of chains


##########
###### 2A. Fit lf thermal responses with uniform priors (Arctic): Briere ----
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
data <- data.lf.arctic

##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
lf.arctic.bri.uni <- jags(
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
# save(lf.arctic.bri.uni, file = "R-scripts/R2jags-objects/lf.arctic.bri.uni.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/lf.arctic.bri.uni.Rdata")


## Diagnostics ----
##### Examine output
lf.arctic.bri.uni$BUGSoutput$summary[1:5,]
mcmcplot(lf.arctic.bri.uni)

# Extract the DIC for future model comparisons
lf.arctic.bri.uni$BUGSoutput$DIC

## Plot data + fit ----
df.lf.arctic.bri.uni <- data.frame(lf.arctic.bri.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.lf.arctic.bri.uni)

##### Plot
plot.lf.arctic.bri.uni <- df.lf.arctic.bri.uni %>%
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.),
              fill = "grey",
              alpha = 0.5) +
  geom_line(aes(y = mean), color = "#868686FF", linewidth = 1) +
  geom_point(data = data,
             aes(x = temp, y = trait, colour = species),
             size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) +
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(x = expression(paste("Temperature (", degree, "C)")), y = "mosquito adult lifespan (days)") +
  # Customize legend
  scale_color_discrete(name = "Species",
                       labels = c("Ae. cinereus",
                                  "Ae. communis",
                                  "Ae. impiger",
                                  "Ae. punctor",
                                  "Ae. vexans")) +
  theme_bw()

plot.lf.arctic.bri.uni

# ggsave("figures/lf.arctic.bri.uni.png", plot.lf.arctic.bri.uni,
#        width = 10.3, height = 5.6)



##########
###### 2B. Fit lf thermal responses for priors (non-Arctic species): Briere ----
##########

##### Temp sequence for derived quantity calculations
# For priors - fewer temps for derived calculations makes it go faster
Temp.xs <- seq(0, 45, 0.5)
N.Temp.xs <-length(Temp.xs)


##### Set data
data <- data.lf.nonarctic

## Create a unique id for each species-study combination
data <- data %>% 
  group_by(species, citation) %>% 
  mutate(unique_id = cur_group_id())


## Set priors
prior <- data.frame(q = c(0, 0.01),
                    T0 = c(0, 20),
                    Tm = c(30, 45),
                    sigma_q = c(0, 0.0001),
                    sigma_T0 = c(0, 10),
                    sigma_Tm = c(0, 10)
)

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
lf.nonarctic.bri.uni.raneff <- jags(
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
# save(lf.nonarctic.bri.uni.raneff, file = "R-scripts/R2jags-objects/lf.nonarctic.bri.uni.raneff.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/lf.nonarctic.bri.uni.raneff.Rdata")


## Diagnostics ----
##### Examine output
lf.nonarctic.bri.uni.raneff$BUGSoutput$summary[1:8,]
mcmcplot(lf.nonarctic.bri.uni.raneff)

# Extract the DIC for future model comparisons
lf.nonarctic.bri.uni.raneff$BUGSoutput$DIC


## Plot data + fit ----
df.lf.nonarctic.bri.uni.raneff <- data.frame(lf.nonarctic.bri.uni.raneff$BUGSoutput$summary)[-(1:8),]

## Extract the model prediction
## Overall curve
df.lf.nonarctic.bri.uni.raneff.pop <- df.lf.nonarctic.bri.uni.raneff %>% 
  filter(grepl("z.trait.mu.pred.pop", rownames(df.lf.nonarctic.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)


## Unique ID 1: Ae. aegypti (Beserra 2009)
df.lf.nonarctic.bri.uni.1 <- df.lf.nonarctic.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[1,*]"), rownames(df.lf.nonarctic.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 1)

## Unique ID 2: Ae. aegypti (Goindin et al. 2015)
df.lf.nonarctic.bri.uni.2 <- df.lf.nonarctic.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[2,*]"), rownames(df.lf.nonarctic.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 2)

## Unique ID 3: Ae. aegypti (Huxley et al. 2021)
df.lf.nonarctic.bri.uni.3 <- df.lf.nonarctic.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[3,*]"), rownames(df.lf.nonarctic.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 3)

## Unique ID 4: Ae. aegypti (Huxley et al. 2022)
df.lf.nonarctic.bri.uni.4 <- df.lf.nonarctic.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[4,*]"), rownames(df.lf.nonarctic.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 4)

## Unique ID 5: Ae. aegypti (Marinho et al. 2016.)
df.lf.nonarctic.bri.uni.5 <- df.lf.nonarctic.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[5,*]"), rownames(df.lf.nonarctic.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 5)

## Unique ID 6: Ae. aegypti (Rocha-Santos et al. 2021)
df.lf.nonarctic.bri.uni.6 <- df.lf.nonarctic.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[6,*]"), rownames(df.lf.nonarctic.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 6)

## Unique ID 7: Ae. aegypti (Yang et al. 2009)
df.lf.nonarctic.bri.uni.7 <- df.lf.nonarctic.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[7,*]"), rownames(df.lf.nonarctic.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 7)

## Unique ID 8: Ae. albopictus (Calado and Navarro-Silva 2002)
df.lf.nonarctic.bri.uni.8 <- df.lf.nonarctic.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[8,*]"), rownames(df.lf.nonarctic.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 8)

## Unique ID 9: Ae. albopictus (Marini et al. 2020)
df.lf.nonarctic.bri.uni.9 <- df.lf.nonarctic.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[9,*]"), rownames(df.lf.nonarctic.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 9)

## Unique ID 10: Ae. albopictus (Tsuda et al. 1994)
df.lf.nonarctic.bri.uni.10 <- df.lf.nonarctic.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[10,*]"), rownames(df.lf.nonarctic.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 10)

## Unique ID 11: Ae. sierrensis
df.lf.nonarctic.bri.uni.11 <- df.lf.nonarctic.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[11,*]"), rownames(df.lf.nonarctic.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(unique_id = 11)

## Combine the model prediciton of all three unique groups into a dataframe
df.lf.nonarctic.bri.uni.raneff.sp <- rbind(df.lf.nonarctic.bri.uni.1,
                                           df.lf.nonarctic.bri.uni.2,
                                           df.lf.nonarctic.bri.uni.3,
                                           df.lf.nonarctic.bri.uni.4,
                                           df.lf.nonarctic.bri.uni.5,
                                           df.lf.nonarctic.bri.uni.6,
                                           df.lf.nonarctic.bri.uni.7,
                                           df.lf.nonarctic.bri.uni.8,
                                           df.lf.nonarctic.bri.uni.9,
                                           df.lf.nonarctic.bri.uni.10,
                                           df.lf.nonarctic.bri.uni.11
                                           ) 

## Change unique_id into factor type
df.lf.nonarctic.bri.uni.raneff.sp$unique_id <- as.factor(df.lf.nonarctic.bri.uni.raneff.sp$unique_id)


##### Plot
plot.lf.nonarctic.bri.uni.raneff <- ggplot(data = df.lf.nonarctic.bri.uni.raneff.pop, 
                                            aes(x = temp)) +
  ## Overall TPC
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.),
              fill = "grey",
              alpha = 0.5) +
  ## a separate TPC (and credible interval) for each unique group
  geom_ribbon(data = df.lf.nonarctic.bri.uni.raneff.sp, aes(ymin = X2.5., ymax = X97.5., fill = unique_id),
              alpha = 0.5) +
  geom_line(aes(y = mean), color = "black", linewidth = 1) +
  geom_line(data = df.lf.nonarctic.bri.uni.raneff.sp, aes(y = mean, color = unique_id)) +
  geom_point(data = data,
             aes(x = temp, y = trait, colour = as.factor(unique_id)),
             size = 2) +
  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), y = "Mosquito adult lifespan (days)") +
  # Customize legend
  # scale_colour_discrete(name = element_blank(),
  #                       labels = c("Ae. aegypti (Goindin et al. 2015)",
  #                                  "Ae. aegypti (Huxley et al. 2021)",
  #                                  "Ae. aegypti (Huxley et al. 2022)",
  #                                  "Ae. aegypti (Marinho et al. 2016)",
  #                                  "Ae. aegypti (Rocha-Santos et al. 2021)",
  #                                  "Ae. aegypti (Yang et al. 2009)",
  #                                  "Ae. albopictus (Calado and Navarro-Silva 2002)",
  #                                  "Ae. albopictus (Marini et al. 2020)",
  #                                  "Ae. albopictus (Tsuda et al. 1994)",
  #                                  "Ae. sierrensis")) +
  theme_bw()


plot.lf.nonarctic.bri.uni.raneff

# ggsave("figures/lf.nonarctic.bri.uni.raneff.png", plot.lf.nonarctic.bri.uni.raneff,
#        width = 10.3, height = 5.6)



##########
###### 2C. Fit gamma distributions to lf prior thermal responses: Briere ----
##########

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
lf.arctic.prior.cf.dists <- data.frame(q = as.vector(lf.nonarctic.bri.uni.raneff$BUGSoutput$sims.list$cf.q),
                                       T0 = as.vector(lf.nonarctic.bri.uni.raneff$BUGSoutput$sims.list$cf.T0),
                                       Tm = as.vector(lf.nonarctic.bri.uni.raneff$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
lf.arctic.prior.gamma.fits = apply(lf.arctic.prior.cf.dists, 2, 
                                    function(df) fitdistr(df, "gamma")$estimate)


lf.hypers <- lf.arctic.prior.gamma.fits
# save(lf.hypers, file = "R-scripts/R2jags-objects/lfhypers.bri.Rsave")



##########
###### 2D. Fit lf thermal responses with data-informed priors (Arctic): Briere ----
##########

load("R-scripts/R2jags-objects/lfhypers.bri.Rsave")
lf.arctic.prior.gamma.fits <- lf.hypers


##### Set data
data <- data.lf.arctic
hypers <- lf.arctic.prior.gamma.fits * 0.1

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
lf.arctic.bri.inf <- jags(data = jag.data,
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
# save(lf.arctic.bri.inf, file = "R-scripts/R2jags-objects/lf.arctic.bri.inf.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/lf.arctic.bri.inf.Rdata")


## Diagnostics ----
##### Examine output
lf.arctic.bri.inf$BUGSoutput$summary[1:5,]
mcmcplot(lf.arctic.bri.inf)

# Extract the DIC for future model comparisons
lf.arctic.bri.inf$BUGSoutput$DIC

## Plot data + fit ----
df.lf.arctic.bri.inf <- data.frame(lf.arctic.bri.inf$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.lf.arctic.bri.inf)

##### Plot
plot.lf.arctic.bri.inf <- df.lf.arctic.bri.inf %>% 
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

plot.lf.arctic.bri.inf

# ggsave("figures/lf.arctic.bri.inf.png", plot.lf.arctic.bri.inf,
#        width = 10.3, height = 5.6)



##########
###### 2E. Plot all TPCs for Arctic species in the same graph (for comparison): Briere ----
##########

# Add an identifying column in each model output dataframe
df.lf.arctic.bri.uni <- df.lf.arctic.bri.uni %>% 
  mutate(type = "Briere uniform")

df.lf.arctic.bri.inf <- df.lf.arctic.bri.inf %>% 
  mutate(type = "Briere informative")


# Combine the three dataframes
df.all <- rbind(df.lf.arctic.bri.uni, df.lf.arctic.bri.inf)

df.all$type <- factor(df.all$type, levels = c( "Briere uniform", "Briere informative"))


# Plot
plot.all <- df.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = mean, color = type), linewidth = 1) +
  geom_point(data = data.lf.arctic, aes(x = temp, y = trait), size = 2) +
  #geom_point(data = data.lf.nonarctic, aes(x = temp, y = trait), size = 2) +
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

# ggsave("figures/lf.arctic.bri.all.png", plot.all, width = 10.3, height = 5.6)

lf.arctic.bri.uni$BUGSoutput$DIC
lf.arctic.bri.inf$BUGSoutput$DIC



##########
###### 3A. Fit lf thermal responses with uniform priors (Arctic): Quadratic ----
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
data <- data.lf.arctic

##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# ##### Run JAGS -----
# lf.arctic.quad.uni <- jags(data = jag.data,
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
# save(lf.arctic.quad.uni, file = "R-scripts/R2jags-objects/lf.arctic.quad.uni.Rdata")

# Read the .Rdata
load("R-scripts/R2jags-objects/lf.arctic.quad.uni.Rdata")


## Diagnostics ----
##### Examine output
lf.arctic.quad.uni$BUGSoutput$summary[1:5,]
mcmcplot(lf.arctic.quad.uni)

# Extract the DIC for future model comparisons
lf.arctic.quad.uni$BUGSoutput$DIC

## Plot data + fit ----
df.lf.arctic.quad.uni <- data.frame(lf.arctic.quad.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.lf.arctic.quad.uni)

##### Plot
plot.lf.arctic.quad.uni <- df.lf.arctic.quad.uni %>% 
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

plot.lf.arctic.quad.uni

# ggsave("figures/lf.arctic.quad.uni.png", plot.lf.arctic.quad.uni, 
#        width = 10.3, height = 5.6)


##########
###### 3B. Fit lf thermal responses for priors (non-Arctic species): Quadratic ----
##########

##### Temp sequence for derived quantity calculations
# For priors - fewer temps for derived calculations makes it go faster
Temp.xs <- seq(0, 45, 0.5)
N.Temp.xs <-length(Temp.xs)


##### Set data
data <- data.lf.nonarctic

##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS -----
# lf.nonarctic.quad.uni <- jags(data = jag.data, 
#                              inits = inits, 
#                              parameters.to.save = parameters, 
#                              model.file = "R-scripts/quad_T.txt",
#                              n.thin = nt, 
#                              n.chains = nc, 
#                              n.burnin = nb, 
#                              n.iter = ni, 
#                              DIC = T, 
#                              working.directory = getwd()
# )


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
plot.df.lf.nonarctic.quad.uni <- df.lf.nonarctic.quad.uni %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "grey", alpha = 0.5) +
  geom_line(aes(y = mean), color = "#868686FF", linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait), size = 2 
             , position = "jitter"
  ) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Mosquito adult lifespan (days)"
  ) +
  theme_bw()

plot.df.lf.nonarctic.quad.uni

# ggsave("figures/lf.nonarctic.quad.uni.png", plot.df.lf.nonarctic.quad.uni, 
#        width = 10.3, height = 5.6)


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
plot.df.lf.nonarctic.quad.uni <- df.lf.nonarctic.quad.uni %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "grey", alpha = 0.5) +
  geom_line(aes(y = mean), color = "#868686FF", linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait), size = 2 
             , position = "jitter"
  ) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Mosquito adult lifespan (days)"
  ) +
  theme_bw()

plot.df.lf.nonarctic.quad.uni

# ggsave("figures/lf.nonarctic.quad.uni.png", plot.df.lf.nonarctic.quad.uni,
#        width = 10.3, height = 5.6)


##########
###### 3C. Fit gamma distributions to lf prior thermal responses: Quadratic ----
##########

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
lf.arctic.prior.cf.dists <- data.frame(q = as.vector(lf.nonarctic.quad.uni$BUGSoutput$sims.list$cf.q),
                                        T0 = as.vector(lf.nonarctic.quad.uni$BUGSoutput$sims.list$cf.T0),
                                        Tm = as.vector(lf.nonarctic.quad.uni$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
lf.arctic.prior.gamma.fits = apply(lf.arctic.prior.cf.dists, 2, 
                                    function(df) fitdistr(df, "gamma")$estimate)


lf.hypers <- lf.arctic.prior.gamma.fits
save(lf.hypers, file = "R-scripts/R2jags-objects/lfhypers.quad.Rsave")


##########
###### 3D. Fit lf thermal responses with data-informed priors (Arctic): Quadratic ----
##########

load("R-scripts/R2jags-objects/lfhypers.quad.Rsave")
lf.arctic.prior.gamma.fits <- lf.hypers

##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)

##### Set data
data <- data.lf.arctic
hypers <- lf.arctic.prior.gamma.fits * 0.1


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
lf.arctic.quad.inf <- jags(data = jag.data,
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
# save(lf.arctic.quad.inf, file = "R-scripts/R2jags-objects/lf.arctic.quad.inf.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/lf.arctic.quad.inf.Rdata")


## Diagnostics ----
##### Examine output
lf.arctic.quad.inf$BUGSoutput$summary[1:5,]
mcmcplot(lf.arctic.quad.inf)

# Extract the DIC for future model comparisons
lf.arctic.quad.inf$BUGSoutput$DIC

## Plot data + fit ----
df.lf.arctic.quad.inf <- data.frame(lf.arctic.quad.inf$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.lf.arctic.quad.inf)

##### Plot
plot.lf.arctic.quad.inf <- df.lf.arctic.quad.inf %>% 
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

plot.lf.arctic.quad.inf

# ggsave("figures/lf.arctic.quad.inf.png", plot.lf.arctic.quad.inf, 
#        width = 10.3, height = 5.6)




##########
###### 3E. Plot all three TPCs in the same graph (for comparison) ----
##########

# Add an identifying column in each model output dataframe
df.lf.arctic.quad.uni <- df.lf.arctic.quad.uni %>% 
  mutate(type = "Quadratic uniform")


df.lf.arctic.quad.inf <- df.lf.arctic.quad.inf %>% 
  mutate(type = "Quadratic informative")


# Combine the three dataframes
df.all <- rbind(df.lf.arctic.quad.uni, df.lf.arctic.quad.inf)

##### Plot
plot.all <- df.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = mean, color = type), linewidth = 1) +
  geom_point(data = data.lf.arctic, aes(x = temp, y = trait), size = 2) +
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

# ggsave("figures/lf.arctic.quad.all.png", plot.all, width = 10.3, height = 5.6)


##### Plot all best fitting TPCs for comparison ----

#### DIC ----
lf.arctic.bri.uni$BUGSoutput$DIC
lf.arctic.bri.inf$BUGSoutput$DIC
lf.arctic.quad.uni$BUGSoutput$DIC
lf.arctic.quad.inf$BUGSoutput$DIC

# Combine the three dataframes
df.all <- rbind(df.lf.arctic.bri.uni, 
                df.lf.arctic.bri.inf, 
                df.lf.arctic.quad.uni)



##### Plot
plot.all <- df.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = mean, color = type), linewidth = 1) +
  geom_point(data = data.lf.arctic, aes(x = temp, y = trait), size = 2) +
  #geom_point(data = data.lf.sierrensis, aes(x = temp, y = trait), size = 2) +
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

# ggsave("figures/lf.arctic.arctic.all.png", plot.all, width = 10.3, height = 5.6)



