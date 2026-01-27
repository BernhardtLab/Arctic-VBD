## Lilian Chan, University of Guelph
## Arctic vector-borne disease transmission suitability model
##
## Purpose: use Bayesian inference (JAGS) to fit TPCs for parasite development 
## rate (PDR) for Arctic nematode using data from Varestrongylus eleguneniensis 
## (Kafle et al. 2018) and from Setaria tundra (Laaksonen et al. 2009).
##     1) with uniform priors; and 
##     2) with data-informed priors from Dirofilaria immitis and Wuchereria bancrofti
##
## Varestrongylus eleguneniensis is a nematode infecting caribou and muskoxen in
## the Canadian Arctic. It is transmitted by gastropod.; Dirofilaria immitis is
## a filarial worm causing dirofilariasis in dogs. Wuchereria bancrofti is a 
## filarial nematode causing of lymphatic filariasis in tropical regions.
## 
##
## Table of content:
##    0. Set-up workspace
##
##    1. MCMC settings for all models
##
##    2. Fitting TPC (Briere)
##        A. Fit PDR thermal responses with uniform priors (Arctic species)
##        B. Fit PDR thermal responses (with random effects) for priors (non-Arctic species)
##        C. Fit gamma distributions to PDR prior thermal responses
##        D. Fit PDR thermal responses with data-informed priors (Arctic)
##        E. Plot all three TPCs in the same graph (for comparison)
##
##    3. Fitting TPC (Quadratic)
##        A. Fit PDR thermal responses with uniform priors (Arctic)
##        B. Fit PDR thermal responses (with random effects) for priors (non-Arctic species)
##        C. Fit gamma distributions to PDR prior thermal responses
##        D. Fit PDR thermal responses with data-informed priors (Arctic)
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
data <- read_csv("data-processed/TraitData_PDR.csv")
unique(data$species)

## Convert development time (1/PDR) to development rate (PDR)
data <- data %>% 
  mutate(trait = 1/trait) %>% 
  mutate(trait_name = "PDR") 
     
# Subset data
## Arctic species
data.PDR.arctic <- subset(data, type == "Arctic")

## Non-Arctic species
data.PDR.nonarctic <- subset(data, type == "non-Arctic")


# Plot the raw data
plot.data.PDR <- data %>% 
  ggplot(aes(x = temp, y = trait)) +
  geom_point(aes(colour = citation)) +
  #geom_line(aes(colour = citation)) +
  labs(y = "Parasite development rate (1/days)", x = "Temperature ºC") +
  scale_colour_discrete(name = "Species", labels = c("D. immitis (in Ae. Trivittatus)",
                                                     "V. eleguneniensis",
                                                     "S. tundra",
                                                     "W. bancrofti (in Ae. polynesiensis)",
                                                     "D. immitis (in Ae. aegypti)"
  )) +
  facet_grid(rows = vars(type)) +
  theme_bw()

plot.data.PDR

# ggsave("figures/raw_data/plot.data.PDR.png", plot.data.PDR, , width = 9.83, height = 6.17)



##########
###### 1. MCMC settings for all models ----
##########

# Number of posterior dist elements = [(ni - nb) / nt] * nc = [(45000 - 5000) / 8] * 3 = 15000
ni <- 45000 # number of iterations in each chain
nb <- 5000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 3 # number of chains


##########
###### 2A. Fit PDR thermal responses with uniform priors (Arctic): Briere ----
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
data <- data.PDR.arctic

##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp


##### Set priors
prior <- data.frame(q = c(0, 1),
                    T0 = c(0, 20),
                    Tm = c(20, 45)
)

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, 
                 N.Temp.xs = N.Temp.xs, prior = prior)

##### Run JAGS
PDR.arctic.bri.uni <- jags(
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
# save(PDR.arctic.bri.uni, file = "R-scripts/R2jagsz-objects/PDR.arctic.bri.uni.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/PDR.arctic.bri.uni.Rdata")


## Diagnostics ----
##### Examine output
PDR.arctic.bri.uni$BUGSoutput$summary[1:5,]
mcmcplot(PDR.arctic.bri.uni)

# Extract the DIC for future model comparisons
PDR.arctic.bri.uni$BUGSoutput$DIC

## Plot data + fit ----
df.PDR.arctic.bri.uni <- data.frame(PDR.arctic.bri.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)

head(df.PDR.arctic.bri.uni)

##### Plot
plot.PDR.arctic.bri.uni <- df.PDR.arctic.bri.uni %>%
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.),
              fill = "#4363d8",
              alpha = 0.5) +
  geom_line(aes(y = X50.), color = "blue", linewidth = 1) +
  geom_point(data = data,
             aes(x = temp, y = trait, colour = species),
             size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) +
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(x = expression(paste("Temperature (", degree, "C)")), y = "Development rate (days-1)") +
  # Customize legend
  scale_color_discrete(name = "Species",
                       labels = c("V. eleguneniensis", "S. tundra")) +
  theme_bw()

plot.PDR.arctic.bri.uni

# ggsave("figures/PDR.arctic.bri.uni.png", plot.PDR.arctic.bri.uni,
#        width = 10.3, height = 5.6)



##########
###### 2B. Fit PDR thermal responses for priors (non-Arctic species): Briere ----
##########

##### Temp sequence for derived quantity calculations
# For priors - fewer temps for derived calculations makes it go faster
Temp.xs <- seq(0, 45, 0.5)
N.Temp.xs <-length(Temp.xs)


##### Set data
data <- data.PDR.nonarctic

##### Create a unique id for each species-study combination
data <- data %>% 
  group_by(species, host.species, citation) %>% 
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
PDR.nonarctic.bri.uni.raneff <- jags(
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
save(PDR.nonarctic.bri.uni.raneff, file = "R-scripts/R2jags-objects/PDR.nonarctic.bri.uni.raneff.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/PDR.nonarctic.bri.uni.raneff.Rdata")


## Diagnostics ----
##### Examine output
PDR.nonarctic.bri.uni.raneff$BUGSoutput$summary[1:8,]
mcmcplot(PDR.nonarctic.bri.uni.raneff)

# Extract the DIC for future model comparisons
PDR.nonarctic.bri.uni.raneff$BUGSoutput$DIC


## Plot data + fit ----
df.PDR.nonarctic.bri.uni.raneff <- data.frame(PDR.nonarctic.bri.uni.raneff$BUGSoutput$summary)[-(1:8),]

## Extract the model prediction
## Overall curve
df.PDR.nonarctic.bri.uni.raneff.pop <- df.PDR.nonarctic.bri.uni.raneff %>% 
  filter(grepl("z.trait.mu.pred.pop", rownames(df.PDR.nonarctic.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)


## Unique ID 1: W. bancrofti in Ae. polynesiensis
df.PDR.nonarctic.bri.uni.1 <- df.PDR.nonarctic.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[1,*]"), rownames(df.PDR.nonarctic.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 1)

## Unique ID 2: D. immitis in Ae. aegypti
df.PDR.nonarctic.bri.uni.2 <- df.PDR.nonarctic.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[2,*]"), rownames(df.PDR.nonarctic.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 2)

## Unique ID 3: D. immitis in Ae. trivittatus
df.PDR.nonarctic.bri.uni.3 <- df.PDR.nonarctic.bri.uni.raneff %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[3,*]"), rownames(df.PDR.nonarctic.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(unique_id = 3)

## Combine the model prediciton of all three unique groups into a dataframe
df.PDR.nonarctic.bri.uni.raneff.sp <- rbind(df.PDR.nonarctic.bri.uni.1,
                                            df.PDR.nonarctic.bri.uni.2,
                                            df.PDR.nonarctic.bri.uni.3) 

## Change unique_id into factor type
df.PDR.nonarctic.bri.uni.raneff.sp$unique_id <- as.factor(df.PDR.nonarctic.bri.uni.raneff.sp$unique_id)


##### Plot
plot.PDR.nonarctic.bri.uni.raneff <- ggplot(data = df.PDR.nonarctic.bri.uni.raneff.pop, 
                                            aes(x = temp)) +
  ## Overall TPC
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.),
              fill = "grey",
              alpha = 0.5) +
  ## a separate TPC (and credible interval) for each unique group
  # geom_ribbon(data = df.PDR.nonarctic.bri.uni.raneff.sp, aes(ymin = X2.5., ymax = X97.5., fill = unique_id),
  #             alpha = 0.5) +
  geom_point(data = data,
             aes(x = temp, y = trait, colour = as.factor(unique_id)),
             size = 2) +
  geom_line(data = df.PDR.nonarctic.bri.uni.raneff.sp, aes(y = X50., color = unique_id)) +
  geom_line(aes(y = X50.), color = "black", linewidth = 1) +
  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), y = "Development rate (days-1)") +
  # Customize legend
  scale_colour_discrete(name = element_blank(),
                        labels = c("W. bancrofti in Ae. polynesiensis",
                                   "D. immitis in Ae. aegypti", 
                                   "D. immitis in Ae. trivittatus")) +
  theme_bw()


plot.PDR.nonarctic.bri.uni.raneff

# ggsave("figures/PDR.nonarctic.bri.uni.raneff.png", plot.PDR.nonarctic.bri.uni.raneff,
#        width = 10.3, height = 5.6)


##### No random effect -----

##### Temp sequence for derived quantity calculations
# For priors - fewer temps for derived calculations makes it go faster
Temp.xs <- seq(0, 45, 0.5)
N.Temp.xs <-length(Temp.xs)

##### Set data
data <- data.PDR.nonarctic


##### Set priors
prior <- data.frame(q = c(0, 1),
                    T0 = c(0, 20),
                    Tm = c(20, 45)
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
                 N.Temp.xs = N.Temp.xs, prior = prior)

##### Run JAGS -----
PDR.nonarctic.bri.uni <- jags(data = jag.data,
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
save(PDR.nonarctic.bri.uni, file = "R-scripts/R2jags-objects/PDR.nonarctic.bri.uni.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/PDR.nonarctic.bri.uni.Rdata")


## Diagnostics ----
##### Examine output
PDR.nonarctic.bri.uni$BUGSoutput$summary[1:5,]
mcmcplot(PDR.nonarctic.bri.uni)

# Extract the DIC for future model comparisons
PDR.nonarctic.bri.uni$BUGSoutput$DIC

## Plot data + fit ----
df.PDR.nonarctic.bri.uni <- data.frame(PDR.nonarctic.bri.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)

head(df.PDR.nonarctic.bri.uni)

##### Plot
plot.df.PDR.nonarctic.bri.uni <- df.PDR.nonarctic.bri.uni %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "grey", alpha = 0.5) +
  geom_line(aes(y = X50.), color = "#868686FF", linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait), size = 2 
             , position = "jitter"
  ) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Development rate (days-1)"
  ) +
  theme_bw()

plot.df.PDR.nonarctic.bri.uni

# ggsave("figures/PDR.nonarctic.bri.uni.png", plot.df.PDR.nonarctic.bri.uni, 
#        width = 10.3, height = 5.6)



##########
###### 2C. Fit gamma distributions to PDR prior thermal responses: Briere ----
##########

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
PDR.arctic.prior.cf.dists <- data.frame(q = as.vector(PDR.nonarctic.bri.uni$BUGSoutput$sims.list$cf.q),
                                        T0 = as.vector(PDR.nonarctic.bri.uni$BUGSoutput$sims.list$cf.T0),
                                          Tm = as.vector(PDR.nonarctic.bri.uni$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
PDR.arctic.prior.gamma.fits = apply(PDR.arctic.prior.cf.dists, 2, 
                                    function(df) fitdistr(df, "gamma")$estimate)


PDR.hypers <- PDR.arctic.prior.gamma.fits
# save(PDR.hypers, file = "R-scripts/R2jags-objects/PDRhypers.bri.Rsave")



##########
###### 2D. Fit PDR thermal responses with data-informed priors (Arctic): Briere ----
##########

load("R-scripts/R2jags-objects/PDRhypers.bri.Rsave")
PDR.arctic.prior.gamma.fits <- PDR.hypers


##### Set data
data <- data.PDR.arctic
hypers <- PDR.arctic.prior.gamma.fits * 0.1


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
# save(PDR.arctic.bri.inf, file = "R-scripts/R2jags-objects/PDR.arctic.bri.inf.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/PDR.arctic.bri.inf.Rdata")


## Diagnostics ----
##### Examine output
PDR.arctic.bri.inf$BUGSoutput$summary[1:5,]
mcmcplot(PDR.arctic.bri.inf)

# Extract the DIC for future model comparisons
PDR.arctic.bri.inf$BUGSoutput$DIC

## Plot data + fit ----
df.PDR.arctic.bri.inf <- data.frame(PDR.arctic.bri.inf$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)

head(df.PDR.arctic.bri.inf)

##### Plot
plot.PDR.arctic.bri.inf <- df.PDR.arctic.bri.inf %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "#4363d8", alpha = 0.5) +
  geom_line(aes(y = X50.), color = "blue", linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait, colour = species), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Parasite development rate (days-1)"
  ) +
  scale_colour_discrete(name = "Species", labels = c("V. eleguneniensis", 
                                                     "S. tundra"
  )) +
  theme_bw()

plot.PDR.arctic.bri.inf

# ggsave("figures/PDR.arctic.bri.inf.png", plot.PDR.arctic.bri.inf, 
#        width = 10.3, height = 5.6)



##########
###### 2E. Plot all TPCs for Arctic species in the same graph (for comparison): Briere ----
##########

# Add an identifying column in each model output dataframe
df.PDR.arctic.bri.uni <- df.PDR.arctic.bri.uni %>% 
  mutate(type = "Briere uniform")

df.PDR.arctic.bri.inf <- df.PDR.arctic.bri.inf %>% 
  mutate(type = "Briere informative")


# Combine the three dataframes
df.all <- rbind(df.PDR.arctic.bri.uni, df.PDR.arctic.bri.inf)

df.all$type <- factor(df.all$type, levels = c( "Briere uniform", "Briere informative"))


# Plot
plot.all <- df.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = X50., color = type), linewidth = 1) +
  geom_point(data = data.PDR.arctic, aes(x = temp, y = trait), size = 2) +
  #geom_point(data = data.PDR.nonarctic, aes(x = temp, y = trait), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Development rate (days-1)"
  ) +
  # Customize the colours
  ## ribbon
  scale_fill_manual(values = c("Briere uniform" = "grey",
                               "Briere informative" = "#4363d8",
                               "Briere informative raneff" = "pink")) +
  
  ## line
  scale_color_manual(values = c("Briere uniform" = "#868686FF",
                                "Briere informative" = "blue",
                                "Briere informative raneff" = "red")) +
  theme_bw()

plot.all

#ggsave("figures/PDR.arctic.bri.all.png", plot.all,  width = 10.3, height = 5.6)

PDR.arctic.bri.uni$BUGSoutput$DIC
PDR.arctic.bri.inf$BUGSoutput$DIC



##########
###### 3A. Fit PDR thermal responses with uniform priors (Arctic): Quadratic ----
##########

##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)

##### Set data
data <- data.PDR.arctic

##### Set priors
prior <- data.frame(q = c(0, 1),
                    T0 = c(0, 20),
                    Tm = c(20, 45)
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
                 N.Temp.xs = N.Temp.xs, prior = prior)

# ##### Run JAGS -----
# PDR.arctic.quad.uni <- jags(data = jag.data,
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
# save(PDR.arctic.quad.uni, file = "R-scripts/R2jags-objects/PDR.arctic.quad.uni.Rdata")

# Read the .Rdata
load("R-scripts/R2jags-objects/PDR.arctic.quad.uni.Rdata")


## Diagnostics ----
##### Examine output
PDR.arctic.quad.uni$BUGSoutput$summary[1:5,]
mcmcplot(PDR.arctic.quad.uni)

# Extract the DIC for future model comparisons
PDR.arctic.quad.uni$BUGSoutput$DIC

## Plot data + fit ----
df.PDR.arctic.quad.uni <- data.frame(PDR.arctic.quad.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)

head(df.PDR.arctic.quad.uni)

##### Plot
plot.PDR.arctic.quad.uni <- df.PDR.arctic.quad.uni %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "#4363d8", alpha = 0.5) +
  geom_line(aes(y = X50.), color = "blue", linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Development rate (days-1)"
  ) +
  theme_bw()

plot.PDR.arctic.quad.uni

# ggsave("figures/PDR.arctic.quad.uni.png", plot.PDR.arctic.quad.uni, 
#        width = 10.3, height = 5.6)


##########
###### 3B. Fit PDR thermal responses (with random effects) for priors (non-Arctic species): Quadratic ----
##########

##### Temp sequence for derived quantity calculations
# For priors - fewer temps for derived calculations makes it go faster
Temp.xs <- seq(0, 45, 0.5)
N.Temp.xs <-length(Temp.xs)

##### Set data
data <- data.PDR.nonarctic


##### Set priors
prior <- data.frame(q = c(0, 1),
                    T0 = c(0, 20),
                    Tm = c(20, 45)
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
                 N.Temp.xs = N.Temp.xs, prior = prior)

##### Run JAGS -----
PDR.nonarctic.quad.uni <- jags(data = jag.data,
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
# save(PDR.nonarctic.quad.uni, file = "R-scripts/R2jags-objects/PDR.nonarctic.quad.uni.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/PDR.nonarctic.quad.uni.Rdata")


## Diagnostics ----
##### Examine output
PDR.nonarctic.quad.uni$BUGSoutput$summary[1:5,]
mcmcplot(PDR.nonarctic.quad.uni)

# Extract the DIC for future model comparisons
PDR.nonarctic.quad.uni$BUGSoutput$DIC

## Plot data + fit ----
df.PDR.nonarctic.quad.uni <- data.frame(PDR.nonarctic.quad.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)

head(df.PDR.nonarctic.quad.uni)

##### Plot
plot.df.PDR.nonarctic.quad.uni <- df.PDR.nonarctic.quad.uni %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "grey", alpha = 0.5) +
  geom_line(aes(y = X50.), color = "#868686FF", linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait), size = 2 
             , position = "jitter"
  ) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Development rate (days-1)"
  ) +
  theme_bw()

plot.df.PDR.nonarctic.quad.uni

# ggsave("figures/PDR.nonarctic.quad.uni.png", plot.df.PDR.nonarctic.quad.uni, 
#        width = 10.3, height = 5.6)


##########
###### 3C. Fit gamma distributions to PDR prior thermal responses: Quadratic ----
##########

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
PDR.arctic.prior.cf.dists <- data.frame(q = as.vector(PDR.nonarctic.quad.uni$BUGSoutput$sims.list$cf.q),
                                        T0 = as.vector(PDR.nonarctic.quad.uni$BUGSoutput$sims.list$cf.T0),
                                        Tm = as.vector(PDR.nonarctic.quad.uni$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
PDR.arctic.prior.gamma.fits = apply(PDR.arctic.prior.cf.dists, 2, 
                                      function(df) fitdistr(df, "gamma")$estimate)


PDR.hypers <- PDR.arctic.prior.gamma.fits
save(PDR.hypers, file = "R-scripts/R2jags-objects/PDRhypers.quad.Rsave")


##########
###### 3D. Fit PDR thermal responses with data-informed priors (Arctic): Quadratic ----
##########

load("R-scripts/R2jags-objects/PDRhypers.quad.Rsave")
PDR.arctic.prior.gamma.fits <- PDR.hypers

##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)

##### Set data
data <- data.PDR.arctic
hypers <- PDR.arctic.prior.gamma.fits * 0.1


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
save(PDR.arctic.quad.inf, file = "R-scripts/R2jags-objects/PDR.arctic.quad.inf.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/PDR.arctic.quad.inf.Rdata")


## Diagnostics ----
##### Examine output
PDR.arctic.quad.inf$BUGSoutput$summary[1:5,]
mcmcplot(PDR.arctic.quad.inf)

# Extract the DIC for future model comparisons
PDR.arctic.quad.inf$BUGSoutput$DIC

## Plot data + fit ----
df.PDR.arctic.quad.inf <- data.frame(PDR.arctic.quad.inf$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)

head(df.PDR.arctic.quad.inf)

##### Plot
plot.PDR.arctic.quad.inf <- df.PDR.arctic.quad.inf %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "pink", alpha = 0.5) +
  geom_line(aes(y = X50.), color = "red", linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Development rate (days-1)"
  ) +
  theme_bw()

plot.PDR.arctic.quad.inf

ggsave("figures/PDR.arctic.quad.inf.png", plot.PDR.arctic.quad.inf,
       width = 10.3, height = 5.6)



##########
###### 3E. Plot all three TPCs in the same graph (for comparison) ----
##########

# Add an identifying column in each model output dataframe
df.PDR.arctic.quad.uni <- df.PDR.arctic.quad.uni %>% 
  mutate(type = "Quadratic uniform")


df.PDR.arctic.quad.inf <- df.PDR.arctic.quad.inf %>% 
  mutate(type = "Quadratic informative")


# Combine the three dataframes
df.all <- rbind(df.PDR.arctic.quad.uni, df.PDR.arctic.quad.inf)

##### Plot
plot.all <- df.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = X50., color = type), linewidth = 1) +
  geom_point(data = data.PDR.arctic, aes(x = temp, y = trait), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Development rate (days-1)"
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

# ggsave("figures/PDR.arctic.quad.all.png", plot.all, width = 10.3, height = 5.6)


##### Plot all best fitting TPCs for comparison ----

#### DIC ----
PDR.arctic.bri.uni$BUGSoutput$DIC
PDR.arctic.bri.inf$BUGSoutput$DIC # This is the best fitting TPC
PDR.arctic.quad.uni$BUGSoutput$DIC
PDR.arctic.quad.inf$BUGSoutput$DIC

# Combine the three dataframes
df.all <- rbind(#df.PDR.arctic.bri.uni, 
                df.PDR.arctic.bri.inf, 
                #df.PDR.arctic.quad.uni,
                df.PDR.arctic.quad.inf)



##### Plot
plot.all <- df.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = X50., color = type), linewidth = 1) +
  geom_point(data = data.PDR.arctic, aes(x = temp, y = trait), size = 2) +
  #geom_point(data = data.PDR.sierrensis, aes(x = temp, y = trait), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Development rate (days-1)"
  ) +
  # Customize the colours
  scale_fill_jco() +
  scale_color_jco() +
  # scale_fill_brewer(palette = "Accent") +
  # scale_color_brewer(palette = "Accent") +
  theme_bw()

plot.all

# ggsave("figures/PDR.arctic.all.png", plot.all, width = 10.3, height = 5.6)


##########
###### 4. Process and save model output for plotting ----
##########

## Analyze TPC model
PDR.TPC.analysis <- extractTPC(PDR.arctic.bri.inf, "PDR", Temp.xs)
PDR.predictions.summary <- PDR.TPC.analysis[[1]]
PDR.params.summary <- PDR.TPC.analysis[[2]]
PDR.params.fullposts <- PDR.TPC.analysis[[3]]

write_csv(PDR.predictions.summary, "data-processed/PDR.predictions.summary.csv")
write_csv(PDR.params.summary, "data-processed/PDR.params.summary.csv")
write_csv(PDR.params.fullposts, "data-processed/PDR.params.fullposts.csv")

