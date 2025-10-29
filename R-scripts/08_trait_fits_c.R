## Lilian Chan, University of Guelph
## Arctic vector-borne disease transmission suitability model
##
## Purpose: use Bayesian inference (JAGS) to fit TPCs for Infection efficiency
## (c) using data from Ae. Trivittatus (transmitting Dirofilaria immitis;
## Christensen and Hollander 1978) with uniform priors
##
## 
## Table of content:
##    0. Set-up workspace
##
##    1. MCMC settings for all models
##
##    2. Fitting TPC (Briere)
##        A. Fit c thermal responses with uniform priors
##
##    3. Fitting TPC (Quadratic)
##        A. Fit c thermal responses with uniform priors


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
data <- read_csv("data-processed/TraitData_c.csv")
unique(data$species)

data.c <- data

## Plot raw data
plot.data.c <- data.c %>% 
  mutate(type = "non-Arctic") %>% 
  ggplot() +
  geom_point(aes(x = temp, y = trait, colour = species)) +
  labs(y = "Infection proportion", x = "Temperature ºC") +
  scale_color_discrete(name = "Species", label = "Ae. trivittatus") +
  facet_grid(rows = vars(type)) +
  theme_bw()

plot.data.c

# ggsave("figures/raw_data/plot.data.c.png", plot.data.c, , width = 9.83, height = 6.17)



##########
###### 1. MCMC settings for all models ----
##########

# Number of posterior dist elements = [(ni - nb) / nt] * nc = [(45000 - 5000) / 8] * 3 = 15000
ni <- 45000 # number of iterations in each chain
nb <- 5000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves cery nt iterations in each chain
nc <- 3 # number of chains


## Model to fit this 
sink("R-scripts/briereprob_c.txt")
cat("
    model{
    
    ## Priors
    cf.q ~ dunif(prior[1,1], prior[2,1])
    cf.T0 ~ dunif(prior[1,2], prior[2,2])
    cf.Tm ~ dunif(prior[1,3], prior[2,3])
    cf.sigma ~ dunif(0, 1000)
    cf.tau <- 1 / (cf.sigma * cf.sigma)
    
    ## Likelihood
    for(i in 1:N.obs){
    trait.mu[i] <- cf.q * temp[i] * (temp[i] - cf.T0) * sqrt((cf.Tm - temp[i]) * (cf.Tm > temp[i])) * (cf.T0 < temp[i])
    trait[i] ~ dnorm(trait.mu[i], cf.tau)
    }
    
    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- cf.q * Temp.xs[i] * (Temp.xs[i] - cf.T0) * sqrt((cf.Tm - Temp.xs[i]) * (cf.Tm > Temp.xs[i])) * (cf.T0 < Temp.xs[i]) * (cf.q * Temp.xs[i] * (Temp.xs[i] - cf.T0) * sqrt((cf.Tm - Temp.xs[i]) * (cf.Tm > Temp.xs[i])) < 1) + (cf.q * Temp.xs[i] * (Temp.xs[i] - cf.T0) * sqrt((cf.Tm - Temp.xs[i]) * (cf.Tm > Temp.xs[i])) > 1)
    }
    
    } # close model
    ",fill=T)
sink()

## Set priors
prior <- data.frame(q = c(0, 1),
                    T0 = c(0, 15),
                    Tm = c(35, 50))

##########
###### 2A. Fit c thermal responses with uniform priors (Arctic): Briere ----
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
data <- data.c

##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, 
                 N.Temp.xs = N.Temp.xs, prior = prior)

##### Run JAGS
c.nonarctic.bri.uni <- jags(
  data = jag.data,
  inits = inits,
  parameters.to.save = parameters,
  model.file = "R-scripts/briereprob_c.txt",
  n.thin = nt,
  n.chains = nc,
  n.burnin = nb,
  n.iter = ni,
  DIC = T,
  working.directory = getwd()
)

## Save the model as Rdata 
# save(c.nonarctic.bri.uni, file = "R-scripts/R2jags-objects/c.nonarctic.bri.uni.Rdata")

# Read the .Rdata
load("R-scripts/R2jags-objects/c.nonarctic.bri.uni.Rdata")


## Diagnostics ----
##### Examine output
c.nonarctic.bri.uni$BUGSoutput$summary[1:5,]
mcmcplot(c.nonarctic.bri.uni)

# Extract the DIC for future model comparisons
c.nonarctic.bri.uni$BUGSoutput$DIC

## Plot data + fit ----
df.c.nonarctic.bri.uni <- data.frame(c.nonarctic.bri.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.c.nonarctic.bri.uni)

##### Plot
plot.c.nonarctic.bri.uni <- df.c.nonarctic.bri.uni %>%
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
  labs(x = expression(paste("Temperature (", degree, "C)")), y = "Infection proportion") +
  scale_colour_manual(name = "Species", labels = "Ae. trivittatus", values = "black") +
  theme_bw()

plot.c.nonarctic.bri.uni

# ggsave("figures/c.nonarctic.bri.uni.png", plot.c.nonarctic.bri.uni,
#        width = 10.3, height = 5.6)



##########
###### 2E. Plot all TPCs for Arctic species in the same graph (for comparison): Briere ----
##########

# Add an identifying column in each model output dataframe
df.c.nonarctic.bri.uni <- df.c.nonarctic.bri.uni %>% 
  mutate(type = "Briere uniform")



##########
###### 3A. Fit c thermal responses with uniform priors (Arctic): Quadratic ----
##########

sink("R-scripts/quadprob_c.txt")
cat("
    model{
    
    ## Priors
    cf.q ~ dunif(prior[1,1], prior[2,1])
    cf.T0 ~ dunif(prior[1,2], prior[2,2])
    cf.Tm ~ dunif(prior[1,3], prior[2,3])
    cf.sigma ~ dunif(0, 1000)
    cf.tau <- 1 / (cf.sigma * cf.sigma)
    
    ## Likelihood
    for(i in 1:N.obs){
    trait.mu[i] <- -1 * cf.q * (temp[i] - cf.T0) * (temp[i] - cf.Tm) * (cf.Tm > temp[i]) * (cf.T0 < temp[i])
    trait[i] ~ dnorm(trait.mu[i], cf.tau)
    }
    
    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- (-1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) * (cf.Tm > Temp.xs[i]) * (cf.T0 < Temp.xs[i])) * (-1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) < 1) + (-1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) > 1)
    }
    
    } # close model
    ",fill=T)
sink()

## Set priors
prior <- data.frame(q = c(0, 1),
                    T0 = c(0, 15),
                    Tm = c(35, 50))

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
data <- data.c

##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp

##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, 
                 N.Temp.xs = N.Temp.xs, prior = prior)

# ##### Run JAGS -----
c.nonarctic.quad.uni <- jags(data = jag.data,
                              inits = inits,
                              parameters.to.save = parameters,
                              model.file = "R-scripts/quadprob_c.txt",
                              n.thin = nt,
                              n.chains = nc,
                              n.burnin = nb,
                              n.iter = ni,
                              DIC = T,
                              working.directory = getwd()
)

## Save the model as Rdata 
# save(c.nonarctic.quad.uni, file = "R-scripts/R2jags-objects/c.nonarctic.quad.uni.Rdata")

# Read the .Rdata
load("R-scripts/R2jags-objects/c.nonarctic.quad.uni.Rdata")


## Diagnostics ----
##### Examine output
c.nonarctic.quad.uni$BUGSoutput$summary[1:5,]
mcmcplot(c.nonarctic.quad.uni)

# Extract the DIC for future model comparisons
c.nonarctic.quad.uni$BUGSoutput$DIC

## Plot data + fit ----
df.c.nonarctic.quad.uni <- data.frame(c.nonarctic.quad.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)

head(df.c.nonarctic.quad.uni)

##### Plot
plot.c.nonarctic.quad.uni <- df.c.nonarctic.quad.uni %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "#4363d8", alpha = 0.5) +
  geom_line(aes(y = mean), color = "blue", linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait, colour = species), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Infection proportion"
  ) +
  scale_colour_manual(name = "Species", labels = "Ae. trivittatus", values = "black") +
  theme_bw()

plot.c.nonarctic.quad.uni

# ggsave("figures/c.nonarctic.quad.uni.png", plot.c.nonarctic.quad.uni, 
#        width = 10.3, height = 5.6)



##########
###### 3E. Plot all three TPCs in the same graph (for comparison) ----
##########

# Add an identifying column in each model output dataframe
df.c.nonarctic.quad.uni <- df.c.nonarctic.quad.uni %>% 
  mutate(type = "Quadratic uniform")




##### Plot all best fitting TPCs for comparison ----

#### DIC ----
c.nonarctic.bri.uni$BUGSoutput$DIC
c.nonarctic.quad.uni$BUGSoutput$DIC

# Combine the three dataframes
df.all <- rbind(df.c.nonarctic.bri.uni, 
                df.c.nonarctic.quad.uni)



##### Plot
plot.all <- df.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = mean, color = type), linewidth = 1) +
  geom_point(data = data.c, aes(x = temp, y = trait), size = 2) +
  #geom_point(data = data.c.sierrensis, aes(x = temp, y = trait), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Infection proportion"
  ) +
  # Customize the colours
  # scale_fill_jco() +
  # scale_color_jco() +
  # scale_fill_brewer(palette = "Accent") +
  # scale_color_brewer(palette = "Accent") +
  theme_bw()

plot.all

# ggsave("figures/c.nonarctic.all.png", plot.all, width = 10.3, height = 5.6)
