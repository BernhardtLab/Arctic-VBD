## Lilian Chan, University of Guelph
## Arctic vector-borne disease transmission suitability model
##
## Purpose: Fit thermal performance curves (TPCs) for vector competence (bc) 
## using Bayesian inference (JAGS).
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
## data-processed/TraitData_bc.csv -
##     Synthesized published trait data for bc
##
## Outputs: 
## R-scripts/R2jags-objects/best-fitting-mods/bc.nonarctic.mod.Rdata -
##     Best-fitting TPC models 
##
## data-processed/bc/bc.nonarctic.predictions.summary.csv - 
##     Posterior summary of TPC predictions across temperatures
##
## data-processed/bc/bc.nonarctic.params.summary.csv -
##     Summary statistics of TPC parameters
##
## data-processed/bc/bc.nonarctic.params.fullposts.csv -
##     Full posterior distributions for TPC parameters



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
data.all <- read_csv("data-processed/TraitData_bc.csv")
unique(data.all$host_species)

## Non-Arctic species
data.bc.nonarctic <- subset(data.all, type == "non-Arctic")

## Plot raw data
plot.data.bc.nonarctic <- data.all %>% 
  ggplot(aes(x = temp, y = trait, colour = host_species)) +
  geom_point() +
  labs(y = "Proportion", x = "Temperature ºC") +
  scale_color_discrete(name = "Species") +
  facet_grid(rows = vars(type)) +
  theme_bw()

plot.data.bc.nonarctic



# 1. MCMC Settings for all models ----------------------------------------------

# Number of posterior dist elements = [(ni - nb) / nt] * nc = [(450000 - 50000) / 100] * 3 = 12000
ni <- 450000 # number of iterations in each chain
nb <- 50000 # number of 'burn in' iterations to discard
nt <- 100 # thinning rate - jags saves every nt iterations in each chain
nc <- 3 # number of chains


# 2. Fitting TPC (Briere) ------------------------------------------------------

##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)


##### Set data
data <- data.bc.nonarctic


## Set priors
prior <- data.frame(q = c(0, 1),
                    T0 = c(0, 20),
                    Tm = c(25, 45))


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


##### Run JAGS
set.seed(123) # for reproducibility
bc.nonarctic.bri.uni <- jags(data = jag.data,
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
save(bc.nonarctic.bri.uni, file = "R-scripts/R2jags-objects/all-mods/bc.nonarctic.bri.uni.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/all-mods/bc.nonarctic.bri.uni.Rdata")


## Diagnostics
##### Examine output
bc.nonarctic.bri.uni$BUGSoutput$summary[c("cf.T0", "cf.Tm", "cf.q", "cf.sigma", "deviance"),]
mcmcplot(bc.nonarctic.bri.uni, parms = c("cf.T0", "cf.Tm", "cf.q", "cf.sigma", "deviance"))

# Extract the DIC for future model comparisons
bc.nonarctic.bri.uni$BUGSoutput$DIC

## Plot data + fit
df.bc.nonarctic.bri.uni <- data.frame(bc.nonarctic.bri.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5.,X50., X97.5.)

head(df.bc.nonarctic.bri.uni)


##### Plot
plot.bc.nonarctic.bri.uni <- df.bc.nonarctic.bri.uni %>%
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "#4363d8", alpha = 0.5) +
  geom_line(aes(y = X50.), color = "blue", linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait), size = 2) +

  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Proportion"
       ) +
  theme_bw()

plot.bc.nonarctic.bri.uni

ggsave("figures/bc.nonarctic.bri.uni.png", plot.bc.nonarctic.bri.uni,
       width = 10.3, height = 5.6)




# 3. Fitting TPC (quadratic) ---------------------------------------------------


##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)


##### Set data
data <- data.bc.nonarctic


## Set priors
prior <- data.frame(q = c(0, 1),
                    T0 = c(0, 20),
                    Tm = c(25, 45))


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


# ##### Run JAGS
set.seed(123) # for reproducibility
bc.nonarctic.quad.uni <- jags(data = jag.data,
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
save(bc.nonarctic.quad.uni, file = "R-scripts/R2jags-objects/all-mods/bc.nonarctic.quad.uni.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/all-mods/bc.nonarctic.quad.uni.Rdata")


## Diagnostics
##### Examine output
bc.nonarctic.quad.uni$BUGSoutput$summary[c("cf.T0", "cf.Tm", "cf.q", "cf.sigma", "deviance"),]
mcmcplot(bc.nonarctic.quad.uni, parms = c("cf.T0", "cf.Tm", "cf.q", "cf.sigma", "deviance"))


# Extract the DIC for future model comparisons
bc.nonarctic.quad.uni$BUGSoutput$DIC

## Plot data + fit
df.bc.nonarctic.quad.uni <- data.frame(bc.nonarctic.quad.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5.,X50., X97.5.)

head(df.bc.nonarctic.quad.uni)

##### Plot
plot.bc.nonarctic.quad.uni <- df.bc.nonarctic.quad.uni %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "#4363d8", alpha = 0.5) +
  geom_line(aes(y = X50.), color = "blue", linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait), size = 2) +

  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Proportion"
    ) +
  theme_bw()

plot.bc.nonarctic.quad.uni

ggsave("figures/bc.nonarctic.quad.uni.png", plot.bc.nonarctic.quad.uni,
       width = 10.3, height = 5.6)



# 4. Compare model fit between Briere and Quadratic models ---------------------

##### Find best fitting model #####
# Add an identifying column in each model output dataframe
df.bc.nonarctic.bri.uni <- df.bc.nonarctic.bri.uni %>% 
  mutate(type = "briere")

df.bc.nonarctic.quad.uni <- df.bc.nonarctic.quad.uni %>% 
  mutate(type = "quadratic")

# Combine the two dataframes
df.all <- bind_rows(df.bc.nonarctic.bri.uni, df.bc.nonarctic.quad.uni)

##### Plot
plot.all <- df.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = X50., color = type), linewidth = 1) +
  geom_point(data = data.bc.nonarctic, aes(x = temp, y = trait), size = 2) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Proportion"
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

ggsave("figures/bc.bri.quad.png", plot.all, width = 10.3, height = 5.6)


#### DIC
bc.nonarctic.bri.uni$BUGSoutput$DIC
bc.nonarctic.quad.uni$BUGSoutput$DIC # This is the best fitting TPC


# Save best-fitting TPC in a separate folder
bc.nonarctic.mod <- bc.nonarctic.quad.uni

## Save the model as Rdata 
save(bc.nonarctic.mod, file = "R-scripts/R2jags-objects/best-fitting-mods/bc.nonarctic.mod.Rdata")


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
bc.TPC.analysis <- extractTPC(bc.nonarctic.quad.uni, "bc", Temp.xs)
bc.nonarctic.predictions.summary <- bc.TPC.analysis[[1]]
bc.nonarctic.params.summary <- bc.TPC.analysis[[2]]
bc.nonarctic.params.fullposts <- bc.TPC.analysis[[3]]

write_csv(bc.nonarctic.predictions.summary, "data-processed/bc/bc.nonarctic.predictions.summary.csv")
write_csv(bc.nonarctic.params.summary, "data-processed/bc/bc.nonarctic.params.summary.csv")
write_csv(bc.nonarctic.params.fullposts, "data-processed/bc/bc.nonarctic.params.fullposts.csv")


##### Briere model #####
Temp.xs <- seq(0, 45, 0.1)
bc.TPC.analysis <- extractTPC(bc.nonarctic.bri.uni, "bc", Temp.xs)
bc.nonarctic.predictions.summary <- bc.TPC.analysis[[1]]
bc.nonarctic.params.summary <- bc.TPC.analysis[[2]]
bc.nonarctic.params.fullposts <- bc.TPC.analysis[[3]]

write_csv(bc.nonarctic.predictions.summary, "data-processed/supplemental-analysis/briere-only/bc.nonarctic.predictions.summary.csv")
write_csv(bc.nonarctic.params.summary, "data-processed/supplemental-analysis/briere-only/bc.nonarctic.params.summary.csv")
write_csv(bc.nonarctic.params.fullposts, "data-processed/supplemental-analysis/briere-only/bc.nonarctic.params.fullposts.csv")

