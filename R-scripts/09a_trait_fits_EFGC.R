## Lilian Chan, University of Guelph
## Arctic vector-borne disease transmission suitability model
##
## Purpose: use Bayesian inference (JAGS) to fit TPCs for Eggs per female per 
## gonotrophic cycle (EFGC). I will use data from Aedes albopictus.
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
#data.EFGC.cpip <- subset(data, species == c("pipiens"))



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
###### 2A. Fit EFGC thermal responses: Briere ----
##########

##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)

##### Set data
data <- data.EFGC.nonarctic


## Set priors
prior <- data.frame(q = c(0, 1),
                    T0 = c(0, 20),
                    Tm = c(30, 45)
)

##### inits Function
inits<-function(){list(
  cf.q = 0.01,
  cf.Tm = 35,
  cf.T0 = 5,
  cf.sigma = rlnorm(1))}


##### Parameters to Estimate
parameters <- c("cf.q", "cf.T0", "cf.Tm", "cf.sigma", "z.trait.mu.pred")


##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp


##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, 
                 N.Temp.xs = N.Temp.xs, prior = prior)

##### Run JAGS
EFGC.nonarctic.bri.uni <- jags(
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
save(EFGC.nonarctic.bri.uni, file = "R-scripts/R2jags-objects/EFGC.nonarctic.bri.uni.Rdata")

# Read the .Rdata
load("R-scripts/R2jags-objects/EFGC.nonarctic.bri.uni.Rdata")


## Diagnostics ----
##### Examine output
EFGC.nonarctic.bri.uni$BUGSoutput$summary[1:5,]
mcmcplot(EFGC.nonarctic.bri.uni)

# Extract the DIC for future model comparisons
EFGC.nonarctic.bri.uni$BUGSoutput$DIC


## Plot data + fit ----
df.EFGC.nonarctic.bri.uni <- data.frame(EFGC.nonarctic.bri.uni$BUGSoutput$summary)[-(1:5),] %>%
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)


head(df.EFGC.nonarctic.bri.uni)


##### Plot
plot.EFGC.nonarctic.bri.uni <- df.EFGC.nonarctic.bri.uni %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "#4363d8", alpha = 0.5) +
  geom_line(aes(y = mean), color = "blue", linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait, colour = citation), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Eggs per female per gonotrophic cycle") +
  # Customize legend
  scale_colour_discrete(name = element_blank(),
                        labels = c("Delatte 2009",
                                   "Ezeakacha 2015",
                                   "Yee 2016")) +
  theme_bw()


plot.EFGC.nonarctic.bri.uni

# ggsave("figures/EFGC.nonarctic.bri.uni.png", plot.EFGC.nonarctic.bri.uni,
#        width = 10.3, height = 5.6)


##########
###### 2B. Fit EFGC thermal responses: Quadratic ----
##########

##### Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)

##### Set data
data <- data.EFGC.nonarctic


## Set priors
prior <- data.frame(q = c(0, 1),
                    T0 = c(0, 20),
                    Tm = c(30, 45)
)

##### inits Function
inits<-function(){list(
  cf.q = 0.01,
  cf.Tm = 35,
  cf.T0 = 5,
  cf.sigma = rlnorm(1))}


##### Parameters to Estimate
parameters <- c("cf.q", "cf.T0", "cf.Tm", "cf.sigma", "z.trait.mu.pred")


##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$temp


##### define data for JAGS in a list object
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, 
                 N.Temp.xs = N.Temp.xs, prior = prior)

##### Run JAGS
EFGC.nonarctic.quad.uni <- jags(
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
save(EFGC.nonarctic.quad.uni, file = "R-scripts/R2jags-objects/EFGC.nonarctic.quad.uni.Rdata")

# Read the .Rdata
load("R-scripts/R2jags-objects/EFGC.nonarctic.quad.uni.Rdata")


## Diagnostics ----
##### Examine output
EFGC.nonarctic.quad.uni$BUGSoutput$summary[1:5,]
mcmcplot(EFGC.nonarctic.quad.uni)

# Extract the DIC for future model comparisons
EFGC.nonarctic.quad.uni$BUGSoutput$DIC


## Plot data + fit ----
df.EFGC.nonarctic.quad.uni <- data.frame(EFGC.nonarctic.quad.uni$BUGSoutput$summary)[-(1:5),] %>%
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)


head(df.EFGC.nonarctic.quad.uni)


##### Plot
plot.EFGC.nonarctic.quad.uni <- df.EFGC.nonarctic.quad.uni %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "#4363d8", alpha = 0.5) +
  geom_line(aes(y = mean), color = "blue", linewidth = 1) +
  geom_point(data = data, aes(x = temp, y = trait, colour = citation), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Eggs per female per gonotrophic cycle") +
  # Customize legend
  scale_colour_discrete(name = element_blank(),
                        labels = c("Delatte 2009",
                                   "Ezeakacha 2015",
                                   "Yee 2016")) +
  theme_bw()


plot.EFGC.nonarctic.quad.uni

# ggsave("figures/EFGC.nonarctic.quad.uni.png", plot.EFGC.nonarctic.quad.uni,
#        width = 10.3, height = 5.6)


##########
###### 2C. Compare the TPC fits ----
##########

## DIC
EFGC.nonarctic.bri.uni$BUGSoutput$DIC
EFGC.nonarctic.quad.uni$BUGSoutput$DIC

df.EFGC.nonarctic.bri.uni <- df.EFGC.nonarctic.bri.uni %>% 
  mutate(type = "Briere")

df.EFGC.nonarctic.quad.uni <- df.EFGC.nonarctic.quad.uni %>% 
  mutate(type = "Quadratic")


# Combine the three dataframes
df.all <- rbind(df.EFGC.nonarctic.bri.uni,
                df.EFGC.nonarctic.quad.uni
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
  # Customize legend
  scale_shape_discrete(name = element_blank(),
                        labels = c("Delatte 2009",
                                   "Ezeakacha 2015",
                                   "Yee 2016")) +
  theme_bw()

plot.all

# ggsave("figures/EFGC.nonarctic.all.png", plot.all, width = 10.3, height = 5.6)


##########
###### 3. Process and save model output for plotting ----
##########

## Analyze TPC model
EFGC.TPC.analysis <- extractTPC(EFGC.nonarctic.quad.uni, "EFGC", Temp.xs)
EFGC.predictions.summary <- EFGC.TPC.analysis[[1]]
EFGC.params.summary <- EFGC.TPC.analysis[[2]]
EFGC.params.fullposts <- EFGC.TPC.analysis[[3]]

write_csv(EFGC.predictions.summary, "data-processed/EFGC.predictions.summary.csv")
write_csv(EFGC.params.summary, "data-processed/EFGC.params.summary.csv")
write_csv(EFGC.params.fullposts, "data-processed/EFGC.params.fullposts.csv")


