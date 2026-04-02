## Lilian Chan, University of Guelph
## Arctic vector-borne disease transmission suitability model
##
## Purpose: use Bayesian inference (JAGS) to fit TPCs for biting rate (a). 
##
## 
## Table of content:
##    0. Set-up workspace
##    1. MCMC settings for all models
##    2. Fitting TPC (Briere)
##    3. Fitting TPC (Quadratic)
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
data.all <- read_csv("data-processed/TraitData_a.csv")
unique(data.all$species)


## Plot raw data
plot.data.a <- data.all %>%
  ggplot(aes(x = temp, y = trait)) +
  geom_point(aes(colour = species)) +
  labs(y = "Rate (days-1)", x = "Temperature ºC") +
  scale_colour_discrete(name = "Species", labels = c("Ae. albopictus",
                                                     "Ae. cinereus",
                                                     "Ae. communis",
                                                     "Ae. impiger",
                                                     "Ae. punctor"
  )) +
  facet_grid(rows = vars(type)) +
  theme_bw()

plot.data.a


## Since the Arctic dataset included fewer than three temperature treatments 
## (the minimum requirement for constructing TPC), we combined the Arctic and 
## non-Arctic species data and fitted a single TPC with uniform priors.


## Put all data into the same graph
plot.data.a.combine <- data.all %>%
  ggplot(aes(x = temp, y = trait)) +
  geom_point(aes(colour = species)) +
  labs(y = "Rate (days-1)", x = "Temperature ºC") +
  scale_colour_discrete(name = "Species", labels = c("Ae. albopictus",
                                                     "Ae. cinereus",
                                                     "Ae. communis",
                                                     "Ae. impiger",
                                                     "Ae. punctor"
  )) +
  theme_bw()

plot.data.a.combine



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


##### Set priors
prior <- data.frame(q = c(0, 1),
                    T0 = c(0, 20),
                    Tm = c(25, 45),
                    sigma_q = c(0, 0.001),
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
a.alldata.bri.uni <- jags(data = jag.data,
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
save(a.alldata.bri.uni, file = "R-scripts/R2jags-objects/all-mods/a.alldata.bri.uni.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/all-mods/a.alldata.bri.uni.Rdata")


## Diagnostics
##### Examine output
a.alldata.bri.uni$BUGSoutput$summary[1:8,]
# mcmcplot(a.alldata.bri.uni)

# Extract the DIC for future model comparisons
a.alldata.bri.uni$BUGSoutput$DIC


## Plot data + fit
df.a.alldata.bri.uni <- data.frame(a.alldata.bri.uni$BUGSoutput$summary)[-(1:8),]

## Extract the model prediction
## Overall curve
df.a.alldata.bri.uni.pop <- df.a.alldata.bri.uni %>% 
  filter(grepl("z.trait.mu.pred.pop", rownames(df.a.alldata.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5.,X50., X97.5.)


## Unique ID 1: Ae. albopictus (Delatte 2009)
df.a.alldata.bri.uni.1 <- df.a.alldata.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[1,*]"), rownames(df.a.alldata.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5.,X50., X97.5.) %>% 
  mutate(unique_id = 1)

## Unique ID 2: Ae. albopictus (Marini 2020)
df.a.alldata.bri.uni.2 <- df.a.alldata.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[2,*]"), rownames(df.a.alldata.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5.,X50., X97.5.) %>% 
  mutate(unique_id = 2)

## Unique ID 3: Ae. cinereus
df.a.alldata.bri.uni.3 <- df.a.alldata.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[3,*]"), rownames(df.a.alldata.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5.,X50., X97.5.) %>% 
  mutate(unique_id = 3)


## Unique ID 4: Ae. communis
df.a.alldata.bri.uni.4 <- df.a.alldata.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[4,*]"), rownames(df.a.alldata.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5.,X50., X97.5.) %>% 
  mutate(unique_id = 4)


## Unique ID 5: Ae. impiger
df.a.alldata.bri.uni.5 <- df.a.alldata.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[5,*]"), rownames(df.a.alldata.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5.,X50., X97.5.) %>% 
  mutate(unique_id = 5)


## Unique ID 6: Ae. punctor
df.a.alldata.bri.uni.6 <- df.a.alldata.bri.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[6,*]"), rownames(df.a.alldata.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5.,X50., X97.5.) %>% 
  mutate(unique_id = 6)


## Combine the model prediciton of all three unique groups into a dataframe
df.a.alldata.bri.uni.sp <- rbind(df.a.alldata.bri.uni.1,
                                 df.a.alldata.bri.uni.2,
                                 df.a.alldata.bri.uni.3,
                                 df.a.alldata.bri.uni.4,
                                 df.a.alldata.bri.uni.5,
                                 df.a.alldata.bri.uni.6) 

## Change unique_id into factor type
df.a.alldata.bri.uni.sp$unique_id <- as.factor(df.a.alldata.bri.uni.sp$unique_id)


##### Plot
plot.a.alldata.bri.uni <- ggplot() +
  ## data
  geom_point(data = data,
             aes(x = temp, y = trait, colour = as.factor(unique_id)),
             size = 2) +
  
  ## a separate TPC for each unique group
  geom_line(data = df.a.alldata.bri.uni.sp, 
            aes(x = temp, y = X50., color = unique_id)) +
  
  ## Overall TPC
  geom_ribbon(data = df.a.alldata.bri.uni.pop,
              aes(x = temp, ymin = X2.5., ymax = X97.5.),
              fill = "grey",
              alpha = 0.5) +
  geom_line(data = df.a.alldata.bri.uni.pop,
            aes(x = temp, y = X50.), color = "black", linewidth = 1) +

  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Rate (days-1)") +
  # Customize legend
  scale_colour_discrete(name = element_blank(),
                        labels = c("Ae. albopictus (Delatte 2009)",
                                   "Ae. albopictus (Marini 2020)",
                                   "Ae. cinereus",
                                   "Ae. communis",
                                   "Ae. impiger",
                                   "Ae. punctor")) +
  theme_bw()


plot.a.alldata.bri.uni

ggsave("figures/a.alldata.bri.uni.png", plot.a.alldata.bri.uni,
       width = 10.3, height = 5.6)




# 3. Fitting TPC (quadratic) ---------------------------------------------------


##### Temp sequence for derived quantity calculations
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


##### Set priors
prior <- data.frame(q = c(0, 1),
                    T0 = c(0, 20),
                    Tm = c(25, 45),
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
a.alldata.quad.uni <- jags(data = jag.data,
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
save(a.alldata.quad.uni, file = "R-scripts/R2jags-objects/all-mods/a.alldata.quad.uni.Rdata")

# Read the .Rdata
# load("R-scripts/R2jags-objects/all-mods/a.alldata.quad.uni.Rdata")


## Diagnostics
##### Examine output
a.alldata.quad.uni$BUGSoutput$summary[1:8,]
# mcmcplot(a.alldata.quad.uni)

# Extract the DIC for future model comparisons
a.alldata.quad.uni$BUGSoutput$DIC


## Plot data + fit
df.a.alldata.quad.uni <- data.frame(a.alldata.quad.uni$BUGSoutput$summary)[-(1:8),]

## Extract the model prediction
## Overall curve
df.a.alldata.quad.uni.pop <- df.a.alldata.quad.uni %>% 
  filter(grepl("z.trait.mu.pred.pop", rownames(df.a.alldata.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5.,X50., X97.5.)


## Unique ID 1: Ae. albopictus (Delatte 2009)
df.a.alldata.quad.uni.1 <- df.a.alldata.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[1,*]"), rownames(df.a.alldata.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5.,X50., X97.5.) %>% 
  mutate(unique_id = 1)

## Unique ID 2: Ae. albopictus (Marini 2020)
df.a.alldata.quad.uni.2 <- df.a.alldata.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[2,*]"), rownames(df.a.alldata.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5.,X50., X97.5.) %>% 
  mutate(unique_id = 2)

## Unique ID 3: Ae. cinereus
df.a.alldata.quad.uni.3 <- df.a.alldata.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[3,*]"), rownames(df.a.alldata.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5.,X50., X97.5.) %>% 
  mutate(unique_id = 3)


## Unique ID 4: Ae. communis
df.a.alldata.quad.uni.4 <- df.a.alldata.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[4,*]"), rownames(df.a.alldata.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5.,X50., X97.5.) %>% 
  mutate(unique_id = 4)


## Unique ID 5: Ae. impiger
df.a.alldata.quad.uni.5 <- df.a.alldata.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[5,*]"), rownames(df.a.alldata.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5.,X50., X97.5.) %>% 
  mutate(unique_id = 5)


## Unique ID 6: Ae. punctor
df.a.alldata.quad.uni.6 <- df.a.alldata.quad.uni %>% 
  filter(grepl(glob2rx("z.trait.mu.pred.id[6,*]"), rownames(df.a.alldata.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5.,X50., X97.5.) %>% 
  mutate(unique_id = 6)


## Combine the model prediciton of all three unique groups into a dataframe
df.a.alldata.quad.uni.sp <- rbind(df.a.alldata.quad.uni.1,
                                  df.a.alldata.quad.uni.2,
                                  df.a.alldata.quad.uni.3,
                                  df.a.alldata.quad.uni.4,
                                  df.a.alldata.quad.uni.5,
                                  df.a.alldata.quad.uni.6) 

## Change unique_id into factor type
df.a.alldata.quad.uni.sp$unique_id <- as.factor(df.a.alldata.quad.uni.sp$unique_id)


##### Plot
plot.a.alldata.quad.uni <- ggplot() +
  ## data
  geom_point(data = data,
             aes(x = temp, y = trait, colour = as.factor(unique_id)),
             size = 2) +
  
  ## a separate TPC for each unique group
  geom_line(data = df.a.alldata.quad.uni.sp, 
            aes(x = temp, y = X50., color = unique_id)) +
  
  ## Overall TPC
  geom_ribbon(data = df.a.alldata.quad.uni.pop,
              aes(x = temp, ymin = X2.5., ymax = X97.5.),
              fill = "grey",
              alpha = 0.5) +
  geom_line(data = df.a.alldata.quad.uni.pop,
            aes(x = temp, y = X50.), color = "black", linewidth = 1) +
  
  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Rate (days-1)") +
  # Customize legend
  scale_colour_discrete(name = element_blank(),
                        labels = c("Ae. albopictus (Delatte 2009)",
                                   "Ae. albopictus (Marini 2020)",
                                   "Ae. cinereus",
                                   "Ae. communis",
                                   "Ae. impiger",
                                   "Ae. punctor")) +
  theme_bw()


plot.a.alldata.quad.uni

ggsave("figures/a.alldata.quad.uni.png", plot.a.alldata.quad.uni,
       width = 10.3, height = 5.6)




# 4. Compare model fit between Quadratic and Briere models ---------------------

##### Find best fitting model #####
# Add an identifying column in each model output dataframe
df.a.alldata.bri.uni.pop <- df.a.alldata.bri.uni.pop %>% 
  mutate(type = "briere")

df.a.alldata.quad.uni.pop <- df.a.alldata.quad.uni.pop %>% 
  mutate(type = "quadratic")

# Combine the two dataframes
df.all <- bind_rows(df.a.alldata.bri.uni.pop, df.a.alldata.quad.uni.pop)

##### Plot
plot.all <- df.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = X50., color = type), linewidth = 1) +
  geom_point(data = data.all, aes(x = temp, y = trait), size = 2) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Rate (days-1)"
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

ggsave("figures/a.bri.quad.png", plot.all, width = 10.3, height = 5.6)


#### DIC
a.alldata.bri.uni$BUGSoutput$DIC 
a.alldata.quad.uni$BUGSoutput$DIC 

# Although the quadratic model produced a slightly lower DIC, this difference 
# was negligible (ΔDIC = 1.8), indicating that both models were similarly 
# supported. Because biting rate is expected to exhibit a nonlinear, 
# right-skewed thermal response, we selected the Brière function to parameterize
# the suitability model, which better reflects the biological form of rate-based
# TPCs and has been widely used in previous vector-borne disease models.


# Save best-fitting TPC in a separate folder
a.alldata.mod <- a.alldata.bri.uni

## Save the model as Rdata 
save(a.alldata.mod, file = "R-scripts/R2jags-objects/best-fitting-mods/a.alldata.mod.Rdata")



# 5. Process and save model output for plotting -------------------------------

## Analyze TPC model
# We will create 3 files: 
# a. predictions.summary: showing the mean, median, and 95% credible interval of
#      the predicted trait value at each temp from 0 to 45ºC at a 0.1ºC
# b. params.summary: showing the showing the mean, median, and 95% credible 
#      interval of TPC parameters, Topt, and Tbreadth
# c. params.fullposts: showing the TPC parameter of each MCMC iteration


Temp.xs <- seq(0, 45, 0.1)
a.TPC.analysis <- extractTPC_raneff(a.alldata.bri.uni, "a", Temp.xs)
a.alldata.predictions.summary <- a.TPC.analysis[[1]]
a.alldata.params.summary <- a.TPC.analysis[[2]]
a.alldata.params.fullposts <- a.TPC.analysis[[3]]

write_csv(a.alldata.predictions.summary, "data-processed/a/a.alldata.predictions.summary.csv")
write_csv(a.alldata.params.summary, "data-processed/a/a.alldata.params.summary.csv")
write_csv(a.alldata.params.fullposts, "data-processed/a/a.alldata.params.fullposts.csv")
