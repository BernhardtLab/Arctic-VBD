## Lilian Chan, University of Guelph
## Arctic vector-borne disease transmission suitability model
##
## Purpose: To apply a horizontal shift approach to "cold-shift" TPCs estimated
## in the absence of Arctic species data. To do so, we quantified the horizontal
## shifts in TPCs between Arctic and non-Arctic species data for traits with 
## data from both temperate and Arctic regions, then calculates the mean 
## difference across these traits. The difference is then applied to traits that
## lack Arctic data , allowing non-Arctic derived TPC to be adjusted to better 
## approximate Arctic thermal performance.

## 
## Traits with both Arctic and non-Arctic data: MR, PDR, pLA, EV, lf
## Traits without Arctic species data: bc
##
## Although we used the combined data approach (combine Arctic and non-Arctic data)
## for biting rate (a) and eggs per female per gonotrophic cycle (EFGC), we will 
## also try the cold-shift approach for non-Arctic TPCs fitted for these traits 
## and see if it will fit the Arctic data better.
##
## Table of content:
##    0. Set-up workspace
##    1. Load data and the TPCs
##    2. Quantify the difference between Arctic and non-Arctic TPCs
##    3. Cold-shift
##
##
## Inputs:
## Summary statistics of TPC parameters from Arctic and non-Arctic TPCs:
## data-processed/MDR/MDR.nonarctic.params.summary.csv
## data-processed/MDR/MDR.arctic.params.summary.csv
## data-processed/PDR/PDR.nonarctic.params.summary.csv
## data-processed/PDR/PDR.arctic.params.summary.csv
## data-processed/pLA/pLA.nonarctic.params.summary.csv
## data-processed/pLA/pLA.arctic.params.summary.csv
## data-processed/EV/EV.nonarctic.params.summary.csv
## data-processed/EV/EV.arctic.params.summary.csv
## data-processed/lf/lf.nonarctic.params.summary.csv
## data-processed/lf/lf.arctic.params.summary.csv
##
##
## Outputs: 
## data-processed/bc/bc.arctic.params.fullposts.csv -
##     Full posterior distributions for TPC parameters
##
## data-processed/bc/bc.arctic.params.summary.csv -
##     Summary statistics of TPC parameters
##
## data-processed/bc/bc.arctic.predictions.fullposts.csv -
##     Full posterior distributions for TPC predictions
##
## data-processed/bc/bc.arctic.predictions.summary.csv -
##     Posterior summary of TPC predictions across temperatures



# 0. Set-up workspace ----------------------------------------------------------

library(tidyverse)
library(readxl)
library(janitor)


##### Load functions
source("R-scripts/00_Functions.R")


##### Function plot TPC using TPC parameters #####
briere = function(T, T0, Tm, q){

  b <- c()

  for (i in 1:length(T)){
    if (T[i]>T0 && T[i]<Tm) # When trait value > 0
    {b[i] <- q * T[i] * (T[i]-T0) * (Tm-T[i])**0.5} # Briere function
    else {b[i] <- 0}
  }

  b # return output

}

quad = function(T, T0, Tm, q){

  b <- c()

  for (i in 1:length(T)){
    if (T[i]>T0 && T[i]<Tm) # When trait value > 0
    {b[i] <- -1 * q * (T[i]-T0) * (T[i] - Tm)} # Quadratic function
    else {b[i] <- 0}
  }

  b # return output

}

##### Temp sequence for derived quantity calculations
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)



# 1. Load data and the TPCs ----------------------------------------------------

##### MDR #####
# Read non-Arctic TPC parameters
MDR.nonarctic.params <- read_csv("data-processed/MDR/MDR.nonarctic.params.summary.csv")

MDR.nonarctic.params <- MDR.nonarctic.params %>% 
  filter(term %in% c("cf.T0", "cf.Tm", "cf.q")) %>% 
  select(term, median)

# Read Arctic TPC parameters
MDR.arctic.params <- read_csv("data-processed/MDR/MDR.arctic.params.summary.csv")

MDR.arctic.params <- MDR.arctic.params %>% 
  filter(term %in% c("cf.T0", "cf.Tm", "cf.q")) %>% 
  select(term, median)

# Combine non-Arctic and Arctic TPC parameters
MDR.params <- left_join(MDR.nonarctic.params, MDR.arctic.params, by = "term")
colnames(MDR.params) <- c("parameter", "nonarctic", "arctic")
MDR.params$trait <- "MDR"


##### PDR #####
# Read non-Arctic TPC parameters
PDR.nonarctic.params <- read_csv("data-processed/PDR/PDR.nonarctic.params.summary.csv")

PDR.nonarctic.params <- PDR.nonarctic.params %>% 
  filter(term %in% c("cf.T0", "cf.Tm", "cf.q")) %>% 
  select(term, median)

# Read Arctic TPC parameters
PDR.arctic.params <- read_csv("data-processed/PDR/PDR.arctic.params.summary.csv")

PDR.arctic.params <- PDR.arctic.params %>% 
  filter(term %in% c("cf.T0", "cf.Tm", "cf.q")) %>% 
  select(term, median)

# Combine non-Arctic and Arctic TPC parameters
PDR.params <- left_join(PDR.nonarctic.params, PDR.arctic.params, by = "term")
colnames(PDR.params) <- c("parameter", "nonarctic", "arctic")
PDR.params$trait <- "PDR"



##### pLA #####
# Read non-Arctic TPC parameters
pLA.nonarctic.params <- read_csv("data-processed/pLA/pLA.nonarctic.params.summary.csv")

pLA.nonarctic.params <- pLA.nonarctic.params %>% 
  filter(term %in% c("cf.T0", "cf.Tm", "cf.q")) %>% 
  select(term, median)

# Read Arctic TPC parameters
pLA.arctic.params <- read_csv("data-processed/pLA/pLA.arctic.params.summary.csv")

pLA.arctic.params <- pLA.arctic.params %>% 
  filter(term %in% c("cf.T0", "cf.Tm", "cf.q")) %>% 
  select(term, median)

# Combine non-Arctic and Arctic TPC parameters
pLA.params <- left_join(pLA.nonarctic.params, pLA.arctic.params, by = "term")
colnames(pLA.params) <- c("parameter", "nonarctic", "arctic")
pLA.params$trait <- "pLA"



##### EV #####
# Read non-Arctic TPC parameters
EV.nonarctic.params <- read_csv("data-processed/EV/EV.nonarctic.params.summary.csv")

EV.nonarctic.params <- EV.nonarctic.params %>% 
  filter(term %in% c("cf.T0", "cf.Tm", "cf.q")) %>% 
  select(term, median)

# Read Arctic TPC parameters
EV.arctic.params <- read_csv("data-processed/EV/EV.arctic.params.summary.csv")

EV.arctic.params <- EV.arctic.params %>% 
  filter(term %in% c("cf.T0", "cf.Tm", "cf.q")) %>% 
  select(term, median)

# Combine non-Arctic and Arctic TPC parameters
EV.params <- left_join(EV.nonarctic.params, EV.arctic.params, by = "term")
colnames(EV.params) <- c("parameter", "nonarctic", "arctic")
EV.params$trait <- "EV"


##### lf #####
# Read non-Arctic TPC parameters
lf.nonarctic.params <- read_csv("data-processed/lf/lf.nonarctic.params.summary.csv")

lf.nonarctic.params <- lf.nonarctic.params %>% 
  filter(term %in% c("cf.T0", "cf.Tm", "cf.q")) %>% 
  select(term, median)

# Read Arctic TPC parameters
lf.arctic.params <- read_csv("data-processed/lf/lf.arctic.params.summary.csv")

lf.arctic.params <- lf.arctic.params %>% 
  filter(term %in% c("cf.T0", "cf.Tm", "cf.q")) %>% 
  select(term, median)

# Combine non-Arctic and Arctic TPC parameters
lf.params <- left_join(lf.nonarctic.params, lf.arctic.params, by = "term")
colnames(lf.params) <- c("parameter", "nonarctic", "arctic")
lf.params$trait <- "lf"

##### bc #####
# Read model
load("R-scripts/R2jags-objects/best-fitting-mods/bc.nonarctic.mod.Rdata")

##### a #####
# Read model
load("R-scripts/R2jags-objects/best-fitting-mods/a.nonarctic.mod.Rdata")

##### EFGC #####
# Read model
load("R-scripts/R2jags-objects/best-fitting-mods/EFGC.nonarctic.mod.Rdata")


# 2. Quantify the difference between Arctic and non-Arctic TPCs -------------------
# Combine the TPC parameters of all the traits into a single dataset
params.list <- bind_rows(MDR.params, PDR.params, pLA.params, EV.params, lf.params)
params.list

# Calculate the difference
params.list <- params.list %>% 
  mutate(diff = nonarctic - arctic)

# Change the name of the parameters
params.list <- params.list %>% 
  mutate(parameter = case_when(parameter == "cf.Tm" ~ "Tmax",
                               parameter == "cf.T0"~ "Tmin",
                               parameter == "cf.q" ~ "q"))

write_csv(params.list, "data-processed/tpc_parameters_diff.csv")

## Because Arctic data were concentrated at low temperatures, we have greater 
## confidence for estimates of Tmin difference than those for Tmax. Hence, the 
## non-Arctic TPCs for traits without Arctic data were shifted toward lower 
## temperatures by the mean Tmin difference.


# Calculate the mean difference in Tmin between non-Arctic and Arctic TPCs
T0.diff.summary <- params.list %>% 
  filter(parameter == "Tmin") %>% 
  summarise(mean_Tmin_diff = mean(diff),
            sd_Tmin_diff = sd(diff))

T0.diff <- T0.diff.summary$mean_Tmin_diff
T0.diff



# 3. Cold-shift ----------------------------------------------------------------

## 3a. vector competence (bc) ------------------------------------------------
bc.iter.param <- data.frame(T0 = bc.nonarctic.mod$BUGSoutput$sims.list$cf.T0,
                            Tm = bc.nonarctic.mod$BUGSoutput$sims.list$cf.Tm,
                            q = bc.nonarctic.mod$BUGSoutput$sims.list$cf.q
                            )



# Perform the hot-old shift
bc.iter.param <- bc.iter.param %>% 
  mutate(new.T0 = T0 - T0.diff,
         new.Tm = Tm - T0.diff)

# We will create 4 files: 
# a. params.fullposts: showing the TPC parameter of each MCMC iteration
#
# b. params.summary: showing the showing the mean, median, and 95% credible 
#      interval of TPC parameters, Topt, and Tbreadth
#
# c. prediction.fullposts: showing the predicted trait values of each MCMC
#      iterations at each temperature from 0 to 45ºC at a 0.1ºC 
#      (rows = MCMC iteration, cols = temp)
#
# d. predictions.summary: showing the mean, median, and 95% credible interval of
#      the predicted trait value at each temp from 0 to 45ºC at a 0.1ºC

##### Create a dataframe showing the TPC parameters for each iteration #####
bc.params.fullposts <- bc.iter.param %>% 
  dplyr::select(new.T0, new.Tm, q) %>% 
  mutate(Tbreadth = new.Tm - new.T0)

bc.params.fullposts$iteration <- seq(1:nrow(bc.iter.param)) # Add a column indicating the number of MCMC iteration
bc.params.fullposts <- relocate(bc.params.fullposts, iteration, .before = new.T0)

bc.params.fullposts$trait <- "bc" # Add a column indicating trait name

colnames(bc.params.fullposts) <- c("iteration", "cf.T0", "cf.Tm", "cf.q", "Tbreadth", "trait")

## Save output
write_csv(bc.params.fullposts, "data-processed/bc/bc.arctic.params.fullposts.csv")



##### Calculate trait values based on new TPC parameters #####
bc.arctic.preds <- data.frame() # Initialize an empty dataframe
                             
for (i in 1:nrow(bc.iter.param)) {
  # Calculate trait values for each MCMC iteration
  iter.preds <- quad(Temp.xs, 
                     T0 = bc.iter.param$new.T0[i],
                     Tm = bc.iter.param$new.Tm[i], 
                     q = bc.iter.param$q[i]) 
  
  # Add the trait values as a column
  if (i == 1) {
    bc.arctic.preds <- iter.preds
    }
  else {
    bc.arctic.preds <- bind_cols(bc.arctic.preds, iter.preds)
    }
}

# Transpose the dataset so that each row is a MCMC iteration (for calculating TPC parameters summary in the next step)
bc.arctic.preds <- as.data.frame(t(bc.arctic.preds))
colnames(bc.arctic.preds) <- seq(1:ncol(bc.arctic.preds))


##### Summary of the TPC parameters (Tmin, Tmax, q, Topt) #####
bc.T0 <- data.frame(term = "cf.T0",
                   mean = mean(bc.params.fullposts$cf.T0),
                   sd = sd(bc.params.fullposts$cf.T0),
                           lowerCI = quantile(bc.params.fullposts$cf.T0, 0.025)[[1]],
                           lowerQ = quantile(bc.params.fullposts$cf.T0, 0.25)[[1]],
                           median =  quantile(bc.params.fullposts$cf.T0, 0.5)[[1]],
                           upperQ = quantile(bc.params.fullposts$cf.T0, 0.75)[[1]],
                           upperCI = quantile(bc.params.fullposts$cf.T0, 0.975)[[1]],
                           trait = "bc")

bc.Tm <- data.frame(term = "cf.Tm",
                   mean = mean(bc.params.fullposts$cf.Tm),
                   sd = sd(bc.params.fullposts$cf.Tm),
                   lowerCI = quantile(bc.params.fullposts$cf.Tm, 0.025)[[1]],
                   lowerQ = quantile(bc.params.fullposts$cf.Tm, 0.25)[[1]],
                   median =  quantile(bc.params.fullposts$cf.Tm, 0.5)[[1]],
                   upperQ = quantile(bc.params.fullposts$cf.Tm, 0.75)[[1]],
                   upperCI = quantile(bc.params.fullposts$cf.Tm, 0.975)[[1]],
                   trait = "bc")

bc.q <- data.frame(term = "cf.q",
                  mean = mean(bc.params.fullposts$cf.q),
                  sd = sd(bc.params.fullposts$cf.q),
                  lowerCI = quantile(bc.params.fullposts$cf.q, 0.025)[[1]],
                  lowerQ = quantile(bc.params.fullposts$cf.q, 0.25)[[1]],
                  median =  quantile(bc.params.fullposts$cf.q, 0.5)[[1]],
                  upperQ = quantile(bc.params.fullposts$cf.q, 0.75)[[1]],
                  upperCI = quantile(bc.params.fullposts$cf.q, 0.975)[[1]],
                  trait = "bc")

# Calculate Topt for each iteration and calculate summary statistics (mean, sd, & quantiles)
bc.Topt <- calcToptQuants(bc.arctic.preds, "bc", Temp.xs)

# Add Topt and Tbreadth to parameters summary data frame
bc.params.summary <- bind_rows(bc.T0, bc.Tm, bc.q, bc.Topt)


## Save output
write_csv(bc.params.summary, "data-processed/bc/bc.arctic.params.summary.csv")


## Since vector competence is a proportion, it cannot be greater than 1
## Replace values greater than 1 to 1
## Do this step after calculating TPC parameters
bc.arctic.preds <- replace(bc.arctic.preds, bc.arctic.preds > 1, 1)

## Save output
write_csv(bc.arctic.preds, "data-processed/bc/bc.arctic.predictions.fullposts.csv")


##### Create a dataframe showing the mean, median and CIs of trait values at each temp
bc.predictions.summary <- calcPostQuants(bc.arctic.preds, "bc", Temp.xs)

## Save output
write_csv(bc.predictions.summary, "data-processed/bc/bc.arctic.predictions.summary.csv")


##### Plot
# bc.predictions.summary <- read.csv("data-processed/bc.arctic.predictions.summary.csv")
# bc.params.summary <- read.csv("data-processed/bc.arctic.params.summary.csv")


plot.bc <- bc.predictions.summary %>% 
  ggplot() +
  geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI), fill = "#009E73", alpha = 0.5) +
  geom_line(aes(x = temperature, y = median), color = "#009E73", linewidth = 1) +
  # Customize the axes and labels
  scale_x_continuous(limits = c(0, 46)) + 
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "vector competence") +
  annotate("text", x = 0, y = 0.95, label = expression(paste(italic("bc"))), size = 5) +
  theme_bw()

plot.bc

ggsave("figures/bc.arctic.quad.coldshift.png", plot.bc, width = 10.3, height = 5.6)


## Compare original and adjusted TPC
df.bc.nonarctic.quad.uni <- data.frame(bc.nonarctic.mod$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5.,X50., X97.5.)

plot.bc.arctic.nonarctic <- bc.predictions.summary %>% 
  ggplot() +
  # Original
  geom_ribbon(data = df.bc.nonarctic.quad.uni, aes(x = temp, ymin = X2.5., ymax = X97.5.),
              fill = "#868686FF", alpha = 0.5) +
  geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI), 
              fill = "#009E73", alpha = 0.5) +
  
  geom_line(data = df.bc.nonarctic.quad.uni, aes(x = temp, y = X50.), 
            color = "#868686FF", linewidth = 1) +
  geom_line(aes(x = temperature, y = median), color = "#009E73", linewidth = 1) +
  
  # Customize the axes and labels
  scale_x_continuous(limits = c(0, 46)) + 
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "vector competence") +
  annotate("text", x = 2, y = 0.9, label = expression(paste(italic("bc"))), size = 5) +
  theme_bw()

plot.bc.arctic.nonarctic

ggsave("figures/bc.arctic.nonarctic.png", plot.bc.arctic.nonarctic, width = 10.3, height = 5.6)


## 3b. biting rate (a) ------------------------------------------------
a.iter.param <- data.frame(T0 = a.nonarctic.mod$BUGSoutput$sims.list$cf.T0,
                           Tm = a.nonarctic.mod$BUGSoutput$sims.list$cf.Tm,
                           q = a.nonarctic.mod$BUGSoutput$sims.list$cf.q
                           )



# Perform the hot-old shift
a.iter.param <- a.iter.param %>% 
  mutate(new.T0 = T0 - T0.diff,
         new.Tm = Tm - T0.diff,
         new.q = q * (4*Tm + 3*T0 + sqrt(16*Tm^2 - 16*Tm*T0 + 9*T0^2))/
           (4*Tm + 3*T0 + sqrt(16*Tm^2 - 16*Tm*T0 + 9*T0^2) - 10*T0.diff)
         )

# We will create 4 files: 
# a. params.fullposts: showing the TPC parameter of each MCMC iteration
#
# b. params.summary: showing the showing the mean, median, and 95% credible 
#      interval of TPC parameters, Topt, and Tbreadth
#
# c. prediction.fullposts: showing the predicted trait values of each MCMC
#      iterations at each temperature from 0 to 45ºC at a 0.1ºC 
#      (rows = MCMC iteration, cols = temp)
#
# d. predictions.summary: showing the mean, median, and 95% credible interval of
#      the predicted trait value at each temp from 0 to 45ºC at a 0.1ºC

##### Create a dataframe showing the TPC parameters for each iteration #####
a.params.fullposts <- a.iter.param %>% 
  dplyr::select(new.T0, new.Tm, new.q) %>% 
  mutate(Tbreadth = new.Tm - new.T0)

a.params.fullposts$iteration <- seq(1:nrow(a.iter.param)) # Add a column indicating the number of MCMC iteration
a.params.fullposts <- relocate(a.params.fullposts, iteration, .before = new.T0)

a.params.fullposts$trait <- "a" # Add a column indicating trait name

colnames(a.params.fullposts) <- c("iteration", "cf.T0", "cf.Tm", "cf.q", "Tbreadth", "trait")

## Save output
write_csv(a.params.fullposts, "data-processed/a/a.arctic.params.fullposts.csv")



##### Calculate trait values based on new TPC parameters #####
a.arctic.preds <- data.frame() # Initialize an empty dataframe

for (i in 1:nrow(a.iter.param)) {
  # Calculate trait values for each MCMC iteration
  iter.preds <- briere(Temp.xs, 
                     T0 = a.iter.param$new.T0[i],
                     Tm = a.iter.param$new.Tm[i], 
                     q = a.iter.param$new.q[i]) 
  
  # Add the trait values as a column
  if (i == 1) {
    a.arctic.preds <- iter.preds
  }
  else {
    a.arctic.preds <- bind_cols(a.arctic.preds, iter.preds)
  }
}

# Transpose the dataset so that each row is a MCMC iteration (for calculating TPC parameters summary in the next step)
a.arctic.preds <- as.data.frame(t(a.arctic.preds))
colnames(a.arctic.preds) <- seq(1:ncol(a.arctic.preds))


##### Summary of the TPC parameters (Tmin, Tmax, q, Topt) #####
a.T0 <- data.frame(term = "cf.T0",
                    mean = mean(a.params.fullposts$cf.T0),
                    sd = sd(a.params.fullposts$cf.T0),
                    lowerCI = quantile(a.params.fullposts$cf.T0, 0.025)[[1]],
                    lowerQ = quantile(a.params.fullposts$cf.T0, 0.25)[[1]],
                    median =  quantile(a.params.fullposts$cf.T0, 0.5)[[1]],
                    upperQ = quantile(a.params.fullposts$cf.T0, 0.75)[[1]],
                    upperCI = quantile(a.params.fullposts$cf.T0, 0.975)[[1]],
                    trait = "a")

a.Tm <- data.frame(term = "cf.Tm",
                    mean = mean(a.params.fullposts$cf.Tm),
                    sd = sd(a.params.fullposts$cf.Tm),
                    lowerCI = quantile(a.params.fullposts$cf.Tm, 0.025)[[1]],
                    lowerQ = quantile(a.params.fullposts$cf.Tm, 0.25)[[1]],
                    median =  quantile(a.params.fullposts$cf.Tm, 0.5)[[1]],
                    upperQ = quantile(a.params.fullposts$cf.Tm, 0.75)[[1]],
                    upperCI = quantile(a.params.fullposts$cf.Tm, 0.975)[[1]],
                    trait = "a")

a.q <- data.frame(term = "cf.q",
                   mean = mean(a.params.fullposts$cf.q),
                   sd = sd(a.params.fullposts$cf.q),
                   lowerCI = quantile(a.params.fullposts$cf.q, 0.025)[[1]],
                   lowerQ = quantile(a.params.fullposts$cf.q, 0.25)[[1]],
                   median =  quantile(a.params.fullposts$cf.q, 0.5)[[1]],
                   upperQ = quantile(a.params.fullposts$cf.q, 0.75)[[1]],
                   upperCI = quantile(a.params.fullposts$cf.q, 0.975)[[1]],
                   trait = "a")

# Calculate Topt for each iteration and calculate summary statistics (mean, sd, & quantiles)
a.Topt <- calcToptQuants(a.arctic.preds, "a", Temp.xs)

# Add Topt and Tbreadth to parameters summary data frame
a.params.summary <- bind_rows(a.T0, a.Tm, a.q, a.Topt)


## Save output
write_csv(a.params.summary, "data-processed/a/a.arctic.params.summary.csv")

## Save output
write_csv(a.arctic.preds, "data-processed/a/a.arctic.predictions.fullposts.csv")


##### Create a dataframe showing the mean, median and CIs of trait values at each temp
a.predictions.summary <- calcPostQuants(a.arctic.preds, "a", Temp.xs)

## Save output
write_csv(a.predictions.summary, "data-processed/a/a.arctic.predictions.summary.csv")


##### Plot
# a.predictions.summary <- read.csv("data-processed/a.arctic.predictions.summary.csv")
# a.params.summary <- read.csv("data-processed/a.arctic.params.summary.csv")

data.a.all <- read_csv("data-processed/TraitData_a.csv")
unique(data.a.all$species)

# Subset data
## Arctic species
data.a.arctic <- subset(data.a.all, type == "Arctic")

## Non-Arctic species
data.a.nonarctic <- subset(data.a.all, type == "non-Arctic")

plot.a <- a.predictions.summary %>% 
  ggplot() +
  geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI), fill = "#E69F00", alpha = 0.5) +
  geom_line(aes(x = temperature, y = median), color = "#E69F00", linewidth = 1) +
  geom_point(data = data.a.arctic, aes(x = temp, y = trait, colour = species), size = 2) +
  
  # Customize the axes and labels
  scale_x_continuous(limits = c(0, 46)) + 
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Rate (days-1)") +
  annotate("text", x = 0, y = 0.4, label = expression(paste(italic("a"))), size = 5) +
  theme_bw()

plot.a

ggsave("figures/a.arctic.bri.coldshift.png", plot.a, width = 10.3, height = 5.6)


## Compare original and adjusted TPC
df.a.nonarctic.bri.uni <- data.frame(a.nonarctic.mod$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5.,X50., X97.5.)

plot.a.arctic.nonarctic <- a.predictions.summary %>% 
  ggplot() +
  # Original
  geom_ribbon(data = df.a.nonarctic.bri.uni, aes(x = temp, ymin = X2.5., ymax = X97.5.),
              fill = "#868686FF", alpha = 0.5) +
  geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI), 
              fill = "#E69F00", alpha = 0.5) +
  
  geom_line(data = df.a.nonarctic.bri.uni, aes(x = temp, y = X50.), 
            color = "#868686FF", linewidth = 1) +
  geom_line(aes(x = temperature, y = median), color = "#E69F00", linewidth = 1) +
  geom_point(data = data.a.all, aes(x = temp, y = trait, colour = type), size = 2) +
  
  # Customize the axes and labels
  scale_x_continuous(limits = c(0, 46)) + 
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Rate (days-1)") +
  annotate("text", x = 2, y = 0.4, label = expression(paste(italic("a"))), size = 5) +
  theme_bw()

plot.a.arctic.nonarctic

ggsave("figures/a.arctic.nonarctic.png", plot.a.arctic.nonarctic, width = 10.3, height = 5.6)



## 3c. eggs per female per gonotrophic cycle (EFGC) ----------------------------
EFGC.iter.param <- data.frame(T0 = EFGC.nonarctic.mod$BUGSoutput$sims.list$cf.T0,
                              Tm = EFGC.nonarctic.mod$BUGSoutput$sims.list$cf.Tm,
                              q = EFGC.nonarctic.mod$BUGSoutput$sims.list$cf.q
                              )



# Perform the hot-old shift
EFGC.iter.param <- EFGC.iter.param %>% 
  mutate(new.T0 = T0 - T0.diff,
         new.Tm = Tm - T0.diff)

# We will create 4 files: 
# a. params.fullposts: showing the TPC parameter of each MCMC iteration
#
# b. params.summary: showing the showing the mean, median, and 95% credible 
#      interval of TPC parameters, Topt, and Tbreadth
#
# c. prediction.fullposts: showing the predicted trait values of each MCMC
#      iterations at each temperature from 0 to 45ºC at a 0.1ºC 
#      (rows = MCMC iteration, cols = temp)
#
# d. predictions.summary: showing the mean, median, and 95% credible interval of
#      the predicted trait value at each temp from 0 to 45ºC at a 0.1ºC

##### Create a dataframe showing the TPC parameters for each iteration #####
EFGC.params.fullposts <- EFGC.iter.param %>% 
  dplyr::select(new.T0, new.Tm, q) %>% 
  mutate(Tbreadth = new.Tm - new.T0)

EFGC.params.fullposts$iteration <- seq(1:nrow(EFGC.iter.param)) # Add a column indicating the number of MCMC iteration
EFGC.params.fullposts <- relocate(EFGC.params.fullposts, iteration, .before = new.T0)

EFGC.params.fullposts$trait <- "EFGC" # Add a column indicating trait name

colnames(EFGC.params.fullposts) <- c("iteration", "cf.T0", "cf.Tm", "cf.q", "Tbreadth", "trait")

## Save output
write_csv(EFGC.params.fullposts, "data-processed/EFGC/EFGC.arctic.params.fullposts.csv")



##### Calculate trait values based on new TPC parameters #####
EFGC.arctic.preds <- data.frame() # Initialize an empty dataframe

for (i in 1:nrow(EFGC.iter.param)) {
  # Calculate trait values for each MCMC iteration
  iter.preds <- quad(Temp.xs, 
                     T0 = EFGC.iter.param$new.T0[i],
                     Tm = EFGC.iter.param$new.Tm[i], 
                     q = EFGC.iter.param$q[i]) 
  
  # Add the trait values as a column
  if (i == 1) {
    EFGC.arctic.preds <- iter.preds
  }
  else {
    EFGC.arctic.preds <- bind_cols(EFGC.arctic.preds, iter.preds)
  }
}

# Transpose the dataset so that each row is a MCMC iteration (for calculating TPC parameters summary in the next step)
EFGC.arctic.preds <- as.data.frame(t(EFGC.arctic.preds))
colnames(EFGC.arctic.preds) <- seq(1:ncol(EFGC.arctic.preds))


##### Summary of the TPC parameters (Tmin, Tmax, q, Topt) #####
EFGC.T0 <- data.frame(term = "cf.T0",
                      mean = mean(EFGC.params.fullposts$cf.T0),
                      sd = sd(EFGC.params.fullposts$cf.T0),
                      lowerCI = quantile(EFGC.params.fullposts$cf.T0, 0.025)[[1]],
                      lowerQ = quantile(EFGC.params.fullposts$cf.T0, 0.25)[[1]],
                      median =  quantile(EFGC.params.fullposts$cf.T0, 0.5)[[1]],
                      upperQ = quantile(EFGC.params.fullposts$cf.T0, 0.75)[[1]],
                      upperCI = quantile(EFGC.params.fullposts$cf.T0, 0.975)[[1]],
                      trait = "EFGC")

EFGC.Tm <- data.frame(term = "cf.Tm",
                      mean = mean(EFGC.params.fullposts$cf.Tm),
                      sd = sd(EFGC.params.fullposts$cf.Tm),
                      lowerCI = quantile(EFGC.params.fullposts$cf.Tm, 0.025)[[1]],
                      lowerQ = quantile(EFGC.params.fullposts$cf.Tm, 0.25)[[1]],
                      median =  quantile(EFGC.params.fullposts$cf.Tm, 0.5)[[1]],
                      upperQ = quantile(EFGC.params.fullposts$cf.Tm, 0.75)[[1]],
                      upperCI = quantile(EFGC.params.fullposts$cf.Tm, 0.975)[[1]],
                      trait = "EFGC")

EFGC.q <- data.frame(term = "cf.q",
                     mean = mean(EFGC.params.fullposts$cf.q),
                     sd = sd(EFGC.params.fullposts$cf.q),
                     lowerCI = quantile(EFGC.params.fullposts$cf.q, 0.025)[[1]],
                     lowerQ = quantile(EFGC.params.fullposts$cf.q, 0.25)[[1]],
                     median =  quantile(EFGC.params.fullposts$cf.q, 0.5)[[1]],
                     upperQ = quantile(EFGC.params.fullposts$cf.q, 0.75)[[1]],
                     upperCI = quantile(EFGC.params.fullposts$cf.q, 0.975)[[1]],
                     trait = "EFGC")

# Calculate Topt for each iteration and calculate summary statistics (mean, sd, & quantiles)
EFGC.Topt <- calcToptQuants(EFGC.arctic.preds, "EFGC", Temp.xs)

# Add Topt and Tbreadth to parameters summary data frame
EFGC.params.summary <- bind_rows(EFGC.T0, EFGC.Tm, EFGC.q, EFGC.Topt)


## Save output
write_csv(EFGC.params.summary, "data-processed/EFGC/EFGC.arctic.params.summary.csv")


## Save output
write_csv(EFGC.arctic.preds, "data-processed/EFGC/EFGC.arctic.predictions.fullposts.csv")


##### Create a dataframe showing the mean, median and CIs of trait values at each temp
EFGC.predictions.summary <- calcPostQuants(EFGC.arctic.preds, "EFGC", Temp.xs)

## Save output
write_csv(EFGC.predictions.summary, "data-processed/EFGC/EFGC.arctic.predictions.summary.csv")


##### Plot
# EFGC.predictions.summary <- read.csv("data-processed/EFGC.arctic.predictions.summary.csv")
# EFGC.params.summary <- read.csv("data-processed/EFGC.arctic.params.summary.csv")

# Load data
data.EFGC.all <- read_csv("data-processed/TraitData_EFGC.csv")
unique(data.EFGC.all$species)

# Subset data
## Arctic species
data.EFGC.arctic <- subset(data.EFGC.all, type == "Arctic")

## Non-Arctic species
data.EFGC.nonarctic <- subset(data.EFGC.all, type == "non-Arctic")

plot.EFGC <- EFGC.predictions.summary %>% 
  ggplot() +
  geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI), fill = "#56B4E9", alpha = 0.5) +
  geom_line(aes(x = temperature, y = median), color = "#56B4E9", linewidth = 1) +
  geom_point(data = data.EFGC.arctic, aes(x = temp, y = trait, colour = species), size = 2) +
  
  # Customize the axes and labels
  scale_x_continuous(limits = c(0, 46)) + 
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Eggs") +
  annotate("text", x = 1, y = 80, label = expression(paste(italic("EFGC"))), size = 5) +
  theme_bw()

plot.EFGC

ggsave("figures/EFGC.arctic.quad.coldshift.png", plot.EFGC, width = 10.3, height = 5.6)


## Compare original and adjusted TPC
df.EFGC.nonarctic.quad.uni <- data.frame(EFGC.nonarctic.mod$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5.,X50., X97.5.)

plot.EFGC.arctic.nonarctic <- EFGC.predictions.summary %>% 
  ggplot() +
  # Original
  geom_ribbon(data = df.EFGC.nonarctic.quad.uni, aes(x = temp, ymin = X2.5., ymax = X97.5.),
              fill = "#868686FF", alpha = 0.5) +
  geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI), 
              fill = "#56B4E9", alpha = 0.5) +
  
  geom_line(data = df.EFGC.nonarctic.quad.uni, aes(x = temp, y = X50.), 
            color = "#868686FF", linewidth = 1) +
  geom_line(aes(x = temperature, y = median), color = "#56B4E9", linewidth = 1) +
  
  geom_point(data = data.EFGC.all, aes(x = temp, y = trait, colour = type), size = 2) +
  # Customize the axes and labels
  scale_x_continuous(limits = c(0, 46)) + 
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Eggs") +
  annotate("text", x = 2, y = 100, label = expression(paste(italic("EFGC"))), size = 5) +
  theme_bw()

plot.EFGC.arctic.nonarctic

ggsave("figures/EFGC.arctic.nonarctic.png", plot.EFGC.arctic.nonarctic, width = 10.3, height = 5.6)


## Compare moderate-case vs cold shift
load("R-scripts/R2jags-objects/all-mods/EFGC.alldata.quad.uni.Rdata")

df.EFGC.alldata.quad.uni <- data.frame(EFGC.alldata.quad.uni$BUGSoutput$summary)

## Extract the Overall curve
df.EFGC.alldata.quad.uni.pop <- df.EFGC.alldata.quad.uni %>% 
  filter(grepl("z.trait.mu.pred.pop", rownames(df.EFGC.alldata.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(type = "moderate")

EFGC.predictions.summary <- EFGC.predictions.summary %>% 
  mutate(type = "cold-shift")

colnames(df.EFGC.alldata.quad.uni.pop)
colnames(EFGC.predictions.summary) <- c("temp", "X2.5.", "X97.5.", "lowerQ", "upperQ", "mean", "X50.", "trait", "type")

df.EFGC.alldata.coldshift <- bind_rows(df.EFGC.alldata.quad.uni.pop, EFGC.predictions.summary)

##### Plot
plot.EFGC.alldata.coldshift <- df.EFGC.alldata.coldshift %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = X50., color = type), linewidth = 1) +
  geom_point(data = data.EFGC.arctic, aes(x = temp, y = trait), size = 2) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Eggs"
  ) +
  # Customize the colours
  ## ribbon
  scale_fill_manual(values = c("cold-shift" = "#4363d8", 
                               "moderate" = "black")) +
  ## line
  scale_color_manual(values = c("cold-shift" = "blue", 
                                "moderate" = "black")) +
  theme_bw()

plot.EFGC.alldata.coldshift

