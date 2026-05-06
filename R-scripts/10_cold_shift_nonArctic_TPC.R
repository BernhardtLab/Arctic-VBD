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
## Traits with both Arctic and non-Arctic data: MR, PDR, pLA, EV
## Traits without Arctic species data: bc
##
## Since we used the combined data approach (combine Arctic and non-Arctic data)
## for biting rate (a), eggs per female per gonotrophic cycle (EFGC), and adult 
## lifespan (lf), so we excluded these traits for the cold-shift approach.
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


##### bc #####
# Read model
load("R-scripts/R2jags-objects/best-fitting-mods/bc.nonarctic.mod.Rdata")




# 2. Quantify the difference between Arctic and non-Arctic TPCs -------------------
# Combine the TPC parameters of all the traits into a single dataset
params.list <- bind_rows(MDR.params, PDR.params, pLA.params, EV.params)
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
T0.diff <- params.list %>% 
  filter(parameter == "Tmin") %>% 
  summarise(mean_Tmin_diff = mean(diff)) # mean Tmin offset is 6.03ºC

T0.diff <- T0.diff$mean_Tmin_diff
T0.diff




# 3. Cold-shift ----------------------------------------------------------------

## 3a. infection proportion (c) ------------------------------------------------
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


## Since infection efficiency is a proportion, it cannot be greater than 1
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
  geom_line(data = df.bc.nonarctic.quad.uni, aes(x = temp, y = X50.), 
            color = "#868686FF", linewidth = 1) +
  
  geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI), 
              fill = "#009E73", alpha = 0.5) +
  geom_line(aes(x = temperature, y = median), color = "#009E73", linewidth = 1) +
  
  # Customize the axes and labels
  scale_x_continuous(limits = c(0, 46)) + 
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "vector competence") +
  annotate("text", x = 0, y = 0.95, label = expression(paste(italic("bc"))), size = 5) +
  theme_bw()

plot.bc.arctic.nonarctic

ggsave("figures/bc.arctic.nonarctic.png", plot.bc.arctic.nonarctic, width = 10.3, height = 5.6)

