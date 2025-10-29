## Lilian Chan, University of Guelph
## Arctic vector-borne disease transmission suitability model
##
## Purpose: Use trait thermal response posterior distributions from JAGS to calculate R0(T)
## 
## Table of content:
##    0. Set-up workspace
##    1. Load R2jags model output
##    2. Specify functions to calculate mean & quantiles, define R0
##    3. Calculate R0


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
library(ggpubr)

##########
###### 1. Load R2jags model output ----
##########

## Mosquito development rate (MDR)
load("R-scripts/R2jags-objects/MDR.arctic.bri.inf.Rdata")

## Parasite development rate (PDR)
load("R-scripts/R2jags-objects/PDR.arctic.bri.inf.Rdata")

## biting rate (a)
load("R-scripts/R2jags-objects/a.alldata.bri.uni.raneff.Rdata")

## Lifetime egg production (B)
load("R-scripts/R2jags-objects/B.nonarctic.bri.uni.Rdata")

## Adult lifespan (lf)
load("R-scripts/R2jags-objects/lf.arctic.bri.inf.Rdata")

## Egg viability
load("R-scripts/R2jags-objects/EV.arctic.quad.inf.Rdata")

## Larval-to-adult survival (pLA)
load("R-scripts/R2jags-objects/pLA.arctic.quad.inf.Rdata")

## c
load("R-scripts/R2jags-objects/c.nonarctic.quad.uni.Rdata")


#####  Pull out the derived/predicted values:
a.preds <- a.alldata.bri.uni.raneff$BUGSoutput$sims.list$z.trait.mu.pred.pop ## Only get the population-level fit
c.preds <- c.nonarctic.quad.uni$BUGSoutput$sims.list$z.trait.mu.pred
EV.preds <- EV.arctic.quad.inf$BUGSoutput$sims.list$z.trait.mu.pred
pLA.preds <- pLA.arctic.quad.inf$BUGSoutput$sims.list$z.trait.mu.pred
lf.preds <- lf.arctic.bri.inf$BUGSoutput$sims.list$z.trait.mu.pred
B.preds <- B.nonarctic.bri.uni$BUGSoutput$sims.list$z.trait.mu.pred
MDR.preds <- MDR.arctic.bri.inf$BUGSoutput$sims.list$z.trait.mu.pred
PDR.preds <- PDR.arctic.bri.inf$BUGSoutput$sims.list$z.trait.mu.pred


##########
###### 2. Specify functions to calculate mean & quantiles, define R0
##########

############# Specify function to calculate mean & quantiles
calcPostQuants = function(input, grad.xs) {
  
  # Get length of gradient
  N.grad.xs <- length(grad.xs)
  
  # Create output dataframe
  output.df <- data.frame("mean" = numeric(N.Temp.xs), "median" = numeric(N.Temp.xs), 
                          "lowerCI" = numeric(N.Temp.xs), "upperCI" = numeric(N.Temp.xs), 
                          "lowerQuartile" = numeric(N.Temp.xs), "upperQuartile" = numeric(N.Temp.xs), temp = grad.xs)
  
  # Calculate mean & quantiles
  for(i in 1:N.grad.xs){
    output.df$mean[i] <- mean(input[ ,i])
    output.df$median[i] <- quantile(input[ ,i], 0.5, na.rm = TRUE)
    output.df$lowerCI[i] <- quantile(input[ ,i], 0.025, na.rm = TRUE)
    output.df$upperCI[i] <- quantile(input[ ,i], 0.975, na.rm = TRUE)
    output.df$lowerQuartile[i] <- quantile(input[ ,i], 0.25, na.rm = TRUE)
    output.df$upperQuartile[i] <- quantile(input[ ,i], 0.75, na.rm = TRUE)
  }
  
  output.df # return output
  
}


############# Specify two different R0 functions and M
# **Both are written to take lifespan as argument instead of mortality rate (mu)**

# Creating a small constant to keep denominators from being zero
ec <- 0.000001

# Define R0 with bc as one value
R0 = function(a, bc, lf, PDR, B, EV, pLA, MDR){
  (a^2 * bc * exp(-(1/(lf+ec))*(1/(PDR+ec))) * B * EV * pLA * MDR * lf^2)^0.5
}


##########
###### 3. Calculate R0 ----
##########

## Columns = temp from 0 to 45ºC at a 0.1ºC interval, Rows = 15000 iterations
R0.calc <- R0(a.preds, c.preds, lf.preds, PDR.preds, B.preds, EV.preds, pLA.preds, MDR.preds)

##### Temp sequence for derived quantity calculations
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)


# Get R0s mean, median, upper + lower CIs
R0.out <- calcPostQuants(R0.calc, Temp.xs)

## Calculate relative R0
R0.out <- R0.out %>% 
  mutate(scaled_mean = R0.out$mean / max(R0.out$mean)) %>% 
  mutate(scaled_median = R0.out$median / max(R0.out$median))
  
## Plot R0
plot.R0 <- ggplot(data = R0.out) +
  geom_line(aes(x = temp, y = scaled_median), colour = "black", linewidth = 1) +
  scale_x_continuous(limits = c(10, 35)) +
  labs(x = expression(paste("Temperature (", degree, "C)")), y = expression(paste("Relative R"[0]))) +
  theme_bw()

plot.R0

# ggsave("figures/R0.png", plot.R0, width = 10.3, height = 5.6)

##########
###### 3. Calculate T0, Tm and peak R0 (and CI) ----
##########

calcT0TmPeak = function(input, temp.list) {
  # Create a dataframe to store the output
  output.df <- data.frame("T0" = numeric(nrow(input)), 
                       "Tmax" = numeric(nrow(input)),
                       "peak" = numeric(nrow(input)))
  
  for (i in 1:nrow(input)) { # loop through each row of the input (MCMC step)
    
    ## Create vector of list of indices where R0 > 0
    index.list <- which(input[i,] > 0)
    length.index.list <- length(index.list)
    
    # Store T0
    ifelse(
      # If R0 > 0 at the lowest temp (should be 0.0ºC, or index 1), then 0ºC will be T0
      index.list[1] == 1, output.df$T0[i] <- temp.list[1], 
      # Otherwise, store the corresponding temp of the index before the first value of index.list
      output.df$T0[i] <- temp.list[index.list[1] - 1])
    
    # Store Tm 
    ifelse(
      # If R0 > 0 at the highest temp (should be 45.0ºC, or index 451), then 45.0ºC will be Tmax
      temp.list[index.list[length.index.list]] == length.index.list, output.df$Tmax[i] <- length.index.list, 
      # Otherwise, store the corresponding temp of the index after the last value of index.list
      output.df$Tmax[i] <- temp.list[index.list[length.index.list] + 1])
    
    # Store Peak
    max.R0 <- which.max(input[i,]) # index with the max R0
    output.df$peak[i] <- temp.list[max.R0] # Store the corresponding temp
  }

  return(output.df)
}


R0.dist <- calcT0TmPeak(R0.calc, Temp.xs)


## Calculate the median and CI of T0, Tmax, and peak

# Create output df
R0.viz.out <- data.frame(parameter = c("Tmin", "Tmax", "peak"),
                         med = numeric(3), 
                         lowerCI = numeric(3), 
                         upperCI = numeric(3), 
                         lowerQ = numeric(3), 
                         upperQ = numeric(3))


# Tmin
R0.viz.out$med[1] <- median(R0.dist$T0)
R0.viz.out$lowerCI[1] <- quantile(R0.dist$T0, 0.025)
R0.viz.out$upperCI[1] <- quantile(R0.dist$T0, 0.975)
R0.viz.out$lowerQ[1] <- quantile(R0.dist$T0, 0.25)
R0.viz.out$upperQ[1] <- quantile(R0.dist$T0, 0.75)

# Tmax
R0.viz.out$med[2] <- median(R0.dist$Tmax)
R0.viz.out$lowerCI[2] <- quantile(R0.dist$Tmax, 0.025)
R0.viz.out$upperCI[2] <- quantile(R0.dist$Tmax, 0.975)
R0.viz.out$lowerQ[2] <- quantile(R0.dist$Tmax, 0.25)
R0.viz.out$upperQ[2] <- quantile(R0.dist$Tmax, 0.75)

# peak
R0.viz.out$med[3] <- median(R0.dist$peak)
R0.viz.out$lowerCI[3] <- quantile(R0.dist$peak, 0.025)
R0.viz.out$upperCI[3] <- quantile(R0.dist$peak, 0.975)
R0.viz.out$lowerQ[3] <- quantile(R0.dist$peak, 0.25)
R0.viz.out$upperQ[3] <- quantile(R0.dist$peak, 0.75)


R0.viz.out$parameter <- factor(R0.viz.out$parameter,
                               levels = c("Tmax", "peak", "Tmin"))

# Save output
# write.csv(R0.viz.out, "output/R0Vizout.csv")

# Plot

plot.R0.viz <- ggplot(data = R0.viz.out) +
  geom_linerange(aes(xmin = lowerQ, xmax = upperQ, y = parameter, colour = parameter), 
                 size = 1) +
  geom_linerange(aes(xmin = lowerCI, xmax = upperCI, y = parameter, colour = parameter), 
                 size = 0.5) + #default size = 0.5
  geom_point(aes(x = med, y = parameter, colour = parameter)) +
  scale_x_continuous(limits = c(10, 35)) +
  labs(x = expression(paste("Temperature (", degree, "C)"))) +
  scale_y_discrete(labels=c("Tmin" = expression(paste("T"[min])), 
                            "peak" = expression(paste("T"[opt])), 
                            "Tmax" = expression(paste("T"[max])))) +
  theme(axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

plot.R0.viz

# ggsave("figures/R0.viz.png", plot.R0.viz, width = 10.3, height = 5.6)


## Put the R0 and R0.viz together
plot.R0.all <- ggarrange(plot.R0, plot.R0.viz, nrow = 2, align = "v", heights = c(2,1))

plot.R0.all

# ggsave("figures/R0.all.png", plot.R0.all, width = 10.3, height = 5.6)
