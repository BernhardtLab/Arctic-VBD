## Lilian Chan, University of Guelph
## Arctic vector-borne disease transmission suitability model
##
## Purpose: Use trait thermal response posterior distributions from JAGS to calculate suitability S(T)
## 
## Table of content:
##    0. Set-up workspace
##    1. Load R2jags model output
##    2. Calculate S(T)
##    3. Sensitivity analysis - partial derivatives


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
library(ggpubr) # For ggarrange
library(grafify)

##### Load functions
source("R-scripts/00_Functions.R")


##########
###### 1. Load R2jags model output ----
##########

## biting rate (a)
load("R-scripts/R2jags-objects/a.alldata.bri.uni.raneff.Rdata")

## c
load("R-scripts/R2jags-objects/c.nonarctic.quad.uni.Rdata")

## Adult lifespan (lf)
load("R-scripts/R2jags-objects/lf.arctic.quad.inf.raneff.Rdata")

## Parasite development rate (PDR)
load("R-scripts/R2jags-objects/PDR.arctic.bri.inf.Rdata")

## Lifetime egg production (B)
load("R-scripts/R2jags-objects/B.alldata.bri.uni.Rdata")

## Egg viability (EV)
load("R-scripts/R2jags-objects/EV.arctic.quad.inf.Rdata")

## Larval-to-adult survival (pLA)
load("R-scripts/R2jags-objects/pLA.arctic.quad.inf.Rdata")

## Mosquito development rate (MDR)
load("R-scripts/R2jags-objects/MDR.arctic.bri.inf.Rdata")


#####  Pull out the derived/predicted values:
a.preds <- a.alldata.bri.uni.raneff$BUGSoutput$sims.list$z.trait.mu.pred.pop ## Only get the population-level fit
c.preds <- c.nonarctic.quad.uni$BUGSoutput$sims.list$z.trait.mu.pred
lf.preds <- lf.arctic.quad.inf.raneff$BUGSoutput$sims.list$z.trait.mu.pred.pop ## Only get the population-level fit
PDR.preds <- PDR.arctic.bri.inf$BUGSoutput$sims.list$z.trait.mu.pred
B.preds <- B.alldata.bri.uni$BUGSoutput$sims.list$z.trait.mu.pred
EV.preds <- EV.arctic.quad.inf$BUGSoutput$sims.list$z.trait.mu.pred
pLA.preds <- pLA.arctic.quad.inf$BUGSoutput$sims.list$z.trait.mu.pred
MDR.preds <- MDR.arctic.bri.inf$BUGSoutput$sims.list$z.trait.mu.pred




##########
###### 2. Calculate S(T) ----
##########

## Columns = temp from 0 to 45ºC at a 0.1ºC interval, Rows = 15000 MCMC iterations
S.calc <- S(a.preds, c.preds, lf.preds, PDR.preds, 1, EV.preds, pLA.preds, MDR.preds)

##### Temp sequence for derived quantity calculations
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)


# Get S mean, median, upper + lower CIs
S.out <- calcPostQuants(S.calc, Temp.xs)

# Save output
# write.csv(S.out, "data-processed/S.noB.output.raw.csv")


## Calculate relative S(T)
S.out.median <- S.out %>% 
  mutate(scaled_median = S.out$median / max(S.out$median))

# write.csv(S.out.median, "data-processed/S.noB.output.median.csv")

S.out.upperCI <- S.out %>% 
  mutate(scaled_median = S.out$median / max(S.out$upperCI)) %>%
  mutate(scaled_lowerCI = S.out$lowerCI / max(S.out$upperCI)) %>%
  mutate(scaled_upperCI = S.out$upperCI / max(S.out$upperCI )) %>%
  mutate(scaled_lowerQuartile = S.out$lowerQuartile / max(S.out$upperCI)) %>%
  mutate(scaled_upperQuartile = S.out$upperQuartile / max(S.out$upperCI ))
  


## Plot S
plot.S <- ggplot(data = S.out.upperCI) +
  geom_ribbon(aes(x = temp, ymin = scaled_lowerCI, ymax = scaled_upperCI),
              fill = "grey",
              alpha = 0.5) +
  geom_ribbon(aes(x = temp, ymin = scaled_lowerQuartile, ymax = scaled_upperQuartile),
              fill = "grey",
              alpha = 0.7) +
  geom_line(aes(x = temp, y = scaled_median), colour = "black", linewidth = 1) +
  # geom_line(aes(x = temp, y = scaled_lowerQuartile), colour = "black", linetype = "dotted", linewidth = 1) +
  # geom_line(aes(x = temp, y = scaled_upperQuartile), colour = "black", linetype = "dotted", linewidth = 1) +
  # scale_x_continuous(limits = c(10, 35)) +
  scale_x_continuous(limits = c(10, 40)) +
  labs(x = expression(paste("Temperature (", degree, "C)")), y = "Suitability (S)") +
  theme_bw()

plot.S

# ggsave("figures/S.CI.noB.png", plot.S, width = 10.3, height = 5.6)

plot.S <- ggplot(data = S.out.median) +
  geom_line(aes(x = temp, y = scaled_median), colour = "black", linewidth = 1) +
  # scale_x_continuous(limits = c(10, 35)) +
  scale_x_continuous(limits = c(10, 40)) +
  labs(x = expression(paste("Temperature (", degree, "C)")), y = "Suitability (S)") +
  theme_bw()

plot.S

# ggsave("figures/S.noB.png", plot.S, width = 10.3, height = 5.6)


##########
###### 3. Calculate T0, Tm and peak S (and CI) ----
##########

calcT0TmPeak = function(input, temp.list) {
  # Create a dataframe to store the output
  output.df <- data.frame("T0" = numeric(nrow(input)), 
                       "Tmax" = numeric(nrow(input)),
                       "peak" = numeric(nrow(input)))
  
  for (i in 1:nrow(input)) { # loop through each row of the input (MCMC step)
    
    ## Create vector of list of indices where S > 0
    index.list <- which(input[i,] > 0)
    length.index.list <- length(index.list)
    
    # Store T0
    ifelse(
      # If S > 0 at the lowest temp (should be 0.0ºC, or index 1), then 0ºC will be T0
      index.list[1] == 1, output.df$T0[i] <- temp.list[1], 
      # Otherwise, store the corresponding temp of the index before the first value of index.list
      output.df$T0[i] <- temp.list[index.list[1] - 1])
    
    # Store Tm 
    ifelse(
      # If S > 0 at the highest temp (should be 45.0ºC, or index 451), then 45.0ºC will be Tmax
      temp.list[index.list[length.index.list]] == length.index.list, output.df$Tmax[i] <- length.index.list, 
      # Otherwise, store the corresponding temp of the index after the last value of index.list
      output.df$Tmax[i] <- temp.list[index.list[length.index.list] + 1])
    
    # Store Peak
    max.S <- which.max(input[i,]) # index with the max S
    output.df$peak[i] <- temp.list[max.S] # Store the corresponding temp
  }

  return(output.df)
}


S.dist <- calcT0TmPeak(S.calc, Temp.xs)


## Calculate the median and CI of T0, Tmax, and peak

# Create output df
S.viz.out <- data.frame(parameter = c("Tmin", "Tmax", "peak"),
                        mean = numeric(3),
                        med = numeric(3), 
                        lowerCI = numeric(3), 
                        upperCI = numeric(3), 
                        lowerQ = numeric(3), 
                        upperQ = numeric(3))


# Tmin
S.viz.out$mean[1] <- mean(S.dist$T0)
S.viz.out$med[1] <- median(S.dist$T0)
S.viz.out$lowerCI[1] <- quantile(S.dist$T0, 0.025)
S.viz.out$upperCI[1] <- quantile(S.dist$T0, 0.975)
S.viz.out$lowerQ[1] <- quantile(S.dist$T0, 0.25)
S.viz.out$upperQ[1] <- quantile(S.dist$T0, 0.75)

# Tmax
S.viz.out$mean[2] <- mean(S.dist$Tmax)
S.viz.out$med[2] <- median(S.dist$Tmax)
S.viz.out$lowerCI[2] <- quantile(S.dist$Tmax, 0.025)
S.viz.out$upperCI[2] <- quantile(S.dist$Tmax, 0.975)
S.viz.out$lowerQ[2] <- quantile(S.dist$Tmax, 0.25)
S.viz.out$upperQ[2] <- quantile(S.dist$Tmax, 0.75)

# peak
S.viz.out$mean[3] <- mean(S.dist$peak)
S.viz.out$med[3] <- median(S.dist$peak)
S.viz.out$lowerCI[3] <- quantile(S.dist$peak, 0.025)
S.viz.out$upperCI[3] <- quantile(S.dist$peak, 0.975)
S.viz.out$lowerQ[3] <- quantile(S.dist$peak, 0.25)
S.viz.out$upperQ[3] <- quantile(S.dist$peak, 0.75)


S.viz.out$parameter <- factor(S.viz.out$parameter,
                               levels = c("Tmax", "peak", "Tmin"))

# Save output
# write.csv(S.viz.out, "data-processed/S.noB.Vizout.csv")

# Plot

plot.S.viz <- ggplot(data = S.viz.out) +
  geom_linerange(aes(xmin = lowerQ, xmax = upperQ, y = parameter, colour = parameter), 
                 linewidth = 1) +
  geom_linerange(aes(xmin = lowerCI, xmax = upperCI, y = parameter, colour = parameter), 
                 linewidth = 0.5) + #default size = 0.5
  geom_point(aes(x = med, y = parameter, colour = parameter)) +
  # scale_x_continuous(limits = c(10, 35)) +
  scale_x_continuous(limits = c(10, 40)) +
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

plot.S.viz

# ggsave("figures/S.viz.noB.png", plot.S.viz, width = 10.3, height = 5.6)


## Put the S and S.viz together
plot.S.all <- ggarrange(plot.S, plot.S.viz, nrow = 2, align = "v", heights = c(2,1))

plot.S.all

# ggsave("figures/S.all.noB.png", plot.S.all, width = 10.3, height = 5.6)


##########
###### 3. Sensitivity Analysis - partial derivatives
##########

# Temperature levels and # MCMC steps
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)
nMCMC <- 15000


##### Calculate trait means
a.m <- colMeans(a.preds)
bc.m <- colMeans(c.preds)
lf.m <- colMeans(lf.preds)
PDR.m <- colMeans(PDR.preds)

EV.m <- colMeans(EV.preds)
pLA.m <- colMeans(pLA.preds)
MDR.m <- colMeans(MDR.preds)


# R0 <- expression((a^2 * bc * exp(-(1/(lf+ec))*(1/(PDR+ec))) * B * EV * pLA * MDR * lf^2)^0.5)
# 
# dS_dlf <- D(R0, "lf")
# dS_dlf



# Calculate sensitivity using partial derivatives
# SA <- SensitivityAnalysis_pd(a.alldata.bri.uni.raneff, c.nonarctic.quad.uni, 
#                              lf.arctic.quad.inf.raneff, PDR.arctic.bri.inf, 
#                              B.alldata.bri.uni, EV.arctic.quad.inf, 
#                              pLA.arctic.quad.inf, MDR.arctic.bri.inf,
#                              a.m, bc.m, lf.m, PDR.m, B.m, EV.m, pLA.m, MDR.m)

SA <- SensitivityAnalysis_pd_noB(a.alldata.bri.uni.raneff, c.nonarctic.quad.uni, 
                                 lf.arctic.quad.inf.raneff, PDR.arctic.bri.inf, 
                                 EV.arctic.quad.inf, pLA.arctic.quad.inf, 
                                 MDR.arctic.bri.inf,
                                 a.m, bc.m, lf.m, PDR.m, 1, EV.m, pLA.m, MDR.m)


# Get sensitivity posteriors for each parameter and summarize them
dS.da		<- calcPostQuants(as.data.frame(SA[[1]]), Temp.xs)
dS.dbc		<- calcPostQuants(as.data.frame(SA[[2]]), Temp.xs)
dS.dlf		<- calcPostQuants(as.data.frame(SA[[3]]), Temp.xs)
dS.dPDR	<- calcPostQuants(as.data.frame(SA[[4]]), Temp.xs)

dS.dEV 	<- calcPostQuants(as.data.frame(SA[[5]]), Temp.xs)
dS.dpLA 	<- calcPostQuants(as.data.frame(SA[[6]]), Temp.xs)
dS.dMDR		<- calcPostQuants(as.data.frame(SA[[7]]), Temp.xs)
dS.dT		<- calcPostQuants(as.data.frame(SA[[8]]), Temp.xs)


##### Plot results
plot.SA <- ggplot() +
  geom_line(data = dS.da, aes(x = temp, y = median, colour = "a")) +
  geom_line(data = dS.dbc, aes(x = temp, y = median, colour = "bc")) +
  geom_line(data = dS.dlf, aes(x = temp, y = median, colour = "lf")) +
  geom_line(data = dS.dPDR, aes(x = temp, y = median, colour = "PDR")) +
  
  geom_line(data = dS.dEV, aes(x = temp, y = median, , colour = "EV")) +
  geom_line(data = dS.dpLA, aes(x = temp, y = median, colour = "pLA")) +
  geom_line(data = dS.dMDR, aes(x = temp, y = median, colour = "MDR")) +
  geom_line(data = dS.dT, aes(x = temp, y = median, , colour = "S"), linewidth = 1) +
  # geom_vline(aes(xintercept = 15.4), linetype = "dashed") +
  # geom_vline(aes(xintercept = 22.1), linetype = "dashed") +
  # geom_vline(aes(xintercept = 20.0), linetype = "dashed") +
  # scale_x_continuous(limits = c(10, 35)) +
  scale_x_continuous(limits = c(10, 40)) +
  labs(x = expression(paste("Temperature (", degree, "C)")),
       y = "Relative sensitivity") +
  # scale_colour_grafify(palette = "okabe_ito",
  #                      name = element_blank(), # No legend title
  #                      breaks = c("S", "a", "bc", "lf", "PDR", "B", "EV", "pLA", "MDR"),
  #                      labels = c("S", "a", "bc", "lf", "PDR", "B", "EV", "pLA", "MDR")) +
  scale_colour_manual(values = c("#000000", "#E69F00", "#009E73", "#0072B2",
                                 "#CC79A7", #"#56B4E9",
                                 "#F5C710", "#999999", "#D55E00"),
                       name = element_blank(), # No legend title
                       breaks = c("S", "a", "bc", "lf", "PDR", "EV", "pLA", "MDR"),
                       labels = c("S", "a", "bc", "lf", "PDR", "EV", "pLA", "MDR")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

plot.SA

# ggsave("figures/SA.noB.png", plot.SA, width = 10.3, height = 5.6)


plot.everything <- ggarrange(plot.S, plot.S.viz, plot.SA, nrow = 3, align = "v", heights = c(2,1))
plot.everything

# ggsave("figures/S.and.SA.noB.png", plot.everything, width = 10.3, height = 5.6)


##### TPC summary ----

a.alldata.bri.uni.raneff$BUGSoutput$summary[1:8,]
c.nonarctic.quad.uni$BUGSoutput$summary[1:5,]
lf.arctic.quad.inf.raneff$BUGSoutput$summary[1:8,]
PDR.arctic.bri.inf$BUGSoutput$summary[1:5,]
B.alldata.bri.uni$BUGSoutput$summary[1:5,]
EV.arctic.quad.inf$BUGSoutput$summary[1:5,]
pLA.arctic.quad.inf$BUGSoutput$summary[1:5,]
MDR.arctic.bri.inf$BUGSoutput$summary[1:5,]



##########
###### 8. Packaging output for mapping
##########

head(S.out)

## Goal: create maps to compare the spatial distribution of days of thermal 
## suitability for transmission.

## We use two thermal suitability thresholds: S(T) > 0.001 and S(T) > 0.5, both 
## with a posterior probability greater than 0.975. This conservative thresholds
## ensure that transmission is almost certainly not excluded by temperature and 
## can minimize type I error (inclusion of inclusion of unsuitable areas and 
## prevent overestimation of potential risk)

## We will do that by scale the lower CI contour to itself get relative R0/S(T) for mapping
S.output.lowerCI <- data.frame(temp = S.out$temp,
                               scaled_lowerCI = S.out$lowerCI/max(S.out$lowerCI))


# Check output
plot(scaled_lowerCI ~ temp, data = S.output.lowerCI, type = "l")

# write.csv(S.output.lowerCI, "data-processed/S.noB.output.lowerCI.csv")
