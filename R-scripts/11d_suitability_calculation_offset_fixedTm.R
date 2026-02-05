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
#load("R-scripts/R2jags-objects/c.nonarctic.quad.uni.Rdata")
# adjusted

## Adult lifespan (lf)
load("R-scripts/R2jags-objects/lf.arctic.quad.inf.raneff.Rdata")

## Parasite development rate (PDR)
load("R-scripts/R2jags-objects/PDR.arctic.bri.inf.Rdata")

## Eggs per female per gonotrophic cycle (EFGC)
#load("R-scripts/R2jags-objects/EFGC.nonarctic.quad.uni.Rdata")

## Egg viability (EV)
load("R-scripts/R2jags-objects/EV.arctic.quad.inf.Rdata")

## Larval-to-adult survival (pLA)
load("R-scripts/R2jags-objects/pLA.arctic.quad.inf.Rdata")

## Mosquito development rate (MDR)
load("R-scripts/R2jags-objects/MDR.arctic.bri.inf.Rdata")


#####  Pull out the derived/predicted values:
a.preds <- a.alldata.bri.uni.raneff$BUGSoutput$sims.list$z.trait.mu.pred.pop ## Only get the population-level fit
c.preds <- read_csv("data-processed/c.arctic.predictions.fullposts.fixedTm.csv")
c.preds <- as.matrix(c.preds)

lf.preds <- lf.arctic.quad.inf.raneff$BUGSoutput$sims.list$z.trait.mu.pred.pop ## Only get the population-level fit
PDR.preds <- PDR.arctic.bri.inf$BUGSoutput$sims.list$z.trait.mu.pred
EFGC.preds <- read_csv("data-processed/EFGC.arctic.predictions.fullposts.fixedTm.csv")
EFGC.preds <- as.matrix(EFGC.preds)

EV.preds <- EV.arctic.quad.inf$BUGSoutput$sims.list$z.trait.mu.pred
pLA.preds <- pLA.arctic.quad.inf$BUGSoutput$sims.list$z.trait.mu.pred
MDR.preds <- MDR.arctic.bri.inf$BUGSoutput$sims.list$z.trait.mu.pred


## Pull out the full posterior distributions of TPC parameters
a.params.fullposts <- read.csv("data-processed/a.params.fullposts.csv")
c.params.fullposts <- read.csv("data-processed/c.params.fullposts.csv")
lf.params.fullposts <- read.csv("data-processed/lf.params.fullposts.csv")
PDR.params.fullposts <- read.csv("data-processed/PDR.params.fullposts.csv")
EFGC.params.fullposts <- read.csv("data-processed/EFGC.params.fullposts.csv")
EV.params.fullposts <- read.csv("data-processed/EV.params.fullposts.csv")
pLA.params.fullposts <- read.csv("data-processed/pLA.params.fullposts.csv")
MDR.params.fullposts <- read.csv("data-processed/MDR.params.fullposts.csv")


##########
###### 2. Calculate S(T) ----
##########

## Columns = temp from 0 to 45ºC at a 0.1ºC interval, Rows = 15000 MCMC iterations
S.calc <- S(a.preds, c.preds, lf.preds, PDR.preds, EFGC.preds, EV.preds, pLA.preds, MDR.preds)

##### Temp sequence for derived quantity calculations
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)

################ Debugging area: suitability all zeros across temp gradient #####################

# colnames(S.calc) <- Temp.xs
# 
# output <- S.calc %>% 
#   mutate(iteration = rownames(.)) %>% # add column with iteration number 
#   pivot_longer(!iteration, names_to = "temperature", values_to = "trait_value") %>% # convert to long format (3 columns: iteration, temp, & trait value)
#   group_by(iteration) %>% 
#   slice_max(order_by = trait_value, n = 1) %>% # for each iteration, select row with highest value for the trait
#   ungroup() %>% 
#   mutate(temperature = as.numeric(temperature)) %>% 
#   mutate(iteration = as.numeric(iteration)) %>% 
#   arrange(iteration)
# 
# n_distinct(output$iteration) == nrow(output)
# 
# output %>% filter(trait_value == 0.1) %>% 
#   arrange(iteration) %>% 
#   distinct(iteration)
# 
# 
# check <- data.frame(temp = Temp.xs, 
#                     a = a.preds[637,], 
#                     c = c.preds[637,], 
#                     lf = lf.preds[637,],
#                     PDR = PDR.preds[637,],
#                     EFGC = EFGC.preds[637,],
#                     EV = EV.preds[637,],
#                     pLA = pLA.preds[637,],
#                     MDR = MDR.preds[637,]) # lf are all zeros, PDR are very low (1e-6 max)
# 
# 
# check <- check %>% 
#   mutate(S = S(a, c, lf, PDR, EFGC, EV, pLA, MDR))
# 
# test.S <- data.frame(matrix(0, nrow = 1, ncol = 451))
# calcDerivedTPCParamPosteriors(as.data.frame(test.S), Temp.xs)

#####################################################################################################

## Because of problems with the priors in some trait TPCs (e.g. T0>Tm or q 
## extremely small), the predicted trait values are basically zeros across temp 
## gradient for some MCMC iteration in those traits, thus the suitability 
## prediction is also zero.
## Now I'll just filter those problematic MCMC iterations out (ask Joey for better solution)
which(rowSums(S.calc[])==0)

# which(rowSums(a.preds[])<1e-7)
# which(rowSums(c.preds[])<1e-7)
# which(rowSums(lf.preds[])<1e-7)
# which(rowSums(PDR.preds[])<1e-7)
# which(rowSums(EFGC.preds[])<1e-7)
# which(rowSums(EV.preds[])<1e-7)
# which(rowSums(pLA.preds[])<1e-7)
# which(rowSums(MDR.preds[])<1e-7)

## Seems like only PDR has problematic MCMC iterations
## I'll remove problematic MCMC iterations from all traits to keep the total number of iteration consistent

a.preds <- a.preds[rowSums(S.calc[])>0,]
c.preds <- c.preds[rowSums(S.calc[])>0,]
lf.preds <- lf.preds[rowSums(S.calc[])>0,]
EFGC.preds <- EFGC.preds[rowSums(S.calc[])>0,]
EV.preds <- EV.preds[rowSums(S.calc[])>0,]
pLA.preds <- pLA.preds[rowSums(S.calc[])>0,]
MDR.preds <- MDR.preds[rowSums(S.calc[])>0,]
PDR.preds <- PDR.preds[rowSums(S.calc[])>0,] 

S.calc <- S.calc[rowSums(S.calc[])>0,]


# save all 15000 MCMC iterations of suitability calculation (451 row (temp), columns = temp and 15000 MCMC iterations)
S.calc.iter <- data.frame(Temp.xs, t(S.calc))
colnames(S.calc.iter) <- c("temp", paste0("iter", seq(1:nrow(S.calc))))

# Save output
write.csv(S.calc.iter, "data-processed/S.offset.fixedTm.calc.iter.csv")

# Get the Tmin, Topt and Tmax for each MCMC iteration
S.calc.iter.summary <- calcDerivedTPCParamPosteriors(S.calc, Temp.xs)
S.calc.iter.summary$iter <- seq(1:nrow(S.calc))
write.csv(S.calc.iter.summary, "data-processed/S.offset.fixedTm.calc.iter.summary.csv")


# Get S mean, median, upper + lower CIs
S.out <- calcPostQuants(as.data.frame(S.calc), "S", Temp.xs)

# Save output
write.csv(S.out, "data-processed/S.offset.fixedTm.output.raw.csv")


## Calculate relative S(T)
S.out.median <- S.out %>% 
  mutate(scaled_median = S.out$median / max(S.out$median))

write.csv(S.out.median, "data-processed/S.offset.fixedTm.output.median.csv")

S.out.upperCI <- S.out %>% 
  mutate(scaled_median = S.out$median / max(S.out$upperCI)) %>%
  mutate(scaled_lowerCI = S.out$lowerCI / max(S.out$upperCI)) %>%
  mutate(scaled_upperCI = S.out$upperCI / max(S.out$upperCI )) %>%
  mutate(scaled_lowerQ = S.out$lowerQ / max(S.out$upperCI)) %>%
  mutate(scaled_upperQ= S.out$upperQ / max(S.out$upperCI ))
  


# Plot S
plot.S <- ggplot(data = S.out.upperCI) +
  geom_ribbon(aes(x = temperature, ymin = scaled_lowerCI, ymax = scaled_upperCI),
              fill = "grey",
              alpha = 0.5) +
  geom_ribbon(aes(x = temperature, ymin = scaled_lowerQ, ymax = scaled_upperQ),
              fill = "grey",
              alpha = 0.7) +
  geom_line(aes(x = temperature, y = scaled_median), colour = "black", linewidth = 1) +
  # geom_line(aes(x = temp, y = scaled_lowerQuartile), colour = "black", linetype = "dotted", linewidth = 1) +
  # geom_line(aes(x = temp, y = scaled_upperQuartile), colour = "black", linetype = "dotted", linewidth = 1) +
  # scale_x_continuous(limits = c(10, 35)) +
  scale_x_continuous(limits = c(10, 40)) +
  labs(title = "A",
       x = expression(paste("Temperature (", degree, "C)")), 
       y = "Suitability (S)") +
  theme_bw()

plot.S

ggsave("figures/S.offset.fixedTm.CI.png", plot.S, width = 10.3, height = 5.6)

# plot.S <- ggplot(data = S.out.median) +
#   geom_line(aes(x = temperature, y = scaled_median), colour = "black", linewidth = 1) +
#   # scale_x_continuous(limits = c(10, 35)) +
#   scale_x_continuous(limits = c(10, 40)) +
#   labs(x = expression(paste("Temperature (", degree, "C)")), y = "Suitability (S)") +
#   theme_bw()
# 
# plot.S
# 
# ggsave("figures/S.offset.fixedTm.png", plot.S, width = 10.3, height = 5.6)


##########
###### 3. Calculate T0, Tm and peak S (and CI) ----
##########

# calcT0TmPeak = function(input, temp.list) {
#   # Create a dataframe to store the output
#   output.df <- data.frame("T0" = numeric(nrow(input)), 
#                        "Tmax" = numeric(nrow(input)),
#                        "peak" = numeric(nrow(input)))
#   
#   for (i in 1:nrow(input)) { # loop through each row of the input (MCMC step)
#     
#     ## Create vector of list of indices where S > 0
#     index.list <- which(input[i,] > 0)
#     length.index.list <- length(index.list) ## should be 451
#     
#     # Store T0
#     ifelse(
#       # If S > 0 at the lowest temp (should be 0.0ºC, or index 1), then 0ºC will be T0
#       index.list[1] == 1, output.df$T0[i] <- temp.list[1], 
#       # Otherwise, store the corresponding temp of the index before the first value of index.list
#       output.df$T0[i] <- temp.list[index.list[1] - 1])
#     
#     # Store Tm 
#     ifelse(
#       # If S > 0 at the highest temp (should be 45.0ºC, or index 451), then 45.0ºC will be Tmax
#       temp.list[index.list[length.index.list]] == length.index.list, output.df$Tmax[i] <- length.index.list, 
#       # Otherwise, store the corresponding temp of the index after the last value of index.list
#       output.df$Tmax[i] <- temp.list[index.list[length.index.list] + 1])
#     
#     # Calculate Tbreadth from Tmin and Tmax
#     output.df$Tbreadth[i] <- output.df$Tmax[i] - output.df$T0[i]
#     
#     # Store Peak
#     max.S <- which.max(input[i,]) # index with the max S
#     output.df$peak[i] <- temp.list[max.S] # Store the corresponding temp
#   }
# 
#   return(output.df)
# }




S.viz.out <- extractDerivedTPC(as.data.frame(S.calc), "S", Temp.xs)

S.viz.out <- S.viz.out %>% 
  mutate(term = case_when(term == "cf.Tm" ~ "Tmax",
                          term == "cf.T0"~ "Tmin",
                          term == "Topt" ~ "Topt",
                          term == "Tbreadth" ~ "Tbreadth"))

# ## Calculate the median and CI of T0, Tmax, and peak
# 
# # Create output df
# S.viz.out <- data.frame(term = c("Tmin", "Tmax", "peak"),
#                         mean = numeric(3),
#                         med = numeric(3), 
#                         lowerCI = numeric(3), 
#                         upperCI = numeric(3), 
#                         lowerQ = numeric(3), 
#                         upperQ = numeric(3))
# 
# 
# # Tmin
# S.viz.out$mean[1] <- mean(S.dist$T0)
# S.viz.out$med[1] <- median(S.dist$T0)
# S.viz.out$lowerCI[1] <- quantile(S.dist$T0, 0.025)
# S.viz.out$upperCI[1] <- quantile(S.dist$T0, 0.975)
# S.viz.out$lowerQ[1] <- quantile(S.dist$T0, 0.25)
# S.viz.out$upperQ[1] <- quantile(S.dist$T0, 0.75)
# 
# # Tmax
# S.viz.out$mean[2] <- mean(S.dist$Tmax)
# S.viz.out$med[2] <- median(S.dist$Tmax)
# S.viz.out$lowerCI[2] <- quantile(S.dist$Tmax, 0.025)
# S.viz.out$upperCI[2] <- quantile(S.dist$Tmax, 0.975)
# S.viz.out$lowerQ[2] <- quantile(S.dist$Tmax, 0.25)
# S.viz.out$upperQ[2] <- quantile(S.dist$Tmax, 0.75)
# 
# # peak
# S.viz.out$mean[3] <- mean(S.dist$peak)
# S.viz.out$med[3] <- median(S.dist$peak)
# S.viz.out$lowerCI[3] <- quantile(S.dist$peak, 0.025)
# S.viz.out$upperCI[3] <- quantile(S.dist$peak, 0.975)
# S.viz.out$lowerQ[3] <- quantile(S.dist$peak, 0.25)
# S.viz.out$upperQ[3] <- quantile(S.dist$peak, 0.75)


S.viz.out$term <- factor(S.viz.out$term,
                         levels = c("Tmax", "Topt", "Tmin", "Tbreadth"))

# Save output
write.csv(S.viz.out, "data-processed/S.offset.fixedTm.vizout.csv")

# Plot
plot.S.viz <- S.viz.out %>% 
  filter(term != "Tbreadth") %>% 
  ggplot() +
  geom_linerange(aes(xmin = lowerQ, xmax = upperQ, y = term, colour = term), 
                 linewidth = 1) +
  geom_linerange(aes(xmin = lowerCI, xmax = upperCI, y = term, colour = term), 
                 linewidth = 0.5) + #default size = 0.5
  geom_point(aes(x = median, y = term, colour = term)) +
  # scale_x_continuous(limits = c(10, 35)) +
  scale_x_continuous(limits = c(10, 40)) +
  labs(title = "B",
       x = expression(paste("Temperature (", degree, "C)"))) +
  scale_y_discrete(labels=c("Tmin" = expression(paste("T"[min])), 
                            "Topt" = expression(paste("T"[opt])), 
                            "Tmax" = expression(paste("T"[max])))) +
  theme(axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

plot.S.viz

ggsave("figures/S.offset.fixedTm.viz.png", plot.S.viz, width = 10.3, height = 5.6)


## Put the S and S.viz together
plot.S.all <- ggarrange(plot.S, plot.S.viz, nrow = 2, align = "v", heights = c(2,1))

plot.S.all

ggsave("figures/S.offset.fixedTm.all.png", plot.S.all, width = 10.3, height = 5.6)


##########
###### 3. Sensitivity Analysis - partial derivatives
##########

# Temperature levels and # MCMC steps
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)
nMCMC <- nrow(S.calc)


##### Calculate trait means at each temperature
a.m <- colMeans(a.preds)
c.m <- colMeans(c.preds)
lf.m <- colMeans(lf.preds)
PDR.m <- colMeans(PDR.preds)
EFGC.m <- colMeans(EFGC.preds)
EV.m <- colMeans(EV.preds)
pLA.m <- colMeans(pLA.preds)
MDR.m <- colMeans(MDR.preds)


# R0 <- expression((a^2 * bc * exp(-(1/(lf+ec))*(1/(PDR+ec))) * B * EV * pLA * MDR * lf^2)^0.5)
R0 <- expression((a^3 * bc * exp(-(1/(lf+ec))*(1/(PDR+ec))) * EFGC * EV * pLA * MDR * lf^3)^0.5)

dS_da <- D(R0, "a")
dS_da


# Calculate sensitivity using partial derivatives
SA <- SensitivityAnalysis_pd_offset(a.preds, c.preds, lf.preds, PDR.preds, 
                                    EFGC.preds, EV.preds, pLA.preds, MDR.preds,
                                    a.params.fullposts, c.params.fullposts, 
                                    lf.params.fullposts, PDR.params.fullposts, 
                                    EFGC.params.fullposts, EV.params.fullposts, 
                                    pLA.params.fullposts, MDR.params.fullposts,
                                    a.m, c.m, lf.m, PDR.m, EFGC.m, EV.m, pLA.m, MDR.m)


# Get sensitivity posteriors for each term and summarize them
dS.da		<- calcPostQuants(as.data.frame(SA[[1]]), "a", Temp.xs)
dS.dc		<- calcPostQuants(as.data.frame(SA[[2]]), "c", Temp.xs)
dS.dlf		<- calcPostQuants(as.data.frame(SA[[3]]), "lf", Temp.xs)
dS.dPDR	<- calcPostQuants(as.data.frame(SA[[4]]), "PDR", Temp.xs)
dS.dEFGC		<- calcPostQuants(as.data.frame(SA[[5]]), "EFGC", Temp.xs)
dS.dEV 	<- calcPostQuants(as.data.frame(SA[[6]]),"EV",  Temp.xs)
dS.dpLA 	<- calcPostQuants(as.data.frame(SA[[7]]), "pLA", Temp.xs)
dS.dMDR		<- calcPostQuants(as.data.frame(SA[[8]]), "MDR", Temp.xs)
dS.dT		<- calcPostQuants(as.data.frame(SA[[9]]), "S", Temp.xs)



##### Plot results
plot.SA <- ggplot() +
  geom_line(data = dS.da, aes(x = temperature, y = median, colour = "a")) +
  geom_line(data = dS.dc, aes(x = temperature, y = median, colour = "c")) +
  geom_line(data = dS.dlf, aes(x = temperature, y = median, colour = "lf")) +
  geom_line(data = dS.dPDR, aes(x = temperature, y = median, colour = "PDR")) +
  geom_line(data = dS.dEFGC, aes(x = temperature, y = median, , colour = "EFGC")) +
  geom_line(data = dS.dEV, aes(x = temperature, y = median, , colour = "EV")) +
  geom_line(data = dS.dpLA, aes(x = temperature, y = median, colour = "pLA")) +
  geom_line(data = dS.dMDR, aes(x = temperature, y = median, colour = "MDR")) +
  geom_line(data = dS.dT, aes(x = temperature, y = median, , colour = "S"), linewidth = 1) +
  # geom_vline(aes(xintercept = 15.4), linetype = "dashed") +
  # geom_vline(aes(xintercept = 22.1), linetype = "dashed") +
  # geom_vline(aes(xintercept = 20.0), linetype = "dashed") +
  # scale_x_continuous(limits = c(10, 35)) +
  scale_x_continuous(limits = c(10, 40)) +
  labs(title = "C",
       x = expression(paste("Temperature (", degree, "C)")),
       y = "Relative sensitivity") +
  # scale_colour_grafify(palette = "okabe_ito",
  #                      name = element_blank(), # No legend title
  #                      breaks = c("S", "a", "bc", "lf", "PDR", "EFGC", "EV", "pLA", "MDR"),
  #                      labels = c("S", "a", "bc", "lf", "PDR", "EFGC", "EV", "pLA", "MDR")) +
  scale_colour_manual(values = c("S" = "#000000", "a" = "#E69F00", "c" = "#009E73",
                                 "lf" = "#0072B2", "PDR" = "#CC79A7", 
                                 "EFGC" = "#56B4E9", "EV" = "#F5C710", 
                                 "pLA" = "#999999", "MDR" = "#D55E00"),
                       name = element_blank(), # No legend title
                       breaks = c("S", "a", "c", "lf", "PDR", "EFGC", "EV", "pLA", "MDR"),
                       labels = c("S", "a", "c", "lf", "PDR", "EFGC",  "EV", "pLA", "MDR")) +
  # theme(panel.grid.major = element_blank(),
  #       panel.grid.minor = element_blank(),
  #       panel.background = element_blank(),
  #       panel.border = element_rect(colour = "black", fill = NA))
  theme_bw()

plot.SA

ggsave("figures/SA.offset.fixedTm.png", plot.SA, width = 10.3, height = 5.6)


plot.everything <- ggarrange(plot.S, plot.S.viz, plot.SA, nrow = 3, 
                             align = "v", heights = c(2,1,2))
plot.everything

ggsave("figures/S.and.SA.offset.fixedTm.png", plot.everything, width = 10.3, height = 5.6)


##### TPC summary ----

a.alldata.bri.uni.raneff$BUGSoutput$summary[1:8,]
c.nonarctic.quad.uni$BUGSoutput$summary[1:5,]
lf.arctic.quad.inf.raneff$BUGSoutput$summary[1:8,]
PDR.arctic.bri.inf$BUGSoutput$summary[1:5,]
EFGC.nonarctic.quad.uni$BUGSoutput$summary[1:5,]
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
S.output.lowerCI <- data.frame(temp = S.out$temperature,
                               scaled_lowerCI = S.out$lowerCI/max(S.out$lowerCI))


# Check output
plot(scaled_lowerCI ~ temp, data = S.output.lowerCI, type = "l")

write.csv(S.output.lowerCI, "data-processed/S.offset.fixedTm.output.lowerCI.csv")
