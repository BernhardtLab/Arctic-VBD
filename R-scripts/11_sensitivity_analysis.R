## Lilian Chan, University of Guelph
## Arctic vector-borne disease transmission suitability model
##
## Purpose: Use trait thermal response posterior distributions from JAGS to calculate suitability S(T)
## 
## Table of content:
##    0. Set-up workspace
##    1. Load R2jags model output
##    2. Extracted predicted values and calculate trait means across temperature gradient
##    3. Sensitivity analysis - partial derivatives
##    4. Plot results


##########
###### 0. Set-up workspace ----
##########

library(tidyverse)
library(janitor)
library(ggsci)
library(ggpubr)
library(grafify) # For colour palette

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
load("R-scripts/R2jags-objects/lf.arctic.bri.inf.Rdata")

## Parasite development rate (PDR)
load("R-scripts/R2jags-objects/PDR.arctic.bri.inf.Rdata")

## Lifetime egg production (B)
load("R-scripts/R2jags-objects/B.alldata.bri.uni.Rdata")

## Egg viability
load("R-scripts/R2jags-objects/EV.arctic.quad.inf.Rdata")

## Larval-to-adult survival (pLA)
load("R-scripts/R2jags-objects/pLA.arctic.quad.inf.Rdata")

## Mosquito development rate (MDR)
load("R-scripts/R2jags-objects/MDR.arctic.bri.inf.Rdata")




# Temperature levels and # MCMC steps
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)
nMCMC <- 15000


##########
###### 2. Extracted predicted values and calculate trait means across temperature gradient
##########

#####  Pull out the derived/predicted values:
a.preds <- a.alldata.bri.uni.raneff$BUGSoutput$sims.list$z.trait.mu.pred.pop ## Only get the population-level fit
c.preds <- c.nonarctic.quad.uni$BUGSoutput$sims.list$z.trait.mu.pred
lf.preds <- lf.arctic.bri.inf$BUGSoutput$sims.list$z.trait.mu.pred
PDR.preds <- PDR.arctic.bri.inf$BUGSoutput$sims.list$z.trait.mu.pred
B.preds <- B.alldata.bri.uni$BUGSoutput$sims.list$z.trait.mu.pred
EV.preds <- EV.arctic.quad.inf$BUGSoutput$sims.list$z.trait.mu.pred
pLA.preds <- pLA.arctic.quad.inf$BUGSoutput$sims.list$z.trait.mu.pred
MDR.preds <- MDR.arctic.bri.inf$BUGSoutput$sims.list$z.trait.mu.pred


##### Calculate trait means
a.m <- colMeans(a.preds)
bc.m <- colMeans(c.preds)
lf.m <- colMeans(lf.preds)
PDR.m <- colMeans(PDR.preds)
B.m <- colMeans(B.preds)
EV.m <- colMeans(EV.preds)
pLA.m <- colMeans(pLA.preds)
MDR.m <- colMeans(MDR.preds)


# R0 <- expression((a^2 * bc * exp(-(1/(lf+ec))*(1/(PDR+ec))) * B * EV * pLA * MDR * lf^2)^0.5)
# 
# dS_dlf <- D(R0, "lf")
# dS_dlf


##########
###### 3. Sensitivity Analysis - partial derivatives
##########

# Calculate sensitivity using partial derivatives
SA <- SensitivityAnalysis_pd(a.alldata.bri.uni.raneff, c.nonarctic.quad.uni, 
                             lf.arctic.bri.inf, PDR.arctic.bri.inf, 
                             B.alldata.bri.uni, EV.arctic.quad.inf, 
                             pLA.arctic.quad.inf, MDR.arctic.bri.inf,
                             a.m, bc.m, lf.m, PDR.m, B.m, EV.m, pLA.m, MDR.m)


# Get sensitivity posteriors for each parameter and summarize them
dS.da		<- calcPostQuants(as.data.frame(SA[[1]]), Temp.xs)
dS.dbc		<- calcPostQuants(as.data.frame(SA[[2]]), Temp.xs)
dS.dlf		<- calcPostQuants(as.data.frame(SA[[3]]), Temp.xs)
dS.dPDR	<- calcPostQuants(as.data.frame(SA[[4]]), Temp.xs)
dS.dB		<- calcPostQuants(as.data.frame(SA[[5]]), Temp.xs)
dS.dEV 	<- calcPostQuants(as.data.frame(SA[[6]]), Temp.xs)
dS.dpLA 	<- calcPostQuants(as.data.frame(SA[[7]]), Temp.xs)
dS.dMDR		<- calcPostQuants(as.data.frame(SA[[8]]), Temp.xs)
dS.dT		<- calcPostQuants(as.data.frame(SA[[9]]), Temp.xs)



##########
###### 4. Plot results
##########

plot.SA <- ggplot() +
  geom_line(data = dS.da, aes(x = temp, y = median, colour = "a")) +
  geom_line(data = dS.dbc, aes(x = temp, y = median, colour = "bc")) +
  geom_line(data = dS.dlf, aes(x = temp, y = median, colour = "lf")) +
  geom_line(data = dS.dPDR, aes(x = temp, y = median, colour = "PDR")) +
  geom_line(data = dS.dB, aes(x = temp, y = median, , colour = "B")) +
  geom_line(data = dS.dEV, aes(x = temp, y = median, , colour = "EV")) +
  geom_line(data = dS.dpLA, aes(x = temp, y = median, colour = "pLA")) +
  geom_line(data = dS.dMDR, aes(x = temp, y = median, colour = "MDR")) +
  geom_line(data = dS.dT, aes(x = temp, y = median, , colour = "S"), linewidth = 1) +
  # geom_vline(aes(xintercept = 15.4), linetype = "dashed") +
  # geom_vline(aes(xintercept = 22.1), linetype = "dashed") +
  # geom_vline(aes(xintercept = 20.0), linetype = "dashed") +
  scale_x_continuous(limits = c(10, 35)) +
  labs(x = expression(paste("Temperature (", degree, "C)")),
       y = "Relative sensitivity") +
  scale_colour_grafify(palette = "okabe_ito", 
                       name = element_blank(), # No legend title
                       breaks = c("S", "a", "bc", "lf", "PDR", "B", "EV", "pLA", "MDR"),
                       labels = c("S", "a", "bc", "lf", "PDR", "B", "EV", "pLA", "MDR")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

plot.SA

# ggsave("figures/sensitivity_analysis.png", plot.SA, width = 10.3, height = 5.6)

plot.everything <- ggarrange(plot.S, plot.S.viz, plot.SA, nrow = 3, align = "v", heights = c(2,1))
plot.everything
