## Lilian Chan, University of Guelph
## Arctic vector-borne disease transmission suitability model
##
## Purpose: Quantifies the difference between TPCs between Arctic and non-Arctic 
## species data for each trait, then calculates the average offset across 
## multiple traits. The estimated offset is then applied to traits that lack 
## Arctic data , allowing non-Arctic derived TPC to be adjusted to better 
## approximate Arctic thermal performance
## 
## Traits with both Arctic and non-Arctic data: MR, PDR, pLA, EV, lf
## Traits without Arctic species data: EFGC, c
##
## Arctic data from biting rate (a) only has 1 temp, so we excluded this trait 
## for this exercise.
##
## Table of content:
##    0. Set-up workspace
##    1. 
##    2. 
##    3. 

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


###### Function plot TPC based on parameters
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

##########
###### 1. Load R2jags model output ----
##########

## c
load("R-scripts/R2jags-objects/c.nonarctic.quad.uni.Rdata")


## Eggs per female per gonotrophic cycle (EFGC)
load("R-scripts/R2jags-objects/EFGC.nonarctic.quad.uni.Rdata")




##############  Parasite development rate (PDR) ---------------------------------------------------------------
## Arctic:
load("R-scripts/R2jags-objects/PDR.arctic.bri.inf.Rdata")

## Non-Arctic:
load("R-scripts/R2jags-objects/PDR.nonarctic.bri.uni.Rdata")


# TPC parameters
PDR.nonarctic.params <- data.frame(PDR.nonarctic.bri.uni$BUGSoutput$summary)[(1:3),] %>% 
  dplyr::select(X50.) #median

PDR.arctic.params <- data.frame(PDR.arctic.bri.inf$BUGSoutput$summary)[(1:3),] %>% 
  dplyr::select(X50.)


Temp.xs <- seq(-5, 50, 0.1)

PDR.arctic.preds <- data.frame(temp = Temp.xs,
                               preds = briere(Temp.xs, 
                                              T0 = PDR.arctic.params[1,],
                                              Tm = PDR.arctic.params[2,],
                                              q = PDR.arctic.params[3,]),
                               type = "Arctic")

PDR.nonarctic.preds <- data.frame(temp = Temp.xs,
                                  preds = briere(Temp.xs, 
                                                 T0 = PDR.nonarctic.params[1,],
                                                 Tm = PDR.nonarctic.params[2,], 
                                                 q = PDR.nonarctic.params[3,]),
                                  type = "non-Arctic")



PDR.preds <- bind_rows(PDR.arctic.preds, PDR.nonarctic.preds)

##### Plot
plot.PDR <- ggplot(data = PDR.preds, aes(x = temp, y = preds)) +
  geom_line(aes(color = type), linewidth = 1) +
  
  scale_colour_manual(values = c("Arctic" = "#CC79A7", "non-Arctic" = "black"),
                      name = element_blank(), # No legend title
  ) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = parse(text = "Parasite~development~rate~(day^-1)")) +
  annotate("text", x = -2, y = 0.1, label = expression(paste(italic("PDR"))), size = 5) +
  theme_bw()


plot.PDR


############## Mosquito development rate (MDR) ---------------------------------------------------------------
## Arctic:
load("R-scripts/R2jags-objects/MDR.arctic.bri.inf.Rdata")

## Non-Arctic:
load("R-scripts/R2jags-objects/MDR.nonarctic.bri.uni.Rdata")

# TPC parameters
MDR.nonarctic.params <- data.frame(MDR.nonarctic.bri.uni$BUGSoutput$summary)[(1:3),] %>% 
  dplyr::select(X50.) 

MDR.arctic.params <- data.frame(MDR.arctic.bri.inf$BUGSoutput$summary)[(1:3),] %>% 
  dplyr::select(X50.)


Temp.xs <- seq(-5, 50, 0.1)

MDR.arctic.preds <- data.frame(temp = Temp.xs,
                               preds = briere(Temp.xs, 
                                              T0 = MDR.arctic.params[1,],
                                              Tm = MDR.arctic.params[2,],
                                              q = MDR.arctic.params[3,]),
                               type = "Arctic")

MDR.nonarctic.preds <- data.frame(temp = Temp.xs,
                                  preds = briere(Temp.xs, 
                                                T0 = MDR.nonarctic.params[1,],
                                                Tm = MDR.nonarctic.params[2,], 
                                                q = MDR.nonarctic.params[3,]),
                                  type = "non-Arctic")



MDR.preds <- bind_rows(MDR.arctic.preds, MDR.nonarctic.preds)

##### Plot
plot.MDR <- ggplot(data = MDR.preds, aes(x = temp, y = preds)) +
  geom_line(aes(color = type), linewidth = 1) +
  
  scale_colour_manual(values = c("Arctic" = "#D55E00", "non-Arctic" = "black"),
                      name = element_blank(), # No legend title
                      ) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = parse(text = "Mosquito~development~rate~(day^-1)")) +
  annotate("text", x = -2, y = 0.13, label = expression(paste(italic("MDR"))), size = 5) +
  theme_bw()

plot.MDR



##############  Egg viability (EV) ---------------------------------------------------------------
## Arctic:
load("R-scripts/R2jags-objects/EV.arctic.quad.inf.Rdata")

## Non-Arctic:
load("R-scripts/R2jags-objects/EV.nonarctic.quad.uni.raneff.Rdata")

# TPC parameters
EV.nonarctic.params <- data.frame(EV.nonarctic.quad.uni.raneff$BUGSoutput$summary)[(1:3),] %>% 
  dplyr::select(X50.) 

EV.arctic.params <- data.frame(EV.arctic.quad.inf$BUGSoutput$summary)[(1:3),] %>% 
  dplyr::select(X50.)


Temp.xs <- seq(-5, 50, 0.1)

EV.arctic.preds <- data.frame(temp = Temp.xs,
                               preds = quad(Temp.xs, 
                                            T0 = EV.arctic.params[1,],
                                            Tm = EV.arctic.params[2,],
                                            q = EV.arctic.params[3,]),
                               type = "Arctic")

EV.nonarctic.preds <- data.frame(temp = Temp.xs,
                                  preds = quad(Temp.xs, 
                                               T0 = EV.nonarctic.params[1,],
                                               Tm = EV.nonarctic.params[2,], 
                                               q = EV.nonarctic.params[3,]),
                                  type = "non-Arctic")



EV.preds <- bind_rows(EV.arctic.preds, EV.nonarctic.preds)

##### Plot
plot.EV <- ggplot(data = EV.preds, aes(x = temp, y = preds)) +
  geom_line(aes(color = type), linewidth = 1) +
  
  scale_colour_manual(values = c("Arctic" = "#F5C710", "non-Arctic" = "black"),
                      name = element_blank(), # No legend title
  ) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Egg viability (%)") +
  annotate("text", x = -2, y = 0.92, label = expression(paste(italic("EV"))), size = 5) +
  theme_bw()

plot.EV



##############  Larval-to-adult survival (pLA) ---------------------------------------------------------------
## Arctic:
load("R-scripts/R2jags-objects/pLA.arctic.quad.inf.Rdata")

## Non-Arctic:
load("R-scripts/R2jags-objects/pLA.nonarctic.quad.uni.raneff.Rdata")

# TPC parameters
pLA.nonarctic.params <- data.frame(pLA.nonarctic.quad.uni.raneff$BUGSoutput$summary)[(1:3),] %>% 
  dplyr::select(X50.) 

pLA.arctic.params <- data.frame(pLA.arctic.quad.inf$BUGSoutput$summary)[(1:3),] %>% 
  dplyr::select(X50.)


pLA.arctic.preds <- data.frame(temp = Temp.xs,
                              preds = quad(Temp.xs, 
                                           T0 = pLA.arctic.params[1,],
                                           Tm = pLA.arctic.params[2,],
                                           q = pLA.arctic.params[3,]),
                              type = "Arctic")

pLA.nonarctic.preds <- data.frame(temp = Temp.xs,
                                 preds = quad(Temp.xs, 
                                              T0 = pLA.nonarctic.params[1,],
                                              Tm = pLA.nonarctic.params[2,], 
                                              q = pLA.nonarctic.params[3,]),
                                 type = "non-Arctic")



pLA.preds <- bind_rows(pLA.arctic.preds, pLA.nonarctic.preds)

##### Plot
plot.pLA <- ggplot(data = pLA.preds, aes(x = temp, y = preds)) +
  geom_line(aes(color = type), linewidth = 1) +
  
  scale_colour_manual(values = c("Arctic" = "#999999", "non-Arctic" = "black"),
                      name = element_blank(), # No legend title
  ) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Larval survival (%)") +
  annotate("text", x = -2, y = 0.8, label = expression(paste(italic("pLA"))), size = 5) +
  theme_bw()

plot.pLA


##############  Adult lifespan (lf) ---------------------------------------------------------------
## Arctic:
load("R-scripts/R2jags-objects/lf.arctic.quad.inf.raneff.Rdata")

## Non-Arctic:
load("R-scripts/R2jags-objects/lf.nonarctic.quad.uni.raneff.Rdata")

# TPC parameters
lf.nonarctic.params <- data.frame(lf.nonarctic.quad.uni.raneff$BUGSoutput$summary)[(1:3),] %>% 
  dplyr::select(X50.) 

lf.arctic.params <- data.frame(lf.arctic.quad.inf.raneff$BUGSoutput$summary)[(1:3),] %>% 
  dplyr::select(X50.)


lf.arctic.preds <- data.frame(temp = Temp.xs,
                               preds = quad(Temp.xs, 
                                            T0 = lf.arctic.params[1,],
                                            Tm = lf.arctic.params[2,],
                                            q = lf.arctic.params[3,]),
                               type = "Arctic")

lf.nonarctic.preds <- data.frame(temp = Temp.xs,
                                  preds = quad(Temp.xs, 
                                               T0 = lf.nonarctic.params[1,],
                                               Tm = lf.nonarctic.params[2,], 
                                               q = lf.nonarctic.params[3,]),
                                  type = "non-Arctic")



lf.preds <- bind_rows(lf.arctic.preds, lf.nonarctic.preds)

##### Plot
plot.lf <- ggplot(data = lf.preds, aes(x = temp, y = preds)) +
  geom_line(aes(color = type), linewidth = 1) +
  
  scale_colour_manual(values = c("Arctic" = "#0072B2", "non-Arctic" = "black"),
                      name = element_blank(), # No legend title
  ) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Mosquito adult lifespan (days)") +
  annotate("text", x = -2, y = 42, label = expression(paste(italic("lf"))), size = 5) +
  theme_bw()

plot.lf



plot.traits <- ggarrange(plot.lf, plot.EV, plot.pLA, plot.PDR, plot.MDR,
                         nrow = 5, ncol = 1, align = "hv") + 
  bgcolor("white")       

plot.traits

# ggsave("figures/trait.arctic.vs.nonarctic.png", plot.traits, width = 9, height = 15)
# 