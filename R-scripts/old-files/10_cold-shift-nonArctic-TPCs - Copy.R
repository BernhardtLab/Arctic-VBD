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
##    1. Visualize the difference
##    2. Calculate difference between Arctic and non-Arctic TPCs
##    3. Adjusting non-Arctic TPCs
##       a) Approach 1: hot-cold shift
##       b) Approach 2: Adjust Tmin and q (constant Tmax)

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
###### 1. Visualize the difference ----
##########

##############  Parasite development rate (PDR) ---------------------------------------------------------------
## Arctic:
load("R-scripts/R2jags-objects/PDR.arctic.bri.inf.Rdata")


## Non-Arctic:
load("R-scripts/R2jags-objects/PDR.nonarctic.bri.uni.raneff.Rdata")


# TPC parameters
PDR.nonarctic.params <- data.frame(PDR.nonarctic.bri.uni.raneff$BUGSoutput$summary)[(1:3),] %>% 
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
  
  scale_colour_manual(values = c("Arctic" = "red", "non-Arctic" = "black"),
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


## Plot fit + CI for Arctic and non-Arctic TPC ----
Temp.xs <- seq(0, 45, 0.5)

df.PDR.nonarctic.bri.uni.raneff <- data.frame(PDR.nonarctic.bri.uni.raneff$BUGSoutput$summary)[-(1:8),]

df.PDR.nonarctic.bri.uni.raneff.pop <- df.PDR.nonarctic.bri.uni.raneff %>% 
  filter(grepl("z.trait.mu.pred.pop", rownames(df.PDR.nonarctic.bri.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(type = "non-Arctic")

Temp.xs <- seq(0, 45, 0.1)

df.PDR.arctic.bri.inf <- data.frame(PDR.arctic.bri.inf$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(type = "Arctic")

df.PDR.all <- rbind(df.PDR.arctic.bri.inf, df.PDR.nonarctic.bri.uni.raneff.pop)

plot.PDR.all <- df.PDR.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = X50., color = type), linewidth = 1) +
  #geom_point(data = data.MDR.arctic, aes(x = temp, y = trait), size = 2) +
  #geom_point(data = data.MDR.nonarctic, aes(x = temp, y = trait), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = parse(text = "Parasite~development~rate~(day^-1)")) +
  annotate("text", x = 1, y = 0.25, label = expression(paste(italic("PDR"))), size = 5) +
  # Customize the colours
  ## ribbon
  scale_fill_manual(values = c("Arctic" = "pink",
                               "non-Arctic" = "#868686FF")) +
  ## line
  scale_color_manual(values = c("Arctic" = "red",
                                "non-Arctic" = "#868686FF")) +
  theme_bw()

plot.PDR.all

# ggsave("tpc/PDR.all.png", plot.PDR.all, width = 10.3, height = 5.6)

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
  
  scale_colour_manual(values = c("Arctic" = "red", "non-Arctic" = "black"),
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

## Plot fit + CI for Arctic and non-Arctic TPC ----
Temp.xs <- seq(0, 45, 0.5)

df.MDR.nonarctic.bri.uni <- data.frame(MDR.nonarctic.bri.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(type = "non-Arctic")

Temp.xs <- seq(0, 45, 0.1)

df.MDR.arctic.bri.inf <- data.frame(MDR.arctic.bri.inf$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)%>% 
  mutate(type = "Arctic")

df.MDR.all <- rbind(df.MDR.arctic.bri.inf, df.MDR.nonarctic.bri.uni)

plot.MDR.all <- df.MDR.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = X50., color = type), linewidth = 1) +
  #geom_point(data = data.MDR.arctic, aes(x = temp, y = trait), size = 2) +
  #geom_point(data = data.MDR.nonarctic, aes(x = temp, y = trait), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = parse(text = "Mosquito~development~rate~(day^-1)")) +
  annotate("text", x = 1, y = 0.15, label = expression(paste(italic("MDR"))), size = 5) +
  # Customize the colours
  ## ribbon
  scale_fill_manual(values = c("Arctic" = "pink",
                               "non-Arctic" = "#868686FF")) +
  ## line
  scale_color_manual(values = c("Arctic" = "red",
                                "non-Arctic" = "#868686FF")) +
  theme_bw()

plot.MDR.all

# ggsave("tpc/MDR.all.png", plot.MDR.all, width = 10.3, height = 5.6)


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
  
  scale_colour_manual(values = c("Arctic" = "blue", "non-Arctic" = "black"),
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


## Plot fit + CI for Arctic and non-Arctic TPCs ----
Temp.xs <- seq(0, 45, 0.5)

df.EV.nonarctic.quad.uni.raneff <- data.frame(EV.nonarctic.quad.uni.raneff$BUGSoutput$summary)[-(1:8),]

df.EV.nonarctic.quad.uni.raneff.pop <- df.EV.nonarctic.quad.uni.raneff %>% 
  filter(grepl("z.trait.mu.pred.pop", rownames(df.EV.nonarctic.quad.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(type = "non-Arctic")

Temp.xs <- seq(0, 45, 0.1)

df.EV.arctic.quad.inf <- data.frame(EV.arctic.quad.inf$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(type = "Arctic")

df.EV.all <- rbind(df.EV.arctic.quad.inf, df.EV.nonarctic.quad.uni.raneff.pop)

plot.EV.all <- df.EV.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = X50., color = type), linewidth = 1) +
  #geom_point(data = data.MDR.arctic, aes(x = temp, y = trait), size = 2) +
  #geom_point(data = data.MDR.nonarctic, aes(x = temp, y = trait), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Egg viability (%)") +
  annotate("text", x = 1, y = 1, label = expression(paste(italic("EV"))), size = 5) +
  # Customize the colours
  ## ribbon
  scale_fill_manual(values = c("Arctic" = "#4363d8",
                               "non-Arctic" = "#868686FF")) +
  ## line
  scale_color_manual(values = c("Arctic" = "blue",
                                "non-Arctic" = "#868686FF")) +
  theme_bw()

plot.EV.all

# ggsave("tpc/EV.all.png", plot.EV.all, width = 10.3, height = 5.6)

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
  
  scale_colour_manual(values = c("Arctic" = "blue", "non-Arctic" = "black"),
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


## Plot fit + CI for Arctic and non-Arctic TPCs ----
Temp.xs <- seq(0, 45, 0.5)

df.pLA.nonarctic.quad.uni.raneff <- data.frame(pLA.nonarctic.quad.uni.raneff$BUGSoutput$summary)[-(1:8),]

df.pLA.nonarctic.quad.uni.raneff.pop <- df.pLA.nonarctic.quad.uni.raneff %>% 
  filter(grepl("z.trait.mu.pred.pop", rownames(df.pLA.nonarctic.quad.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(type = "non-Arctic")

Temp.xs <- seq(0, 45, 0.1)

df.pLA.arctic.quad.inf <- data.frame(pLA.arctic.quad.inf$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(type = "Arctic")

df.pLA.all <- rbind(df.pLA.arctic.quad.inf, df.pLA.nonarctic.quad.uni.raneff.pop)

plot.pLA.all <- df.pLA.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = X50., color = type), linewidth = 1) +
  #geom_point(data = data.MDR.arctic, aes(x = temp, y = trait), size = 2) +
  #geom_point(data = data.MDR.nonarctic, aes(x = temp, y = trait), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Egg viability (%)") +
  annotate("text", x = 1, y = 1, label = expression(paste(italic("pLA"))), size = 5) +
  # Customize the colours
  ## ribbon
  scale_fill_manual(values = c("Arctic" = "#4363d8",
                               "non-Arctic" = "#868686FF")) +
  ## line
  scale_color_manual(values = c("Arctic" = "blue",
                                "non-Arctic" = "#868686FF")) +
  theme_bw()

plot.pLA.all

# ggsave("tpc/pLA.all.png", plot.pLA.all, width = 10.3, height = 5.6)

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
  
  scale_colour_manual(values = c("Arctic" = "blue", "non-Arctic" = "black"),
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


## Plot fit + CI for Arctic and non-Arctic TPCs ----
Temp.xs <- seq(0, 45, 0.5)

df.lf.nonarctic.quad.uni.raneff <- data.frame(lf.nonarctic.quad.uni.raneff$BUGSoutput$summary)[-(1:8),]

df.lf.nonarctic.quad.uni.raneff.pop <- df.lf.nonarctic.quad.uni.raneff %>% 
  filter(grepl("z.trait.mu.pred.pop", rownames(df.lf.nonarctic.quad.uni.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(type = "non-Arctic")

Temp.xs <- seq(0, 45, 0.1)

df.lf.arctic.quad.inf.raneff <- data.frame(lf.arctic.quad.inf.raneff$BUGSoutput$summary)[-(1:8),]

df.lf.arctic.quad.inf.raneff.pop <- df.lf.arctic.quad.inf.raneff %>% 
  filter(grepl("z.trait.mu.pred.pop", rownames(df.lf.arctic.quad.inf.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(type = "Arctic")

df.lf.all <- rbind(df.lf.arctic.quad.inf.raneff.pop, df.lf.nonarctic.quad.uni.raneff.pop)

plot.lf.all <- df.lf.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = X50., color = type), linewidth = 1) +
  #geom_point(data = data.MDR.arctic, aes(x = temp, y = trait), size = 2) +
  #geom_point(data = data.MDR.nonarctic, aes(x = temp, y = trait), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Mosquito adult lifespan (days)") +
  annotate("text", x = 1, y = 75, label = expression(paste(italic("lf"))), size = 5) +
  # Customize the colours
  ## ribbon
  scale_fill_manual(values = c("Arctic" = "#4363d8",
                               "non-Arctic" = "#868686FF")) +
  ## line
  scale_color_manual(values = c("Arctic" = "blue",
                                "non-Arctic" = "#868686FF")) +
  theme_bw()

plot.lf.all

# ggsave("tpc/lf.all.png", plot.lf.all, width = 10.3, height = 5.6)


##############  Infection proportion (c) ---------------------------------------------------------------
load("R-scripts/R2jags-objects/c.nonarctic.quad.uni.Rdata")

# TPC parameters
c.nonarctic.params <- data.frame(c.nonarctic.quad.uni$BUGSoutput$summary)[(1:3),] %>% 
  dplyr::select(X50.) #median

Temp.xs <- seq(-5, 50, 0.1)


c.nonarctic.preds <- data.frame(temp = Temp.xs,
                                preds = quad(Temp.xs, 
                                             T0 = c.nonarctic.params[1,],
                                             Tm = c.nonarctic.params[2,], 
                                             q = c.nonarctic.params[3,]),
                                type = "non-Arctic")

c.preds <- bind_rows(c.nonarctic.preds)

##### Plot
plot.c <- ggplot(data = c.preds, aes(x = temp, y = preds)) +
  geom_line(aes(color = type), linewidth = 1) +
  
  scale_colour_manual(values = c("Arctic" = "red", "non-Arctic" = "black"),
                      name = element_blank(), # No legend title
  ) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Infection proportion") +
  annotate("text", x = -2, y = 0.8, label = expression(paste(italic("c"))), size = 5) +
  theme_bw()


plot.c


## Plot fit + CI for non-Arctic TPC ----

Temp.xs <- seq(0, 45, 0.1)

df.c.nonarctic.quad.uni <- data.frame(c.nonarctic.quad.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50.,X97.5.) %>% 
  mutate(type = "non-Arctic")

df.c.all <- rbind(df.c.nonarctic.quad.uni)

plot.c.all <- df.c.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = X50., color = type), linewidth = 1) +
  #geom_point(data = data.MDR.arctic, aes(x = temp, y = trait), size = 2) +
  #geom_point(data = data.MDR.nonarctic, aes(x = temp, y = trait), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Infection proportion") +
  annotate("text", x = 1, y = 1, label = expression(paste(italic("c"))), size = 5) +
  # Customize the colours
  ## ribbon
  scale_fill_manual(values = c("Arctic" = "pink",
                               "non-Arctic" = "#868686FF")) +
  ## line
  scale_color_manual(values = c("Arctic" = "red",
                                "non-Arctic" = "#868686FF")) +
  theme_bw()

plot.c.all

# ggsave("tpc/c.all.png", plot.c.all, width = 10.3, height = 5.6)

##############  Eggs per female per gonotrophic cycle (EFGC) ---------------------------------------------------------------

load("R-scripts/R2jags-objects/EFGC.nonarctic.quad.uni.Rdata")


# TPC parameters
EFGC.nonarctic.params <- data.frame(EFGC.nonarctic.quad.uni$BUGSoutput$summary)[(1:3),] %>% 
  dplyr::select(X50.) #median

Temp.xs <- seq(-5, 50, 0.1)


EFGC.nonarctic.preds <- data.frame(temp = Temp.xs,
                                   preds = quad(Temp.xs, 
                                                T0 = EFGC.nonarctic.params[1,],
                                                Tm = EFGC.nonarctic.params[2,], 
                                                q = EFGC.nonarctic.params[3,]),
                                   type = "non-Arctic")

EFGC.preds <- bind_rows(EFGC.nonarctic.preds)

##### Plot
plot.EFGC <- ggplot(data = EFGC.preds, aes(x = temp, y = preds)) +
  geom_line(aes(color = type), linewidth = 1) +
  
  scale_colour_manual(values = c("Arctic" = "red", "non-Arctic" = "black"),
                      name = element_blank(), # No legend title
  ) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Eggs per female per gonotrophic cycle") +
  annotate("text", x = -2, y = 60, label = expression(paste(italic("EFGC"))), size = 5) +
  theme_bw()


plot.EFGC


## Plot fit + CI for non-Arctic TPC ----

Temp.xs <- seq(0, 45, 0.1)

df.EFGC.nonarctic.quad.uni <- data.frame(EFGC.nonarctic.quad.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.) %>% 
  mutate(type = "non-Arctic")

df.all <- rbind(df.EFGC.nonarctic.quad.uni)

plot.EFGC.all <- df.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = X50., color = type), linewidth = 1) +
  #geom_point(data = data.MDR.arctic, aes(x = temp, y = trait), size = 2) +
  #geom_point(data = data.MDR.nonarctic, aes(x = temp, y = trait), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(0, 41)) + 
  #scale_y_continuous(limits = c(-0.005, 0.19)) +
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Eggs per female per gonotrophic cycle") +
  annotate("text", x = 1, y = 70, label = expression(paste(italic("EFGC"))), size = 5) +
  # Customize the colours
  ## ribbon
  scale_fill_manual(values = c("Arctic" = "pink",
                               "non-Arctic" = "#868686FF")) +
  ## line
  scale_color_manual(values = c("Arctic" = "red",
                                "non-Arctic" = "#868686FF")) +
  theme_bw()

plot.EFGC.all

# ggsave("tpc/EFGC.all.png", plot.EFGC.all, width = 10.3, height = 5.6)

#--------------------------------------------------------------------------------

plot.traits <- ggarrange(plot.lf, plot.c,
                         plot.EV, plot.EFGC,
                         plot.pLA, NULL,
                         plot.PDR, NULL,
                         plot.MDR, NULL,
                         nrow = 5, ncol = 2, align = "hv") + 
  bgcolor("white")       

plot.traits

# ggsave("figures/trait.arctic.vs.nonarctic.png", plot.traits, width = 12, height = 15)
                 
plot.traits.all <- ggarrange(plot.lf.all, plot.c.all, 
                             plot.EV.all, plot.EFGC.all,
                             plot.pLA.all, NULL,
                             plot.PDR.all, NULL,
                             plot.MDR.all,
                             nrow = 5, ncol = 2, align = "hv") + 
  bgcolor("white")       

plot.traits.all

# ggsave("figures/trait.arctic.vs.nonarctic.all.png", plot.traits.all, width = 18, height = 15)
 

#--------------------------------------------------------------------------------

##########
###### 2. Calculate difference between Arctic and non-Arctic TPCs ----
##########

EV.params <- cbind(EV.nonarctic.params, EV.arctic.params)
colnames(EV.params) <- c("nonarctic", "arctic")
EV.params$trait <- "EV"
EV.params <- rownames_to_column(EV.params, var = "parameter")

lf.params <- cbind(lf.nonarctic.params, lf.arctic.params)
colnames(lf.params) <- c("nonarctic", "arctic")
lf.params$trait <- "lf"
lf.params <- rownames_to_column(lf.params, var = "parameter")

pLA.params <- cbind(pLA.nonarctic.params, pLA.arctic.params)
colnames(pLA.params) <- c("nonarctic", "arctic")
pLA.params$trait <- "pLA"
pLA.params <- rownames_to_column(pLA.params, var = "parameter")

PDR.params <- cbind(PDR.nonarctic.params, PDR.arctic.params)
colnames(PDR.params) <- c("nonarctic", "arctic")
PDR.params$trait <- "PDR"
PDR.params <- rownames_to_column(PDR.params, var = "parameter")

MDR.params <- cbind(MDR.nonarctic.params, MDR.arctic.params)
colnames(MDR.params) <- c("nonarctic", "arctic")
MDR.params$trait <- "MDR"
MDR.params <- rownames_to_column(MDR.params, var = "parameter")

params.list <- rbind(lf.params, EV.params, pLA.params, PDR.params, MDR.params)
params.list

params.list <- params.list %>% 
  mutate(diff = nonarctic - arctic)

params.list <- params.list %>% 
  mutate(parameter = case_when(parameter == "cf.Tm" ~ "Tmax",
                               parameter == "cf.T0"~ "Tmin",
                               parameter == "cf.q" ~ "q"))

# Calculate the mean difference in Tmin between non-Arctic and Arctic TPCs
T0.diff <- params.list %>% 
  filter(parameter == "Tmin") %>% 
  summarise(mean_Tmin_diff = mean(diff)) # mean Tmin offset is 3.68ºC

T0.diff <- T0.diff$mean_Tmin_diff
T0.diff

##########
###### 3. Adjusting non-Arctic TPCs ----
##########

## Now we have two approach to adjust the non-Arctic TPCs for c and EFGC:
## 1. shift the whole TPC left by the offset (hot-cold shift)
##
##    As quadratic model fits better for both traits, this can be done by 
##    shifting both Tmin and Tmax left by the offset
##
## 2. Keep Tmax the same, shift Tmin by the offset and adjust q such that the 
##    maximum trait value remains the same



##########
###### 3a. Approach 1: hot-cold shift ----
##########

#### 3ai. infection proportion (c) ----
c.iter.param <- data.frame(T0 = c.nonarctic.quad.uni$BUGSoutput$sims.list$cf.T0,
                           Tm = c.nonarctic.quad.uni$BUGSoutput$sims.list$cf.Tm,
                           q = c.nonarctic.quad.uni$BUGSoutput$sims.list$cf.q
)


# Test if plugging params into the quadratic model will get the same pred 
# c.iter.param[1,]
#
# c.iter.pred <- data.frame(pred = c.nonarctic.quad.uni$BUGSoutput$sims.list$z.trait.mu.pred)
#
# test.pred <- data.frame(temp = Temp.xs,
#                         test = quad(T = Temp.xs, T0 = 9.310141, Tm = 36.78479, q = 0.004806817),
#                         actual = t(c.iter.pred[1,]))
# 
# colnames(test.pred) <- c("temp", "test", "actual")
# 
# test.pred$equal <- ifelse(test.pred$test == test.pred2, T, F)
# Conclusion: yes

##### Perform the hot-old shift
c.iter.param <- c.iter.param %>% 
  mutate(new.T0 = T0 - T0.diff,
         new.Tm = Tm - T0.diff)

## Create a dataframe showing the TPC parameters for each iteration
c.params.fullposts <- c.iter.param %>% 
  dplyr::select(new.T0, new.Tm, q)

c.params.fullposts$iteration <- seq(1:nrow(c.iter.param)) # Add a column indicating the number of MCMC iteration
c.params.fullposts <- relocate(c.params.fullposts, iteration, .before = new.T0)

c.params.fullposts$trait <- "c" # Add a column indicating trait name

colnames(c.params.fullposts) <- c("iteration", "cf.T0", "cf.Tm", "cf.q", "trait")

## Save output
# write_csv(c.params.fullposts, "data-processed/c.arctic.params.fullposts.csv")



##### Calculate trait values based on new TPC parameters
c.arctic.preds <- data.frame() # Initialize an empty dataframe
                             
for (i in 1:nrow(c.iter.param)) {
  # Calculate trait values for each MCMC iteration
  iter.preds <- quad(Temp.xs, 
                    T0 = c.iter.param$new.T0[i],
                    Tm = c.iter.param$new.Tm[i], 
                    q = c.iter.param$q[i]) 
  
  # Add the trait values as a column
  if (i == 1) {
    c.arctic.preds <- iter.preds
    }
  else {
    c.arctic.preds <- bind_cols(c.arctic.preds, iter.preds)
    }
}

# Transpose the dataset so that each row is a MCMC iteration (for calculating TPC parameters summary in the next step)
c.arctic.preds <- as.data.frame(t(c.arctic.preds))
colnames(c.arctic.preds) <- seq(1:ncol(c.arctic.preds))


##### Calculate the mean, median, and CIs of the TPC parameters (Tmin, Tmax, q, Topt)
c.T0 <- data.frame(term = "cf.T0",
                   mean = mean(c.params.fullposts$cf.T0),
                   sd = sd(c.params.fullposts$cf.T0),
                           lowerCI = quantile(c.params.fullposts$cf.T0, 0.025)[[1]],
                           lowerQ = quantile(c.params.fullposts$cf.T0, 0.25)[[1]],
                           median =  quantile(c.params.fullposts$cf.T0, 0.5)[[1]],
                           upperQ = quantile(c.params.fullposts$cf.T0, 0.75)[[1]],
                           upperCI = quantile(c.params.fullposts$cf.T0, 0.975)[[1]],
                           trait = "c")

c.Tm <- data.frame(term = "cf.Tm",
                   mean = mean(c.params.fullposts$cf.Tm),
                   sd = sd(c.params.fullposts$cf.Tm),
                   lowerCI = quantile(c.params.fullposts$cf.Tm, 0.025)[[1]],
                   lowerQ = quantile(c.params.fullposts$cf.Tm, 0.25)[[1]],
                   median =  quantile(c.params.fullposts$cf.Tm, 0.5)[[1]],
                   upperQ = quantile(c.params.fullposts$cf.Tm, 0.75)[[1]],
                   upperCI = quantile(c.params.fullposts$cf.Tm, 0.975)[[1]],
                   trait = "c")

c.q <- data.frame(term = "cf.q",
                  mean = mean(c.params.fullposts$cf.q),
                  sd = sd(c.params.fullposts$cf.q),
                  lowerCI = quantile(c.params.fullposts$cf.q, 0.025)[[1]],
                  lowerQ = quantile(c.params.fullposts$cf.q, 0.25)[[1]],
                  median =  quantile(c.params.fullposts$cf.q, 0.5)[[1]],
                  upperQ = quantile(c.params.fullposts$cf.q, 0.75)[[1]],
                  upperCI = quantile(c.params.fullposts$cf.q, 0.975)[[1]],
                  trait = "c")

# Calculate Topt for each iteration and calculate summary statistics (mean, sd, & quantiles)
c.Topt <- calcToptQuants(c.arctic.preds, "c", Temp.xs)

# Add Topt and Tbreadth to parameters summary data frame
c.params.summary <- bind_rows(c.T0, c.Tm, c.q, c.Topt)


## Save output
# write_csv(c.params.summary, "data-processed/c.arctic.params.summary.csv")


## Since infection efficiency is a proportion, it cannot be greater than 1
## Replace values greater than 1 to 1
## Do this step after calculating TPC parameters
c.arctic.preds <- replace(c.arctic.preds, c.arctic.preds > 1, 1)

## Save output
# write_csv(c.arctic.preds, "data-processed/c.arctic.predictions.fullposts.csv")


##### Create a dataframe showing the mean, median and CIs of trait values at each temp
c.predictions.summary <- calcPostQuants(c.arctic.preds, "c", Temp.xs)

## Save output
# write_csv(c.predictions.summary, "data-processed/c.arctic.predictions.summary.csv")


##### Plot
c.predictions.summary <- read.csv("data-processed/c.arctic.predictions.summary.csv")
c.params.summary <- read.csv("data-processed/c.arctic.params.summary.csv")


plot.c <- c.predictions.summary %>% 
  ggplot() +
  geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI), fill = "#009E73", alpha = 0.5) +
  geom_line(aes(x = temperature, y = median), color = "#009E73", linewidth = 1) +
  # Customize the axes and labels
  scale_x_continuous(limits = c(0, 46)) + 
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Infection proportion") +
  annotate("text", x = 0, y = 0.95, label = expression(paste(italic("c"))), size = 5) +
  theme_bw()

plot.c

# ggsave("figures/c.arctic.quad.adjust.png", plot.c, width = 10.3, height = 5.6)

## Compare original and adjusted TPC
plot.c.all <- c.predictions.summary %>% 
  ggplot() +
  # Original
  geom_ribbon(data = df.c.nonarctic.quad.uni, aes(x = temp, ymin = X2.5., ymax = X97.5.), fill = "#868686FF", alpha = 0.5) +
  geom_line(data = df.c.nonarctic.quad.uni, aes(x = temp, y = X50.), color = "#868686FF", linewidth = 1) +
  
  geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI), fill = "#009E73", alpha = 0.5) +
  geom_line(aes(x = temperature, y = median), color = "#009E73", linewidth = 1) +
  # Customize the axes and labels
  scale_x_continuous(limits = c(0, 46)) + 
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Infection proportion") +
  annotate("text", x = 0, y = 0.95, label = expression(paste(italic("c"))), size = 5) +
  theme_bw()

plot.c.all

# ggsave("figures/c.arctic.vs.nonarctic.png", plot.c.all, width = 10.3, height = 5.6)





#### 3aii. Eggs per female per gonotrophic cycle (EFGC) ----
EFGC.iter.param <- data.frame(T0 = EFGC.nonarctic.quad.uni$BUGSoutput$sims.list$cf.T0,
                           Tm = EFGC.nonarctic.quad.uni$BUGSoutput$sims.list$cf.Tm,
                           q = EFGC.nonarctic.quad.uni$BUGSoutput$sims.list$cf.q
)


##### Perform the hot-old shift
EFGC.iter.param <- EFGC.iter.param %>% 
  mutate(new.T0 = T0 - T0.diff,
         new.Tm = Tm - T0.diff)

## Create a dataframe showing the TPC parameters for each iteration
EFGC.params.fullposts <- EFGC.iter.param %>% 
  dplyr::select(new.T0, new.Tm, q)

EFGC.params.fullposts$iteration <- seq(1:nrow(EFGC.iter.param)) # Add a column indicating the number of MCMC iteration
EFGC.params.fullposts <- relocate(EFGC.params.fullposts, iteration, .before = new.T0)

EFGC.params.fullposts$trait <- "EFGC" # Add a column indicating trait name

colnames(EFGC.params.fullposts) <- c("iteration", "cf.T0", "cf.Tm", "cf.q", "trait")

## Save output
# write_csv(EFGC.params.fullposts, "data-processed/EFGC.arctic.params.fullposts.csv")



##### Calculate trait values based on new TPC parameters
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

## Save output
# write_csv(EFGC.arctic.preds, "data-processed/EFGC.arctic.predictions.fullposts.csv")


##### Calculate the mean, median, and CIs of the TPC parameters (Tmin, Tmax, q, Topt)
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
# write_csv(EFGC.params.summary, "data-processed/EFGC.arctic.params.summary.csv")


##### Create a dataframe showing the mean, median and CIs of trait values at each temp
EFGC.predictions.summary <- calcPostQuants(EFGC.arctic.preds, "EFGC", Temp.xs)

## Save output
# write_csv(EFGC.predictions.summary, "data-processed/EFGC.arctic.predictions.summary.csv")


##### Plot
EFGC.predictions.summary <- read.csv("data-processed/EFGC.arctic.predictions.summary.csv")
EFGC.params.summary <- read.csv("data-processed/EFGC.arctic.params.summary.csv")


plot.EFGC <- EFGC.predictions.summary %>% 
  ggplot() +
  geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI), fill = "#0072B2", alpha = 0.5) +
  geom_line(aes(x = temperature, y = median), color = "#0072B2", linewidth = 1) +
  # Customize the axes and labels
  scale_x_continuous(limits = c(0, 46)) + 
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Eggs per female per gonotrophic cycle") +
  annotate("text", x = 2, y = 70, label = expression(paste(italic("EFGC"))), size = 5) +
  theme_bw()

plot.EFGC

# ggsave("figures/EFGC.arctic.quad.adjust.png", plot.EFGC, width = 10.3, height = 5.6)

## Compare original and adjusted TPC
plot.EFGC.all <- EFGC.predictions.summary %>% 
  ggplot() +
  # Original
  geom_ribbon(data = df.EFGC.nonarctic.quad.uni, aes(x = temp, ymin = X2.5., ymax = X97.5.), fill = "#868686FF", alpha = 0.5) +
  geom_line(data = df.EFGC.nonarctic.quad.uni, aes(x = temp, y = X50.), color = "#868686FF", linewidth = 1) +
  
  geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI), fill = "#0072B2", alpha = 0.5) +
  geom_line(aes(x = temperature, y = median), color = "#0072B2", linewidth = 1) +
  # Customize the axes and labels
  scale_x_continuous(limits = c(0, 46)) + 
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Eggs per female per gonotrophic cycle") +
  annotate("text", x = 2, y = 70, label = expression(paste(italic("EFGC"))), size = 5) +
  theme_bw()

plot.EFGC.all

# ggsave("figures/EFGC.arctic.vs.nonarctic.png", plot.EFGC.all, width = 10.3, height = 5.6)




##########
###### 3b. Approach 2: Fixed Tm and maximum trait values ----
##########


## Because the maximum trait values are always at T = (Tmax + Tmin)/2 (i.e. mid-point),
## We can substitute T = (Tmax + Tmin)/2 into the quadratic equation.
## Pmax = -q((Tm+T0)/2 - T0)((Tm+T0)/2 - Tm) --> Pmax = q(Tm-T0)^2/4
## Therefore, if Pmax is the same for two sets of q and T0, the relationships 
## between q and T0 are:
## q(Tm-T0)^2/4 = q'(Tm-T0')^2/4 (where q' and T0' are new q and Tmin)
## --> q(Tm-T0)^2 = q'(Tm-T0')^2
## --> q' = q(Tm-T0)^2/(Tm-T0')^2


#### 3bi. infection proportion (c) ----

## Get the original TPC parameters
c.iter.param <- data.frame(T0 = c.nonarctic.quad.uni$BUGSoutput$sims.list$cf.T0,
                           Tm = c.nonarctic.quad.uni$BUGSoutput$sims.list$cf.Tm,
                           q = c.nonarctic.quad.uni$BUGSoutput$sims.list$cf.q
                           )


c.iter.param <- c.iter.param %>% 
  mutate(new.T0 = T0 - T0.diff,
         new.q = q * (Tm - T0)**2/ (Tm - new.T0)**2)


## Create a dataframe showing the TPC parameters for each iteration
c.params.fullposts <- c.iter.param %>% 
  dplyr::select(new.T0, Tm, new.q)

c.params.fullposts$iteration <- seq(1:nrow(c.iter.param)) # Add a column indicating the number of MCMC iteration
c.params.fullposts <- relocate(c.params.fullposts, iteration, .before = new.T0)

c.params.fullposts$trait <- "c" # Add a column indicating trait name

colnames(c.params.fullposts) <- c("iteration", "cf.T0", "cf.Tm", "cf.q", "trait")

## Save output
write_csv(c.params.fullposts, "data-processed/c.arctic.params.fullposts.fixedTm.csv")



##### Calculate trait values based on new TPC parameters
c.arctic.preds <- data.frame() # Initialize an empty dataframe

for (i in 1:nrow(c.iter.param)) {
  # Calculate trait values for each MCMC iteration
  iter.preds <- quad(Temp.xs, 
                     T0 = c.iter.param$new.T0[i],
                     Tm = c.iter.param$Tm[i], 
                     q = c.iter.param$new.q[i]) 
  
  # Add the trait values as a column
  if (i == 1) {
    c.arctic.preds <- iter.preds
  }
  else {
    c.arctic.preds <- bind_cols(c.arctic.preds, iter.preds)
  }
}

# Transpose the dataset so that each row is a MCMC iteration (for calculating TPC parameters summary in the next step)
c.arctic.preds <- as.data.frame(t(c.arctic.preds))
colnames(c.arctic.preds) <- seq(1:ncol(c.arctic.preds))


##### Calculate the mean, median, and CIs of the TPC parameters (Tmin, Tmax, q, Topt)
c.T0 <- data.frame(term = "cf.T0",
                   mean = mean(c.params.fullposts$cf.T0),
                   sd = sd(c.params.fullposts$cf.T0),
                   lowerCI = quantile(c.params.fullposts$cf.T0, 0.025)[[1]],
                   lowerQ = quantile(c.params.fullposts$cf.T0, 0.25)[[1]],
                   median =  quantile(c.params.fullposts$cf.T0, 0.5)[[1]],
                   upperQ = quantile(c.params.fullposts$cf.T0, 0.75)[[1]],
                   upperCI = quantile(c.params.fullposts$cf.T0, 0.975)[[1]],
                   trait = "c")

c.Tm <- data.frame(term = "cf.Tm",
                   mean = mean(c.params.fullposts$cf.Tm),
                   sd = sd(c.params.fullposts$cf.Tm),
                   lowerCI = quantile(c.params.fullposts$cf.Tm, 0.025)[[1]],
                   lowerQ = quantile(c.params.fullposts$cf.Tm, 0.25)[[1]],
                   median =  quantile(c.params.fullposts$cf.Tm, 0.5)[[1]],
                   upperQ = quantile(c.params.fullposts$cf.Tm, 0.75)[[1]],
                   upperCI = quantile(c.params.fullposts$cf.Tm, 0.975)[[1]],
                   trait = "c")

c.q <- data.frame(term = "cf.q",
                  mean = mean(c.params.fullposts$cf.q),
                  sd = sd(c.params.fullposts$cf.q),
                  lowerCI = quantile(c.params.fullposts$cf.q, 0.025)[[1]],
                  lowerQ = quantile(c.params.fullposts$cf.q, 0.25)[[1]],
                  median =  quantile(c.params.fullposts$cf.q, 0.5)[[1]],
                  upperQ = quantile(c.params.fullposts$cf.q, 0.75)[[1]],
                  upperCI = quantile(c.params.fullposts$cf.q, 0.975)[[1]],
                  trait = "c")

# Calculate Topt for each iteration and calculate summary statistics (mean, sd, & quantiles)
c.Topt <- calcToptQuants(c.arctic.preds, "c", Temp.xs)

# Add Topt and Tbreadth to parameters summary data frame
c.params.summary <- bind_rows(c.T0, c.Tm, c.q, c.Topt)


## Save output
write_csv(c.params.summary, "data-processed/c.arctic.params.summary.fixedTm.csv")


## Since infection efficiency is a proportion, it cannot be greater than 1
## Replace values greater than 1 to 1
## Do this step after calculating TPC parameters
c.arctic.preds <- replace(c.arctic.preds, c.arctic.preds > 1, 1)

## Save output
write_csv(c.arctic.preds, "data-processed/c.arctic.predictions.fullposts.fixedTm.csv")


##### Create a dataframe showing the mean, median and CIs of trait values at each temp
c.predictions.summary <- calcPostQuants(c.arctic.preds, "c", Temp.xs)

## Save output
write_csv(c.predictions.summary, "data-processed/c.arctic.predictions.summary.fixedTm.csv")


##### Plot
c.predictions.summary <- read.csv("data-processed/c.arctic.predictions.summary.fixedTm.csv")
c.params.summary <- read.csv("data-processed/c.arctic.params.summary.fixedTm.csv")


plot.c <- c.predictions.summary %>% 
  ggplot() +
  geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI), fill = "#009E73", alpha = 0.5) +
  geom_line(aes(x = temperature, y = median), color = "#009E73", linewidth = 1) +
  # Customize the axes and labels
  scale_x_continuous(limits = c(0, 46)) + 
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Infection proportion") +
  annotate("text", x = 0, y = 0.95, label = expression(paste(italic("c"))), size = 5) +
  theme_bw()

plot.c

ggsave("figures/c.arctic.quad.adjust.fixedTm.png", plot.c, width = 10.3, height = 5.6)

## Compare original and adjusted TPC
plot.c.all <- c.predictions.summary %>% 
  ggplot() +
  # Original
  geom_ribbon(data = df.c.nonarctic.quad.uni, aes(x = temp, ymin = X2.5., ymax = X97.5.), fill = "#868686FF", alpha = 0.5) +
  geom_line(data = df.c.nonarctic.quad.uni, aes(x = temp, y = X50.), color = "#868686FF", linewidth = 1) +
  
  geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI), fill = "#009E73", alpha = 0.5) +
  geom_line(aes(x = temperature, y = median), color = "#009E73", linewidth = 1) +
  # Customize the axes and labels
  scale_x_continuous(limits = c(0, 46)) + 
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Infection proportion") +
  annotate("text", x = 0, y = 0.95, label = expression(paste(italic("c"))), size = 5) +
  theme_bw()

plot.c.all

ggsave("figures/c.arctic.vs.nonarctic.fixedTm.png", plot.c.all, width = 10.3, height = 5.6)



#### 3bii. Eggs per female per gonotrophic cycle (EFGC) ----

## Get the original TPC parameters
EFGC.iter.param <- data.frame(T0 = EFGC.nonarctic.quad.uni$BUGSoutput$sims.list$cf.T0,
                              Tm = EFGC.nonarctic.quad.uni$BUGSoutput$sims.list$cf.Tm,
                              q = EFGC.nonarctic.quad.uni$BUGSoutput$sims.list$cf.q
)


##### Perform the hot-old shift
EFGC.iter.param <- EFGC.iter.param %>% 
  mutate(new.T0 = T0 - T0.diff,
         new.q = q * (Tm - T0)**2/ (Tm - new.T0)**2)


## Create a dataframe showing the TPC parameters for each iteration
EFGC.params.fullposts <- EFGC.iter.param %>% 
  dplyr::select(new.T0, Tm, new.q)

EFGC.params.fullposts$iteration <- seq(1:nrow(EFGC.iter.param)) # Add a column indicating the number of MCMC iteration
EFGC.params.fullposts <- relocate(EFGC.params.fullposts, iteration, .before = new.T0)

EFGC.params.fullposts$trait <- "EFGC" # Add a column indicating trait name

colnames(EFGC.params.fullposts) <- c("iteration", "cf.T0", "cf.Tm", "cf.q", "trait")

## Save output
write_csv(EFGC.params.fullposts, "data-processed/EFGC.arctic.params.fullposts.fixedTm.csv")



##### Calculate trait values based on new TPC parameters
EFGC.arctic.preds <- data.frame() # Initialize an empty dataframe

for (i in 1:nrow(EFGC.iter.param)) {
  # Calculate trait values for each MCMC iteration
  iter.preds <- quad(Temp.xs, 
                     T0 = EFGC.iter.param$new.T0[i],
                     Tm = EFGC.iter.param$Tm[i], 
                     q = EFGC.iter.param$new.q[i]) 
  
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

## Save output
write_csv(EFGC.arctic.preds, "data-processed/EFGC.arctic.predictions.fullposts.fixedTm.csv")


##### Calculate the mean, median, and CIs of the TPC parameters (Tmin, Tmax, q, Topt)
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
write_csv(EFGC.params.summary, "data-processed/EFGC.arctic.params.summary.fixedTm.csv")


##### Create a dataframe showing the mean, median and CIs of trait values at each temp
EFGC.predictions.summary <- calcPostQuants(EFGC.arctic.preds, "EFGC", Temp.xs)

## Save output
write_csv(EFGC.predictions.summary, "data-processed/EFGC.arctic.predictions.summary.fixedTm.csv")


##### Plot
EFGC.predictions.summary <- read.csv("data-processed/EFGC.arctic.predictions.summary.fixedTm.csv")
EFGC.params.summary <- read.csv("data-processed/EFGC.arctic.params.summary.fixedTm.csv")


plot.EFGC <- EFGC.predictions.summary %>% 
  ggplot() +
  geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI), fill = "#0072B2", alpha = 0.5) +
  geom_line(aes(x = temperature, y = median), color = "#0072B2", linewidth = 1) +
  # Customize the axes and labels
  scale_x_continuous(limits = c(0, 46)) + 
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Eggs per female per gonotrophic cycle") +
  annotate("text", x = 2, y = 70, label = expression(paste(italic("EFGC"))), size = 5) +
  theme_bw()

plot.EFGC

ggsave("figures/EFGC.arctic.quad.adjust.fixedTm.png", plot.EFGC, width = 10.3, height = 5.6)

## Compare original and adjusted TPC
plot.EFGC.all <- EFGC.predictions.summary %>% 
  ggplot() +
  # Original
  geom_ribbon(data = df.EFGC.nonarctic.quad.uni, aes(x = temp, ymin = X2.5., ymax = X97.5.), fill = "#868686FF", alpha = 0.5) +
  geom_line(data = df.EFGC.nonarctic.quad.uni, aes(x = temp, y = X50.), color = "#868686FF", linewidth = 1) +
  
  geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI), fill = "#0072B2", alpha = 0.5) +
  geom_line(aes(x = temperature, y = median), color = "#0072B2", linewidth = 1) +
  # Customize the axes and labels
  scale_x_continuous(limits = c(0, 46)) + 
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Eggs per female per gonotrophic cycle") +
  annotate("text", x = 2, y = 70, label = expression(paste(italic("EFGC"))), size = 5) +
  theme_bw()

plot.EFGC.all

ggsave("figures/EFGC.arctic.vs.nonarctic.fixedTm.png", plot.EFGC.all, width = 10.3, height = 5.6)

