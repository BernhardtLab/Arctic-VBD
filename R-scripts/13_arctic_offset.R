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
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(type = "non-Arctic")

Temp.xs <- seq(0, 45, 0.1)

df.PDR.arctic.bri.inf <- data.frame(PDR.arctic.bri.inf$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(type = "Arctic")

df.PDR.all <- rbind(df.PDR.arctic.bri.inf, df.PDR.nonarctic.bri.uni.raneff.pop)

plot.PDR.all <- df.PDR.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = mean, color = type), linewidth = 1) +
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
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(type = "non-Arctic")

Temp.xs <- seq(0, 45, 0.1)

df.MDR.arctic.bri.inf <- data.frame(MDR.arctic.bri.inf$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.)%>% 
  mutate(type = "Arctic")

df.MDR.all <- rbind(df.MDR.arctic.bri.inf, df.MDR.nonarctic.bri.uni)

plot.MDR.all <- df.MDR.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = mean, color = type), linewidth = 1) +
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

ggsave("tpc/MDR.all.png", plot.MDR.all, width = 10.3, height = 5.6)


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
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(type = "non-Arctic")

Temp.xs <- seq(0, 45, 0.1)

df.EV.arctic.quad.inf <- data.frame(EV.arctic.quad.inf$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(type = "Arctic")

df.EV.all <- rbind(df.EV.arctic.quad.inf, df.EV.nonarctic.quad.uni.raneff.pop)

plot.EV.all <- df.EV.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = mean, color = type), linewidth = 1) +
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

ggsave("tpc/EV.all.png", plot.EV.all, width = 10.3, height = 5.6)

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
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(type = "non-Arctic")

Temp.xs <- seq(0, 45, 0.1)

df.pLA.arctic.quad.inf <- data.frame(pLA.arctic.quad.inf$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(type = "Arctic")

df.pLA.all <- rbind(df.pLA.arctic.quad.inf, df.pLA.nonarctic.quad.uni.raneff.pop)

plot.pLA.all <- df.pLA.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = mean, color = type), linewidth = 1) +
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

ggsave("tpc/pLA.all.png", plot.pLA.all, width = 10.3, height = 5.6)

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
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(type = "non-Arctic")

Temp.xs <- seq(0, 45, 0.1)

df.lf.arctic.quad.inf.raneff <- data.frame(lf.arctic.quad.inf.raneff$BUGSoutput$summary)[-(1:8),]

df.lf.arctic.quad.inf.raneff.pop <- df.lf.arctic.quad.inf.raneff %>% 
  filter(grepl("z.trait.mu.pred.pop", rownames(df.lf.arctic.quad.inf.raneff))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(type = "Arctic")

df.lf.all <- rbind(df.lf.arctic.quad.inf.raneff.pop, df.lf.nonarctic.quad.uni.raneff.pop)

plot.lf.all <- df.lf.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = mean, color = type), linewidth = 1) +
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

ggsave("tpc/lf.all.png", plot.lf.all, width = 10.3, height = 5.6)


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
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(type = "non-Arctic")

df.c.all <- rbind(df.c.nonarctic.quad.uni)

plot.c.all <- df.c.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = mean, color = type), linewidth = 1) +
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

ggsave("tpc/c.all.png", plot.c.all, width = 10.3, height = 5.6)

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
  dplyr::select(temp, mean, sd, X2.5., X97.5.) %>% 
  mutate(type = "non-Arctic")

df.all <- rbind(df.EFGC.nonarctic.quad.uni)

plot.EFGC.all <- df.all %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = type), alpha = 0.5) +
  geom_line(aes(y = mean, color = type), linewidth = 1) +
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

ggsave("tpc/EFGC.all.png", plot.EFGC.all, width = 10.3, height = 5.6)

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
##    
##    No idea how to do this yet. Will come back to it


##########
###### 3a. Approach 1: hot-cold shift ----
##########

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

### Plot

c.nonarctic.preds.new <- data.frame(temp = Temp.xs,
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



