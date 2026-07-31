## Lilian Chan, University of Guelph
## Arctic vector-borne disease transmission suitability model
##
## Purpose: Investigate the relationship between thermal limits (Tmin, Tmax) and 
## latitude
## 
## Table of content:
##    0. Set-up workspace
##    1. Latitudinal regression
##    2. Summary table


# 0. Set-up workspace -----------------------------------------------------
library(tidyverse)
library(readxl)
library(janitor)
library(R2jags)
library(MASS)
library(cowplot)
library(RColorBrewer)
library(metafor)

# Load functions
source("R-scripts/00_Functions.R")

#  1. Latitudinal regression ---------------------------------------------------


##### Mosquito adult lifespan (lf) #####

# Load R2jags model output
load("R-scripts/R2jags-objects/all-mods/lf.nonarctic.quad.uni.Rdata") # Non-arctic
load("R-scripts/R2jags-objects/all-mods/lf.arctic.quad.uni.Rdata") # Arctic

# Load data
data.lf.all <- read_csv("data-processed/TraitData_lf.csv")
data.lf.all <- data.lf.all %>%  # absolute latitude
  mutate(latitude = abs(latitude))

# Subset data
data.lf.arctic <- subset(data.lf.all, type == "Arctic") # Arctic species
data.lf.nonarctic <- subset(data.lf.all, type == "non-Arctic") # Non-Arctic species

###### Non-Arctic ######
# Get Tmin, Tmax, and q from each random effect
lf.sims <- lf.nonarctic.quad.uni$BUGSoutput$sims.list

# parameter values for each random effect
lf.T0.nonarctic.fullpost <- sweep(lf.sims$T0, 1, lf.sims$cf.T0, "+")
lf.Tm.nonarctic.fullpost <- sweep(lf.sims$Tm, 1, lf.sims$cf.Tm, "+")
lf.q.nonarctic.fullpost  <- sweep(lf.sims$q,  1, lf.sims$cf.q,  "+")

# summarize posterior distributions
lf.T0.nonarctic <- t(apply(lf.T0.nonarctic.fullpost, 2, quantile, c(0.025,0.5,0.975)))
lf.Tm.nonarctic <- t(apply(lf.Tm.nonarctic.fullpost, 2, quantile, c(0.025,0.5,0.975)))
lf.q.nonarctic <- t(apply(lf.q.nonarctic.fullpost, 2, quantile, c(0.025,0.5,0.975)))

# change column names
colnames(lf.T0.nonarctic) <- c("T0_2.5", "T0_50", "T0_97.5")
colnames(lf.Tm.nonarctic) <- c("Tm_2.5", "Tm_50", "Tm_97.5")
colnames(lf.q.nonarctic)  <- c("q_2.5",  "q_50",  "q_97.5")


lf.nonarctic.id.info <- data.lf.nonarctic %>%
  group_by(species, citation) %>% 
  mutate(unique_id = cur_group_id()) %>% 
  group_by(unique_id) %>% 
  dplyr::select(unique_id, genus, species, citation, latitude, longitude, type) %>% 
  unique() %>% # drop duplicate
  arrange(unique_id)

# Data from Tsuda et al. 1994 were collected from two locations. We will use the mean latitude
lf.nonarctic.id.info <- lf.nonarctic.id.info %>% 
  mutate(latitude = ifelse(
    str_detect(citation, "Tsuda"),
    mean(latitude[str_detect(citation, "Tsuda")]),
    latitude)
    ) %>% 
  mutate(longitude = ifelse(
    str_detect(citation, "Tsuda"),
    NA,
    longitude)
    ) %>% 
  unique() # drop duplicate
  

lf.TPC.pars.nonarctic <- bind_cols(lf.nonarctic.id.info, lf.T0.nonarctic, lf.Tm.nonarctic, lf.q.nonarctic)
lf.TPC.pars.nonarctic

###### Arctic ######

# Get Tmin, Tmax, and q from each random effect
lf.sims.arctic <- lf.arctic.quad.uni$BUGSoutput$sims.list

# parameter values for each random effect
lf.T0.arctic.fullpost <- sweep(lf.sims.arctic$T0, 1, lf.sims.arctic$cf.T0, "+")
lf.Tm.arctic.fullpost <- sweep(lf.sims.arctic$Tm, 1, lf.sims.arctic$cf.Tm, "+")
lf.q.arctic.fullpost  <- sweep(lf.sims.arctic$q,  1, lf.sims.arctic$cf.q,  "+")

# summarize posterior distributions
lf.T0.arctic <- t(apply(lf.T0.arctic.fullpost, 2, quantile, c(0.025,0.5,0.975)))
lf.Tm.arctic <- t(apply(lf.Tm.arctic.fullpost, 2, quantile, c(0.025,0.5,0.975)))
lf.q.arctic <- t(apply(lf.q.arctic.fullpost, 2, quantile, c(0.025,0.5,0.975)))

# change column names
colnames(lf.T0.arctic) <- c("T0_2.5", "T0_50", "T0_97.5")
colnames(lf.Tm.arctic) <- c("Tm_2.5", "Tm_50", "Tm_97.5")
colnames(lf.q.arctic)  <- c("q_2.5",  "q_50",  "q_97.5")

lf.arctic.id.info <- data.lf.arctic %>%
  group_by(species, citation) %>% 
  mutate(unique_id = cur_group_id()) %>% 
  group_by(unique_id) %>% 
  dplyr::select(unique_id, genus, species, citation, latitude, longitude, type) %>% 
  unique() %>% # drop duplicate
  arrange(unique_id)


lf.TPC.pars.arctic <- bind_cols(lf.arctic.id.info, lf.T0.arctic, lf.Tm.arctic, lf.q.arctic)
lf.TPC.pars <- bind_rows(lf.TPC.pars.arctic, lf.TPC.pars.nonarctic)


###### TPC parameter-latitudinal analysis ######

# Tmin
## Calculate sampling variances
lf.T0.vars <- bind_cols(lf.T0.arctic.fullpost, lf.T0.nonarctic.fullpost)
lf.T0.vars <- apply(lf.T0.vars, 2, var)
lf.TPC.pars$T0_var <- lf.T0.vars

## Regression
lf.T0.fit <- rma.mv(yi = T0_50, 
                    V = T0_var,
                    mods = ~latitude,
                    data = lf.TPC.pars,
                    method = "REML")
  
summary(lf.T0.fit)


## Get predicted values
lf.T0.pred <- data.frame(latitude = seq(min(lf.TPC.pars$latitude),
                                        max(lf.TPC.pars$latitude),
                                        length.out = 100))
lf.T0.newdata <- predict(lf.T0.fit,
                         newmods = lf.T0.pred$latitude)

lf.T0.pred$pred <- lf.T0.newdata$pred
lf.T0.pred$lower <- lf.T0.newdata$ci.lb
lf.T0.pred$upper <- lf.T0.newdata$ci.ub


# Tmax
## Calculate sampling variances
lf.Tm.vars <- bind_cols(lf.Tm.arctic.fullpost, lf.Tm.nonarctic.fullpost)
lf.Tm.vars <- apply(lf.Tm.vars, 2, var)
lf.TPC.pars$Tm_var <- lf.Tm.vars

## Regression
lf.Tm.fit <- rma.mv(yi = Tm_50, 
                    V = Tm_var,
                    mods = ~latitude,
                    data = lf.TPC.pars,
                    method = "REML")

summary(lf.Tm.fit)

## Get predicted values
lf.Tm.pred <- data.frame(latitude = seq(min(lf.TPC.pars$latitude),
                                        max(lf.TPC.pars$latitude),
                                        length.out = 100))
lf.Tm.newdata <- predict(lf.Tm.fit,
                         newmods = lf.Tm.pred$latitude)

lf.Tm.pred$pred <- lf.Tm.newdata$pred
lf.Tm.pred$lower <- lf.Tm.newdata$ci.lb
lf.Tm.pred$upper <- lf.Tm.newdata$ci.ub


# q
## Calculate sampling variances
lf.q.vars <- bind_cols(lf.q.arctic.fullpost, lf.q.nonarctic.fullpost)
lf.q.vars <- apply(lf.q.vars, 2, var)
lf.TPC.pars$q_var <- lf.q.vars

## Regression
lf.q.fit <- rma.mv(yi = q_50, 
                   V = q_var,
                   mods = ~latitude,
                   data = lf.TPC.pars,
                   method = "REML")

summary(lf.q.fit)

## Get predicted values
lf.q.pred <- data.frame(latitude = seq(min(lf.TPC.pars$latitude),
                                        max(lf.TPC.pars$latitude),
                                        length.out = 100))
lf.q.newdata <- predict(lf.q.fit,
                         newmods = lf.q.pred$latitude)

lf.q.pred$pred <- lf.q.newdata$pred
lf.q.pred$lower <- lf.q.newdata$ci.lb
lf.q.pred$upper <- lf.q.newdata$ci.ub


# Plot
lf.Tmin.lat <- lf.TPC.pars %>% 
  ggplot(aes(x = latitude)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  # data
  geom_errorbar(aes(ymin = T0_2.5, ymax = T0_97.5), width = 1) +
  geom_point(aes(y = T0_50, colour = species)) +
  
  # fitted rma.mv model
  geom_line(data = lf.T0.pred, aes(y = pred)) +
  geom_ribbon(data = lf.T0.pred, aes(ymin = lower, ymax = upper), alpha = 0.2) + 
  
  labs(
    x = "absolute latitude",
    y = "Tmin",
    title = expression(paste("Adult Lifespan (",italic(lf),")"))
  ) +
  
  scale_colour_manual(values = c("albopictus" = "#67000D",
                                 "vexans" = "#9ECAE1",
                                 "punctor" = "#4292C6",
                                 "cinereus" = "#2171B5",
                                 "communis" = "#08519C",
                                 "impiger" = "#08306B"),
                      name = element_blank(), # No legend title
                      breaks = c("albopictus",
                                 "vexans",
                                 "punctor",
                                 "cinereus",
                                 "communis",
                                 "Ae. impiger"),
                      labels = c("Ae. albopictus",
                                 "Ae. vexans",
                                 "Ae. punctor",
                                 "Ae. cinereus",
                                 "Ae. communis",
                                 "Ae. impiger")) +
  
  theme_bw()

lf.Tmin.lat


lf.Tmax.lat <- lf.TPC.pars %>% 
  ggplot(aes(x = latitude)) +
  geom_vline(xintercept = 0, linetype = "dashed") +

  # data
  geom_errorbar(aes(ymin = Tm_2.5, ymax = Tm_97.5), width = 1) +
  geom_point(aes(y = Tm_50, colour = species)) +
  
  # fitted rma.mv model
  geom_line(data = lf.Tm.pred, aes(y = pred)) +
  geom_ribbon(data = lf.Tm.pred, aes(ymin = lower, ymax = upper), alpha = 0.2) + 
  
  labs(
    x = "absolute latitude",
    y = "Tmax",
    title = expression(paste("Adult Lifespan (",italic(lf),")"))
  ) +
  
  scale_colour_manual(values = c("albopictus" = "#67000D",
                                 "vexans" = "#9ECAE1",
                                 "punctor" = "#4292C6",
                                 "cinereus" = "#2171B5",
                                 "communis" = "#08519C",
                                 "impiger" = "#08306B"),
                      name = element_blank(), # No legend title
                      breaks = c("albopictus",
                                 "vexans",
                                 "punctor",
                                 "cinereus",
                                 "communis",
                                 "Ae. impiger"),
                      labels = c("Ae. albopictus",
                                 "Ae. vexans",
                                 "Ae. punctor",
                                 "Ae. cinereus",
                                 "Ae. communis",
                                 "Ae. impiger")) +
  
  theme_bw()

lf.Tmax.lat

lf.q.lat <- lf.TPC.pars %>% 
  ggplot(aes(x = latitude)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  # data
  geom_errorbar(aes(ymin = q_2.5, ymax = q_97.5), width = 1) +
  geom_point(aes(y = q_50, colour = species)) +
  
  # fitted rma.mv model
  geom_line(data = lf.q.pred, aes(y = pred)) +
  geom_ribbon(data = lf.q.pred, aes(ymin = lower, ymax = upper), alpha = 0.2) + 
  
  labs(
    x = "absolute latitude",
    y = "q",
    title = expression(paste("Adult Lifespan (",italic(lf),")"))
  ) +
  
  scale_colour_manual(values = c("albopictus" = "#67000D",
                                 "vexans" = "#9ECAE1",
                                 "punctor" = "#4292C6",
                                 "cinereus" = "#2171B5",
                                 "communis" = "#08519C",
                                 "impiger" = "#08306B"),
                      name = element_blank(), # No legend title
                      breaks = c("albopictus",
                                 "vexans",
                                 "punctor",
                                 "cinereus",
                                 "communis",
                                 "Ae. impiger"),
                      labels = c("Ae. albopictus",
                                 "Ae. vexans",
                                 "Ae. punctor",
                                 "Ae. cinereus",
                                 "Ae. communis",
                                 "Ae. impiger")) +
  
  theme_bw()

lf.q.lat

TPC.params.lat <- plot_grid(lf.Tmin.lat, lf.Tmax.lat, lf.q.lat, align = "v", ncol = 1)
TPC.params.lat

ggsave("figures/lf.quad.TPC.params.lat.png", TPC.params.lat,
       width = 10.3, height = 10)



##### Egg  viability (EV) #####

# Load R2jags model output
load("R-scripts/R2jags-objects/all-mods/EV.nonarctic.quad.uni.Rdata") # Non-arctic
load("R-scripts/R2jags-objects/all-mods/EV.arctic.quad.uni.Rdata") # Arctic

# Load data
data.EV.all <- read_csv("data-processed/TraitData_EV.csv")
data.EV.all <- data.EV.all %>% # absolute latitude
mutate(latitude = abs(latitude)) 

# Subset data
data.EV.arctic <- subset(data.EV.all, type == "Arctic") # Arctic species
data.EV.nonarctic <- subset(data.EV.all, type == "non-Arctic") # Non-Arctic species

###### Non-Arctic ######
# Get Tmin, Tmax, and q from each random effect
EV.sims <- EV.nonarctic.quad.uni$BUGSoutput$sims.list

# parameter values for each unique_id
EV.T0.nonarctic.fullpost <- sweep(EV.sims$T0, 1, EV.sims$cf.T0, "+")
EV.Tm.nonarctic.fullpost <- sweep(EV.sims$Tm, 1, EV.sims$cf.Tm, "+")
EV.q.nonarctic.fullpost  <- sweep(EV.sims$q,  1, EV.sims$cf.q,  "+")

# summarize posterior distributions
EV.T0.nonarctic <- t(apply(EV.T0.nonarctic.fullpost, 2, quantile, c(0.025,0.5,0.975)))
EV.Tm.nonarctic <- t(apply(EV.Tm.nonarctic.fullpost, 2, quantile, c(0.025,0.5,0.975)))
EV.q.nonarctic <- t(apply(EV.q.nonarctic.fullpost, 2, quantile, c(0.025,0.5,0.975)))

# change column names
colnames(EV.T0.nonarctic) <- c("T0_2.5", "T0_50", "T0_97.5")
colnames(EV.Tm.nonarctic) <- c("Tm_2.5", "Tm_50", "Tm_97.5")
colnames(EV.q.nonarctic)  <- c("q_2.5",  "q_50",  "q_97.5")


EV.nonarctic.id.info <- data.EV.nonarctic %>%
  group_by(species, citation) %>% 
  mutate(unique_id = cur_group_id()) %>% 
  group_by(unique_id) %>% 
  dplyr::select(unique_id, genus, species, citation, latitude, longitude, type) %>% 
  unique() %>% # drop duplicate
  arrange(unique_id)

EV.TPC.pars.nonarctic <- bind_cols(EV.nonarctic.id.info, EV.T0.nonarctic, EV.Tm.nonarctic, EV.q.nonarctic)
EV.TPC.pars.nonarctic

###### Arctic ######
EV.T0.arctic.fullpost <- EV.arctic.quad.uni$BUGSoutput$sims.list$cf.T0
EV.Tm.arctic.fullpost <- EV.arctic.quad.uni$BUGSoutput$sims.list$cf.Tm
EV.q.arctic.fullpost <- EV.arctic.quad.uni$BUGSoutput$sims.list$cf.q

# summarize posterior distributions
EV.T0.arctic <- t(quantile(EV.T0.arctic.fullpost, c(0.025, 0.5, 0.975)))
EV.Tm.arctic <- t(quantile(EV.Tm.arctic.fullpost, c(0.025, 0.5, 0.975)))
EV.q.arctic <- t(quantile(EV.q.arctic.fullpost, c(0.025, 0.5, 0.975)))

# change column names
colnames(EV.T0.arctic) <- c("T0_2.5", "T0_50", "T0_97.5")
colnames(EV.Tm.arctic) <- c("Tm_2.5", "Tm_50", "Tm_97.5")
colnames(EV.q.arctic)  <- c("q_2.5",  "q_50",  "q_97.5")

EV.arctic.id.info <- data.EV.arctic %>%
  dplyr::select(genus, species, citation, latitude, longitude, type) %>% 
  unique() # drop duplicate

EV.TPC.pars.arctic <- bind_cols(EV.arctic.id.info, EV.T0.arctic, EV.Tm.arctic, EV.q.arctic)
EV.TPC.pars <- bind_rows(EV.TPC.pars.arctic, EV.TPC.pars.nonarctic)


###### TPC parameter-latitudinal analysis ######

# Tmin
## Calculate sampling variances
EV.T0.vars <- bind_cols(EV.T0.arctic.fullpost, EV.T0.nonarctic.fullpost)
EV.T0.vars <- apply(EV.T0.vars, 2, var)
EV.TPC.pars$T0_var <- EV.T0.vars

## Regression
EV.T0.fit <- rma.mv(yi = T0_50, 
                    V = T0_var,
                    mods = ~latitude,
                    data = EV.TPC.pars,
                    method = "REML")

summary(EV.T0.fit)


## Get predicted values
EV.T0.pred <- data.frame(latitude = seq(min(EV.TPC.pars$latitude),
                                        max(EV.TPC.pars$latitude),
                                        length.out = 100))
EV.T0.newdata <- predict(EV.T0.fit,
                         newmods = EV.T0.pred$latitude)

EV.T0.pred$pred <- EV.T0.newdata$pred
EV.T0.pred$lower <- EV.T0.newdata$ci.lb
EV.T0.pred$upper <- EV.T0.newdata$ci.ub


# Tmax
## Calculate sampling variances
EV.Tm.vars <- bind_cols(EV.Tm.arctic.fullpost, EV.Tm.nonarctic.fullpost)
EV.Tm.vars <- apply(EV.Tm.vars, 2, var)
EV.TPC.pars$Tm_var <- EV.Tm.vars

## Regression
EV.Tm.fit <- rma.mv(yi = Tm_50, 
                    V = Tm_var,
                    mods = ~latitude,
                    data = EV.TPC.pars,
                    method = "REML")

summary(EV.Tm.fit)

## Get predicted values
EV.Tm.pred <- data.frame(latitude = seq(min(EV.TPC.pars$latitude),
                                        max(EV.TPC.pars$latitude),
                                        length.out = 100))
EV.Tm.newdata <- predict(EV.Tm.fit,
                         newmods = EV.Tm.pred$latitude)

EV.Tm.pred$pred <- EV.Tm.newdata$pred
EV.Tm.pred$lower <- EV.Tm.newdata$ci.lb
EV.Tm.pred$upper <- EV.Tm.newdata$ci.ub


# q
## Calculate sampling variances
EV.q.vars <- bind_cols(EV.q.arctic.fullpost, EV.q.nonarctic.fullpost)
EV.q.vars <- apply(EV.q.vars, 2, var)
EV.TPC.pars$q_var <- EV.q.vars

## Regression
EV.q.fit <- rma.mv(yi = q_50, 
                   V = q_var,
                   mods = ~latitude,
                   data = EV.TPC.pars,
                   method = "REML")

summary(EV.q.fit)

## Get predicted values
EV.q.pred <- data.frame(latitude = seq(min(EV.TPC.pars$latitude),
                                       max(EV.TPC.pars$latitude),
                                       length.out = 100))
EV.q.newdata <- predict(EV.q.fit,
                        newmods = EV.q.pred$latitude)

EV.q.pred$pred <- EV.q.newdata$pred
EV.q.pred$lower <- EV.q.newdata$ci.lb
EV.q.pred$upper <- EV.q.newdata$ci.ub


# Plot
EV.Tmin.lat <- EV.TPC.pars %>% 
  ggplot(aes(x = latitude)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  # data
  geom_errorbar(aes(ymin = T0_2.5, ymax = T0_97.5), width = 1) +
  geom_point(aes(y = T0_50, colour = species)) +
  
  # fitted rma.mv model
  geom_line(data = EV.T0.pred, aes(y = pred)) +
  geom_ribbon(data = EV.T0.pred, aes(ymin = lower, ymax = upper), alpha = 0.2) + 
  
  labs(
    x = "absolute latitude",
    y = "Tmin",
    title = expression(paste("Egg Viability (",italic(EV),")"))
  ) +
  
  scale_colour_manual(values = c("albopictus" = "#67000D",
                                 "triseriatus" = "#CB181D",
                                 "dorsalis" = "#FB6A4A", 
                                 "nigromaculis" = "#FC9272",
                                 "vexans" = "#08519C"),
                      name = element_blank(), # No legend title
                      breaks = c("albopictus",
                                 "triseriatus",
                                 "dorsalis", 
                                 "nigromaculis",
                                 "vexans"),
                      labels = c("Ae. albopictus", "Ae. triseriatus",
                                 "Ae. dorsalis", "Ae. nigromaculis", 
                                 "Ae. vexans")) +
  
  theme_bw()

EV.Tmin.lat

EV.Tmax.lat <- EV.TPC.pars %>% 
  ggplot(aes(x = latitude)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  # data
  geom_errorbar(aes(ymin = Tm_2.5, ymax = Tm_97.5), width = 1) +
  geom_point(aes(y = Tm_50, colour = species)) +
  
  # fitted rma.mv model
  geom_line(data = EV.Tm.pred, aes(y = pred)) +
  geom_ribbon(data = EV.Tm.pred, aes(ymin = lower, ymax = upper), alpha = 0.2) + 
  
  labs(
    x = "absolute latitude",
    y = "Tmax",
    title = expression(paste("Egg Viability (",italic(EV),")"))
  ) +
  
  scale_colour_manual(values = c("albopictus" = "#67000D",
                                 "triseriatus" = "#CB181D",
                                 "dorsalis" = "#FB6A4A", 
                                 "nigromaculis" = "#FC9272",
                                 "vexans" = "#08519C"),
                      name = element_blank(), # No legend title
                      breaks = c("albopictus",
                                 "triseriatus",
                                 "dorsalis", 
                                 "nigromaculis",
                                 "vexans"),
                      labels = c("Ae. albopictus", "Ae. triseriatus",
                                 "Ae. dorsalis", "Ae. nigromaculis", 
                                 "Ae. vexans")) +
  
  theme_bw()

EV.Tmax.lat

EV.q.lat <- EV.TPC.pars %>% 
  ggplot(aes(x = latitude)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  # data
  geom_errorbar(aes(x = latitude, ymin = q_2.5, ymax = q_97.5), width = 1) +
  geom_point(aes(x = latitude, y = q_50, colour = species)) +
  
  # fitted rma.mv model
  geom_line(data = EV.q.pred, aes(y = pred)) +
  geom_ribbon(data = EV.q.pred, aes(ymin = lower, ymax = upper), alpha = 0.2) + 
  
  labs(
    x = "absolute latitude",
    y = "q",
    title = expression(paste("Egg Viability (",italic(EV),")"))
  ) +
  
  scale_colour_manual(values = c("albopictus" = "#67000D",
                                 "triseriatus" = "#CB181D",
                                 "dorsalis" = "#FB6A4A", 
                                 "nigromaculis" = "#FC9272",
                                 "vexans" = "#08519C"),
                      name = element_blank(), # No legend title
                      breaks = c("albopictus",
                                 "triseriatus",
                                 "dorsalis", 
                                 "nigromaculis",
                                 "vexans"),
                      labels = c("Ae. albopictus", "Ae. triseriatus",
                                 "Ae. dorsalis", "Ae. nigromaculis", 
                                 "Ae. vexans")) +
  
  theme_bw()

EV.q.lat

TPC.params.lat <- plot_grid(EV.Tmin.lat, EV.Tmax.lat, EV.q.lat, align = "v", ncol = 1)
TPC.params.lat

ggsave("figures/EV.quad.TPC.params.lat.png", TPC.params.lat,
       width = 10.3, height = 10)



##### Larval-to-adult survival (pLA) #####

# Load R2jags model output
load("R-scripts/R2jags-objects/all-mods/pLA.nonarctic.quad.uni.Rdata") # Non-arctic
load("R-scripts/R2jags-objects/all-mods/pLA.arctic.quad.uni.Rdata") # Arctic

# Load data
data.pLA.all <- read_csv("data-processed/TraitData_pLA.csv")
data.pLA.all <- data.pLA.all %>% # absolute latitude
  mutate(latitude = abs(latitude)) 

# Subset data
data.pLA.arctic <- subset(data.pLA.all, type == "Arctic") # Arctic species
data.pLA.nonarctic <- subset(data.pLA.all, type == "non-Arctic") # Non-Arctic species

###### Non-Arctic ######
# Get Tmin, Tmax, and q from each random effect
pLA.sims <- pLA.nonarctic.quad.uni$BUGSoutput$sims.list

# parameter values for each unique_id
pLA.T0.nonarctic.fullpost <- sweep(pLA.sims$T0, 1, pLA.sims$cf.T0, "+")
pLA.Tm.nonarctic.fullpost <- sweep(pLA.sims$Tm, 1, pLA.sims$cf.Tm, "+")
pLA.q.nonarctic.fullpost  <- sweep(pLA.sims$q,  1, pLA.sims$cf.q,  "+")

# summarize posterior distributions
pLA.T0.nonarctic <- t(apply(pLA.T0.nonarctic.fullpost, 2, quantile, c(0.025,0.5,0.975)))
pLA.Tm.nonarctic <- t(apply(pLA.Tm.nonarctic.fullpost, 2, quantile, c(0.025,0.5,0.975)))
pLA.q.nonarctic <- t(apply(pLA.q.nonarctic.fullpost, 2, quantile, c(0.025,0.5,0.975)))

# change column names
colnames(pLA.T0.nonarctic) <- c("T0_2.5", "T0_50", "T0_97.5")
colnames(pLA.Tm.nonarctic) <- c("Tm_2.5", "Tm_50", "Tm_97.5")
colnames(pLA.q.nonarctic)  <- c("q_2.5",  "q_50",  "q_97.5")


pLA.nonarctic.id.info <- data.pLA.nonarctic %>%
  group_by(species, citation) %>% 
  mutate(unique_id = cur_group_id()) %>% 
  group_by(unique_id) %>% 
  dplyr::select(unique_id, genus, species, citation, latitude, longitude, type) %>% 
  unique() %>% # drop duplicate
  arrange(unique_id)

pLA.TPC.pars.nonarctic <- bind_cols(pLA.nonarctic.id.info, pLA.T0.nonarctic, pLA.Tm.nonarctic, pLA.q.nonarctic)
pLA.TPC.pars.nonarctic

###### Arctic ######
pLA.T0.arctic.fullpost <- pLA.arctic.quad.uni$BUGSoutput$sims.list$cf.T0
pLA.Tm.arctic.fullpost <- pLA.arctic.quad.uni$BUGSoutput$sims.list$cf.Tm
pLA.q.arctic.fullpost <- pLA.arctic.quad.uni$BUGSoutput$sims.list$cf.q

# summarize posterior distributions
pLA.T0.arctic <- t(quantile(pLA.T0.arctic.fullpost, c(0.025, 0.5, 0.975)))
pLA.Tm.arctic <- t(quantile(pLA.Tm.arctic.fullpost, c(0.025, 0.5, 0.975)))
pLA.q.arctic <- t(quantile(pLA.q.arctic.fullpost, c(0.025, 0.5, 0.975)))

# change column names
colnames(pLA.T0.arctic) <- c("T0_2.5", "T0_50", "T0_97.5")
colnames(pLA.Tm.arctic) <- c("Tm_2.5", "Tm_50", "Tm_97.5")
colnames(pLA.q.arctic)  <- c("q_2.5",  "q_50",  "q_97.5")

# Since Arctic data do not have sufficient groups to include random effects, we will use the mean latitude
pLA.arctic.id.info <- data.frame(
  unique_id = 0,
  genus = "Aedes",
  species = "Arctic spp.",
  citation = NA,
  latitude = mean(data.pLA.arctic$latitude),
  longitude = NA,
  type = "Arctic"
  )
  

pLA.TPC.pars.arctic <- bind_cols(pLA.arctic.id.info, pLA.T0.arctic, pLA.Tm.arctic, pLA.q.arctic)

pLA.TPC.pars <- bind_rows(pLA.TPC.pars.arctic, pLA.TPC.pars.nonarctic)


###### TPC parameter-latitudinal analysis ######

# Tmin
## Calculate sampling variances
pLA.T0.vars <- bind_cols(pLA.T0.arctic.fullpost, pLA.T0.nonarctic.fullpost)
pLA.T0.vars <- apply(pLA.T0.vars, 2, var)
pLA.TPC.pars$T0_var <- pLA.T0.vars


## Regression
pLA.T0.fit <- rma.mv(yi = T0_50, 
                    V = T0_var,
                    mods = ~latitude,
                    data = pLA.TPC.pars,
                    method = "REML")

summary(pLA.T0.fit)


## Get predicted values
pLA.T0.pred <- data.frame(latitude = seq(min(pLA.TPC.pars$latitude),
                                        max(pLA.TPC.pars$latitude),
                                        length.out = 100))
pLA.T0.newdata <- predict(pLA.T0.fit,
                         newmods = pLA.T0.pred$latitude)

pLA.T0.pred$pred <- pLA.T0.newdata$pred
pLA.T0.pred$lower <- pLA.T0.newdata$ci.lb
pLA.T0.pred$upper <- pLA.T0.newdata$ci.ub


# Tmax
## Calculate sampling variances
pLA.Tm.vars <- bind_cols(pLA.Tm.arctic.fullpost, pLA.Tm.nonarctic.fullpost)
pLA.Tm.vars <- apply(pLA.Tm.vars, 2, var)
pLA.TPC.pars$Tm_var <- pLA.Tm.vars

## Regression
pLA.Tm.fit <- rma.mv(yi = Tm_50, 
                    V = Tm_var,
                    mods = ~latitude,
                    data = pLA.TPC.pars,
                    method = "REML")

summary(pLA.Tm.fit)

## Get predicted values
pLA.Tm.pred <- data.frame(latitude = seq(min(pLA.TPC.pars$latitude),
                                        max(pLA.TPC.pars$latitude),
                                        length.out = 100))
pLA.Tm.newdata <- predict(pLA.Tm.fit,
                         newmods = pLA.Tm.pred$latitude)

pLA.Tm.pred$pred <- pLA.Tm.newdata$pred
pLA.Tm.pred$lower <- pLA.Tm.newdata$ci.lb
pLA.Tm.pred$upper <- pLA.Tm.newdata$ci.ub


# q
## Calculate sampling variances
pLA.q.vars <- bind_cols(pLA.q.arctic.fullpost, pLA.q.nonarctic.fullpost)
pLA.q.vars <- apply(pLA.q.vars, 2, var)
pLA.TPC.pars$q_var <- pLA.q.vars

## Regression
pLA.q.fit <- rma.mv(yi = q_50, 
                   V = q_var,
                   mods = ~latitude,
                   data = pLA.TPC.pars,
                   method = "REML")

summary(pLA.q.fit)

## Get predicted values
pLA.q.pred <- data.frame(latitude = seq(min(pLA.TPC.pars$latitude),
                                       max(pLA.TPC.pars$latitude),
                                       length.out = 100))
pLA.q.newdata <- predict(pLA.q.fit,
                        newmods = pLA.q.pred$latitude)

pLA.q.pred$pred <- pLA.q.newdata$pred
pLA.q.pred$lower <- pLA.q.newdata$ci.lb
pLA.q.pred$upper <- pLA.q.newdata$ci.ub


# Plot
pLA.Tmin.lat <- pLA.TPC.pars %>% 
  ggplot(aes(x = latitude)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  # data
  geom_errorbar(aes(ymin = T0_2.5, ymax = T0_97.5), width = 1) +
  geom_point(aes(y = T0_50, colour = species)) +
  
  # fitted rma.mv model
  geom_line(data = pLA.T0.pred, aes(y = pred)) +
  geom_ribbon(data = pLA.T0.pred, aes(ymin = lower, ymax = upper), alpha = 0.2) + 
  
  labs(
    x = "absolute latitude",
    y = "Tmin",
    title = expression(paste("Larval-to-Adult Survival (",italic(pLA),")"))
  ) +
  
  scale_colour_manual(values = c("sollicitans" = "#67000D",
                                 "triseriatus" = "#CB181D",
                                 "albopictus" = "#FB6A4A", 
                                 "nigromaculis" = "#FC9272",
                                 "Arctic spp." = "#08519C"),
                      name = element_blank(), # No legend title
                      breaks = c("sollicitans",
                                 "triseriatus",
                                 "albopictus", 
                                 "nigromaculis",
                                 "Arctic spp."),
                      labels = c("Ae. sollicitans", "Ae. triseriatus",
                                 "Ae. albopictus", "Ae. nigromaculis", 
                                 "Arctic spp.")) +
  
  theme_bw()

pLA.Tmin.lat

pLA.Tmax.lat <- pLA.TPC.pars %>% 
  ggplot(aes(x = latitude)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  # data
  geom_errorbar(aes(ymin = Tm_2.5, ymax = Tm_97.5), width = 1) +
  geom_point(aes(y = Tm_50, colour = species)) +
  
  # fitted rma.mv model
  geom_line(data = pLA.Tm.pred, aes(y = pred)) +
  geom_ribbon(data = pLA.Tm.pred, aes(ymin = lower, ymax = upper), alpha = 0.2) + 
  
  labs(
    x = "absolute latitude",
    y = "Tmax",
    title = expression(paste("Larval-to-Adult Survival (",italic(pLA),")"))
  ) +
  
  scale_colour_manual(values = c("sollicitans" = "#67000D",
                                 "triseriatus" = "#CB181D",
                                 "albopictus" = "#FB6A4A", 
                                 "nigromaculis" = "#FC9272",
                                 "Arctic spp." = "#08519C"),
                      name = element_blank(), # No legend title
                      breaks = c("sollicitans",
                                 "triseriatus",
                                 "albopictus", 
                                 "nigromaculis",
                                 "Arctic spp."),
                      labels = c("Ae. sollicitans", "Ae. triseriatus",
                                 "Ae. albopictus", "Ae. nigromaculis", 
                                 "Arctic spp.")) +
  
  theme_bw()

pLA.Tmax.lat

pLA.q.lat <- pLA.TPC.pars %>% 
  ggplot(aes(x = latitude)) +
  geom_vline(xintercept = 0, linetype = "dashed") +

  # data
  geom_errorbar(aes(ymin = q_2.5, ymax = q_97.5), width = 1) +
  geom_point(aes(y = q_50, colour = species)) +
  
  # fitted rma.mv model
  geom_line(data = pLA.q.pred, aes(y = pred)) +
  geom_ribbon(data = pLA.q.pred, aes(ymin = lower, ymax = upper), alpha = 0.2) + 
  
  labs(
    x = "absolute latitude",
    y = "q",
    title = expression(paste("Larval-to-Adult Survival (",italic(pLA),")"))
  ) +
  
  scale_colour_manual(values = c("sollicitans" = "#67000D",
                                 "triseriatus" = "#CB181D",
                                 "albopictus" = "#FB6A4A", 
                                 "nigromaculis" = "#FC9272",
                                 "Arctic spp." = "#08519C"),
                      name = element_blank(), # No legend title
                      breaks = c("sollicitans",
                                 "triseriatus",
                                 "albopictus", 
                                 "nigromaculis",
                                 "Arctic spp."),
                      labels = c("Ae. sollicitans", "Ae. triseriatus",
                                 "Ae. albopictus", "Ae. nigromaculis", 
                                 "Arctic spp.")) +
  
  theme_bw()

pLA.q.lat

TPC.params.lat <- plot_grid(pLA.Tmin.lat, pLA.Tmax.lat, pLA.q.lat, align = "v", ncol = 1)
TPC.params.lat

ggsave("figures/pLA.quad.TPC.params.lat.png", TPC.params.lat,
       width = 10.3, height = 10)


##### Mosquito development rate (MDR) #####

# Load R2jags model output
load("R-scripts/R2jags-objects/all-mods/MDR.nonarctic.quad.uni.Rdata") # Non-arctic
load("R-scripts/R2jags-objects/all-mods/MDR.arctic.quad.uni.Rdata") # Arctic

# Load data
data.MDR.all <- read_csv("data-processed/TraitData_MDR.csv")
data.MDR.all <- data.MDR.all %>%  # absolute latitude
  mutate(latitude = abs(latitude))

# Subset data
data.MDR.arctic <- subset(data.MDR.all, type == "Arctic") # Arctic species
data.MDR.nonarctic <- subset(data.MDR.all, type == "non-Arctic") # Non-Arctic species

###### Non-Arctic ######
# Get Tmin, Tmax, and q from each random effect
MDR.sims <- MDR.nonarctic.quad.uni$BUGSoutput$sims.list

# parameter values for each random effect
MDR.T0.nonarctic.fullpost <- sweep(MDR.sims$T0, 1, MDR.sims$cf.T0, "+")
MDR.Tm.nonarctic.fullpost <- sweep(MDR.sims$Tm, 1, MDR.sims$cf.Tm, "+")
MDR.q.nonarctic.fullpost <- sweep(MDR.sims$q,  1, MDR.sims$cf.q,  "+")

# summarize posterior distributions
MDR.T0.nonarctic <- t(apply(MDR.T0.nonarctic.fullpost, 2, quantile, c(0.025,0.5,0.975)))
MDR.Tm.nonarctic <- t(apply(MDR.Tm.nonarctic.fullpost, 2, quantile, c(0.025,0.5,0.975)))
MDR.q.nonarctic <- t(apply(MDR.q.nonarctic.fullpost, 2, quantile, c(0.025,0.5,0.975)))

# change column names
colnames(MDR.T0.nonarctic) <- c("T0_2.5", "T0_50", "T0_97.5")
colnames(MDR.Tm.nonarctic) <- c("Tm_2.5", "Tm_50", "Tm_97.5")
colnames(MDR.q.nonarctic)  <- c("q_2.5",  "q_50",  "q_97.5")


MDR.nonarctic.id.info <- data.MDR.nonarctic %>%
  group_by(species, citation) %>% 
  mutate(unique_id = cur_group_id()) %>% 
  group_by(unique_id) %>% 
  dplyr::select(unique_id, genus, species, citation, latitude, longitude, type) %>% 
  unique() %>% # drop duplicate
  arrange(unique_id)

MDR.TPC.pars.nonarctic <- bind_cols(MDR.nonarctic.id.info, MDR.T0.nonarctic, MDR.Tm.nonarctic, MDR.q.nonarctic)
MDR.TPC.pars.nonarctic

###### Arctic ######
MDR.T0.arctic.fullpost <- MDR.arctic.quad.uni$BUGSoutput$sims.list$cf.T0
MDR.Tm.arctic.fullpost <- MDR.arctic.quad.uni$BUGSoutput$sims.list$cf.Tm
MDR.q.arctic.fullpost <- MDR.arctic.quad.uni$BUGSoutput$sims.list$cf.q

# summarize posterior distributions
MDR.T0.arctic <- t(quantile(MDR.T0.arctic.fullpost, c(0.025, 0.5, 0.975)))
MDR.Tm.arctic <- t(quantile(MDR.Tm.arctic.fullpost, c(0.025, 0.5, 0.975)))
MDR.q.arctic <- t(quantile(MDR.q.arctic.fullpost, c(0.025, 0.5, 0.975)))

# change column names
colnames(MDR.T0.arctic) <- c("T0_2.5", "T0_50", "T0_97.5")
colnames(MDR.Tm.arctic) <- c("Tm_2.5", "Tm_50", "Tm_97.5")
colnames(MDR.q.arctic)  <- c("q_2.5",  "q_50",  "q_97.5")

# Since Arctic data do not have sufficient groups to include random effects, we will use the mean latitude
MDR.arctic.id.info <- data.frame(
  unique_id = 0,
  genus = "Aedes",
  species = "Arctic spp.",
  citation = NA,
  latitude = mean(data.MDR.arctic$latitude),
  longitude = NA,
  type = "Arctic"
  )


MDR.TPC.pars.arctic <- bind_cols(MDR.arctic.id.info, MDR.T0.arctic, MDR.Tm.arctic, MDR.q.arctic)
MDR.TPC.pars <- bind_rows(MDR.TPC.pars.arctic, MDR.TPC.pars.nonarctic)

###### TPC parameter-latitudinal analysis ######

# Tmin
## Calculate sampling variances
MDR.T0.vars <- bind_cols(MDR.T0.arctic.fullpost, MDR.T0.nonarctic.fullpost)
MDR.T0.vars <- apply(MDR.T0.vars, 2, var)
MDR.TPC.pars$T0_var <- MDR.T0.vars

## Regression
MDR.T0.fit <- rma.mv(yi = T0_50, 
                     V = T0_var,
                     mods = ~latitude,
                     data = MDR.TPC.pars,
                     method = "REML")

summary(MDR.T0.fit)


## Get predicted values
MDR.T0.pred <- data.frame(latitude = seq(min(MDR.TPC.pars$latitude),
                                         max(MDR.TPC.pars$latitude),
                                         length.out = 100))
MDR.T0.newdata <- predict(MDR.T0.fit,
                          newmods = MDR.T0.pred$latitude)

MDR.T0.pred$pred <- MDR.T0.newdata$pred
MDR.T0.pred$lower <- MDR.T0.newdata$ci.lb
MDR.T0.pred$upper <- MDR.T0.newdata$ci.ub


# Tmax
## Calculate sampling variances
MDR.Tm.vars <- bind_cols(MDR.Tm.arctic.fullpost, MDR.Tm.nonarctic.fullpost)
MDR.Tm.vars <- apply(MDR.Tm.vars, 2, var)
MDR.TPC.pars$Tm_var <- MDR.Tm.vars

## Regression
MDR.Tm.fit <- rma.mv(yi = Tm_50, 
                     V = Tm_var,
                     mods = ~latitude,
                     data = MDR.TPC.pars,
                     method = "REML")

summary(MDR.Tm.fit)

## Get predicted values
MDR.Tm.pred <- data.frame(latitude = seq(min(MDR.TPC.pars$latitude),
                                         max(MDR.TPC.pars$latitude),
                                         length.out = 100))
MDR.Tm.newdata <- predict(MDR.Tm.fit,
                          newmods = MDR.Tm.pred$latitude)

MDR.Tm.pred$pred <- MDR.Tm.newdata$pred
MDR.Tm.pred$lower <- MDR.Tm.newdata$ci.lb
MDR.Tm.pred$upper <- MDR.Tm.newdata$ci.ub


# q
## Calculate sampling variances
MDR.q.vars <- bind_cols(MDR.q.arctic.fullpost, MDR.q.nonarctic.fullpost)
MDR.q.vars <- apply(MDR.q.vars, 2, var)
MDR.TPC.pars$q_var <- MDR.q.vars

## Regression
MDR.q.fit <- rma.mv(yi = q_50, 
                    V = q_var,
                    mods = ~latitude,
                    data = MDR.TPC.pars,
                    method = "REML")

summary(MDR.q.fit)

## Get predicted values
MDR.q.pred <- data.frame(latitude = seq(min(MDR.TPC.pars$latitude),
                                        max(MDR.TPC.pars$latitude),
                                        length.out = 100))
MDR.q.newdata <- predict(MDR.q.fit,
                         newmods = MDR.q.pred$latitude)

MDR.q.pred$pred <- MDR.q.newdata$pred
MDR.q.pred$lower <- MDR.q.newdata$ci.lb
MDR.q.pred$upper <- MDR.q.newdata$ci.ub


# Plot
MDR.Tmin.lat <- MDR.TPC.pars %>% 
  ggplot(aes(x = latitude)) +
  geom_vline(xintercept = 0, linetype = "dashed") +

  # data
  geom_errorbar(aes( ymin = T0_2.5, ymax = T0_97.5), width = 1) +
  geom_point(aes(y = T0_50, colour = species)) +
  
  # fitted rma.mv model
  geom_line(data = MDR.T0.pred, aes(y = pred)) +
  geom_ribbon(data = MDR.T0.pred, aes(ymin = lower, ymax = upper), alpha = 0.2) + 
  
  labs(
    x = "absolute latitude",
    y = "Tmin",
    title = expression(paste("Mosquito Development Rate (",italic(MDR),")"))
  ) +
  
  scale_colour_manual(values = c("albopictus" = "#67000D",
                                 "sollicitans" = "#CB181D",
                                 "triseriatus" = "#FB6A4A", 
                                 "nigromaculis" = "#FC9272",
                                 "Arctic spp." = "#08519C"),
                      name = element_blank(), # No legend title
                      breaks = c("albopictus", 
                                 "sollicitans",
                                 "triseriatus",
                                 "nigromaculis",
                                 "Arctic spp."),
                      labels = c("Ae. albopictus", 
                                 "Ae. sollicitans", 
                                 "Ae. triseriatus",
                                 "Ae. nigromaculis", 
                                 "Arctic spp.")) +
  
  theme_bw()

MDR.Tmin.lat

MDR.Tmax.lat <- MDR.TPC.pars %>% 
  ggplot(aes(x = latitude)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  # data
  geom_errorbar(aes(ymin = Tm_2.5, ymax = Tm_97.5), width = 1) +
  geom_point(aes(y = Tm_50, colour = species)) +
  
  # fitted rma.mv model
  geom_line(data = MDR.Tm.pred, aes(y = pred)) +
  geom_ribbon(data = MDR.Tm.pred, aes(ymin = lower, ymax = upper), alpha = 0.2) + 
  
  labs(
    x = "Absolute latitude",
    y = "Tmax",
    title = expression(paste("Mosquito Development Rate (",italic(MDR),")"))
  ) +
  
  scale_colour_manual(values = c("albopictus" = "#67000D",
                                 "sollicitans" = "#CB181D",
                                 "triseriatus" = "#FB6A4A", 
                                 "nigromaculis" = "#FC9272",
                                 "Arctic spp." = "#08519C"),
                      name = element_blank(), # No legend title
                      breaks = c("albopictus", 
                                 "sollicitans",
                                 "triseriatus",
                                 "nigromaculis",
                                 "Arctic spp."),
                      labels = c("Ae. albopictus", 
                                 "Ae. sollicitans", 
                                 "Ae. triseriatus",
                                 "Ae. nigromaculis", 
                                 "Arctic spp.")) +
  
  theme_bw()

MDR.Tmax.lat

MDR.q.lat <- MDR.TPC.pars %>% 
  ggplot(aes(x = latitude)) +
  geom_vline(xintercept = 0, linetype = "dashed") +

  # data
  geom_errorbar(aes(ymin = q_2.5, ymax = q_97.5), width = 1) +
  geom_point(aes(y = q_50, colour = species)) +
  
  # fitted rma.mv model
  geom_line(data = MDR.q.pred, aes(y = pred)) +
  geom_ribbon(data = MDR.q.pred, aes(ymin = lower, ymax = upper), alpha = 0.2) + 
  
  labs(
    x = "absolute latitude",
    y = "q",
    title = expression(paste("Mosquito Development Rate (",italic(MDR),")"))
  ) +
  
  scale_colour_manual(values = c("albopictus" = "#67000D",
                                 "sollicitans" = "#CB181D",
                                 "triseriatus" = "#FB6A4A", 
                                 "nigromaculis" = "#FC9272",
                                 "Arctic spp." = "#08519C"),
                      name = element_blank(), # No legend title
                      breaks = c("albopictus", 
                                 "sollicitans",
                                 "triseriatus",
                                 "nigromaculis",
                                 "Arctic spp."),
                      labels = c("Ae. albopictus", 
                                 "Ae. sollicitans", 
                                 "Ae. triseriatus",
                                 "Ae. nigromaculis", 
                                 "Arctic spp.")) +
  
  theme_bw()

MDR.q.lat

TPC.params.lat <- plot_grid(MDR.Tmin.lat, MDR.Tmax.lat, MDR.q.lat, align = "v", ncol = 1)
TPC.params.lat

ggsave("figures/MDR.quad.TPC.params.lat.png", TPC.params.lat,
       width = 10.3, height = 10)


##### Pathogen development rate (PDR) #####

# Load R2jags model output
load("R-scripts/R2jags-objects/all-mods/PDR.nonarctic.bri.uni.Rdata") # Non-arctic
load("R-scripts/R2jags-objects/all-mods/PDR.arctic.bri.uni.Rdata") # Arctic

# Load data
data.PDR.all <- read_csv("data-processed/TraitData_PDR.csv")
data.PDR.all <- data.PDR.all %>%# absolute latitude
  mutate(latitude = abs(latitude))

# Subset data
data.PDR.arctic <- subset(data.PDR.all, type == "Arctic") # Arctic species
data.PDR.nonarctic <- subset(data.PDR.all, type == "non-Arctic") # Non-Arctic species

###### Non-Arctic ######
# Get Tmin, Tmax, and q from each random effect
PDR.sims <- PDR.nonarctic.bri.uni$BUGSoutput$sims.list

# parameter values for each unique_id
PDR.T0.nonarctic.fullpost <- sweep(PDR.sims$T0, 1, PDR.sims$cf.T0, "+")
PDR.Tm.nonarctic.fullpost <- sweep(PDR.sims$Tm, 1, PDR.sims$cf.Tm, "+")
PDR.q.nonarctic.fullpost  <- sweep(PDR.sims$q,  1, PDR.sims$cf.q,  "+")

# summarize posterior distributions
PDR.T0.nonarctic <- t(apply(PDR.T0.nonarctic.fullpost, 2, quantile, c(0.025,0.5,0.975)))
PDR.Tm.nonarctic <- t(apply(PDR.Tm.nonarctic.fullpost, 2, quantile, c(0.025,0.5,0.975)))
PDR.q.nonarctic <- t(apply(PDR.q.nonarctic.fullpost, 2, quantile, c(0.025,0.5,0.975)))

# change column names
colnames(PDR.T0.nonarctic) <- c("T0_2.5", "T0_50", "T0_97.5")
colnames(PDR.Tm.nonarctic) <- c("Tm_2.5", "Tm_50", "Tm_97.5")
colnames(PDR.q.nonarctic)  <- c("q_2.5",  "q_50",  "q_97.5")


PDR.nonarctic.id.info <- data.PDR.nonarctic %>%
  group_by(paras_species, host_species, citation) %>% 
  mutate(unique_id = cur_group_id()) %>% 
  group_by(unique_id) %>% 
  dplyr::select(unique_id, paras_genus, paras_species, host_genus, host_species,
                citation, latitude, longitude, type) %>% 
  unique() %>% # drop duplicate
  arrange(unique_id)

PDR.TPC.pars.nonarctic <- bind_cols(PDR.nonarctic.id.info, PDR.T0.nonarctic, PDR.Tm.nonarctic, PDR.q.nonarctic)
PDR.TPC.pars.nonarctic


###### Arctic ######
PDR.T0.arctic.fullpost <- PDR.arctic.bri.uni$BUGSoutput$sims.list$cf.T0
PDR.Tm.arctic.fullpost <- PDR.arctic.bri.uni$BUGSoutput$sims.list$cf.Tm
PDR.q.arctic.fullpost <- PDR.arctic.bri.uni$BUGSoutput$sims.list$cf.q

# summarize posterior distributions
PDR.T0.arctic <- t(quantile(PDR.T0.arctic.fullpost, c(0.025, 0.5, 0.975)))
PDR.Tm.arctic <- t(quantile(PDR.Tm.arctic.fullpost, c(0.025, 0.5, 0.975)))
PDR.q.arctic <- t(quantile(PDR.q.arctic.fullpost, c(0.025, 0.5, 0.975)))

# change column names
colnames(PDR.T0.arctic) <- c("T0_2.5", "T0_50", "T0_97.5")
colnames(PDR.Tm.arctic) <- c("Tm_2.5", "Tm_50", "Tm_97.5")
colnames(PDR.q.arctic)  <- c("q_2.5",  "q_50",  "q_97.5")


# Since Arctic data do not have sufficient groups to include random effects, we will use the mean latitude
PDR.arctic.id.info <- data.frame(
  unique_id = 0,
  paras_genus = NA,
  paras_species = "Arctic spp.",
  host_genus = "Aedes",
  host_species = "Arctic spp.",
  citation = NA,
  latitude = mean(data.PDR.arctic$latitude),
  longitude = NA,
  type = "Arctic"
)


PDR.TPC.pars.arctic <- bind_cols(PDR.arctic.id.info, PDR.T0.arctic, PDR.Tm.arctic, PDR.q.arctic)
PDR.TPC.pars <- bind_rows(PDR.TPC.pars.arctic, PDR.TPC.pars.nonarctic)


###### TPC parameter-latitudinal analysis ######

# Tmin
## Calculate sampling variances
PDR.T0.vars <- bind_cols(PDR.T0.arctic.fullpost, PDR.T0.nonarctic.fullpost)
PDR.T0.vars <- apply(PDR.T0.vars, 2, var)
PDR.TPC.pars$T0_var <- PDR.T0.vars

## Regression
PDR.T0.fit <- rma.mv(yi = T0_50, 
                    V = T0_var,
                    mods = ~latitude,
                    data = PDR.TPC.pars,
                    method = "REML")

summary(PDR.T0.fit)


## Get predicted values
PDR.T0.pred <- data.frame(latitude = seq(min(PDR.TPC.pars$latitude),
                                        max(PDR.TPC.pars$latitude),
                                        length.out = 100))
PDR.T0.newdata <- predict(PDR.T0.fit,
                         newmods = PDR.T0.pred$latitude)

PDR.T0.pred$pred <- PDR.T0.newdata$pred
PDR.T0.pred$lower <- PDR.T0.newdata$ci.lb
PDR.T0.pred$upper <- PDR.T0.newdata$ci.ub


# Tmax
## Calculate sampling variances
PDR.Tm.vars <- bind_cols(PDR.Tm.arctic.fullpost, PDR.Tm.nonarctic.fullpost)
PDR.Tm.vars <- apply(PDR.Tm.vars, 2, var)
PDR.TPC.pars$Tm_var <- PDR.Tm.vars

## Regression
PDR.Tm.fit <- rma.mv(yi = Tm_50, 
                     V = Tm_var,
                     mods = ~latitude,
                     data = PDR.TPC.pars,
                     method = "REML")

summary(PDR.Tm.fit)

## Get predicted values
PDR.Tm.pred <- data.frame(latitude = seq(min(PDR.TPC.pars$latitude),
                                        max(PDR.TPC.pars$latitude),
                                        length.out = 100))
PDR.Tm.newdata <- predict(PDR.Tm.fit,
                         newmods = PDR.Tm.pred$latitude)

PDR.Tm.pred$pred <- PDR.Tm.newdata$pred
PDR.Tm.pred$lower <- PDR.Tm.newdata$ci.lb
PDR.Tm.pred$upper <- PDR.Tm.newdata$ci.ub


# q
## Calculate sampling variances
PDR.q.vars <- bind_cols(PDR.q.arctic.fullpost, PDR.q.nonarctic.fullpost)
PDR.q.vars <- apply(PDR.q.vars, 2, var)
PDR.TPC.pars$q_var <- PDR.q.vars

## Regression
PDR.q.fit <- rma.mv(yi = q_50, 
                    V = q_var,
                    mods = ~latitude,
                    data = PDR.TPC.pars,
                    method = "REML")

summary(PDR.q.fit)

## Get predicted values
PDR.q.pred <- data.frame(latitude = seq(min(PDR.TPC.pars$latitude),
                                       max(PDR.TPC.pars$latitude),
                                       length.out = 100))
PDR.q.newdata <- predict(PDR.q.fit,
                        newmods = PDR.q.pred$latitude)

PDR.q.pred$pred <- PDR.q.newdata$pred
PDR.q.pred$lower <- PDR.q.newdata$ci.lb
PDR.q.pred$upper <- PDR.q.newdata$ci.ub


# Plot
PDR.Tmin.lat <- PDR.TPC.pars %>% 
  ggplot(aes(x = latitude)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  # data
  geom_errorbar(aes(ymin = T0_2.5, ymax = T0_97.5), width = 1) +
  geom_point(aes(y = T0_50, colour = paras_species)) +
  
  # fitted rma.mv model
  geom_line(data = PDR.T0.pred, aes(y = pred)) +
  geom_ribbon(data = PDR.T0.pred, aes(ymin = lower, ymax = upper), alpha = 0.2) + 
  
  
  labs(
    x = "Absolute latitude",
    y = "Tmin",
    title = expression(paste("Pathogen Development Rate (",italic(PDR),")"))
  ) +
  
  scale_colour_manual(values = c("bancrofti" = "#67000D",
                                 "immitis" = "#FB6A4A",
                                 "Arctic spp." = "#08519C"),
                      name = element_blank(), # No legend title
                      breaks = c("bancrofti",
                                 "immitis",
                                 "Arctic spp."),
                      labels = c("Wuchereria bancrofti",
                                 "Dirofilaria immitis",
                                 "Arctic spp.")) +
  
  theme_bw()

PDR.Tmin.lat

PDR.Tmax.lat <- PDR.TPC.pars %>% 
  ggplot(aes(x = latitude)) +
  geom_vline(xintercept = 0, linetype = "dashed") +

  # data
  geom_errorbar(aes(ymin = Tm_2.5, ymax = Tm_97.5), width = 1) +
  geom_point(aes(y = Tm_50, colour = paras_species)) +
  
  # fitted rma.mv model
  geom_line(data = PDR.Tm.pred, aes(y = pred)) +
  geom_ribbon(data = PDR.Tm.pred, aes(ymin = lower, ymax = upper), alpha = 0.2) + 
  
  labs(
    x = "Absolute latitude",
    y = "Tmax",
    title = expression(paste("Pathogen Development Rate (",italic(PDR),")"))
  ) +
  
  scale_colour_manual(values = c("bancrofti" = "#67000D",
                                 "immitis" = "#FB6A4A",
                                 "Arctic spp." = "#08519C"),
                      name = element_blank(), # No legend title
                      breaks = c("bancrofti",
                                 "immitis",
                                 "Arctic spp."),
                      labels = c("Wuchereria bancrofti",
                                 "Dirofilaria immitis",
                                 "Arctic spp.")) +
  
  theme_bw()

PDR.Tmax.lat

PDR.q.lat <- PDR.TPC.pars %>% 
  ggplot(aes(x = latitude)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  # data
  geom_errorbar(aes(ymin = q_2.5, ymax = q_97.5), width = 1) +
  geom_point(aes(y = q_50, colour = paras_species)) +
  
  # fitted rma.mv model
  geom_line(data = PDR.q.pred, aes(y = pred)) +
  geom_ribbon(data = PDR.q.pred, aes(ymin = lower, ymax = upper), alpha = 0.2) + 
  labs(
    x = "Absolute latitude",
    y = "q",
    title = expression(paste("Pathogen Development Rate (",italic(PDR),")"))
  ) +
  
  scale_colour_manual(values = c("bancrofti" = "#67000D",
                                 "immitis" = "#FB6A4A",
                                 "Arctic spp." = "#08519C"),
                      name = element_blank(), # No legend title
                      breaks = c("bancrofti",
                                 "immitis",
                                 "Arctic spp."),
                      labels = c("Wuchereria bancrofti",
                                 "Dirofilaria immitis",
                                 "Arctic spp.")) +
  
  theme_bw()

PDR.q.lat

TPC.params.lat <- plot_grid(PDR.Tmin.lat, PDR.Tmax.lat, PDR.q.lat, align = "v", ncol = 1)
TPC.params.lat

ggsave("figures/PDR.bri.TPC.params.lat.png", TPC.params.lat,
       width = 10.3, height = 10)



##### Eggs per female per gonotrophic cycle (EFGC) #####

# Load R2jags model output
load("R-scripts/R2jags-objects/all-mods/EFGC.alldata.quad.uni.Rdata") 

# Load data
data.EFGC.all <- read_csv("data-processed/TraitData_EFGC.csv")
data.EFGC.all <- data.EFGC.all %>%  # absolute latitude
  mutate(latitude = abs(latitude))


# Get Tmin, Tmax, and q from each random effect
EFGC.sims <- EFGC.alldata.quad.uni$BUGSoutput$sims.list

# parameter values for each random effect
EFGC.T0.alldata.fullpost <- sweep(EFGC.sims$T0, 1, EFGC.sims$cf.T0, "+")
EFGC.Tm.alldata.fullpost <- sweep(EFGC.sims$Tm, 1, EFGC.sims$cf.Tm, "+")
EFGC.q.alldata.fullpost  <- sweep(EFGC.sims$q,  1, EFGC.sims$cf.q,  "+")

# summarize posterior distributions
EFGC.T0.alldata <- t(apply(EFGC.T0.alldata.fullpost, 2, quantile, c(0.025,0.5,0.975)))
EFGC.Tm.alldata <- t(apply(EFGC.Tm.alldata.fullpost, 2, quantile, c(0.025,0.5,0.975)))
EFGC.q.alldata <- t(apply(EFGC.q.alldata.fullpost, 2, quantile, c(0.025,0.5,0.975)))

# change column names
colnames(EFGC.T0.alldata) <- c("T0_2.5", "T0_50", "T0_97.5")
colnames(EFGC.Tm.alldata) <- c("Tm_2.5", "Tm_50", "Tm_97.5")
colnames(EFGC.q.alldata)  <- c("q_2.5",  "q_50",  "q_97.5")


EFGC.alldata.id.info <- data.EFGC.all %>%
  group_by(species, citation) %>% 
  mutate(unique_id = cur_group_id()) %>% 
  group_by(unique_id) %>% 
  dplyr::select(unique_id, genus, species, citation, latitude, longitude, type) %>% 
  unique() %>% # drop duplicate
  arrange(unique_id)


EFGC.TPC.pars<- bind_cols(EFGC.alldata.id.info, EFGC.T0.alldata, EFGC.Tm.alldata, EFGC.q.alldata)
EFGC.TPC.pars


###### TPC parameter-latitudinal analysis ######

# Tmin
## Calculate sampling variances
EFGC.T0.vars <- EFGC.T0.alldata.fullpost
EFGC.T0.vars <- apply(EFGC.T0.vars, 2, var)
EFGC.TPC.pars$T0_var <- EFGC.T0.vars

## Regression
EFGC.T0.fit <- rma.mv(yi = T0_50, 
                    V = T0_var,
                    mods = ~latitude,
                    data = EFGC.TPC.pars,
                    method = "REML")

summary(EFGC.T0.fit)


## Get predicted values
EFGC.T0.pred <- data.frame(latitude = seq(min(EFGC.TPC.pars$latitude),
                                        max(EFGC.TPC.pars$latitude),
                                        length.out = 100))
EFGC.T0.newdata <- predict(EFGC.T0.fit,
                         newmods = EFGC.T0.pred$latitude)

EFGC.T0.pred$pred <- EFGC.T0.newdata$pred
EFGC.T0.pred$lower <- EFGC.T0.newdata$ci.lb
EFGC.T0.pred$upper <- EFGC.T0.newdata$ci.ub


# Tmax
## Calculate sampling variances
EFGC.Tm.vars <- EFGC.Tm.alldata.fullpost
EFGC.Tm.vars <- apply(EFGC.Tm.vars, 2, var)
EFGC.TPC.pars$Tm_var <- EFGC.Tm.vars

## Regression
EFGC.Tm.fit <- rma.mv(yi = Tm_50, 
                    V = Tm_var,
                    mods = ~latitude,
                    data = EFGC.TPC.pars,
                    method = "REML")

summary(EFGC.Tm.fit)

## Get predicted values
EFGC.Tm.pred <- data.frame(latitude = seq(min(EFGC.TPC.pars$latitude),
                                        max(EFGC.TPC.pars$latitude),
                                        length.out = 100))
EFGC.Tm.newdata <- predict(EFGC.Tm.fit,
                         newmods = EFGC.Tm.pred$latitude)

EFGC.Tm.pred$pred <- EFGC.Tm.newdata$pred
EFGC.Tm.pred$lower <- EFGC.Tm.newdata$ci.lb
EFGC.Tm.pred$upper <- EFGC.Tm.newdata$ci.ub


# q
## Calculate sampling variances
EFGC.q.vars <- EFGC.q.alldata.fullpost
EFGC.q.vars <- apply(EFGC.q.vars, 2, var)
EFGC.TPC.pars$q_var <- EFGC.q.vars

## Regression
EFGC.q.fit <- rma.mv(yi = q_50, 
                   V = q_var,
                   mods = ~latitude,
                   data = EFGC.TPC.pars,
                   method = "REML")

summary(EFGC.q.fit)

## Get predicted values
EFGC.q.pred <- data.frame(latitude = seq(min(EFGC.TPC.pars$latitude),
                                       max(EFGC.TPC.pars$latitude),
                                       length.out = 100))
EFGC.q.newdata <- predict(EFGC.q.fit,
                        newmods = EFGC.q.pred$latitude)

EFGC.q.pred$pred <- EFGC.q.newdata$pred
EFGC.q.pred$lower <- EFGC.q.newdata$ci.lb
EFGC.q.pred$upper <- EFGC.q.newdata$ci.ub


# Plot
EFGC.Tmin.lat <- EFGC.TPC.pars %>% 
  ggplot(aes(x = latitude)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  # data
  geom_errorbar(aes(ymin = T0_2.5, ymax = T0_97.5), width = 1) +
  geom_point(aes(y = T0_50, colour = species)) +
  
  # fitted rma.mv model
  geom_line(data = EFGC.T0.pred, aes(y = pred)) +
  geom_ribbon(data = EFGC.T0.pred, aes(ymin = lower, ymax = upper), alpha = 0.2) + 
  
  labs(
    x = "absolute latitude",
    y = "Tmin",
    title = expression(paste("Eggs per Female \nper Gonotrophic Cycle (",italic(EFGC),")"))
  ) +
  
  scale_colour_manual(values = c("albopictus" = "#67000D",
                                 "hexodontus" = "#9ECAE1",
                                 "punctor" = "#4292C6",
                                 "cinereus" = "#2171B5",
                                 "communis" = "#08519C",
                                 "impiger" = "#08306B"),
                      name = element_blank(), # No legend title
                      breaks = c("albopictus",
                                 "hexodontus",
                                 "punctor",
                                 "cinereus",
                                 "communis",
                                 "impiger"),
                      labels = c("Ae. albopictus",
                                 "Ae. hexodontus",
                                 "Ae. punctor",
                                 "Ae. cinereus",
                                 "Ae. communis",
                                 "Ae. impiger")) +
  
  theme_bw()

EFGC.Tmin.lat


EFGC.Tmax.lat <- EFGC.TPC.pars %>% 
  ggplot(aes(x = latitude)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  # data
  geom_errorbar(aes(ymin = Tm_2.5, ymax = Tm_97.5), width = 1) +
  geom_point(aes(y = Tm_50, colour = species)) +
  
  # fitted rma.mv model
  geom_line(data = EFGC.Tm.pred, aes(y = pred)) +
  geom_ribbon(data = EFGC.Tm.pred, aes(ymin = lower, ymax = upper), alpha = 0.2) + 
  
  labs(
    x = "absolute latitude",
    y = "Tmax",
    title = expression(paste("Eggs per Female \nper Gonotrophic Cycle (",italic(EFGC),")"))
  ) +
  
  scale_colour_manual(values = c("albopictus" = "#67000D",
                                 "hexodontus" = "#9ECAE1",
                                 "punctor" = "#4292C6",
                                 "cinereus" = "#2171B5",
                                 "communis" = "#08519C",
                                 "impiger" = "#08306B"),
                      name = element_blank(), # No legend title
                      breaks = c("albopictus",
                                 "hexodontus",
                                 "punctor",
                                 "cinereus",
                                 "communis",
                                 "impiger"),
                      labels = c("Ae. albopictus",
                                 "Ae. hexodontus",
                                 "Ae. punctor",
                                 "Ae. cinereus",
                                 "Ae. communis",
                                 "Ae. impiger")) +
  
  theme_bw()

EFGC.Tmax.lat

EFGC.q.lat <- EFGC.TPC.pars %>% 
  ggplot(aes(x = latitude)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  # data
  geom_errorbar(aes(ymin = q_2.5, ymax = q_97.5), width = 1) +
  geom_point(aes(y = q_50, colour = species)) +
  
  # fitted rma.mv model
  geom_line(data = EFGC.q.pred, aes(y = pred)) +
  geom_ribbon(data = EFGC.q.pred, aes(ymin = lower, ymax = upper), alpha = 0.2) + 
  
  labs(
    x = "absolute latitude",
    y = "q",
    title = expression(paste("Eggs per Female \nper Gonotrophic Cycle (",italic(EFGC),")"))
  ) +
  
  scale_colour_manual(values = c("albopictus" = "#67000D",
                                 "hexodontus" = "#9ECAE1",
                                 "punctor" = "#4292C6",
                                 "cinereus" = "#2171B5",
                                 "communis" = "#08519C",
                                 "impiger" = "#08306B"),
                      name = element_blank(), # No legend title
                      breaks = c("albopictus",
                                 "hexodontus",
                                 "punctor",
                                 "cinereus",
                                 "communis",
                                 "impiger"),
                      labels = c("Ae. albopictus",
                                 "Ae. hexodontus",
                                 "Ae. punctor",
                                 "Ae. cinereus",
                                 "Ae. communis",
                                 "Ae. impiger")) +
  
  theme_bw()

EFGC.q.lat

TPC.params.lat <- plot_grid(EFGC.Tmin.lat, EFGC.Tmax.lat, EFGC.q.lat, align = "v", ncol = 1)
TPC.params.lat

ggsave("figures/EFGC.quad.TPC.params.lat.png", TPC.params.lat,
       width = 10.3, height = 10)



##### biting rate (a) #####

# Load R2jags model output
load("R-scripts/R2jags-objects/all-mods/a.alldata.bri.uni.Rdata") 

# Load data
data.a.all <- read_csv("data-processed/TraitData_a.csv")
data.a.all <- data.a.all %>%  # absolute latitude
  mutate(latitude = abs(latitude))


# Get Tmin, Tmax, and q from each random effect
a.sims <- a.alldata.bri.uni$BUGSoutput$sims.list

# parameter values for each random effect
a.T0.alldata.fullpost <- sweep(a.sims$T0, 1, a.sims$cf.T0, "+")
a.Tm.alldata.fullpost <- sweep(a.sims$Tm, 1, a.sims$cf.Tm, "+")
a.q.alldata.fullpost  <- sweep(a.sims$q,  1, a.sims$cf.q,  "+")

# summarize posterior distributions
a.T0.alldata <- t(apply(a.T0.alldata.fullpost, 2, quantile, c(0.025,0.5,0.975)))
a.Tm.alldata <- t(apply(a.Tm.alldata.fullpost, 2, quantile, c(0.025,0.5,0.975)))
a.q.alldata <- t(apply(a.q.alldata.fullpost, 2, quantile, c(0.025,0.5,0.975)))

# change column names
colnames(a.T0.alldata) <- c("T0_2.5", "T0_50", "T0_97.5")
colnames(a.Tm.alldata) <- c("Tm_2.5", "Tm_50", "Tm_97.5")
colnames(a.q.alldata)  <- c("q_2.5",  "q_50",  "q_97.5")


a.alldata.id.info <- data.a.all %>%
  group_by(species, citation) %>% 
  mutate(unique_id = cur_group_id()) %>% 
  group_by(unique_id) %>% 
  dplyr::select(unique_id, genus, species, citation, latitude, longitude, type) %>% 
  unique() %>% # drop duplicate
  arrange(unique_id)


a.TPC.pars<- bind_cols(a.alldata.id.info, a.T0.alldata, a.Tm.alldata, a.q.alldata)
a.TPC.pars


###### TPC parameter-latitudinal analysis ######

# Tmin
## Calculate sampling variances
a.T0.vars <- a.T0.alldata.fullpost
a.T0.vars <- apply(a.T0.vars, 2, var)
a.TPC.pars$T0_var <- a.T0.vars

## Regression
a.T0.fit <- rma.mv(yi = T0_50, 
                      V = T0_var,
                      mods = ~latitude,
                      data = a.TPC.pars,
                      method = "REML")

summary(a.T0.fit)


## Get predicted values
a.T0.pred <- data.frame(latitude = seq(min(a.TPC.pars$latitude),
                                          max(a.TPC.pars$latitude),
                                          length.out = 100))
a.T0.newdata <- predict(a.T0.fit,
                           newmods = a.T0.pred$latitude)

a.T0.pred$pred <- a.T0.newdata$pred
a.T0.pred$lower <- a.T0.newdata$ci.lb
a.T0.pred$upper <- a.T0.newdata$ci.ub


# Tmax
## Calculate sampling variances
a.Tm.vars <- a.Tm.alldata.fullpost
a.Tm.vars <- apply(a.Tm.vars, 2, var)
a.TPC.pars$Tm_var <- a.Tm.vars

## Regression
a.Tm.fit <- rma.mv(yi = Tm_50, 
                      V = Tm_var,
                      mods = ~latitude,
                      data = a.TPC.pars,
                      method = "REML")

summary(a.Tm.fit)

## Get predicted values
a.Tm.pred <- data.frame(latitude = seq(min(a.TPC.pars$latitude),
                                          max(a.TPC.pars$latitude),
                                          length.out = 100))
a.Tm.newdata <- predict(a.Tm.fit,
                           newmods = a.Tm.pred$latitude)

a.Tm.pred$pred <- a.Tm.newdata$pred
a.Tm.pred$lower <- a.Tm.newdata$ci.lb
a.Tm.pred$upper <- a.Tm.newdata$ci.ub


# q
## Calculate sampling variances
a.q.vars <- a.q.alldata.fullpost
a.q.vars <- apply(a.q.vars, 2, var)
a.TPC.pars$q_var <- a.q.vars

## Regression
a.q.fit <- rma.mv(yi = q_50, 
                     V = q_var,
                     mods = ~latitude,
                     data = a.TPC.pars,
                     method = "REML")

summary(a.q.fit)

## Get predicted values
a.q.pred <- data.frame(latitude = seq(min(a.TPC.pars$latitude),
                                         max(a.TPC.pars$latitude),
                                         length.out = 100))
a.q.newdata <- predict(a.q.fit,
                          newmods = a.q.pred$latitude)

a.q.pred$pred <- a.q.newdata$pred
a.q.pred$lower <- a.q.newdata$ci.lb
a.q.pred$upper <- a.q.newdata$ci.ub


# Plot
a.Tmin.lat <- a.TPC.pars %>% 
  ggplot(aes(x = latitude)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  # data
  geom_errorbar(aes(ymin = T0_2.5, ymax = T0_97.5), width = 1) +
  geom_point(aes(y = T0_50, colour = species)) +
  
  # fitted rma.mv model
  geom_line(data = a.T0.pred, aes(y = pred)) +
  geom_ribbon(data = a.T0.pred, aes(ymin = lower, ymax = upper), alpha = 0.2) + 
  
  labs(
    x = "absolute latitude",
    y = "Tmin",
    title = expression(paste("Biting Rate (",italic(a),")"))
  ) +
  
  scale_colour_manual(values = c("albopictus" = "#67000D",
                                 "punctor" = "#4292C6",
                                 "cinereus" = "#2171B5",
                                 "communis" = "#08519C",
                                 "impiger" = "#08306B"),
                      name = element_blank(), # No legend title
                      breaks = c("albopictus",
                                 "punctor",
                                 "cinereus",
                                 "communis",
                                 "impiger"),
                      labels = c("Ae. albopictus",
                                 "Ae. punctor",
                                 "Ae. cinereus",
                                 "Ae. communis",
                                 "Ae. impiger")) +
  
  theme_bw()

a.Tmin.lat


a.Tmax.lat <- a.TPC.pars %>% 
  ggplot(aes(x = latitude)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  # data
  geom_errorbar(aes(ymin = Tm_2.5, ymax = Tm_97.5), width = 1) +
  geom_point(aes(y = Tm_50, colour = species)) +
  
  # fitted rma.mv model
  geom_line(data = a.Tm.pred, aes(y = pred)) +
  geom_ribbon(data = a.Tm.pred, aes(ymin = lower, ymax = upper), alpha = 0.2) + 
  
  labs(
    x = "absolute latitude",
    y = "Tmax",
    title = expression(paste("Biting Rate (",italic(a),")"))
  ) +
  
  scale_colour_manual(values = c("albopictus" = "#67000D",
                                 "punctor" = "#4292C6",
                                 "cinereus" = "#2171B5",
                                 "communis" = "#08519C",
                                 "impiger" = "#08306B"),
                      name = element_blank(), # No legend title
                      breaks = c("albopictus",
                                 "punctor",
                                 "cinereus",
                                 "communis",
                                 "impiger"),
                      labels = c("Ae. albopictus",
                                 "Ae. punctor",
                                 "Ae. cinereus",
                                 "Ae. communis",
                                 "Ae. impiger")) +
  
  theme_bw()

a.Tmax.lat

a.q.lat <- a.TPC.pars %>% 
  ggplot(aes(x = latitude)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  # data
  geom_errorbar(aes(ymin = q_2.5, ymax = q_97.5), width = 1) +
  geom_point(aes(y = q_50, colour = species)) +
  
  # fitted rma.mv model
  geom_line(data = a.q.pred, aes(y = pred)) +
  geom_ribbon(data = a.q.pred, aes(ymin = lower, ymax = upper), alpha = 0.2) + 
  
  labs(
    x = "absolute latitude",
    y = "q",
    title = expression(paste("Biting Rate (",italic(a),")"))
  ) +
  
  scale_colour_manual(values = c("albopictus" = "#67000D",
                                 "punctor" = "#4292C6",
                                 "cinereus" = "#2171B5",
                                 "communis" = "#08519C",
                                 "impiger" = "#08306B"),
                      name = element_blank(), # No legend title
                      breaks = c("albopictus",
                                 "punctor",
                                 "cinereus",
                                 "communis",
                                 "impiger"),
                      labels = c("Ae. albopictus",
                                 "Ae. punctor",
                                 "Ae. cinereus",
                                 "Ae. communis",
                                 "Ae. impiger")) +
  
  theme_bw()

a.q.lat

TPC.params.lat <- plot_grid(a.Tmin.lat, a.Tmax.lat, a.q.lat, align = "v", ncol = 1)
TPC.params.lat

ggsave("figures/a.bri.TPC.params.lat.png", TPC.params.lat,
       width = 10.3, height = 10)



##### Plot everything together #####
Tmin.lat <- plot_grid(pLA.Tmin.lat, MDR.Tmin.lat, lf.Tmin.lat, 
                      PDR.Tmin.lat, EV.Tmin.lat, EFGC.Tmin.lat, a.Tmin.lat,
                      align = "v", ncol = 2, labels = "AUTO")
Tmin.lat

ggsave("figures/FigS7-lat.Tmin.png", Tmin.lat,
       width = 12, height = 6)


Tmax.lat <- plot_grid(pLA.Tmax.lat, MDR.Tmax.lat, lf.Tmax.lat, 
                      PDR.Tmax.lat, EV.Tmax.lat, EFGC.Tmax.lat, a.Tmax.lat,
                      align = "v", ncol = 2, labels = "AUTO")
Tmax.lat

ggsave("figures/FigS8-lat.Tmax.png", Tmax.lat,
       width = 12, height = 6)

q.lat <- plot_grid(pLA.q.lat, MDR.q.lat, lf.q.lat, 
                   PDR.q.lat, EV.q.lat, EFGC.q.lat, a.q.lat,
                   align = "v", ncol = 2, labels = "AUTO")
q.lat

ggsave("figures/FigS9-lat.q.png", q.lat,
       width = 12, height = 6)


# 2. Summary table -------------------------------------------------------------

