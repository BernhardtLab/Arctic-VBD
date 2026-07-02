## Lilian Chan, University of Guelph
## Arctic vector-borne disease transmission suitability model
##
## Purpose: Investigate the relationship between thermal limits (Tmin, Tmax) and 
## latitude
## 
## Table of content:
##    0. Set-up workspace
##    1. Load data and model output


# 0. Set-up workspace -----------------------------------------------------
library(tidyverse)
library(readxl)
library(janitor)
library(R2jags)
library(MASS)
library(cowplot)
library(RColorBrewer)
library(metafor)


#  1. Load data and model output -----------------------------------------------

# ##### Biting rate (a) #####
# 
# # Load R2jags model output
# load("R-scripts/R2jags-objects/all-mods/a.alldata.bri.uni.Rdata")
# 
# # Load data
# data.a.all <- read_csv("data-processed/TraitData_a.csv")
# 
# # Get Tmin, Tmax, and q from each random effect
# a.sims <- a.alldata.bri.uni$BUGSoutput$sims.list
# 
# # parameter values for each unique_id
# a.T0.alldata <- sweep(a.sims$T0, 1, a.sims$cf.T0, "+")
# a.Tm.alldata <- sweep(a.sims$Tm, 1, a.sims$cf.Tm, "+")
# a.q.alldata  <- sweep(a.sims$q,  1, a.sims$cf.q,  "+")
# 
# # summarize posterior distributions
# a.T0.alldata <- t(apply(a.T0.alldata, 2, quantile, c(0.025,0.5,0.975)))
# a.Tm.alldata <- t(apply(a.Tm.alldata, 2, quantile, c(0.025,0.5,0.975)))
# a.q.alldata <- t(apply(a.q.alldata, 2, quantile, c(0.025,0.5,0.975)))
# 
# # change column names
# colnames(a.T0.alldata) <- c("T0_2.5", "T0_50", "T0_97.5")
# colnames(a.Tm.alldata) <- c("Tm_2.5", "Tm_50", "Tm_97.5")
# colnames(a.q.alldata)  <- c("q_2.5",  "q_50",  "q_97.5")
# 
# 
# a.alldata.id.info <- data.a.all %>%
#   group_by(species, citation) %>% 
#   mutate(unique_id = cur_group_id()) %>% 
#   group_by(unique_id) %>% 
#   dplyr::select(unique_id, genus, species, citation, latitude, longitude, type) %>% 
#   unique() %>% # drop duplicate
#   arrange(unique_id)
# 
# a.TPC.pars.alldata <- bind_cols(a.alldata.id.info, a.T0.alldata, a.Tm.alldata, a.q.alldata)
# a.TPC.pars.alldata
# 
# 
# ###### TPC parameter-latitudinal analysis ######
# 
# a.Tmin.lat <- ggplot() +
#   geom_vline(xintercept = 0, linetype = "dashed") +
#   
#   geom_errorbar(data = a.TPC.pars.alldata, 
#                 aes(x = latitude, ymin = T0_2.5, ymax = T0_97.5), width = 1) +
#   geom_point(data = a.TPC.pars.alldata, 
#              aes(x = latitude, y = T0_50, colour = species)) +
#   
#   labs(
#     x = expression(paste("Latitude")),
#     y = "Tmin",
#     title = "A) a Tmin, median & 95% CI"
#   ) +
#   
#   scale_colour_manual(values = c("albopictus" = "#CB181D",
#                                  "punctor" = "#4292C6",
#                                  "cinereus" = "#2171B5",
#                                  "communis" = "#08519C",
#                                  "impiger" = "#08306B"),
#                       name = element_blank(), # No legend title
#                       breaks = c("albopictus",
#                                  "punctor",
#                                  "cinereus",
#                                  "communis",
#                                  "impiger"),
#                       labels = c("Ae. albopictus", "Ae. punctor", "Ae. cinereus",
#                                  "Ae. communis", "Ae. impiger")) +
#   # 
#   theme_bw()
# 
# a.Tmin.lat
# 
# a.Tmax.lat <- ggplot() +
#   geom_vline(xintercept = 0, linetype = "dashed") +
#   
#   # Non-Arctic
#   geom_errorbar(data = a.TPC.pars.alldata, 
#                 aes(x = latitude, ymin = Tm_2.5, ymax = Tm_97.5), width = 1) +
#   geom_point(data = a.TPC.pars.alldata, 
#              aes(x = latitude, y = Tm_50, colour = species)) +
#   
#   labs(
#     x = expression(paste("Latitude")),
#     y = "Tmax",
#     title = "A) a Tmax, median & 95% CI"
#   ) +
#   
#   scale_colour_manual(values = c("albopictus" = "#CB181D",
#                                  "punctor" = "#4292C6",
#                                  "cinereus" = "#2171B5",
#                                  "communis" = "#08519C",
#                                  "impiger" = "#08306B"),
#                       name = element_blank(), # No legend title
#                       breaks = c("albopictus",
#                                  "punctor",
#                                  "cinereus",
#                                  "communis",
#                                  "impiger"),
#                       labels = c("Ae. albopictus", "Ae. punctor", "Ae. cinereus",
#                                  "Ae. communis", "Ae. impiger")) +
#   
#   theme_bw()
# 
# a.Tmax.lat
# 
# a.q.lat <- ggplot() +
#   geom_vline(xintercept = 0, linetype = "dashed") +
#   
#   # Non-Arctic
#   geom_errorbar(data = a.TPC.pars.alldata, 
#                 aes(x = latitude, ymin = q_2.5, ymax = q_97.5), width = 1) +
#   geom_point(data = a.TPC.pars.alldata, 
#              aes(x = latitude, y = q_50, colour = species)) +
#   
#   labs(
#     x = expression(paste("Latitude")),
#     y = "q",
#     title = "A) a q, median & 95% CI"
#   ) +
#   
#   scale_colour_manual(values = c("albopictus" = "#CB181D",
#                                  "punctor" = "#4292C6",
#                                  "cinereus" = "#2171B5",
#                                  "communis" = "#08519C",
#                                  "impiger" = "#08306B"),
#                       name = element_blank(), # No legend title
#                       breaks = c("albopictus",
#                                  "punctor",
#                                  "cinereus",
#                                  "communis",
#                                  "impiger"),
#                       labels = c("Ae. albopictus", "Ae. punctor", "Ae. cinereus",
#                                  "Ae. communis", "Ae. impiger")) +
#   
#   theme_bw()
# 
# a.q.lat
# 
# TPC.params.lat <- plot_grid(a.Tmin.lat, a.Tmax.lat, a.q.lat, align = "v", ncol = 1)
# TPC.params.lat
# 
# ggsave("figures/a.bri.TPC.params.lat.png", TPC.params.lat,
#        width = 10.3, height = 10)
# 
# ##### Eggs per female per gonotrophic cycle (EFGC) #####
# 
# # Load R2jags model output
# load("R-scripts/R2jags-objects/all-mods/EFGC.alldata.quad.uni.Rdata")
# 
# # Load data
# data.EFGC.all <- read_csv("data-processed/TraitData_EFGC.csv")
# 
# # Get Tmin, Tmax, and q from each random effect
# EFGC.sims <- EFGC.alldata.quad.uni$BUGSoutput$sims.list
# 
# # parameter values for each unique_id
# EFGC.T0.alldata <- sweep(EFGC.sims$T0, 1, EFGC.sims$cf.T0, "+")
# EFGC.Tm.alldata <- sweep(EFGC.sims$Tm, 1, EFGC.sims$cf.Tm, "+")
# EFGC.q.alldata  <- sweep(EFGC.sims$q,  1, EFGC.sims$cf.q,  "+")
# 
# # summarize posterior distributions
# EFGC.T0.alldata <- t(apply(EFGC.T0.alldata, 2, quantile, c(0.025,0.5,0.975)))
# EFGC.Tm.alldata <- t(apply(EFGC.Tm.alldata, 2, quantile, c(0.025,0.5,0.975)))
# EFGC.q.alldata <- t(apply(EFGC.q.alldata, 2, quantile, c(0.025,0.5,0.975)))
# 
# # change column names
# colnames(EFGC.T0.alldata) <- c("T0_2.5", "T0_50", "T0_97.5")
# colnames(EFGC.Tm.alldata) <- c("Tm_2.5", "Tm_50", "Tm_97.5")
# colnames(EFGC.q.alldata)  <- c("q_2.5",  "q_50",  "q_97.5")
# 
# 
# EFGC.alldata.id.info <- data.EFGC.all %>%
#   group_by(species, citation) %>% 
#   mutate(unique_id = cur_group_id()) %>% 
#   group_by(unique_id) %>% 
#   dplyr::select(unique_id, genus, species, citation, latitude, longitude, type) %>% 
#   unique() %>% # drop duplicate
#   arrange(unique_id)
# 
# EFGC.TPC.pars.alldata <- bind_cols(EFGC.alldata.id.info, EFGC.T0.alldata, EFGC.Tm.alldata, EFGC.q.alldata)
# EFGC.TPC.pars.alldata
# 
# 
# ###### TPC parameter-latitudinal analysis ######
# 
# EFGC.Tmin.lat <- ggplot() +
#   geom_vline(xintercept = 0, linetype = "dashed") +
#   
#   geom_errorbar(data = EFGC.TPC.pars.alldata, 
#                 aes(x = latitude, ymin = T0_2.5, ymax = T0_97.5), width = 1) +
#   geom_point(data = EFGC.TPC.pars.alldata, 
#              aes(x = latitude, y = T0_50, colour = species)) +
#   
#   labs(
#     x = expression(paste("Latitude")),
#     y = "Tmin",
#     title = "B) EFGC Tmin, median & 95% CI"
#   ) +
#   
#   scale_colour_manual(values = c("albopictus" = "#CB181D",
#                                  "hexodontus" = "#9ECAE1",
#                                  "punctor" = "#4292C6",
#                                  "cinereus" = "#2171B5",
#                                  "communis" = "#08519C",
#                                  "impiger" = "#08306B"),
#                       name = element_blank(), # No legend title
#                       breaks = c("albopictus",
#                                  "hexodontus",
#                                  "punctor",
#                                  "cinereus",
#                                  "communis",
#                                  "impiger"),
#                       labels = c("Ae. albopictus", "Ae. hexodontus",
#                                  "Ae. punctor", "Ae. cinereus",
#                                  "Ae. communis", "Ae. impiger")) +
#   # 
#   theme_bw()
# 
# EFGC.Tmin.lat
# 
# EFGC.Tmax.lat <- ggplot() +
#   geom_vline(xintercept = 0, linetype = "dashed") +
#   
#   # Non-Arctic
#   geom_errorbar(data = EFGC.TPC.pars.alldata, 
#                 aes(x = latitude, ymin = Tm_2.5, ymax = Tm_97.5), width = 1) +
#   geom_point(data = EFGC.TPC.pars.alldata, 
#              aes(x = latitude, y = Tm_50, colour = species)) +
#   
#   labs(
#     x = expression(paste("Latitude")),
#     y = "Tmax",
#     title = "B) EFGC Tmax, median & 95% CI"
#   ) +
#   
#   scale_colour_manual(values = c("albopictus" = "#CB181D",
#                                  "hexodontus" = "#9ECAE1",
#                                  "punctor" = "#4292C6",
#                                  "cinereus" = "#2171B5",
#                                  "communis" = "#08519C",
#                                  "impiger" = "#08306B"),
#                       name = element_blank(), # No legend title
#                       breaks = c("albopictus",
#                                  "hexodontus",
#                                  "punctor",
#                                  "cinereus",
#                                  "communis",
#                                  "impiger"),
#                       labels = c("Ae. albopictus", "Ae. hexodontus",
#                                  "Ae. punctor", "Ae. cinereus",
#                                  "Ae. communis", "Ae. impiger")) +
#   
#   theme_bw()
# 
# EFGC.Tmax.lat
# 
# EFGC.q.lat <- ggplot() +
#   geom_vline(xintercept = 0, linetype = "dashed") +
#   
#   # Non-Arctic
#   geom_errorbar(data = EFGC.TPC.pars.alldata, 
#                 aes(x = latitude, ymin = q_2.5, ymax = q_97.5), width = 1) +
#   geom_point(data = EFGC.TPC.pars.alldata, 
#              aes(x = latitude, y = q_50, colour = species)) +
#   
#   labs(
#     x = expression(paste("Latitude")),
#     y = "q",
#     title = "B) EFGC q, median & 95% CI"
#   ) +
#   
#   scale_colour_manual(values = c("albopictus" = "#CB181D",
#                                  "hexodontus" = "#9ECAE1",
#                                  "punctor" = "#4292C6",
#                                  "cinereus" = "#2171B5",
#                                  "communis" = "#08519C",
#                                  "impiger" = "#08306B"),
#                       name = element_blank(), # No legend title
#                       breaks = c("albopictus",
#                                  "hexodontus",
#                                  "punctor",
#                                  "cinereus",
#                                  "communis",
#                                  "impiger"),
#                       labels = c("Ae. albopictus", "Ae. hexodontus",
#                                  "Ae. punctor", "Ae. cinereus",
#                                  "Ae. communis", "Ae. impiger")) +
#   
#   theme_bw()
# 
# EFGC.q.lat
# 
# TPC.params.lat <- plot_grid(EFGC.Tmin.lat, EFGC.Tmax.lat, EFGC.q.lat, align = "v", ncol = 1)
# TPC.params.lat
# 
# ggsave("figures/EFGC.quad.TPC.params.lat.png", TPC.params.lat,
#        width = 10.3, height = 10)
# 
# 
# ##### Mosquito adult lifespan (lf) #####
# 
# # Load R2jags model output
# load("R-scripts/R2jags-objects/all-mods/lf.alldata.quad.uni.Rdata")
# 
# # Load data
# data.lf.all <- read_csv("data-processed/TraitData_lf.csv")
# 
# # Get Tmin, Tmax, and q from each random effect
# lf.sims <- lf.alldata.quad.uni$BUGSoutput$sims.list
# 
# # parameter values for each unique_id
# lf.T0.alldata <- sweep(lf.sims$T0, 1, lf.sims$cf.T0, "+")
# lf.Tm.alldata <- sweep(lf.sims$Tm, 1, lf.sims$cf.Tm, "+")
# lf.q.alldata  <- sweep(lf.sims$q,  1, lf.sims$cf.q,  "+")
# 
# # summarize posterior distributions
# lf.T0.alldata <- t(apply(lf.T0.alldata, 2, quantile, c(0.025,0.5,0.975)))
# lf.Tm.alldata <- t(apply(lf.Tm.alldata, 2, quantile, c(0.025,0.5,0.975)))
# lf.q.alldata <- t(apply(lf.q.alldata, 2, quantile, c(0.025,0.5,0.975)))
# 
# # change column names
# colnames(lf.T0.alldata) <- c("T0_2.5", "T0_50", "T0_97.5")
# colnames(lf.Tm.alldata) <- c("Tm_2.5", "Tm_50", "Tm_97.5")
# colnames(lf.q.alldata)  <- c("q_2.5",  "q_50",  "q_97.5")
# 
# 
# lf.alldata.id.info <- data.lf.all %>%
#   group_by(species, citation) %>% 
#   mutate(unique_id = cur_group_id()) %>% 
#   group_by(unique_id) %>% 
#   dplyr::select(unique_id, genus, species, citation, latitude, longitude, type) %>% 
#   unique() %>% # drop duplicate
#   arrange(unique_id)
# 
# # The two data points from Tsuda et al. 1994 were collected from two different locations. 
# # But they have the same unique_id. Change their latitude and longitude to NAs
# lf.alldata.id.info$latitude[lf.alldata.id.info$unique_id == 5] <- NA
# lf.alldata.id.info$longitude[lf.alldata.id.info$unique_id == 5] <- NA
# lf.alldata.id.info <- unique(lf.alldata.id.info)
# 
# lf.TPC.pars.alldata <- bind_cols(lf.alldata.id.info, lf.T0.alldata, lf.Tm.alldata, lf.q.alldata)
# lf.TPC.pars.alldata
# 
# 
# ###### TPC parameter-latitudinal analysis ######
# 
# lf.Tmin.lat <- ggplot() +
#   geom_vline(xintercept = 0, linetype = "dashed") +
#   
#   geom_errorbar(data = lf.TPC.pars.alldata, 
#                 aes(x = latitude, ymin = T0_2.5, ymax = T0_97.5), width = 1) +
#   geom_point(data = lf.TPC.pars.alldata, 
#              aes(x = latitude, y = T0_50, colour = species)) +
#   
#   labs(
#     x = expression(paste("Latitude")),
#     y = "Tmin",
#     title = "C) lf Tmin, median & 95% CI"
#   ) +
#   
#   scale_colour_manual(values = c("albopictus" = "#CB181D",
#                                  "vexans" = "#9ECAE1",
#                                  "punctor" = "#4292C6",
#                                  "cinereus" = "#2171B5",
#                                  "communis" = "#08519C",
#                                  "impiger" = "#08306B"),
#                       name = element_blank(), # No legend title
#                       breaks = c("albopictus",
#                                  "vexans",
#                                  "punctor",
#                                  "cinereus",
#                                  "communis",
#                                  "impiger"),
#                       labels = c("Ae. albopictus", "Ae. vexans",
#                                  "Ae. punctor", "Ae. cinereus",
#                                  "Ae. communis", "Ae. impiger")) +
#   
#   theme_bw()
# 
# lf.Tmin.lat
# 
# lf.Tmax.lat <- ggplot() +
#   geom_vline(xintercept = 0, linetype = "dashed") +
#   
#   # Non-Arctic
#   geom_errorbar(data = lf.TPC.pars.alldata, 
#                 aes(x = latitude, ymin = Tm_2.5, ymax = Tm_97.5), width = 1) +
#   geom_point(data = lf.TPC.pars.alldata, 
#              aes(x = latitude, y = Tm_50, colour = species)) +
#   
#   labs(
#     x = expression(paste("Latitude")),
#     y = "Tmax",
#     title = "C) lf Tmax, median & 95% CI"
#   ) +
#   
#   scale_colour_manual(values = c("albopictus" = "#CB181D",
#                                  "vexans" = "#9ECAE1",
#                                  "punctor" = "#4292C6",
#                                  "cinereus" = "#2171B5",
#                                  "communis" = "#08519C",
#                                  "impiger" = "#08306B"),
#                       name = element_blank(), # No legend title
#                       breaks = c("albopictus",
#                                  "vexans",
#                                  "punctor",
#                                  "cinereus",
#                                  "communis",
#                                  "impiger"),
#                       labels = c("Ae. albopictus", "Ae. vexans",
#                                  "Ae. punctor", "Ae. cinereus",
#                                  "Ae. communis", "Ae. impiger")) +
#   
#   theme_bw()
# 
# lf.Tmax.lat
# 
# lf.q.lat <- ggplot() +
#   geom_vline(xintercept = 0, linetype = "dashed") +
#   
#   # Non-Arctic
#   geom_errorbar(data = lf.TPC.pars.alldata, 
#                 aes(x = latitude, ymin = q_2.5, ymax = q_97.5), width = 1) +
#   geom_point(data = lf.TPC.pars.alldata, 
#              aes(x = latitude, y = q_50, colour = species)) +
#   
#   labs(
#     x = expression(paste("Latitude")),
#     y = "q",
#     title = "C) lf q, median & 95% CI"
#   ) +
#   
#   scale_colour_manual(values = c("albopictus" = "#CB181D",
#                                  "vexans" = "#9ECAE1",
#                                  "punctor" = "#4292C6",
#                                  "cinereus" = "#2171B5",
#                                  "communis" = "#08519C",
#                                  "impiger" = "#08306B"),
#                       name = element_blank(), # No legend title
#                       breaks = c("albopictus",
#                                  "vexans",
#                                  "punctor",
#                                  "cinereus",
#                                  "communis",
#                                  "impiger"),
#                       labels = c("Ae. albopictus", "Ae. vexans",
#                                  "Ae. punctor", "Ae. cinereus",
#                                  "Ae. communis", "Ae. impiger")) +
#   
#   theme_bw()
# 
# lf.q.lat
# 
# TPC.params.lat <- plot_grid(lf.Tmin.lat, lf.Tmax.lat, lf.q.lat, align = "v", ncol = 1)
# TPC.params.lat
# 
# ggsave("figures/lf.quad.TPC.params.lat.png", TPC.params.lat,
#        width = 10.3, height = 10)



##### Egg  viability (EV) #####

# Load R2jags model output
load("R-scripts/R2jags-objects/all-mods/EV.nonarctic.quad.uni.Rdata") # Non-arctic
load("R-scripts/R2jags-objects/all-mods/EV.arctic.quad.uni.Rdata") # Arctic

# Load data
data.EV.all <- read_csv("data-processed/TraitData_EV.csv")
data.EV.all <- data.EV.all %>% 
mutate(latitude = abs(latitude)) 

# Subset data
data.EV.arctic <- subset(data.EV.all, type == "Arctic") # Arctic species
data.EV.nonarctic <- subset(data.EV.all, type == "non-Arctic") # Non-Arctic species

###### Non-Arctic ######
# Get Tmin, Tmax, and q from each random effect
EV.sims <- EV.nonarctic.quad.uni$BUGSoutput$sims.list

# parameter values for each unique_id
EV.T0.nonarctic <- sweep(EV.sims$T0, 1, EV.sims$cf.T0, "+")
EV.Tm.nonarctic <- sweep(EV.sims$Tm, 1, EV.sims$cf.Tm, "+")
EV.q.nonarctic  <- sweep(EV.sims$q,  1, EV.sims$cf.q,  "+")

# summarize posterior distributions
EV.T0.nonarctic <- t(apply(EV.T0.nonarctic, 2, quantile, c(0.025,0.5,0.975)))
EV.Tm.nonarctic <- t(apply(EV.Tm.nonarctic, 2, quantile, c(0.025,0.5,0.975)))
EV.q.nonarctic <- t(apply(EV.q.nonarctic, 2, quantile, c(0.025,0.5,0.975)))

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
EV.T0.arctic <- EV.arctic.quad.uni$BUGSoutput$sims.list$cf.T0
EV.Tm.arctic <- EV.arctic.quad.uni$BUGSoutput$sims.list$cf.Tm
EV.q.arctic <- EV.arctic.quad.uni$BUGSoutput$sims.list$cf.q

EV.T0.arctic <- t(quantile(EV.T0.arctic, c(0.025, 0.5, 0.975)))
EV.Tm.arctic <- t(quantile(EV.Tm.arctic, c(0.025, 0.5, 0.975)))
EV.q.arctic <- t(quantile(EV.q.arctic, c(0.025, 0.5, 0.975)))

colnames(EV.T0.arctic) <- c("T0_2.5", "T0_50", "T0_97.5")
colnames(EV.Tm.arctic) <- c("Tm_2.5", "Tm_50", "Tm_97.5")
colnames(EV.q.arctic)  <- c("q_2.5",  "q_50",  "q_97.5")

EV.arctic.id.info <- data.EV.arctic %>%
  dplyr::select(genus, species, citation, latitude, longitude, type) %>% 
  unique() # drop duplicate

EV.TPC.pars.arctic <- bind_cols(EV.arctic.id.info, EV.T0.arctic, EV.Tm.arctic, EV.q.arctic)
EV.TPC.pars <- bind_rows(EV.TPC.pars.arctic, EV.TPC.pars.nonarctic)

###### TPC parameter-latitudinal analysis ######

EV.Tmin.lat <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  geom_errorbar(data = EV.TPC.pars, 
                aes(x = latitude, ymin = T0_2.5, ymax = T0_97.5), width = 1) +
  geom_point(data = EV.TPC.pars, 
             aes(x = latitude, y = T0_50, colour = species)) +
  
  labs(
    x = "absolute latitude",
    y = "Tmin",
    title = "D) EV Tmin, median & 95% CI"
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

EV.Tmax.lat <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  geom_errorbar(data = EV.TPC.pars, 
                aes(x = latitude, ymin = Tm_2.5, ymax = Tm_97.5), width = 1) +
  geom_point(data = EV.TPC.pars, 
             aes(x = latitude, y = Tm_50, colour = species)) +

  
  labs(
    x = "absolute latitude",
    y = "Tmax",
    title = "D) EV Tmax, median & 95% CI"
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

EV.q.lat <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  geom_errorbar(data = EV.TPC.pars, 
                aes(x = latitude, ymin = q_2.5, ymax = q_97.5), width = 1) +
  geom_point(data = EV.TPC.pars, 
             aes(x = latitude, y = q_50, colour = species)) +
  
  labs(
    x = "absolute latitude",
    y = "q",
    title = "D) EV q, median & 95% CI"
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
data.pLA.all <- data.pLA.all %>% 
  mutate(latitude = abs(latitude)) 

# Subset data
data.pLA.arctic <- subset(data.pLA.all, type == "Arctic") # Arctic species
data.pLA.nonarctic <- subset(data.pLA.all, type == "non-Arctic") # Non-Arctic species

###### Non-Arctic ######
# Get Tmin, Tmax, and q from each random effect
pLA.sims <- pLA.nonarctic.quad.uni$BUGSoutput$sims.list

# parameter values for each unique_id
pLA.T0.nonarctic <- sweep(pLA.sims$T0, 1, pLA.sims$cf.T0, "+")
pLA.Tm.nonarctic <- sweep(pLA.sims$Tm, 1, pLA.sims$cf.Tm, "+")
pLA.q.nonarctic  <- sweep(pLA.sims$q,  1, pLA.sims$cf.q,  "+")

# summarize posterior distributions
pLA.T0.nonarctic <- t(apply(pLA.T0.nonarctic, 2, quantile, c(0.025,0.5,0.975)))
pLA.Tm.nonarctic <- t(apply(pLA.Tm.nonarctic, 2, quantile, c(0.025,0.5,0.975)))
pLA.q.nonarctic <- t(apply(pLA.q.nonarctic, 2, quantile, c(0.025,0.5,0.975)))

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
pLA.T0.arctic <- pLA.arctic.quad.uni$BUGSoutput$sims.list$cf.T0
pLA.Tm.arctic <- pLA.arctic.quad.uni$BUGSoutput$sims.list$cf.Tm
pLA.q.arctic <- pLA.arctic.quad.uni$BUGSoutput$sims.list$cf.q

pLA.T0.arctic <- t(quantile(pLA.T0.arctic, c(0.025, 0.5, 0.975)))
pLA.Tm.arctic <- t(quantile(pLA.Tm.arctic, c(0.025, 0.5, 0.975)))
pLA.q.arctic <- t(quantile(pLA.q.arctic, c(0.025, 0.5, 0.975)))

colnames(pLA.T0.arctic) <- c("T0_2.5", "T0_50", "T0_97.5")
colnames(pLA.Tm.arctic) <- c("Tm_2.5", "Tm_50", "Tm_97.5")
colnames(pLA.q.arctic)  <- c("q_2.5",  "q_50",  "q_97.5")

pLA.arctic.id.info <- data.pLA.arctic %>%
  dplyr::select(genus, species, citation, latitude, longitude, type) %>% 
  unique() # drop duplicate

pLA.TPC.pars.arctic <- bind_cols(pLA.T0.arctic, pLA.Tm.arctic, pLA.q.arctic)
pLA.TPC.pars.arctic$species <- "Arctic spp."


###### TPC parameter-latitudinal analysis ######

pLA.Tmin.lat <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  # Non-Arctic
  geom_errorbar(data = pLA.TPC.pars.nonarctic, 
                aes(x = latitude, ymin = T0_2.5, ymax = T0_97.5), width = 1) +
  geom_point(data = pLA.TPC.pars.nonarctic, 
             aes(x = latitude, y = T0_50, colour = species)) +
  
  # Arctic
  geom_errorbar(data = pLA.TPC.pars.arctic, 
                aes(x = 60, ymin = T0_2.5, ymax = T0_97.5), width = 1) +
  geom_point(data = pLA.TPC.pars.arctic, 
             aes(x = 60, y = T0_50, colour = species)) +
  
  labs(
    x = expression(paste("Latitude")),
    y = "Tmin",
    title = "E) pLA Tmin, median & 95% CI"
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

pLA.Tmax.lat <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  # Non-Arctic
  geom_errorbar(data = pLA.TPC.pars.nonarctic, 
                aes(x = latitude, ymin = Tm_2.5, ymax = Tm_97.5), width = 1) +
  geom_point(data = pLA.TPC.pars.nonarctic, 
             aes(x = latitude, y = Tm_50, colour = species)) +
  
  # Arctic
  geom_errorbar(data = pLA.TPC.pars.arctic, 
                aes(x = 60, ymin = Tm_2.5, ymax = Tm_97.5), width = 1) +
  geom_point(data = pLA.TPC.pars.arctic, 
             aes(x = 60, y = Tm_50, colour = species)) +
  
  labs(
    x = expression(paste("Latitude")),
    y = "Tmax",
    title = "E) pLA Tmax, median & 95% CI"
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

pLA.q.lat <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  # Non-Arctic
  geom_errorbar(data = pLA.TPC.pars.nonarctic, 
                aes(x = latitude, ymin = q_2.5, ymax = q_97.5), width = 1) +
  geom_point(data = pLA.TPC.pars.nonarctic, 
             aes(x = latitude, y = q_50, colour = species)) +
  
  # Arctic
  geom_errorbar(data = pLA.TPC.pars.arctic, 
                aes(x = 60, ymin = q_2.5, ymax = q_97.5), width = 1) +
  geom_point(data = pLA.TPC.pars.arctic, 
             aes(x = 60, y = q_50, colour = species)) +
  
  labs(
    x = expression(paste("Latitude")),
    y = "q",
    title = "E) pLA q, median & 95% CI"
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
data.MDR.all <- data.MDR.all %>% 
  mutate(latitude = abs(latitude))

# Subset data
data.MDR.arctic <- subset(data.MDR.all, type == "Arctic") # Arctic species
data.MDR.nonarctic <- subset(data.MDR.all, type == "non-Arctic") # Non-Arctic species

###### Non-Arctic ######
# Get Tmin, Tmax, and q from each random effect
MDR.sims <- MDR.nonarctic.quad.uni$BUGSoutput$sims.list

# parameter values for each unique_id
MDR.T0.nonarctic <- sweep(MDR.sims$T0, 1, MDR.sims$cf.T0, "+")
MDR.Tm.nonarctic <- sweep(MDR.sims$Tm, 1, MDR.sims$cf.Tm, "+")
MDR.q.nonarctic  <- sweep(MDR.sims$q,  1, MDR.sims$cf.q,  "+")

# summarize posterior distributions
MDR.T0.nonarctic <- t(apply(MDR.T0.nonarctic, 2, quantile, c(0.025,0.5,0.975)))
MDR.Tm.nonarctic <- t(apply(MDR.Tm.nonarctic, 2, quantile, c(0.025,0.5,0.975)))
MDR.q.nonarctic <- t(apply(MDR.q.nonarctic, 2, quantile, c(0.025,0.5,0.975)))

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
MDR.T0.arctic <- MDR.arctic.quad.uni$BUGSoutput$sims.list$cf.T0
MDR.Tm.arctic <- MDR.arctic.quad.uni$BUGSoutput$sims.list$cf.Tm
MDR.q.arctic <- MDR.arctic.quad.uni$BUGSoutput$sims.list$cf.q

MDR.T0.arctic <- t(quantile(MDR.T0.arctic, c(0.025, 0.5, 0.975)))
MDR.Tm.arctic <- t(quantile(MDR.Tm.arctic, c(0.025, 0.5, 0.975)))
MDR.q.arctic <- t(quantile(MDR.q.arctic, c(0.025, 0.5, 0.975)))

colnames(MDR.T0.arctic) <- c("T0_2.5", "T0_50", "T0_97.5")
colnames(MDR.Tm.arctic) <- c("Tm_2.5", "Tm_50", "Tm_97.5")
colnames(MDR.q.arctic)  <- c("q_2.5",  "q_50",  "q_97.5")

MDR.arctic.id.info <- data.MDR.arctic %>%
  dplyr::select(genus, species, citation, latitude, longitude, type) %>% 
  unique() # drop duplicate

MDR.TPC.pars.arctic <- bind_cols(MDR.T0.arctic, MDR.Tm.arctic, MDR.q.arctic)
MDR.TPC.pars.arctic$species <- "Arctic spp."


###### TPC parameter-latitudinal analysis ######

MDR.Tmin.lat <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  # Non-Arctic
  geom_errorbar(data = MDR.TPC.pars.nonarctic, 
                aes(x = latitude, ymin = T0_2.5, ymax = T0_97.5), width = 1) +
  geom_point(data = MDR.TPC.pars.nonarctic, 
             aes(x = latitude, y = T0_50, colour = species)) +
  
  # Arctic
  geom_errorbar(data = MDR.TPC.pars.arctic, 
                aes(x = 60, ymin = T0_2.5, ymax = T0_97.5), width = 1) +
  geom_point(data = MDR.TPC.pars.arctic, 
             aes(x = 60, y = T0_50, colour = species)) +
  
  labs(
    x = expression(paste("Latitude")),
    y = "Tmin",
    title = "F) MDR Tmin, median & 95% CI"
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

MDR.Tmax.lat <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  # Non-Arctic
  geom_errorbar(data = MDR.TPC.pars.nonarctic, 
                aes(x = latitude, ymin = Tm_2.5, ymax = Tm_97.5), width = 1) +
  geom_point(data = MDR.TPC.pars.nonarctic, 
             aes(x = latitude, y = Tm_50, colour = species)) +
  
  # Arctic
  geom_errorbar(data = MDR.TPC.pars.arctic, 
                aes(x = 60, ymin = Tm_2.5, ymax = Tm_97.5), width = 1) +
  geom_point(data = MDR.TPC.pars.arctic, 
             aes(x = 60, y = Tm_50, colour = species)) +
  
  labs(
    x = expression(paste("Latitude")),
    y = "Tmax",
    title = "F) MDR Tmax, median & 95% CI"
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

MDR.q.lat <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  # Non-Arctic
  geom_errorbar(data = MDR.TPC.pars.nonarctic, 
                aes(x = latitude, ymin = q_2.5, ymax = q_97.5), width = 1) +
  geom_point(data = MDR.TPC.pars.nonarctic, 
             aes(x = latitude, y = q_50, colour = species)) +
  
  # Arctic
  geom_errorbar(data = MDR.TPC.pars.arctic, 
                aes(x = 60, ymin = q_2.5, ymax = q_97.5), width = 1) +
  geom_point(data = MDR.TPC.pars.arctic, 
             aes(x = 60, y = q_50, colour = species)) +
  
  labs(
    x = expression(paste("Latitude")),
    y = "q",
    title = "F) MDR q, median & 95% CI"
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
data.PDR.all <- data.PDR.all %>% 
  mutate(latitude = abs(latitude))

# Subset data
data.PDR.arctic <- subset(data.PDR.all, type == "Arctic") # Arctic species
data.PDR.nonarctic <- subset(data.PDR.all, type == "non-Arctic") # Non-Arctic species

###### Non-Arctic ######
# Get Tmin, Tmax, and q from each random effect
PDR.sims <- PDR.nonarctic.bri.uni$BUGSoutput$sims.list

# parameter values for each unique_id
PDR.T0.nonarctic <- sweep(PDR.sims$T0, 1, PDR.sims$cf.T0, "+")
PDR.Tm.nonarctic <- sweep(PDR.sims$Tm, 1, PDR.sims$cf.Tm, "+")
PDR.q.nonarctic  <- sweep(PDR.sims$q,  1, PDR.sims$cf.q,  "+")

# summarize posterior distributions
PDR.T0.nonarctic <- t(apply(PDR.T0.nonarctic, 2, quantile, c(0.025,0.5,0.975)))
PDR.Tm.nonarctic <- t(apply(PDR.Tm.nonarctic, 2, quantile, c(0.025,0.5,0.975)))
PDR.q.nonarctic <- t(apply(PDR.q.nonarctic, 2, quantile, c(0.025,0.5,0.975)))

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
PDR.T0.arctic <- PDR.arctic.bri.uni$BUGSoutput$sims.list$cf.T0
PDR.Tm.arctic <- PDR.arctic.bri.uni$BUGSoutput$sims.list$cf.Tm
PDR.q.arctic <- PDR.arctic.bri.uni$BUGSoutput$sims.list$cf.q

PDR.T0.arctic <- t(quantile(PDR.T0.arctic, c(0.025, 0.5, 0.975)))
PDR.Tm.arctic <- t(quantile(PDR.Tm.arctic, c(0.025, 0.5, 0.975)))
PDR.q.arctic <- t(quantile(PDR.q.arctic, c(0.025, 0.5, 0.975)))

colnames(PDR.T0.arctic) <- c("T0_2.5", "T0_50", "T0_97.5")
colnames(PDR.Tm.arctic) <- c("Tm_2.5", "Tm_50", "Tm_97.5")
colnames(PDR.q.arctic)  <- c("q_2.5",  "q_50",  "q_97.5")

PDR.arctic.id.info <- data.PDR.arctic %>%
  dplyr::select(paras_genus, paras_species, host_genus, host_species, 
                citation, latitude, longitude, type) %>% 
  unique() # drop duplicate

PDR.TPC.pars.arctic <- bind_cols(PDR.T0.arctic, PDR.Tm.arctic, PDR.q.arctic)
PDR.TPC.pars.arctic$species <- "Arctic spp."


###### TPC parameter-latitudinal analysis ######

PDR.Tmin.lat <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  # Non-Arctic
  geom_errorbar(data = PDR.TPC.pars.nonarctic, 
                aes(x = latitude, ymin = T0_2.5, ymax = T0_97.5), width = 1) +
  geom_point(data = PDR.TPC.pars.nonarctic, 
             aes(x = latitude, y = T0_50, colour = paras_species)) +
  
  # Arctic
  geom_errorbar(data = PDR.TPC.pars.arctic, 
                aes(x = 60, ymin = T0_2.5, ymax = T0_97.5), width = 1) +
  geom_point(data = PDR.TPC.pars.arctic, 
             aes(x = 60, y = T0_50, colour = species)) +
  
  labs(
    x = expression(paste("Latitude")),
    y = "Tmin",
    title = "G) PDR Tmin, median & 95% CI"
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

PDR.Tmax.lat <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  # Non-Arctic
  geom_errorbar(data = PDR.TPC.pars.nonarctic, 
                aes(x = latitude, ymin = Tm_2.5, ymax = Tm_97.5), width = 1) +
  geom_point(data = PDR.TPC.pars.nonarctic, 
             aes(x = latitude, y = Tm_50, colour = paras_species)) +
  
  # Arctic
  geom_errorbar(data = PDR.TPC.pars.arctic, 
                aes(x = 60, ymin = Tm_2.5, ymax = Tm_97.5), width = 1) +
  geom_point(data = PDR.TPC.pars.arctic, 
             aes(x = 60, y = Tm_50, colour = species)) +
  
  labs(
    x = expression(paste("Latitude")),
    y = "Tmax",
    title = "G) PDR Tmax, median & 95% CI"
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

PDR.q.lat <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  # Non-Arctic
  geom_errorbar(data = PDR.TPC.pars.nonarctic, 
                aes(x = latitude, ymin = q_2.5, ymax = q_97.5), width = 1) +
  geom_point(data = PDR.TPC.pars.nonarctic, 
             aes(x = latitude, y = q_50, colour = paras_species)) +
  
  # Arctic
  geom_errorbar(data = PDR.TPC.pars.arctic, 
                aes(x = 60, ymin = q_2.5, ymax = q_97.5), width = 1) +
  geom_point(data = PDR.TPC.pars.arctic, 
             aes(x = 60, y = q_50, colour = species)) +
  
  labs(
    x = expression(paste("Latitude")),
    y = "q",
    title = "G) PDR q, median & 95% CI"
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




Tmin.lat <- plot_grid(a.Tmin.lat, EFGC.Tmin.lat, lf.Tmin.lat, EV.Tmin.lat, 
                      pLA.Tmin.lat, MDR.Tmin.lat, PDR.Tmin.lat,
                      align = "v", ncol = 1)
Tmin.lat

ggsave("figures/Tmin.lat.png", Tmin.lat,
       width = 10.3, height = 24)


Tmax.lat <- plot_grid(a.Tmax.lat, EFGC.Tmax.lat, lf.Tmax.lat, EV.Tmax.lat, 
                      pLA.Tmax.lat, MDR.Tmax.lat, PDR.Tmax.lat,
                      align = "v", ncol = 1)
Tmax.lat

ggsave("figures/Tmax.lat.png", Tmax.lat,
       width = 10.3, height = 24)

q.lat <- plot_grid(a.q.lat, EFGC.q.lat, lf.q.lat, EV.q.lat, 
                   pLA.q.lat, MDR.q.lat, PDR.q.lat,
                   align = "v", ncol = 1)
q.lat

ggsave("figures/q.lat.png", q.lat,
       width = 10.3, height = 24)
