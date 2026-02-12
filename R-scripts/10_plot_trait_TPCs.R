## Lilian Chan, University of Guelph
## Arctic vector-borne disease transmission suitability model
##
## Purpose: plotting trait data and TPC fits
## 
## Table of content:
##    0. Set-up workspace
##    1. Load data and model output
##  	2. Plot panels for each trait
##  	4. Manuscript Figure 2
##  	5. Load, process, and plot data for other traits (bc, EIP50, pEA, MDR, gamma)
##	  6. Manuscript Figure S1


##########
###### 0. Set-up workspace ----
##########

library(tidyverse)
library(readxl)
library(janitor)
library(R2jags)
library(MASS)
library(ggsci)
library(ggpubr) # For ggarrange
library(grafify)

##### Load functions
source("R-scripts/00_Functions.R")


##########
###### 1. Load data and model output ----
##########

############## biting rate (a) ---------------------------------------------------------------
## Load data
data.a <- read_csv("data-processed/TraitData_a.csv")

## Convert genotrophic cycle duration (1/a) to biting rate (a)
data.a <- data.a %>% 
  mutate(trait = ifelse(trait_name == "1/a", 1/trait, trait)) %>% 
  mutate(trait_name = ifelse(trait_name == "1/a", "a", trait_name))

# Process trait data for plotting
data.a.summary <- processTraitData(data.a, "a")

a.predictions.summary <- read.csv("data-processed/a.predictions.summary.csv")
a.params.summary <- read.csv("data-processed/a.params.summary.csv")

## Load R2jag model output
# load("R-scripts/R2jags-objects/a.alldata.bri.uni.raneff.Rdata")


############## Infection efficiency (c) ---------------------------------------------------------------
## Load data
data.c <- read_csv("data-processed/TraitData_c.csv")
c.predictions.summary <- read.csv("data-processed/c.arctic.predictions.summary.csv")
c.params.summary <- read.csv("data-processed/c.arctic.params.summary.csv")

# Process trait data for plotting
data.c.summary <- processTraitData(data.c, "c")

## Load R2jags model output
# load("R-scripts/R2jags-objects/c.nonarctic.quad.uni.Rdata")


############## Adult lifespan (lf) ---------------------------------------------------------------
## Load data
data.lf <- read_csv("data-processed/TraitData_lf.csv")
lf.predictions.summary <- read.csv("data-processed/lf.predictions.summary.csv")
lf.params.summary <- read.csv("data-processed/lf.params.summary.csv")

# Process trait data for plotting
data.lf.arctic <- data.lf %>% 
  filter(type == "Arctic")
  
data.lf.summary <- processTraitData(data.lf.arctic, "lf")

## Load R2jags model output
# load("R-scripts/R2jags-objects/lf.arctic.quad.inf.raneff.Rdata")


############## Parasite development rate (PDR) ---------------------------------------------------------------
## Load data
data.PDR <- read_csv("data-processed/TraitData_PDR.csv")

## Convert development time (1/PDR) to development rate (PDR)
data.PDR <- data.PDR %>% 
  mutate(trait = 1/trait) %>% 
  mutate(trait_name = "PDR") 

PDR.predictions.summary <- read.csv("data-processed/PDR.predictions.summary.csv")
PDR.params.summary <- read.csv("data-processed/PDR.params.summary.csv")

# Process trait data for plotting
data.PDR.arctic <- data.PDR %>% 
  filter(type == "Arctic")

data.PDR.summary <- processTraitData(data.PDR.arctic, "PDR")

## Load R2jags model output
# load("R-scripts/R2jags-objects/PDR.arctic.bri.inf.Rdata")


############## Eggs per female per gonotrophic cycle (EFGC) ---------------------------------------------------------------
## Load data
data.EFGC <- read_csv("data-processed/TraitData_EFGC.csv")
EFGC.predictions.summary <- read.csv("data-processed/EFGC.arctic.predictions.summary.csv")
EFGC.params.summary <- read.csv("data-processed/EFGC.arctic.params.summary.csv")

data.EFGC.summary <- processTraitData(data.EFGC, "EFGC")

## Load R2jags model output
# load("R-scripts/R2jags-objects/EFGC.nonarctic.quad.uni.raneff.Rdata")


############## Egg viability (EV) ---------------------------------------------------------------
## Load data
data.EV <- read_csv("data-processed/TraitData_EV.csv")

EV.predictions.summary <- read.csv("data-processed/EV.predictions.summary.csv")
EV.params.summary <- read.csv("data-processed/EV.params.summary.csv")

# Process trait data for plotting
data.EV.arctic <- data.EV %>% 
  filter(type == "Arctic")

data.EV.summary <- processTraitData(data.EV.arctic, "EV")

## Load R2jags model output
# load("R-scripts/R2jags-objects/EV.arctic.quad.inf.Rdata")


############## Larval-to-adult survival (pLA) ---------------------------------------------------------------
## Load data
data.pLA <- read_csv("data-processed/TraitData_pLA.csv")

pLA.predictions.summary <- read.csv("data-processed/pLA.predictions.summary.csv")
pLA.params.summary <- read.csv("data-processed/pLA.params.summary.csv")

# Process trait data for plotting
data.pLA.arctic <- data.pLA %>% 
  filter(type == "Arctic")

data.pLA.summary <- processTraitData(data.pLA.arctic, "pLA")

## Load R2jags model output
# load("R-scripts/R2jags-objects/pLA.arctic.quad.inf.Rdata")


############## Mosquito development rate (MDR) ---------------------------------------------------------------
## Load data
data.MDR <- read_csv("data-processed/TraitData_MDR.csv")

## Convert development time (1/MDR) to development rate (MDR)
data.MDR <- data.MDR %>% 
  mutate(trait = ifelse(trait_name == "1/MDR", 1/trait, trait)) %>% 
  mutate(trait_name = "MDR") 

MDR.predictions.summary <- read.csv("data-processed/MDR.predictions.summary.csv")
MDR.params.summary <- read.csv("data-processed/MDR.params.summary.csv")


# Process trait data for plotting
data.MDR.arctic <- data.MDR %>% 
  filter(type == "Arctic")

data.MDR.summary <- processTraitData(data.MDR.arctic, "MDR")

## Load R2jags model output
# load("R-scripts/R2jags-objects/MDR.arctic.bri.inf.Rdata")




##########
###### 2. Plot panels for each trait ----
##########

############## biting rate (a) ---------------------------------------------------------------
plot.a <- a.predictions.summary %>% 
  ggplot() +
  geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI), fill = "#E69F00", alpha = 0.5) +
  geom_line(aes(x = temperature, y = median), color = "#E69F00", linewidth = 1) +
  #geom_point(data = data.a, aes(x = temp, y = trait), size = 2) +
  geom_pointrange(data = data.a.summary, aes(x = temp, ymin = mean - std_error, ymax = mean + std_error, y = mean), size = 0.5) +
  # Customize the axes and labels
  scale_x_continuous(limits = c(0, 46)) + 
  labs(title = "A", x = expression(paste("Temperature (", degree, "C)")), 
       y = parse(text = "Bite~rate~(day^-1)")) +
  annotate("text", x = 1, y = 0.44, label = expression(paste(italic("a"))), size = 5) +
  theme_bw()

plot.a



############## Infection efficiency (c) ---------------------------------------------------------------
plot.c <- c.predictions.summary %>% 
  ggplot() +
  geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI), fill = "#009E73", alpha = 0.5) +
  geom_line(aes(x = temperature, y = median), color = "#009E73", linewidth = 1) +
  #geom_point(data = data.c, aes(x = temp, y = trait), size = 2) +
  #geom_pointrange(data = data.c.summary, aes(x = temp, ymin = mean - std_error, ymax = mean + std_error, y = mean), size = 0.5) +
  # Customize the axes and labels
  scale_x_continuous(limits = c(0, 46)) + 
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Infection proportion") +
  annotate("text", x = 1, y = 0.97, label = expression(paste(italic("c"))), size = 5) +
  theme_bw()

plot.c


##############  Adult lifespan (lf) ---------------------------------------------------------------
plot.lf <- lf.predictions.summary %>% 
  ggplot() +
  geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI), fill = "#0072B2", alpha = 0.5) +
  geom_line(aes(x = temperature, y = median), colour = "#0072B2", linewidth = 1) +
  #geom_point(data = data.lf.arctic, aes(x = temp, y = trait, color = species), size = 2) +
  geom_pointrange(data = data.lf.summary, aes(x = temp, ymin = mean - std_error, ymax = mean + std_error, y = mean), size = 0.5) +
  # Customize the axes and labels
  scale_x_continuous(limits = c(0, 46)) + 
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Mosquito adult lifespan (days)") +
  annotate("text", x = 1, y = 72, label = expression(paste(italic("lf"))), size = 5) +
  theme_bw()

plot.lf


##############  Parasite development rate (PDR) ---------------------------------------------------------------
plot.PDR <- PDR.predictions.summary %>% 
  ggplot() +
  geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI), fill = "#CC79A7", alpha = 0.5) +
  geom_line(aes(x = temperature, y = median), colour = "#CC79A7", linewidth = 1) +
  #geom_point(data = data.PDR.arctic, aes(x = temp, y = trait, color = species), size = 2) +
  geom_pointrange(data = data.PDR.summary, aes(x = temp, ymin = mean - std_error, ymax = mean + std_error, y = mean), size = 0.5) +
  # Customize the axes and labels
  scale_x_continuous(limits = c(0, 46)) + 
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = parse(text = "Parasite~development~rate~(day^-1)")) +
  annotate("text", x = 1, y = 0.22, label = expression(paste(italic("PDR"))), size = 5) +
  theme_bw()

plot.PDR


##############  Eggs per female per gonotrophic cycle (EFGC) ---------------------------------------------------------------
plot.EFGC <- EFGC.predictions.summary %>% 
  ggplot() +
  geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI), fill = "#56B4E9", alpha = 0.5) +
  geom_line(aes(x = temperature, y = median), colour = "#56B4E9", linewidth = 1) +
  #geom_point(data = data.EFGC.arctic, aes(x = temp, y = trait, color = species), size = 2) +
  #geom_pointrange(data = data.EFGC.summary, aes(x = temp, ymin = mean - std_error, ymax = mean + std_error, y = mean), size = 0.5) +
  # Customize the axes and labels
  scale_x_continuous(limits = c(0, 46)) + 
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Eggs per female per gonotrophic cycle") +
  annotate("text", x = 2, y = 68, label = expression(paste(italic("EFGC"))), size = 5) +
  theme_bw()

plot.EFGC


##############  Egg viability (EV) ---------------------------------------------------------------
plot.EV <- EV.predictions.summary %>% 
  ggplot() +
  geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI), fill = "#F5C710", alpha = 0.5) +
  geom_line(aes(x = temperature, y = median), colour = "#F5C710", linewidth = 1) +
  #geom_point(data = data.EV.arctic, aes(x = temp, y = trait, color = species), size = 2) +
  geom_pointrange(data = data.EV.summary, aes(x = temp, ymin = mean - std_error, ymax = mean + std_error, y = mean), size = 0.5) +
  # Customize the axes and labels
  scale_x_continuous(limits = c(0, 46)) + 
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Egg viability (%)") +
  annotate("text", x = 1, y = 0.98, label = expression(paste(italic("EV"))), size = 5) +
  theme_bw()

plot.EV


##############  Larval-to-adult survival (pLA) ---------------------------------------------------------------
plot.pLA <- pLA.predictions.summary %>% 
  ggplot() +
  geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI), fill = "#999999", alpha = 0.5) +
  geom_line(aes(x = temperature, y = median), colour = "#999999", linewidth = 1) +
  #geom_point(data = data.pLA.arctic, aes(x = temp, y = trait, color = species), size = 2) +
  geom_pointrange(data = data.pLA.summary, aes(x = temp, ymin = mean - std_error, ymax = mean + std_error, y = mean), size = 0.5) +
  # Customize the axes and labels
  scale_x_continuous(limits = c(0, 46)) + 
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Larval survival (%)") +
  annotate("text", x = 1, y = 0.85, label = expression(paste(italic("pLA"))), size = 5) +
  theme_bw()

plot.pLA


##############  Mosquito development rate (MDR) ---------------------------------------------------------------
plot.MDR <- MDR.predictions.summary %>% 
  ggplot() +
  geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI), fill = "#D55E00", alpha = 0.5) +
  geom_line(aes(x = temperature, y = median), colour = "#D55E00", linewidth = 1) +
  #geom_point(data = data.MDR.arctic, aes(x = temp, y = trait, color = species), size = 2) +
  geom_pointrange(data = data.MDR.summary, aes(x = temp, ymin = mean - std_error, ymax = mean + std_error, y = mean), size = 0.5) +
  # Customize the axes and labels
  scale_x_continuous(limits = c(0, 46)) + 
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = parse(text = "Mosquito~development~rate~(day^-1)")) +
  annotate("text", x = 1, y = 0.15, label = expression(paste(italic("MDR"))), size = 5) +
  theme_bw()

plot.MDR


plot.traits <- ggarrange(plot.a, plot.lf, plot.MDR,
                         plot.EFGC, plot.EV, plot.pLA, plot.c, plot.PDR,
                         nrow = 3, ncol = 3, align = "hv", heights = c(5,5,5)) + 
  bgcolor("white")       

plot.traits

ggsave("figures/trait.TPCs.png", plot.traits, width = 15, height = 9)


############## TPC parameters ---------------------------------------------------------------
params.summary <- bind_rows(a.params.summary, c.params.summary, lf.params.summary, PDR.params.summary,
                            EFGC.params.summary, EV.params.summary, pLA.params.summary, MDR.params.summary)

# write_csv(params.summary, "data-processed/params.summary.csv")

params.summary <- params.summary %>% 
  filter(term %in% c("cf.T0", "Topt", "cf.Tm")) %>% 
  mutate(trait = factor(trait, levels = c("a", "PDR", "pLA", "EFGC",
                                                  "lf", "c", "EV", "MDR"))) %>% 
  mutate(term = factor(term, levels = c("cf.Tm", "Topt", "cf.T0"))) 
  
plot.params <- params.summary %>%
  ggplot() +
  geom_linerange(aes(xmin = lowerCI, xmax = upperCI, y = trait, 
                     colour = trait),
                 linewidth = 0.8, position = position_dodge2(width = 0.6)) +
  geom_point(aes(x = median, y = trait, colour = trait), size = 2.5, 
             position = position_dodge2(width = 0.6)) +
  labs(title = "B",
       x = expression(paste("Temperature (", degree, "C)"))) +
  scale_y_discrete(labels=c("cf.T0" = expression(paste("T"[min])),
                            "Topt" = expression(paste("T"[opt])),
                            "cf.Tm" = expression(paste("T"[max]))),
                   breaks = c("cf.T0", "Topt", "cf.Tm")) +
  scale_x_continuous(breaks = seq(0, 60, 10)) +
  #scale_color_manual(values = c("#619CFF", "#F8766D", "#00BA38")) +
  scale_colour_manual(values = c("a" = "#E69F00", "c" = "#009E73", "lf" = "#0072B2",
                                 "PDR" = "#CC79A7", "EFGC" = "#56B4E9", 
                                 "EV" = "#F5C710", "pLA" = "#999999", "MDR" = "#D55E00"),
                      name = element_blank(), # No legend title
                      #breaks = c("a", "c", "lf", "PDR", "EFGC", "EV", "pLA", "MDR"),
                      #labels = c("a", "c", "lf", "PDR", "EFGC", "EV", "pLA", "MDR")
                      # Sort by Tmin:
                      breaks = c("MDR", "EV", "c", "lf","EFGC", "pLA", "PDR", "a"),
                      labels = c("MDR", "EV", "c", "lf","EFG", "CpLA", "PDR", "a")) +
  theme(axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

plot.params




#### Check the symmetry of the TPCs curves ----
prediction.summary <- bind_rows(a.predictions.summary, c.predictions.summary, 
                                lf.predictions.summary, PDR.predictions.summary,
                                EFGC.predictions.summary, EV.predictions.summary,
                                pLA.predictions.summary, MDR.predictions.summary)


prediction.summary <- prediction.summary %>% 
  group_by(trait) %>% 
  mutate(scaled_mean = mean / max(mean)) %>% 
  mutate(scaled_median = median / max(median)) %>% 
  ungroup()


## Load Suitabiliy predictions
prediction.S <- read_csv("data-processed/S.output.median.csv")


plot.traits.scaled <- prediction.summary %>% 
  ggplot(aes(x = temperature, y = scaled_median)) +
  geom_line((aes(colour = trait)), linewidth = 0.8) +
  geom_line(data = prediction.S, aes(x = temperature, y = scaled_median, colour = "S"),
            linewidth = 1.5) +
  labs(title = "C",
       x = expression(paste("Temperature (", degree, "C)")), 
       y = "Trait value (scaled)") +
  scale_colour_manual(values = c("S" = "#000000", "a" = "#E69F00", "c" = "#009E73",
                                 "lf" = "#0072B2", "PDR" = "#CC79A7", "EFGC" = "#56B4E9", 
                                 "EV" = "#F5C710", "pLA" = "#999999", "MDR" = "#D55E00"),
                      name = element_blank(), # No legend title
                      breaks = c("S", "EV", "lf", "c", "MDR", "EFGC", "pLA", "PDR", "a"),
                      labels = c("S", "EV", "lf", "c", "MDR", "EFGC", "pLA", "PDR", "a")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

plot.traits.scaled

ggsave("figures/trait.TPCs.scaled.png", plot.traits.scaled, width = 10.3, height = 5.6)


plot.summary <- ggarrange(plot.params, plot.traits.scaled, align = "hv",
                      nrow = 2, heights = c(2,2)) + bgcolor("white")

plot.summary

plot.all <- ggarrange(plot.traits, plot.summary,
                      nrow = 2, heights = c(3,4)) + bgcolor("white")
plot.all

ggsave("figures/trait.TPCs.summary.png", plot.all, width = 15, height = 18)
