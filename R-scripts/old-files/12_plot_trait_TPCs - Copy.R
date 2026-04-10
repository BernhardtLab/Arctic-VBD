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


# 0. Set-up workspace ----------------------------------------------------------#

library(tidyverse)
library(readxl)
library(janitor)
library(ggsci)
library(ggpubr) # For ggarrange
library(grafify)

##### Load functions
source("R-scripts/00_Functions.R")



#  1. Load data and model output -----------------------------------------------

##### biting rate (a) #####
## Load data
data.a <- read_csv("data-processed/TraitData_a.csv")

## Convert gonotrophic cycle duration (1/a) to biting rate (a)
data.a <- data.a %>% 
  mutate(trait = ifelse(trait_name == "1/a", 1/trait, trait)) %>% 
  mutate(trait_name = ifelse(trait_name == "1/a", "a", trait_name))

# Process trait data for plotting
data.a.alldata.summary <- processTraitData(data.a, "a")

a.alldata.predictions.summary <- read.csv("data-processed/a/a.alldata.predictions.summary.csv")
a.alldata.params.summary <- read.csv("data-processed/a/a.alldata.params.summary.csv")


##### Infection efficiency (c) #####
## Load data
data.c <- read_csv("data-processed/TraitData_c.csv")

# Process trait data for plotting
data.c.nonarctic.summary <- processTraitData(data.c, "c")

## Arctic
c.arctic.predictions.summary <- read.csv("data-processed/c/c.arctic.predictions.summary.csv")
c.arctic.params.summary <- read.csv("data-processed/c/c.arctic.params.summary.csv")

# Non-Arctic
c.nonarctic.predictions.summary <- read.csv("data-processed/c/c.nonarctic.predictions.summary.csv")
c.nonarctic.params.summary <- read.csv("data-processed/c/c.nonarctic.params.summary.csv")


##### Adult lifespan (lf) #####
## Load data
data.lf <- read_csv("data-processed/TraitData_lf.csv")

## Convert mortality rate (1/lf) to lifespan (lf)
data.lf <- data.lf %>% 
  mutate(trait = ifelse(trait_name == "1/lf", 1/trait, trait)) %>% 
  mutate(trait_name = "lf") 

# Subset data
## Arctic species
data.lf.arctic <- subset(data.lf, type == "Arctic")

## Non-Arctic species
data.lf.nonarctic <- subset(data.lf, type == "non-Arctic")

data.lf.arctic.summary <- processTraitData(data.lf.arctic, "lf")
data.lf.nonarctic.summary <- processTraitData(data.lf.nonarctic, "lf")

## Arctic
lf.arctic.predictions.summary <- read.csv("data-processed/lf/lf.arctic.predictions.summary.csv")
lf.arctic.params.summary <- read.csv("data-processed/lf/lf.arctic.params.summary.csv")

# Non-Arctic
lf.nonarctic.predictions.summary <- read.csv("data-processed/lf/lf.nonarctic.predictions.summary.csv")
lf.nonarctic.params.summary <- read.csv("data-processed/lf/lf.nonarctic.params.summary.csv")



##### Parasite development rate (PDR) ##### 
## Load data
data.PDR <- read_csv("data-processed/TraitData_PDR.csv")

## Convert development time (1/PDR) to development rate (PDR)
data.PDR <- data.PDR %>% 
  mutate(trait = 1/trait) %>% 
  mutate(trait_name = "PDR") 

# Subset data
## Arctic species
data.PDR.arctic <- subset(data.PDR, type == "Arctic")

## Non-Arctic species
data.PDR.nonarctic <- subset(data.PDR, type == "non-Arctic")

data.PDR.arctic.summary <- processTraitData(data.PDR.arctic, "PDR")
data.PDR.nonarctic.summary <- processTraitData(data.PDR.nonarctic, "PDR")

## Arctic
PDR.arctic.predictions.summary <- read.csv("data-processed/PDR/PDR.arctic.predictions.summary.csv")
PDR.arctic.params.summary <- read.csv("data-processed/PDR/PDR.arctic.params.summary.csv")

# Non-Arctic
PDR.nonarctic.predictions.summary <- read.csv("data-processed/PDR/PDR.nonarctic.predictions.summary.csv")
PDR.nonarctic.params.summary <- read.csv("data-processed/PDR/PDR.nonarctic.params.summary.csv")



##### Eggs per female per gonotrophic cycle (EFGC) #####
## Load data
data.EFGC <- read_csv("data-processed/TraitData_EFGC.csv")


# Subset data
## Non-Arctic species
data.EFGC.nonarctic <- subset(data.EFGC, type == "non-Arctic")

data.EFGC.nonarctic.summary <- processTraitData(data.EFGC.nonarctic, "EFGC")

## Arctic
EFGC.arctic.predictions.summary <- read.csv("data-processed/EFGC/EFGC.arctic.predictions.summary.csv")
EFGC.arctic.params.summary <- read.csv("data-processed/EFGC/EFGC.arctic.params.summary.csv")

# Non-Arctic
EFGC.nonarctic.predictions.summary <- read.csv("data-processed/EFGC/EFGC.nonarctic.predictions.summary.csv")
EFGC.nonarctic.params.summary <- read.csv("data-processed/EFGC/EFGC.nonarctic.params.summary.csv")


##### Egg viability (EV) #####
## Load data
data.EV <- read_csv("data-processed/TraitData_EV.csv")


# Subset data
## Arctic species
data.EV.arctic <- subset(data.EV, type == "Arctic")

## Non-Arctic species
data.EV.nonarctic <- subset(data.EV, type == "non-Arctic")

data.EV.arctic.summary <- processTraitData(data.EV.arctic, "EV")
data.EV.nonarctic.summary <- processTraitData(data.EV.nonarctic, "EV")

## Arctic
EV.arctic.predictions.summary <- read.csv("data-processed/EV/EV.arctic.predictions.summary.csv")
EV.arctic.params.summary <- read.csv("data-processed/EV/EV.arctic.params.summary.csv")

# Non-Arctic
EV.nonarctic.predictions.summary <- read.csv("data-processed/EV/EV.nonarctic.predictions.summary.csv")
EV.nonarctic.params.summary <- read.csv("data-processed/EV/EV.nonarctic.params.summary.csv")


##### Larval-to-adult survival (pLA) #####
data.pLA <- read_csv("data-processed/TraitData_pLA.csv")

# Subset data
## Arctic species
data.pLA.arctic <- subset(data.pLA, type == "Arctic")

## Non-Arctic species
data.pLA.nonarctic <- subset(data.pLA, type == "non-Arctic")

data.pLA.arctic.summary <- processTraitData(data.pLA.arctic, "pLA")
data.pLA.nonarctic.summary <- processTraitData(data.pLA.nonarctic, "pLA")

## Arctic
pLA.arctic.predictions.summary <- read.csv("data-processed/pLA/pLA.arctic.predictions.summary.csv")
pLA.arctic.params.summary <- read.csv("data-processed/pLA/pLA.arctic.params.summary.csv")

# Non-Arctic
pLA.nonarctic.predictions.summary <- read.csv("data-processed/pLA/pLA.nonarctic.predictions.summary.csv")
pLA.nonarctic.params.summary <- read.csv("data-processed/pLA/pLA.nonarctic.params.summary.csv")


##### Mosquito development rate (MDR) #####
## Load data
data.MDR <- read_csv("data-processed/TraitData_MDR.csv")

## Convert development time (1/MDR) to development rate (MDR)
data.MDR <- data.MDR %>% 
  mutate(trait = ifelse(trait_name == "1/MDR", 1/trait, trait)) %>% 
  mutate(trait_name = "MDR") 

# Subset data
## Arctic species
data.MDR.arctic <- subset(data.MDR, type == "Arctic")

## Non-Arctic species
data.MDR.nonarctic <- subset(data.MDR, type == "non-Arctic")

data.MDR.arctic.summary <- processTraitData(data.MDR.arctic, "MDR")
data.MDR.nonarctic.summary <- processTraitData(data.MDR.nonarctic, "MDR")

## Arctic
MDR.arctic.predictions.summary <- read.csv("data-processed/MDR/MDR.arctic.predictions.summary.csv")
MDR.arctic.params.summary <- read.csv("data-processed/MDR/MDR.arctic.params.summary.csv")

# Non-Arctic
MDR.nonarctic.predictions.summary <- read.csv("data-processed/MDR/MDR.nonarctic.predictions.summary.csv")
MDR.nonarctic.params.summary <- read.csv("data-processed/MDR/MDR.nonarctic.params.summary.csv")




#  2. Plot panels for each trait -----------------------------------------------

##### biting rate (a) #####
plot.a <- a.alldata.predictions.summary %>% 
  ggplot() +
  geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI), fill = "#E69F00", alpha = 0.5) +
  geom_line(aes(x = temperature, y = median), color = "#E69F00", linewidth = 1) +
  
  # Arctic data
  geom_pointrange(data = data.a.alldata.summary, 
                  aes(x = temp, ymin = mean - std_error, ymax = mean + std_error, y = mean),
                  size = 0.5) +
  # Customize the axes and labels
  scale_x_continuous(limits = c(0, 46)) + 
  labs(title = "A", x = expression(paste("Temperature (", degree, "C)")), 
       y = parse(text = "Rate~(day^-1)")) +
  annotate("text", x = 1, y = 0.44, label = expression(paste(italic("a"))), size = 5) +
  theme_bw()

plot.a



##### Infection efficiency (c) #####
plot.c <- c.arctic.predictions.summary %>% 
  ggplot() +
  geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI), fill = "#009E73", alpha = 0.5) +
  geom_line(aes(x = temperature, y = median), color = "#009E73", linewidth = 1) +
  
  # Customize the axes and labels
  scale_x_continuous(limits = c(0, 46)) + 
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Infection proportion") +
  annotate("text", x = 1, y = 0.97, label = expression(paste(italic("c"))), size = 5) +
  theme_bw()

plot.c


##### Adult lifespan (lf) #####
plot.lf <- lf.arctic.predictions.summary %>% 
  ggplot() +
  geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI), fill = "#0072B2", alpha = 0.5) +
  geom_line(aes(x = temperature, y = median), colour = "#0072B2", linewidth = 1) +

  # Arctic data
  geom_pointrange(data = data.lf.arctic.summary, aes(x = temp, ymin = mean - std_error, ymax = mean + std_error, y = mean), size = 0.5) +
  # Customize the axes and labels
  scale_x_continuous(limits = c(0, 46)) + 
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Time (days)") +
  annotate("text", x = 1, y = 145, label = expression(paste(italic("lf"))), size = 5) +
  theme_bw()

plot.lf


##### Parasite development rate (PDR) #####
plot.PDR <- PDR.arctic.predictions.summary %>% 
  ggplot() +
  geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI), fill = "#CC79A7", alpha = 0.5) +
  geom_line(aes(x = temperature, y = median), colour = "#CC79A7", linewidth = 1) +

  # Arctic data
  geom_pointrange(data = data.PDR.arctic.summary, aes(x = temp, ymin = mean - std_error, ymax = mean + std_error, y = mean), size = 0.5) +
  # Customize the axes and labels
  scale_x_continuous(limits = c(0, 46)) + 
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = parse(text = "Development~rate~(day^-1)")) +
  annotate("text", x = 1, y = 0.22, label = expression(paste(italic("PDR"))), size = 5) +
  theme_bw()

plot.PDR


##### Eggs per female per gonotrophic cycle (EFGC) #####
plot.EFGC <- EFGC.arctic.predictions.summary %>% 
  ggplot() +
  geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI), fill = "#56B4E9", alpha = 0.5) +
  geom_line(aes(x = temperature, y = median), colour = "#56B4E9", linewidth = 1) +
  
  # Customize the axes and labels
  scale_x_continuous(limits = c(0, 46)) + 
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Eggs per female per gonotrophic cycle") +
  annotate("text", x = 2, y = 68, label = expression(paste(italic("EFGC"))), size = 5) +
  theme_bw()

plot.EFGC


#####  Egg viability (EV) ##### 
plot.EV <- EV.arctic.predictions.summary %>% 
  ggplot() +
  geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI), fill = "#F5C710", alpha = 0.5) +
  geom_line(aes(x = temperature, y = median), colour = "#F5C710", linewidth = 1) +

  # Arctic data
  geom_pointrange(data = data.EV.arctic.summary, aes(x = temp, ymin = mean - std_error, ymax = mean + std_error, y = mean), size = 0.5) +
  # Customize the axes and labels
  scale_x_continuous(limits = c(0, 46)) + 
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Proportion hatching") +
  annotate("text", x = 1, y = 0.98, label = expression(paste(italic("EV"))), size = 5) +
  theme_bw()

plot.EV


##### Larval-to-adult survival (pLA) ##### 
plot.pLA <- pLA.arctic.predictions.summary %>% 
  ggplot() +
  geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI), fill = "#999999", alpha = 0.5) +
  geom_line(aes(x = temperature, y = median), colour = "#999999", linewidth = 1) +
  
  # Arctic data
  geom_pointrange(data = data.pLA.arctic.summary, aes(x = temp, ymin = mean - std_error, ymax = mean + std_error, y = mean), size = 0.5) +
  # Customize the axes and labels
  scale_x_continuous(limits = c(0, 46)) + 
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Survival probability)") +
  annotate("text", x = 1, y = 0.85, label = expression(paste(italic("pLA"))), size = 5) +
  theme_bw()

plot.pLA


##### Mosquito development rate (MDR) #####
plot.MDR <- MDR.arctic.predictions.summary %>% 
  ggplot() +
  geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI), fill = "#D55E00", alpha = 0.5) +
  geom_line(aes(x = temperature, y = median), colour = "#D55E00", linewidth = 1) +
  
  # Arctic data
  geom_pointrange(data = data.MDR.arctic.summary, aes(x = temp, ymin = mean - std_error, ymax = mean + std_error, y = mean), size = 0.5) +
  # Customize the axes and labels
  scale_x_continuous(limits = c(0, 46)) + 
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = parse(text = "Development~rate~(day^-1)")) +
  annotate("text", x = 1, y = 0.16, label = expression(paste(italic("MDR"))), size = 5) +
  theme_bw()

plot.MDR


##### Plot all traits #####
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
  labs(title = "B",
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

plot.all <- ggarrange(plot.traits, plot.traits.scaled,
                      nrow = 2, heights = c(3,2)) + bgcolor("white")
plot.all

ggsave("figures/trait.TPCs.summary.png", plot.all, width = 15, height = 18)
