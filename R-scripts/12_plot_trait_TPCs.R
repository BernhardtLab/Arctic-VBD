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
library(cowplot)
library(grafify)

##### Load functions
source("R-scripts/00_Functions.R")



#  1. Load data and model output -----------------------------------------------

##### biting rate (a) #####
## Load data
data.a <- read_csv("data-processed/TraitData_a.csv")

# Process trait data for plotting
data.a.alldata.summary <- data.a %>% 
  group_by(temp, type) %>% 
  summarise(mean = mean(trait),
            std_error = sd(trait)/sqrt(n())) %>% 
  mutate(trait = "a")

a.alldata.predictions.summary <- read.csv("data-processed/a/a.alldata.predictions.summary.csv")
a.alldata.params.summary <- read.csv("data-processed/a/a.alldata.params.summary.csv")


##### Vector competence (bc) #####
## Load data
data.bc <- read_csv("data-processed/TraitData_bc.csv")

# Process trait data for plotting
data.bc.nonarctic.summary <- processTraitData(data.bc, "bc")

## Arctic
bc.arctic.predictions.summary <- read.csv("data-processed/bc/bc.arctic.predictions.summary.csv")
bc.arctic.params.summary <- read.csv("data-processed/bc/bc.arctic.params.summary.csv")

# Non-Arctic
bc.nonarctic.predictions.summary <- read.csv("data-processed/bc/bc.nonarctic.predictions.summary.csv")
bc.nonarctic.params.summary <- read.csv("data-processed/bc/bc.nonarctic.params.summary.csv")


##### Adult lifespan (lf) #####
## Load data
data.lf <- read_csv("data-processed/TraitData_lf.csv")

# Process trait data for plotting
data.lf.alldata.summary <- data.lf %>% 
  group_by(temp, type) %>% 
  summarise(mean = mean(trait),
            std_error = sd(trait)/sqrt(n())) %>% 
  mutate(trait = "lf")

lf.alldata.predictions.summary <- read.csv("data-processed/lf/lf.alldata.predictions.summary.csv")
lf.alldata.params.summary <- read.csv("data-processed/lf/lf.alldata.params.summary.csv")



##### Parasite development rate (PDR) ##### 
## Load data
data.PDR <- read_csv("data-processed/TraitData_PDR.csv")

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

# Process trait data for plotting
data.EFGC.alldata.summary <- data.EFGC %>% 
  group_by(temp, type) %>% 
  summarise(mean = mean(trait),
            std_error = sd(trait)/sqrt(n())) %>% 
  mutate(trait = "EFGC")

EFGC.alldata.predictions.summary <- read.csv("data-processed/EFGC/EFGC.alldata.predictions.summary.csv")
EFGC.alldata.params.summary <- read.csv("data-processed/EFGC/EFGC.alldata.params.summary.csv")


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
                  aes(x = temp, ymin = mean - std_error, ymax = mean + std_error,
                      y = mean, colour = type),
                  size = 0.5) +
  # Customize the axes and labels
  scale_x_continuous(limits = c(0, 46)) + 
  labs(title = expression(paste("Biting Rate (",italic(a),")")),
       x = expression(paste("Temperature (", degree, "C)")), 
       y = parse(text = "Rate~(day^-1)")) +
  scale_colour_manual(values = c("Arctic" = "black", "non-Arctic" = "azure4"),
                      name = "Dataset"
                      ) +
  theme_bw() +
  theme(title = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))

plot.a



##### Infection efficiency (c) #####
plot.bc <- bc.arctic.predictions.summary %>% 
  ggplot() +
  geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI), fill = "#009E73", alpha = 0.5) +
  geom_line(aes(x = temperature, y = median), color = "#009E73", linewidth = 1) +
  
  # Customize the axes and labels
  scale_x_continuous(limits = c(0, 46)) + 
  labs(title = expression(paste("Vector Competence (",italic(bc),")")),
       x = expression(paste("Temperature (", degree, "C)")), 
       y = "Proportion") +
  theme_bw() +
  theme(title = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 18))

plot.bc


##### Adult lifespan (lf) #####
plot.lf <- lf.alldata.predictions.summary %>% 
  ggplot() +
  geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI), fill = "#0072B2", alpha = 0.5) +
  geom_line(aes(x = temperature, y = median), colour = "#0072B2", linewidth = 1) +

  # Arctic data
  geom_pointrange(data = data.lf.alldata.summary, 
                  aes(x = temp, ymin = mean - std_error, ymax = mean + std_error,
                      y = mean, colour = type), size = 0.5) +
  # Customize the axes and labels
  scale_x_continuous(limits = c(0, 46)) + 
  labs(title = expression(paste("Adult Lifespan (",italic(lf),")")),
       x = expression(paste("Temperature (", degree, "C)")), 
       y = "Time (days)") +
  scale_colour_manual(values = c("Arctic" = "black", "non-Arctic" = "azure4")) +
  theme_bw() +
  theme(title = element_text(size = 14),
        legend.position="none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 18))

plot.lf


##### pathogen development rate (PDR) #####
plot.PDR <- PDR.arctic.predictions.summary %>% 
  ggplot() +
  geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI), fill = "#CC79A7", alpha = 0.5) +
  geom_line(aes(x = temperature, y = median), colour = "#CC79A7", linewidth = 1) +

  # Arctic data
  geom_pointrange(data = data.PDR.arctic.summary, aes(x = temp, ymin = mean - std_error, ymax = mean + std_error, y = mean), size = 0.5) +
  # Customize the axes and labels
  scale_x_continuous(limits = c(0, 46)) + 
  labs(title = expression(paste("Pathogen Development Rate (",italic(PDR),")")),
       x = expression(paste("Temperature (", degree, "C)")), 
       y = parse(text = "Development~rate~(day^-1)")) +
  theme_bw() +
  theme(title = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 18))

plot.PDR


##### Eggs per female per gonotrophic cycle (EFGC) #####
plot.EFGC <- EFGC.alldata.predictions.summary %>% 
  ggplot() +
  geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI), fill = "#56B4E9", alpha = 0.5) +
  geom_line(aes(x = temperature, y = median), colour = "#56B4E9", linewidth = 1) +
  
  # Arctic data
  geom_pointrange(data = data.EFGC.alldata.summary, 
                  aes(x = temp, ymin = mean - std_error, ymax = mean + std_error,
                      y = mean, colour = type), size = 0.5) +
  # Customize the axes and labels
  scale_x_continuous(limits = c(0, 46)) + 
  labs(title = expression(paste("Eggs per Female \nper Gonotrophic Cycle (",italic(EFGC),")")),
       x = expression(paste("Temperature (", degree, "C)")), 
       y = "Eggs") +
  scale_colour_manual(values = c("Arctic" = "black", "non-Arctic" = "azure4")) +
  theme_bw() +
  theme(title = element_text(size = 14),
        legend.position="none",
        plot.margin = margin(20,5.5,5.5,5.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 18)
        )

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
  labs(title = expression(paste("Egg Viability (",italic(EV),")")),
       x = expression(paste("Temperature (", degree, "C)")), 
       y = "Proportion hatching") +
  theme_bw() +
  theme(title = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 18))

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
  labs(title = expression(paste("Larval-to-Adult Survival (",italic(pLA),")")),
       x = expression(paste("Temperature (", degree, "C)")), 
       y = "Survival probability") +
  theme_bw() +
  theme(title = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 18))

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
  labs(title = expression(paste("Mosquito Development Rate (",italic(MDR),")")),
       x = expression(paste("Temperature (", degree, "C)")), 
       y = parse(text = "Development~rate~(day^-1)")) +
  theme_bw() +
  theme(title = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 18))

plot.MDR


##### Plot all traits #####
# extract the legend from biting rate subplot
legend <- get_legend(plot.a)

# create an empty panel and place legend in top-left corner of the panel
legend_panel <- ggdraw() +
  draw_grob(legend, x = -0.2, y = 0)

legend_panel

plot.traits <- plot_grid(plot.a + theme(legend.position="none"), plot.lf, plot.EFGC,
                         plot.EV, plot.pLA, plot.MDR, 
                         plot.PDR, plot.bc, legend_panel, 
                         ncol = 3,
                         labels = LETTERS[1:8])
plot.traits





#### Compare the position of TPCs curves and suitability along temperature gradient #####
prediction.summary <- bind_rows(a.alldata.predictions.summary, 
                                bc.arctic.predictions.summary, 
                                lf.alldata.predictions.summary, 
                                PDR.arctic.predictions.summary,
                                EFGC.alldata.predictions.summary, 
                                EV.arctic.predictions.summary,
                                pLA.arctic.predictions.summary, 
                                MDR.arctic.predictions.summary)


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
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Trait values (scaled)") +
  scale_colour_manual(values = c("S" = "#000000", "a" = "#E69F00", 
                                 "bc" = "#009E73","lf" = "#0072B2", 
                                 "PDR" = "#CC79A7", "EFGC" = "#56B4E9", 
                                 "EV" = "#F5C710", "pLA" = "#999999", 
                                 "MDR" = "#D55E00"),
                      name = element_blank(), # No legend title
                      breaks = c("S", "a", "bc", "lf", "PDR", "EFGC", "EV", "pLA", "MDR"),
                      labels = c("Suitability", "a", "bc", "lf", "PDR", "EFGC",  "EV", "pLA", "MDR")) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 18),,
        legend.text = element_text(size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

plot.traits.scaled


plot.all <- plot_grid(plot.traits, plot.traits.scaled,
                      ncol = 1,
                      labels = c("", LETTERS[9]), # Label only second plot
                      rel_heights = c(5,2)
                      ) +
  theme(panel.background = element_rect(fill = "white", color = NA))

plot.all

ggsave("figures/trait.TPCs.summary.png", plot.all, width = 12, height = 12)
