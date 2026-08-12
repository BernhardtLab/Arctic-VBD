## Lilian Chan, University of Guelph
## Arctic vector-borne disease transmission suitability model
##
## Purpose: plotting trait data and TPC fits
## 
## Table of content:
##    0. Set-up workspace
##    1. Load data and model output
##  	2. Plot panels for each trait
##    3. Summary table for TPC parameters
##
##
## Outputs: 
## figures/Fig3-trait.TPCs.png -
##     Main text figure 4

# 0. Set-up workspace ----------------------------------------------------------

library(tidyverse)
library(readxl)
library(janitor)
library(ggsci)
library(cowplot)
library(grafify)
library(flextable)

##### Load functions
source("R-scripts/00_Functions.R")



#  1. Load data and model output -----------------------------------------------

##### biting rate (a) #####
## Load data
data.a <- read_csv("data-processed/TraitData_a.csv")

a.alldata.predictions.summary <- read.csv("data-processed/a/a.alldata.predictions.summary.csv")
a.alldata.params.summary <- read.csv("data-processed/a/a.alldata.params.summary.csv")


##### Vector competence (bc) #####
## Load data
data.bc <- read_csv("data-processed/TraitData_bc.csv")

## Arctic
bc.arctic.predictions.summary <- read.csv("data-processed/supplemental-analysis/briere-only/bc.arctic.predictions.summary.csv")
bc.arctic.params.summary <- read.csv("data-processed/supplemental-analysis/briere-only/bc.arctic.params.summary.csv")

# Non-Arctic
bc.nonarctic.predictions.summary <- read.csv("data-processed/supplemental-analysis/briere-only/bc.nonarctic.predictions.summary.csv")
bc.nonarctic.params.summary <- read.csv("data-processed/supplemental-analysis/briere-only/bc.nonarctic.params.summary.csv")
  

##### Adult lifespan (lf) #####
## Load data
data.lf <- read_csv("data-processed/TraitData_lf.csv")

# Subset data
## Arctic species
data.lf.arctic <- subset(data.lf, type == "Arctic")

## Non-Arctic species
data.lf.nonarctic <- subset(data.lf, type == "non-Arctic")

## Arctic
lf.arctic.predictions.summary <- read.csv("data-processed/supplemental-analysis/briere-only/lf.arctic.predictions.summary.csv")
lf.arctic.params.summary <- read.csv("data-processed/supplemental-analysis/briere-only/lf.arctic.params.summary.csv")

# Non-Arctic
lf.nonarctic.predictions.summary <- read.csv("data-processed/supplemental-analysis/briere-only/lf.nonarctic.predictions.summary.csv")
lf.nonarctic.params.summary <- read.csv("data-processed/supplemental-analysis/briere-only/lf.nonarctic.params.summary.csv")


##### Pathogen development rate (PDR) ##### 
## Load data
data.PDR <- read_csv("data-processed/TraitData_PDR.csv")

# Subset data
## Arctic species
data.PDR.arctic <- subset(data.PDR, type == "Arctic")

## Non-Arctic species
data.PDR.nonarctic <- subset(data.PDR, type == "non-Arctic")


## Arctic
PDR.arctic.predictions.summary <- read.csv("data-processed/PDR/PDR.arctic.predictions.summary.csv")
PDR.arctic.params.summary <- read.csv("data-processed/PDR/PDR.arctic.params.summary.csv")

# Non-Arctic
PDR.nonarctic.predictions.summary <- read.csv("data-processed/PDR/PDR.nonarctic.predictions.summary.csv")
PDR.nonarctic.params.summary <- read.csv("data-processed/PDR/PDR.nonarctic.params.summary.csv")


##### Eggs per female per gonotrophic cycle (EFGC) #####
## Load data
data.EFGC <- read_csv("data-processed/TraitData_EFGC.csv")

EFGC.alldata.predictions.summary <- read.csv("data-processed/supplemental-analysis/briere-only/EFGC.alldata.predictions.summary.csv")
EFGC.alldata.params.summary <- read.csv("data-processed/supplemental-analysis/briere-only/EFGC.alldata.params.summary.csv")



##### Egg viability (EV) #####
## Load data
data.EV <- read_csv("data-processed/TraitData_EV.csv")

# Subset data
## Arctic species
data.EV.arctic <- subset(data.EV, type == "Arctic")

## Non-Arctic species
data.EV.nonarctic <- subset(data.EV, type == "non-Arctic")

## Arctic
EV.arctic.predictions.summary <- read.csv("data-processed/supplemental-analysis/briere-only/EV.arctic.predictions.summary.csv")
EV.arctic.params.summary <- read.csv("data-processed/supplemental-analysis/briere-only/EV.arctic.params.summary.csv")

# Non-Arctic
EV.nonarctic.predictions.summary <- read.csv("data-processed/supplemental-analysis/briere-only/EV.nonarctic.predictions.summary.csv")
EV.nonarctic.params.summary <- read.csv("data-processed/supplemental-analysis/briere-only/EV.nonarctic.params.summary.csv")


##### Larval-to-adult survival (pLA) #####
data.pLA <- read_csv("data-processed/TraitData_pLA.csv")

# Subset data
## Arctic species
data.pLA.arctic <- subset(data.pLA, type == "Arctic")

## Non-Arctic species
data.pLA.nonarctic <- subset(data.pLA, type == "non-Arctic")


## Arctic
pLA.arctic.predictions.summary <- read.csv("data-processed/supplemental-analysis/briere-only/pLA.arctic.predictions.summary.csv")
pLA.arctic.params.summary <- read.csv("data-processed/supplemental-analysis/briere-only/pLA.arctic.params.summary.csv")

# Non-Arctic
pLA.nonarctic.predictions.summary <- read.csv("data-processed/supplemental-analysis/briere-only/pLA.nonarctic.predictions.summary.csv")
pLA.nonarctic.params.summary <- read.csv("data-processed/supplemental-analysis/briere-only/pLA.nonarctic.params.summary.csv")


##### Mosquito development rate (MDR) #####
## Load data
data.MDR <- read_csv("data-processed/TraitData_MDR.csv")

# Subset data
## Arctic species
data.MDR.arctic <- subset(data.MDR, type == "Arctic")

## Non-Arctic species
data.MDR.nonarctic <- subset(data.MDR, type == "non-Arctic")


## Arctic
MDR.arctic.predictions.summary <- read.csv("data-processed/supplemental-analysis/briere-only/MDR.arctic.predictions.summary.csv")
MDR.arctic.params.summary <- read.csv("data-processed/supplemental-analysis/briere-only/MDR.arctic.params.summary.csv")

# Non-Arctic
MDR.nonarctic.predictions.summary <- read.csv("data-processed/supplemental-analysis/briere-only/MDR.nonarctic.predictions.summary.csv")
MDR.nonarctic.params.summary <- read.csv("data-processed/supplemental-analysis/briere-only/MDR.nonarctic.params.summary.csv")




#  2. Plot panels for each trait -----------------------------------------------

##### biting rate (a) #####
plot.a <- a.alldata.predictions.summary %>% 
  ggplot() +
  geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI), fill = "#E69F00", alpha = 0.5) +
  geom_line(aes(x = temperature, y = median), color = "#E69F00", linewidth = 1) +
  
  # Arctic data
  geom_point(data = data.a, aes(x = temp, y = trait, colour = type), size = 2) +
  # Customize the axes and labels
  scale_x_continuous(limits = c(0, 46)) + 
  labs(title = expression(paste("Biting Rate (",italic(a),")")),
       x = expression(paste("Temperature (", degree, "C)")), 
       y = parse(text = "Rate~(day^-1)")) +
  scale_colour_manual(values = c("Arctic" = "black", "non-Arctic" = "azure4"),
                      name = "Dataset"
                      ) +
  theme_bw() +
  theme(title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))

plot.a



##### Vector competence (bc) #####
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
  theme(title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16))

plot.bc


##### Adult lifespan (lf) #####
plot.lf <- lf.arctic.predictions.summary %>% 
  ggplot() +
  geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI), fill = "#0072B2", alpha = 0.5) +
  geom_line(aes(x = temperature, y = median), colour = "#0072B2", linewidth = 1) +

  # Arctic data
  geom_point(data = data.lf.arctic, aes(x = temp, y = trait), size = 2) +  # Customize the axes and labels
  scale_x_continuous(limits = c(0, 46)) + 
  labs(title = expression(paste("Adult Lifespan (",italic(lf),")")),
       x = expression(paste("Temperature (", degree, "C)")), 
       y = "Time (days)") +
  theme_bw() +
  theme(title = element_text(size = 12),
        legend.position="none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16))

plot.lf


##### pathogen development rate (PDR) #####
plot.PDR <- PDR.arctic.predictions.summary %>% 
  ggplot() +
  geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI), fill = "#CC79A7", alpha = 0.5) +
  geom_line(aes(x = temperature, y = median), colour = "#CC79A7", linewidth = 1) +

  # Arctic data
  geom_point(data = data.PDR.arctic, aes(x = temp, y = trait), size = 2) +  # Customize the axes and labels  # Customize the axes and labels
  scale_x_continuous(limits = c(0, 46)) + 
  labs(title = expression(paste("Pathogen Development Rate (",italic(PDR),")")),
       x = expression(paste("Temperature (", degree, "C)")), 
       y = parse(text = "Rate~(day^-1)")) +
  theme_bw() +
  theme(title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16))

plot.PDR


##### Eggs per female per gonotrophic cycle (EFGC) #####
plot.EFGC <- EFGC.alldata.predictions.summary %>% 
  ggplot() +
  geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI), fill = "#56B4E9", alpha = 0.5) +
  geom_line(aes(x = temperature, y = median), colour = "#56B4E9", linewidth = 1) +
  
  # Arctic data
  geom_point(data = data.EFGC, 
                  aes(x = temp, y = trait, colour = type), size = 2) +
  # Customize the axes and labels
  scale_x_continuous(limits = c(0, 46)) + 
  labs(title = expression(paste("Eggs per Female \nper Gonotrophic Cycle (",italic(EFGC),")")),
       x = expression(paste("Temperature (", degree, "C)")), 
       y = "Eggs") +
  scale_colour_manual(values = c("Arctic" = "black", "non-Arctic" = "azure4")) +
  theme_bw() +
  theme(title = element_text(size = 12),
        legend.position="none",
        plot.margin = margin(20,5.5,5.5,5.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16)
        )

plot.EFGC


#####  Egg viability (EV) ##### 
plot.EV <- EV.arctic.predictions.summary %>% 
  ggplot() +
  geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI), fill = "#F5C710", alpha = 0.5) +
  geom_line(aes(x = temperature, y = median), colour = "#F5C710", linewidth = 1) +

  # Arctic data
  geom_point(data = data.EV.arctic, aes(x = temp, y = trait), size = 2) +  # Customize the axes and labels  # Customize the axes and labels
  # Customize the axes and labels
  scale_x_continuous(limits = c(0, 46)) + 
  labs(title = expression(paste("Egg Viability (",italic(EV),")")),
       x = expression(paste("Temperature (", degree, "C)")), 
       y = "Proportion hatching") +
  theme_bw() +
  theme(title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16))

plot.EV


##### Larval-to-adult survival (pLA) ##### 
plot.pLA <- pLA.arctic.predictions.summary %>% 
  ggplot() +
  geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI), fill = "#999999", alpha = 0.5) +
  geom_line(aes(x = temperature, y = median), colour = "#999999", linewidth = 1) +
  
  # Arctic data
  geom_point(data = data.pLA.arctic, aes(x = temp, y = trait), size = 2) +  # Customize the axes and labels  # Customize the axes and labels
  # Customize the axes and labels
  scale_x_continuous(limits = c(0, 46)) + 
  labs(title = expression(paste("Larval-to-Adult Survival (",italic(pLA),")")),
       x = expression(paste("Temperature (", degree, "C)")), 
       y = "Survival probability") +
  theme_bw() +
  theme(title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16))

plot.pLA


##### Mosquito development rate (MDR) #####
plot.MDR <- MDR.arctic.predictions.summary %>% 
  ggplot() +
  geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI), fill = "#D55E00", alpha = 0.5) +
  geom_line(aes(x = temperature, y = median), colour = "#D55E00", linewidth = 1) +
  
  # Arctic data
  geom_point(data = data.MDR.arctic, aes(x = temp, y = trait), size = 2) +  # Customize the axes and labels  # Customize the axes and labels
  # Customize the axes and labels
  scale_x_continuous(limits = c(0, 46)) + 
  labs(title = expression(paste("Mosquito Development Rate (",italic(MDR),")")),
       x = expression(paste("Temperature (", degree, "C)")), 
       y = parse(text = "Rate~(day^-1)")) +
  theme_bw() +
  theme(title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16))

plot.MDR


##### Plot all traits #####
# extract the legend from biting rate subplot
legend <- get_legend(plot.a)

# create an empty panel and place legend in top-left corner of the panel
legend_panel <- ggdraw() +
  draw_grob(legend, x = -0.2, y = 0)

legend_panel

plot.traits <- plot_grid(plot.pLA, plot.MDR, plot.lf, 
                         plot.PDR, plot.EV, plot.bc,
                         plot.a + theme(legend.position="none"), plot.EFGC, legend_panel, 
                         ncol = 3,
                         align = "hv",
                         labels = LETTERS[1:8])
plot.traits





#### Compare the position of TPCs curves and suitability along temperature gradient #####
prediction.summary <- bind_rows(a.alldata.predictions.summary, 
                                bc.arctic.predictions.summary, 
                                lf.arctic.predictions.summary, 
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
prediction.S <- read_csv("data-processed/supplemental-analysis/briere-only/S.predictions.summary.csv")


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
                      labels = c("Suitability (S)",
                                 "Biting rate (a)", 
                                 "Vector competence (bc)", 
                                 "Adult lifespan (lf)", 
                                 "Pathogen development\nrate (PDR)", 
                                 "Eggs per gonotrophic\ncycle (EFGC)",  
                                 "Egg viability (EV)", 
                                 "Larval-to-adult\nsurvival (pLA)", 
                                 "Mosquito development\nrate (MDR)")) + 
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 18),,
        legend.text = element_text(size = 12),
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

ggsave("figures/supplemental-analysis/briere-only/Fig3-trait.TPCs.png", plot.all, width = 14, height = 12)


