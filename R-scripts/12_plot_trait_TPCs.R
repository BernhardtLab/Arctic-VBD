## Lilian Chan, University of Guelph
## Arctic vector-borne disease transmission suitability model
##
## Purpose: plotting trait data and TPC fits
## 
## Table of content:
##    0. Set-up workspace
##    1. Load data and model output
##  	2. Plot TPCs
##			A. Arctic TPCs
##			B. Non-Arctic TPCs
##			C. Uniform priors vs data-informed priors
##    3. Summary table for TPC parameters
##    4. Summary table for DIC
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

# Non-Arctic
bc.nonarctic.predictions.summary <- read.csv("data-processed/bc/bc.nonarctic.predictions.summary.csv")
bc.nonarctic.params.summary <- read.csv("data-processed/bc/bc.nonarctic.params.summary.csv")


##### Adult lifespan (lf) #####
## Load data
data.lf <- read_csv("data-processed/TraitData_lf.csv")

# Subset data
## Arctic species
data.lf.arctic <- subset(data.lf, type == "Arctic")

## Non-Arctic species
data.lf.nonarctic <- subset(data.lf, type == "non-Arctic")

## Arctic
lf.arctic.predictions.summary <- read.csv("data-processed/lf/lf.arctic.predictions.summary.csv")
lf.arctic.params.summary <- read.csv("data-processed/lf/lf.arctic.params.summary.csv")

# Non-Arctic
lf.nonarctic.predictions.summary <- read.csv("data-processed/lf/lf.nonarctic.predictions.summary.csv")
lf.nonarctic.params.summary <- read.csv("data-processed/lf/lf.nonarctic.params.summary.csv")


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


## Arctic
MDR.arctic.predictions.summary <- read.csv("data-processed/MDR/MDR.arctic.predictions.summary.csv")
MDR.arctic.params.summary <- read.csv("data-processed/MDR/MDR.arctic.params.summary.csv")

# Non-Arctic
MDR.nonarctic.predictions.summary <- read.csv("data-processed/MDR/MDR.nonarctic.predictions.summary.csv")
MDR.nonarctic.params.summary <- read.csv("data-processed/MDR/MDR.nonarctic.params.summary.csv")




#  2. Plot panels for each trait -----------------------------------------------

## 2A. Arctic TPCs -------------------------------------------------------------

##### biting rate (a) #####
plot.a <- a.alldata.predictions.summary %>% 
  ggplot() +
  geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI), fill = "#E69F00", alpha = 0.5) +
  geom_line(aes(x = temperature, y = median), color = "#E69F00", linewidth = 1) +
  
  # data
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
plot.bc <- bc.nonarctic.predictions.summary %>% 
  ggplot() +
  geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI), fill = "#009E73", alpha = 0.5) +
  geom_line(aes(x = temperature, y = median), color = "#009E73", linewidth = 1) +
  
  # data
  geom_point(data = data.bc, aes(x = temp, y = trait), colour = "azure4", size = 2) +
  
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
                         plot.PDR, plot.EV, legend_panel,
                         plot.a + theme(legend.position="none"), plot.EFGC, plot.bc, 
                         ncol = 3,
                         align = "hv",
                         labels = c(LETTERS[1:5], NA, LETTERS[6:8]))
plot.traits





#### Compare the position of TPCs curves and suitability along temperature gradient #####
prediction.summary <- bind_rows(a.alldata.predictions.summary, 
                                bc.nonarctic.predictions.summary, 
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
prediction.S <- read_csv("data-processed/suitability/S.predictions.summary.csv")


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

ggsave("figures/Fig3-trait.TPCs.png", plot.all, width = 14, height = 12)


## 2B. non-Arctic TPCs -------------------------------------------------------------

##### Adult lifespan (lf) #####
plot.lf.nonarctic <- lf.nonarctic.predictions.summary %>% 
  ggplot() +
  geom_line(aes(x = temperature, y = lowerCI), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
  geom_line(aes(x = temperature, y = upperCI), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
  geom_line(aes(x = temperature, y = median), colour = "grey4", linewidth = 1) +
  
  # Arctic data
  geom_point(data = data.lf.nonarctic, aes(x = temp, y = trait), size = 2, shape = 1) +  # Customize the axes and labels
  scale_x_continuous(limits = c(0, 46)) + 
  labs(title = expression(paste("Adult Lifespan (",italic(lf),")")),
       x = expression(paste("Temperature (", degree, "C)")), 
       y = "Time (days)") +
  theme_bw() +
  theme(title = element_text(size = 12),
        legend.position="none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16))

plot.lf.nonarctic 


##### pathogen development rate (PDR) #####
plot.PDR.nonarctic <- PDR.nonarctic.predictions.summary %>% 
  ggplot() +
  geom_line(aes(x = temperature, y = lowerCI), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
  geom_line(aes(x = temperature, y = upperCI), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
  geom_line(aes(x = temperature, y = median), colour = "grey4", linewidth = 1) +
  
  # Arctic data
  geom_point(data = data.PDR.nonarctic, aes(x = temp, y = trait), size = 2, shape = 1) +  # Customize the axes and labels  # Customize the axes and labels
  scale_x_continuous(limits = c(0, 46)) + 
  labs(title = expression(paste("Pathogen Development Rate (",italic(PDR),")")),
       x = expression(paste("Temperature (", degree, "C)")), 
       y = parse(text = "Rate~(day^-1)")) +
  theme_bw() +
  theme(title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16))

plot.PDR.nonarctic



#####  Egg viability (EV) ##### 
plot.EV.nonarctic <- EV.nonarctic.predictions.summary %>% 
  ggplot() +
  geom_line(aes(x = temperature, y = lowerCI), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
  geom_line(aes(x = temperature, y = upperCI), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
  geom_line(aes(x = temperature, y = median), colour = "grey4", linewidth = 1) +
  
  # Arctic data
  geom_point(data = data.EV.nonarctic, aes(x = temp, y = trait), size = 2, shape = 1) +  # Customize the axes and labels  # Customize the axes and labels
  # Customize the axes and labels
  scale_x_continuous(limits = c(0, 46)) + 
  labs(title = expression(paste("Egg Viability (",italic(EV),")")),
       x = expression(paste("Temperature (", degree, "C)")), 
       y = "Proportion hatching") +
  theme_bw() +
  theme(title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16))

plot.EV.nonarctic


##### Larval-to-adult survival (pLA) ##### 
plot.pLA.nonarctic <- pLA.nonarctic.predictions.summary %>% 
  ggplot() +
  geom_line(aes(x = temperature, y = lowerCI), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
  geom_line(aes(x = temperature, y = upperCI), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
  geom_line(aes(x = temperature, y = median), colour = "grey4", linewidth = 1) +
  
  # Arctic data
  geom_point(data = data.pLA.nonarctic, aes(x = temp, y = trait), size = 2, shape = 1) +  # Customize the axes and labels  # Customize the axes and labels
  # Customize the axes and labels
  scale_x_continuous(limits = c(0, 46)) + 
  labs(title = expression(paste("Larval-to-Adult Survival (",italic(pLA),")")),
       x = expression(paste("Temperature (", degree, "C)")), 
       y = "Survival probability") +
  theme_bw() +
  theme(title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16))

plot.pLA.nonarctic


##### Mosquito development rate (MDR) #####
plot.MDR.nonarctic <- MDR.nonarctic.predictions.summary %>% 
  ggplot() +
  geom_line(aes(x = temperature, y = lowerCI), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
  geom_line(aes(x = temperature, y = upperCI), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
  geom_line(aes(x = temperature, y = median), colour = "grey4", linewidth = 1) +
  
  # Arctic data
  geom_point(data = data.MDR.nonarctic, aes(x = temp, y = trait), size = 2, shape = 1) +  # Customize the axes and labels  # Customize the axes and labels
  # Customize the axes and labels
  scale_x_continuous(limits = c(0, 46)) + 
  labs(title = expression(paste("Mosquito Development Rate (",italic(MDR),")")),
       x = expression(paste("Temperature (", degree, "C)")), 
       y = parse(text = "Rate~(day^-1)")) +
  theme_bw() +
  theme(title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16))

plot.MDR.nonarctic


plot.traits.nonarctic <- plot_grid(plot.pLA.nonarctic, plot.MDR.nonarctic, plot.lf.nonarctic, 
                         plot.PDR.nonarctic, plot.EV.nonarctic,
                         ncol = 3,
                         align = "hv",
                         labels = LETTERS[1:5])
plot.traits.nonarctic


ggsave("figures/FigS1-trait.TPCs.nonarctic.png", plot.traits.nonarctic, width = 14, height = 12)


## 2C. Uniform priors vs data-informed priors ----------------------------------
Temp.xs <- seq(0, 45, 0.1)


##### Larval-to-adult survival (pLA) ##### 

# Briere; uniform priors
# Load TPCs fitted using uniform priors
load("R-scripts/R2jags-objects/all-mods/pLA.arctic.bri.uni.Rdata")

# Extract model prediction
df.pLA.arctic.bri.uni <- data.frame(pLA.arctic.bri.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)


# Plot
plot.pLA.arctic.bri.uni <- df.pLA.arctic.bri.uni %>% 
  ggplot(aes(x = temp)) +
  geom_line(aes(y = X2.5.), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = X97.5.), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = X50.), colour = "grey4", linewidth = 1) +
  
  ## data
  geom_point(data = data.pLA.arctic, aes(x = temp, y = trait), size = 2, shape = 1) +  # Customize the axes and labels  # Customize the axes and labels
  
  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Survival probability",
       title = "A) Briere model; uniform priors"
  ) +
  theme_bw()


plot.pLA.arctic.bri.uni


# Quadratic; uniform priors
# Load TPCs fitted using uniform priors
load("R-scripts/R2jags-objects/all-mods/pLA.arctic.quad.uni.Rdata")

# Extract model prediction
df.pLA.arctic.quad.uni <- data.frame(pLA.arctic.quad.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)

# Plot
plot.pLA.arctic.quad.uni <- df.pLA.arctic.quad.uni %>% 
  ggplot(aes(x = temp)) +
  geom_line(aes(y = X2.5.), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = X97.5.), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = X50.), colour = "grey4", linewidth = 1) +
  
  ## data
  geom_point(data = data.pLA.arctic, aes(x = temp, y = trait), size = 2, shape = 1) +  # Customize the axes and labels  # Customize the axes and labels
  
  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Survival probability",
       title = "B) Quadratic model; uniform priors"
  ) +
  theme_bw()


plot.pLA.arctic.quad.uni


# Briere; data-informed priors
# Load TPCs fitted using infform priors
load("R-scripts/R2jags-objects/all-mods/pLA.arctic.bri.inf.Rdata")

# Extract model prediction
df.pLA.arctic.bri.inf <- data.frame(pLA.arctic.bri.inf$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)

# Plot
plot.pLA.arctic.bri.inf <- df.pLA.arctic.bri.inf %>% 
  ggplot(aes(x = temp)) +
  geom_line(aes(y = X2.5.), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = X97.5.), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = X50.), colour = "grey4", linewidth = 1) +
  
  ## data
  geom_point(data = data.pLA.arctic, aes(x = temp, y = trait), size = 2, shape = 1) +  # Customize the axes and labels  # Customize the axes and labels
  
  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Survival probability",
       title = "C) Briere model; data-informed priors"
  ) +
  theme_bw()


plot.pLA.arctic.bri.inf


# Quadratic; data-informed priors
# Load TPCs fitted using infform priors
load("R-scripts/R2jags-objects/all-mods/pLA.arctic.quad.inf.Rdata")

# Extract model prediction
df.pLA.arctic.quad.inf <- data.frame(pLA.arctic.quad.inf$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)

# Plot
plot.pLA.arctic.quad.inf <- df.pLA.arctic.quad.inf %>% 
  ggplot(aes(x = temp)) +
  geom_line(aes(y = X2.5.), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = X97.5.), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = X50.), colour = "grey4", linewidth = 1) +
  
  ## data
  geom_point(data = data.pLA.arctic, aes(x = temp, y = trait), size = 2, shape = 1) +  # Customize the axes and labels  # Customize the axes and labels
  
  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Survival probability",
       title = "D) Quadratic model; data-informed priors"
  ) +
  theme_bw()


plot.pLA.arctic.quad.inf


plot.pLA.arctic <- plot_grid(plot.pLA.arctic.bri.uni, plot.pLA.arctic.quad.uni,
                            plot.pLA.arctic.bri.inf, plot.pLA.arctic.quad.inf, ncol = 2)
plot.pLA.arctic

ggsave("figures/FigS2-pLA.uni.vs.inf.png", plot.pLA.arctic, width = 12, height = 8)



##### Mosquito development rate (MDR) #####
# Briere; uniform priors
# Load TPCs fitted using uniform priors
load("R-scripts/R2jags-objects/all-mods/MDR.arctic.bri.uni.Rdata")

# Extract model prediction
df.MDR.arctic.bri.uni <- data.frame(MDR.arctic.bri.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)


# Plot
plot.MDR.arctic.bri.uni <- df.MDR.arctic.bri.uni %>% 
  ggplot(aes(x = temp)) +
  geom_line(aes(y = X2.5.), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = X97.5.), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = X50.), colour = "grey4", linewidth = 1) +
  
  ## data
  geom_point(data = data.MDR.arctic, aes(x = temp, y = trait), size = 2, shape = 1) +  # Customize the axes and labels  # Customize the axes and labels
  
  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Development rate (days-1)",
       title = "A) Briere model; uniform priors"
  ) +
  theme_bw()


plot.MDR.arctic.bri.uni


# Quadratic; uniform priors
# Load TPCs fitted using uniform priors
load("R-scripts/R2jags-objects/all-mods/MDR.arctic.quad.uni.Rdata")

# Extract model prediction
df.MDR.arctic.quad.uni <- data.frame(MDR.arctic.quad.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)

# Plot
plot.MDR.arctic.quad.uni <- df.MDR.arctic.quad.uni %>% 
  ggplot(aes(x = temp)) +
  geom_line(aes(y = X2.5.), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = X97.5.), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = X50.), colour = "grey4", linewidth = 1) +
  
  ## data
  geom_point(data = data.MDR.arctic, aes(x = temp, y = trait), size = 2, shape = 1) +  # Customize the axes and labels  # Customize the axes and labels
  
  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Development rate (days-1)",
       title = "B) Quadratic model; uniform priors"
  ) +
  theme_bw()


plot.MDR.arctic.quad.uni


# Briere; data-informed priors
# Load TPCs fitted using infform priors
load("R-scripts/R2jags-objects/all-mods/MDR.arctic.bri.inf.Rdata")

# Extract model prediction
df.MDR.arctic.bri.inf <- data.frame(MDR.arctic.bri.inf$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)

# Plot
plot.MDR.arctic.bri.inf <- df.MDR.arctic.bri.inf %>% 
  ggplot(aes(x = temp)) +
  geom_line(aes(y = X2.5.), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = X97.5.), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = X50.), colour = "grey4", linewidth = 1) +
  
  ## data
  geom_point(data = data.MDR.arctic, aes(x = temp, y = trait), size = 2, shape = 1) +  # Customize the axes and labels  # Customize the axes and labels
  
  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Development rate (days-1)",
       title = "C) Briere model; data-informed priors"
  ) +
  theme_bw()


plot.MDR.arctic.bri.inf


# Quadratic; data-informed priors
# Load TPCs fitted using infform priors
load("R-scripts/R2jags-objects/all-mods/MDR.arctic.quad.inf.Rdata")

# Extract model prediction
df.MDR.arctic.quad.inf <- data.frame(MDR.arctic.quad.inf$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)

# Plot
plot.MDR.arctic.quad.inf <- df.MDR.arctic.quad.inf %>% 
  ggplot(aes(x = temp)) +
  geom_line(aes(y = X2.5.), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = X97.5.), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = X50.), colour = "grey4", linewidth = 1) +
  
  ## data
  geom_point(data = data.MDR.arctic, aes(x = temp, y = trait), size = 2, shape = 1) +  # Customize the axes and labels  # Customize the axes and labels
  
  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Development rate (days-1)",
       title = "D) Quadratic model; data-informed priors"
  ) +
  theme_bw()


plot.MDR.arctic.quad.inf


plot.MDR.arctic <- plot_grid(plot.MDR.arctic.bri.uni, plot.MDR.arctic.quad.uni,
                             plot.MDR.arctic.bri.inf, plot.MDR.arctic.quad.inf, ncol = 2)
plot.MDR.arctic

ggsave("figures/FigS3-MDR.uni.vs.inf.png", plot.MDR.arctic, width = 12, height = 8)



##### Adult lifespan (lf) #####

# Briere; uniform priors
# Load TPCs fitted using uniform priors
load("R-scripts/R2jags-objects/all-mods/lf.arctic.bri.uni.Rdata")

# Extract model prediction
df.lf.arctic.bri.uni <- data.frame(lf.arctic.bri.uni$BUGSoutput$summary)

df.lf.arctic.bri.uni <- df.lf.arctic.bri.uni %>% 
  filter(grepl("z.trait.mu.pred.pop", rownames(df.lf.arctic.bri.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)

# Plot
plot.lf.arctic.bri.uni <- df.lf.arctic.bri.uni %>% 
  ggplot(aes(x = temp)) +
  geom_line(aes(y = X2.5.), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = X97.5.), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = X50.), colour = "grey4", linewidth = 1) +
  
  ## data
  geom_point(data = data.lf.arctic, aes(x = temp, y = trait), size = 2, shape = 1) +  # Customize the axes and labels  # Customize the axes and labels
  
  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Time (days)",
       title = "A) Briere model; uniform priors"
  ) +
  theme_bw()


plot.lf.arctic.bri.uni


# Quadratic; uniform priors
# Load TPCs fitted using uniform priors
load("R-scripts/R2jags-objects/all-mods/lf.arctic.quad.uni.Rdata")

# Extract model prediction
df.lf.arctic.quad.uni <- data.frame(lf.arctic.quad.uni$BUGSoutput$summary)

df.lf.arctic.quad.uni <- df.lf.arctic.quad.uni %>% 
  filter(grepl("z.trait.mu.pred.pop", rownames(df.lf.arctic.quad.uni))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)

# Plot
plot.lf.arctic.quad.uni <- df.lf.arctic.quad.uni %>% 
  ggplot(aes(x = temp)) +
  geom_line(aes(y = X2.5.), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = X97.5.), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = X50.), colour = "grey4", linewidth = 1) +
  
  ## data
  geom_point(data = data.lf.arctic, aes(x = temp, y = trait), size = 2, shape = 1) +  # Customize the axes and labels  # Customize the axes and labels

  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Time (days)",
       title = "B) Quadratic model; uniform priors"
  ) +
  theme_bw()


plot.lf.arctic.quad.uni


# Briere; data-informed priors
# Load TPCs fitted using infform priors
load("R-scripts/R2jags-objects/all-mods/lf.arctic.bri.inf.Rdata")

# Extract model prediction
df.lf.arctic.bri.inf <- data.frame(lf.arctic.bri.inf$BUGSoutput$summary)

df.lf.arctic.bri.inf <- df.lf.arctic.bri.inf %>% 
  filter(grepl("z.trait.mu.pred.pop", rownames(df.lf.arctic.bri.inf))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)

# Plot
plot.lf.arctic.bri.inf <- df.lf.arctic.bri.inf %>% 
  ggplot(aes(x = temp)) +
  geom_line(aes(y = X2.5.), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = X97.5.), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = X50.), colour = "grey4", linewidth = 1) +
  
  ## data
  geom_point(data = data.lf.arctic, aes(x = temp, y = trait), size = 2, shape = 1) +  # Customize the axes and labels  # Customize the axes and labels
  
  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Time (days)",
       title = "C) Briere model; data-informed priors"
  ) +
  theme_bw()


plot.lf.arctic.bri.inf


# Quadratic; data-informed priors
# Load TPCs fitted using infform priors
load("R-scripts/R2jags-objects/all-mods/lf.arctic.quad.inf.Rdata")

# Extract model prediction
df.lf.arctic.quad.inf <- data.frame(lf.arctic.quad.inf$BUGSoutput$summary)

df.lf.arctic.quad.inf <- df.lf.arctic.quad.inf %>% 
  filter(grepl("z.trait.mu.pred.pop", rownames(df.lf.arctic.quad.inf))) %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)

# Plot
plot.lf.arctic.quad.inf <- df.lf.arctic.quad.inf %>% 
  ggplot(aes(x = temp)) +
  geom_line(aes(y = X2.5.), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = X97.5.), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = X50.), colour = "grey4", linewidth = 1) +
  
  ## data
  geom_point(data = data.lf.arctic, aes(x = temp, y = trait), size = 2, shape = 1) +  # Customize the axes and labels  # Customize the axes and labels
  
  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Time (days)",
       title = "D) Quadratic model; data-informed priors"
  ) +
  theme_bw()


plot.lf.arctic.quad.inf


plot.lf.arctic <- plot_grid(plot.lf.arctic.bri.uni, plot.lf.arctic.quad.uni,
                            plot.lf.arctic.bri.inf, plot.lf.arctic.quad.inf, ncol = 2)
plot.lf.arctic

ggsave("figures/FigS4-lf.uni.vs.inf.png", plot.lf.arctic, width = 12, height = 8)


##### pathogen development rate (PDR) #####

# Briere; uniform priors
# Load TPCs fitted using uniform priors
load("R-scripts/R2jags-objects/all-mods/PDR.arctic.bri.uni.Rdata")

# Extract model prediction
df.PDR.arctic.bri.uni <- data.frame(PDR.arctic.bri.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)


# Plot
plot.PDR.arctic.bri.uni <- df.PDR.arctic.bri.uni %>% 
  ggplot(aes(x = temp)) +
  geom_line(aes(y = X2.5.), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = X97.5.), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = X50.), colour = "grey4", linewidth = 1) +
  
  ## data
  geom_point(data = data.PDR.arctic, aes(x = temp, y = trait), size = 2, shape = 1) +  # Customize the axes and labels  # Customize the axes and labels
  
  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Development rate (days-1)",
       title = "A) Briere model; uniform priors"
  ) +
  theme_bw()


plot.PDR.arctic.bri.uni


# Quadratic; uniform priors
# Load TPCs fitted using uniform priors
load("R-scripts/R2jags-objects/all-mods/PDR.arctic.quad.uni.Rdata")

# Extract model prediction
df.PDR.arctic.quad.uni <- data.frame(PDR.arctic.quad.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)

# Plot
plot.PDR.arctic.quad.uni <- df.PDR.arctic.quad.uni %>% 
  ggplot(aes(x = temp)) +
  geom_line(aes(y = X2.5.), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = X97.5.), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = X50.), colour = "grey4", linewidth = 1) +
  
  ## data
  geom_point(data = data.PDR.arctic, aes(x = temp, y = trait), size = 2, shape = 1) +  # Customize the axes and labels  # Customize the axes and labels
  
  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Development rate (days-1)",
       title = "B) Quadratic model; uniform priors"
  ) +
  theme_bw()


plot.PDR.arctic.quad.uni


# Briere; data-informed priors
# Load TPCs fitted using infform priors
load("R-scripts/R2jags-objects/all-mods/PDR.arctic.bri.inf.Rdata")

# Extract model prediction
df.PDR.arctic.bri.inf <- data.frame(PDR.arctic.bri.inf$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)

# Plot
plot.PDR.arctic.bri.inf <- df.PDR.arctic.bri.inf %>% 
  ggplot(aes(x = temp)) +
  geom_line(aes(y = X2.5.), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = X97.5.), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = X50.), colour = "grey4", linewidth = 1) +
  
  ## data
  geom_point(data = data.PDR.arctic, aes(x = temp, y = trait), size = 2, shape = 1) +  # Customize the axes and labels  # Customize the axes and labels
  
  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Development rate (days-1)",
       title = "C) Briere model; data-informed priors"
  ) +
  theme_bw()


plot.PDR.arctic.bri.inf


# Quadratic; data-informed priors
# Load TPCs fitted using infform priors
load("R-scripts/R2jags-objects/all-mods/PDR.arctic.quad.inf.Rdata")

# Extract model prediction
df.PDR.arctic.quad.inf <- data.frame(PDR.arctic.quad.inf$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)

# Plot
plot.PDR.arctic.quad.inf <- df.PDR.arctic.quad.inf %>% 
  ggplot(aes(x = temp)) +
  geom_line(aes(y = X2.5.), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = X97.5.), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = X50.), colour = "grey4", linewidth = 1) +
  
  ## data
  geom_point(data = data.PDR.arctic, aes(x = temp, y = trait), size = 2, shape = 1) +  # Customize the axes and labels  # Customize the axes and labels
  
  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Development rate (days-1)",
       title = "D) Quadratic model; data-informed priors"
  ) +
  theme_bw()


plot.PDR.arctic.quad.inf


plot.PDR.arctic <- plot_grid(plot.PDR.arctic.bri.uni, plot.PDR.arctic.quad.uni,
                             plot.PDR.arctic.bri.inf, plot.PDR.arctic.quad.inf, ncol = 2)
plot.PDR.arctic

ggsave("figures/FigS5-PDR.uni.vs.inf.png", plot.PDR.arctic, width = 12, height = 8)




#####  Egg viability (EV) ##### 

# Briere; uniform priors
# Load TPCs fitted using uniform priors
load("R-scripts/R2jags-objects/all-mods/EV.arctic.bri.uni.Rdata")

# Extract model prediction
df.EV.arctic.bri.uni <- data.frame(EV.arctic.bri.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)


# Plot
plot.EV.arctic.bri.uni <- df.EV.arctic.bri.uni %>% 
  ggplot(aes(x = temp)) +
  geom_line(aes(y = X2.5.), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = X97.5.), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = X50.), colour = "grey4", linewidth = 1) +
  
  ## data
  geom_point(data = data.EV.arctic, aes(x = temp, y = trait), size = 2, shape = 1) +  # Customize the axes and labels  # Customize the axes and labels
  
  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Proportion hatching",
       title = "A) Briere model; uniform priors"
  ) +
  theme_bw()


plot.EV.arctic.bri.uni


# Quadratic; uniform priors
# Load TPCs fitted using uniform priors
load("R-scripts/R2jags-objects/all-mods/EV.arctic.quad.uni.Rdata")

# Extract model prediction
df.EV.arctic.quad.uni <- data.frame(EV.arctic.quad.uni$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)

# Plot
plot.EV.arctic.quad.uni <- df.EV.arctic.quad.uni %>% 
  ggplot(aes(x = temp)) +
  geom_line(aes(y = X2.5.), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = X97.5.), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = X50.), colour = "grey4", linewidth = 1) +
  
  ## data
  geom_point(data = data.EV.arctic, aes(x = temp, y = trait), size = 2, shape = 1) +  # Customize the axes and labels  # Customize the axes and labels
  
  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Proportion hatching",
       title = "B) Quadratic model; uniform priors"
  ) +
  theme_bw()


plot.EV.arctic.quad.uni


# Briere; data-informed priors
# Load TPCs fitted using infform priors
load("R-scripts/R2jags-objects/all-mods/EV.arctic.bri.inf.Rdata")

# Extract model prediction
df.EV.arctic.bri.inf <- data.frame(EV.arctic.bri.inf$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)

# Plot
plot.EV.arctic.bri.inf <- df.EV.arctic.bri.inf %>% 
  ggplot(aes(x = temp)) +
  geom_line(aes(y = X2.5.), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = X97.5.), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = X50.), colour = "grey4", linewidth = 1) +
  
  ## data
  geom_point(data = data.EV.arctic, aes(x = temp, y = trait), size = 2, shape = 1) +  # Customize the axes and labels  # Customize the axes and labels
  
  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Proportion hatching",
       title = "C) Briere model; data-informed priors"
  ) +
  theme_bw()


plot.EV.arctic.bri.inf


# Quadratic; data-informed priors
# Load TPCs fitted using infform priors
load("R-scripts/R2jags-objects/all-mods/EV.arctic.quad.inf.Rdata")

# Extract model prediction
df.EV.arctic.quad.inf <- data.frame(EV.arctic.quad.inf$BUGSoutput$summary)[-(1:5),] %>% 
  mutate(temp = Temp.xs) %>% # Add the corresponding temp to the dataframe
  dplyr::select(temp, mean, sd, X2.5., X50., X97.5.)

# Plot
plot.EV.arctic.quad.inf <- df.EV.arctic.quad.inf %>% 
  ggplot(aes(x = temp)) +
  geom_line(aes(y = X2.5.), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = X97.5.), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = X50.), colour = "grey4", linewidth = 1) +
  
  ## data
  geom_point(data = data.EV.arctic, aes(x = temp, y = trait), size = 2, shape = 1) +  # Customize the axes and labels  # Customize the axes and labels
  
  # Customize the axes and labels
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Proportion hatching",
       title = "D) Quadratic model; data-informed priors"
  ) +
  theme_bw()


plot.EV.arctic.quad.inf


plot.EV.arctic <- plot_grid(plot.EV.arctic.bri.uni, plot.EV.arctic.quad.uni,
                             plot.EV.arctic.bri.inf, plot.EV.arctic.quad.inf, ncol = 2)
plot.EV.arctic

ggsave("figures/FigS6-EV.uni.vs.inf.png", plot.EV.arctic, width = 12, height = 8)




# 3. Summary table for TPC parameters ------------------------------------------
##### TPC summary
# load("R-scripts/R2jags-objects/best-fitting-mods/a.alldata.mod.Rdata")
# load("R-scripts/R2jags-objects/best-fitting-mods/bc.nonarctic.mod.Rdata")
# load("R-scripts/R2jags-objects/best-fitting-mods/lf.arctic.mod.Rdata")
# load("R-scripts/R2jags-objects/best-fitting-mods/PDR.arctic.mod.Rdata")
# load("R-scripts/R2jags-objects/best-fitting-mods/EFGC.alldata.mod.Rdata")
# load("R-scripts/R2jags-objects/best-fitting-mods/EV.arctic.mod.Rdata")
# load("R-scripts/R2jags-objects/best-fitting-mods/pLA.arctic.mod.Rdata")
# load("R-scripts/R2jags-objects/best-fitting-mods/MDR.arctic.mod.Rdata")
# 
# a.alldata.mod$BUGSoutput$summary[c("cf.T0", "cf.Tm", "cf.q", "cf.sigma", "deviance"),]
# bc.nonarctic.mod$BUGSoutput$summary[c("cf.T0", "cf.Tm", "cf.q", "cf.sigma", "deviance"),]
# lf.arctic.mod$BUGSoutput$summary[c("cf.T0", "cf.Tm", "cf.q", "cf.sigma", "deviance"),]
# PDR.arctic.mod$BUGSoutput$summary[c("cf.T0", "cf.Tm", "cf.q", "cf.sigma", "deviance"),]
# EFGC.alldata.mod$BUGSoutput$summary[c("cf.T0", "cf.Tm", "cf.q", "cf.sigma", "deviance"),]
# EV.arctic.mod$BUGSoutput$summary[c("cf.T0", "cf.Tm", "cf.q", "cf.sigma", "deviance"),]
# pLA.arctic.mod$BUGSoutput$summary[c("cf.T0", "cf.Tm", "cf.q", "cf.sigma", "deviance"),]
# MDR.arctic.mod$BUGSoutput$summary[c("cf.T0", "cf.Tm", "cf.q", "cf.sigma", "deviance"),]


# Create metadata for each trait
trait_info <- tribble(~trait, ~func, ~raneff, ~params_summary, ~case,
                      "Biting rate (a)", "B", "Yes", a.alldata.params.summary, "moderate-case",
                      "Vector competence (bc)", "Q", "No", bc.nonarctic.params.summary, "worst-case",
                      "Mosquito adult lifespan (lf)", "Q", "Yes", lf.arctic.params.summary, "best-case",
                      "Pathogen development rate (PDR)", "B", "No", PDR.arctic.params.summary, "best-case",
                      "Eggs per female per gonotrophic cycle (EFGC)", "Q", "Yes", EFGC.alldata.params.summary, "moderate-case",
                      "Egg viability (EV)", "Q", "No", EV.arctic.params.summary, "best-case",
                      "Larval-to-adult survival (pLA)", "Q", "No", pLA.arctic.params.summary, "best-case",
                      "Mosquito development rate (MDR)", "Q", "No", MDR.arctic.params.summary, "best-case"
                      )

# Function to format TPC parameter estimates into "median (lowerCI - upperCI)" format
# Use scientific notation for small values and fixed decimals for 
# larger values
format_estimate <- function(median, lower, upper){
  
  fmt <- function(x){
    
    if(abs(x) < 0.01){
      formatC(x, format = "e", digits = 1)
    } else if (abs(x) < 0.1){
      formatC(x, format = "f", digits = 2)
    }
    else {
      formatC(x, format = "f", digits = 1)
    }
    
  }
  
  paste0(
    fmt(median), " (",
    fmt(lower), "–",
    fmt(upper), ")"
  )
}


make_row <- function(trait, func, raneff, params_summary, case){
  tibble(Trait = trait,
         `F(x)` = func,
         `Random effects` = raneff,
         q = format_estimate(
           median = params_summary$median[params_summary$term == "cf.q"],
           lower = params_summary$lowerCI[params_summary$term == "cf.q"],
           upper = params_summary$upperCI[params_summary$term == "cf.q"]
           ),
         `Tmin (C)` = format_estimate(
           median = params_summary$median[params_summary$term == "cf.T0"],
           lower = params_summary$lowerCI[params_summary$term == "cf.T0"],
           upper = params_summary$upperCI[params_summary$term == "cf.T0"]
         ),
         `Topt (C)` = format_estimate(
           median = params_summary$median[params_summary$term == "Topt"],
           lower = params_summary$lowerCI[params_summary$term == "Topt"],
           upper = params_summary$upperCI[params_summary$term == "Topt"]
         ),
         `Tmax (C)` = format_estimate(
           median = params_summary$median[params_summary$term == "cf.Tm"],
           lower = params_summary$lowerCI[params_summary$term == "cf.Tm"],
           upper = params_summary$upperCI[params_summary$term == "cf.Tm"]
         ),
         Case = case
  )
  }

table <- pmap_dfr(trait_info,
                  make_row)


table1 <- flextable(table)


save_as_docx(
  "Table 1" = table1,
  path = "figures/table1.docx"
)



# 4. Summary table for DIC ------------------------------------------
load("R-scripts/R2jags-objects/all-mods/EV.arctic.quad.inf.Rdata")

# Helper: load one .Rdata file and pull out the DIC, returns NA if file doesn't exist
get_dic <- function(filepath) {
  if (!file.exists(filepath)) return(NA)
  obj_name <- load(filepath)        # load() returns the name(s) of loaded objects
  model    <- get(obj_name[1])      # retrieve object by that name
  model$BUGSoutput$DIC
}

# Define traits and their file "group" (which determines suffix + which models exist)
trait_info <- data.frame(
  trait = c("a", "bc", "lf", "PDR", "EFGC", "EV", "pLA", "MDR"),
  group = c("alldata", "nonarctic", rep("arctic", 2), "alldata", rep("arctic", 3)),
  stringsAsFactors = FALSE
)

# Function to build filepath, or return NA path if that model wasn't fit
make_path <- function(trait, group, model, prior) {
  # informed priors only exist for the "arctic" group
  if (prior == "inf" && group != "arctic") return(NA_character_)
  paste0("R-scripts/R2jags-objects/all-mods/", trait, ".", group, ".", model, ".", prior, ".Rdata")
}

dic_table <- trait_info %>%
  rowwise() %>%
  mutate(
    Briere_Uniform     = get_dic(make_path(trait, group, "bri", "uni")),
    Briere_Informed    = { p <- make_path(trait, group, "bri", "inf"); if (is.na(p)) NA else get_dic(p) },
    Quadratic_Uniform  = get_dic(make_path(trait, group, "quad", "uni")),
    Quadratic_Informed = { p <- make_path(trait, group, "quad", "inf"); if (is.na(p)) NA else get_dic(p) }
  ) %>%
  ungroup() %>%
  select(-group)

# Round for display
dic_table[,-1] <- round(dic_table[,-1], 1)

dic_table[1] <- c("Biting rate (a)", "Vector competence (bc)", 
                  "Adult lifespan (lf)", "Pathogen development rate (PDR)",
                  "Eggs per female per gonotrophic cycle (EFGC)", 
                  "Egg viability (EV)", "Larval-to-adult survival (pLA)",
                  "Mosquito development rate (MDR)")
print(dic_table)

dic_table <- flextable(dic_table)

save_as_docx(
  "Table S4" = dic_table,
  path = "figures/tableS4.docx"
)

