## Lilian Chan, University of Guelph
## Arctic vector-borne disease transmission suitability model
##
## Purpose: Use trait thermal response posterior distributions from JAGS to 
## calculate suitability S(T)
## 
## Table of content:
##    0. Set-up workspace
##    1. Load R2jags model output
##    2. Calculate S(T)
##    3. Sensitivity analysis - partial derivatives
##
##
## Inputs:
## Best-fitting TPC models for Arctic species:
## R-scripts/R2jags-objects/best-fitting-mods/a.alldata.mod.Rdata
## R-scripts/R2jags-objects/best-fitting-mods/lf.arctic.mod.Rdata
## R-scripts/R2jags-objects/best-fitting-mods/PDR.arctic.mod.Rdata
## R-scripts/R2jags-objects/best-fitting-mods/EFGC.alldata.mod.Rdata
## R-scripts/R2jags-objects/best-fitting-mods/EV.arctic.mod.Rdata
## R-scripts/R2jags-objects/best-fitting-mods/pLA.arctic.mod.Rdata
## R-scripts/R2jags-objects/best-fitting-mods/MDR.arctic.mod.Rdata
##
## data-processed/bc/bc.arctic.predictions.fullposts.csv - 
##     Full posterior distributions for bc TPC predictions after cold-shifting
## 
## Full posterior distributions for TPC parameters:
## data-processed/a/a.alldata.params.fullposts.csv
## data-processed/bc/bc.arctic.params.fullposts.csv
## data-processed/lf/lf.arctic.params.fullposts.csv
## data-processed/PDR/PDR.arctic.params.fullposts.csv
## data-processed/EFGC/EFGC.alldata.params.fullposts.csv
## data-processed/EV/EV.arctic.params.fullposts.csv
## data-processed/pLA/pLA.arctic.params.fullposts.csv
## data-processed/MDR/MDR.arctic.params.fullposts.csv
## 
##
## Outputs: 
## figures/Fig4-suitability.sensitivity.png -
##     Main text figure 4
##
## data-processed/suitability/S.predictions.fullposts.csv -
##     Full posterior distributions for suitability calculations
##
## data-processed/suitability/S.params.fullposts.csv - 
##     Full posterior distributions for suitability
##
## data-processed/suitability/S.predictions.summary.csv - 
##     Posterior summary of suitability across temperatures
##
## data-processed/suitability/S.params.summary.csv - 
##     Summary statistics of Tmin, Topt, and Tmax for suitability
## 
## ## data-processed/suitability/S.predictions.scaled.csv -
##     Suitability for mapping



# 0. Set-up workspace ----------------------------------------------------------

library(tidyverse)
library(readxl)
library(janitor)
library(R2jags)
library(ggsci)
library(cowplot)
library(grafify)

##### Load functions
source("R-scripts/00_Functions.R")


# 1. Load R2jags model output ----------------------------------------------------------

## biting rate (a)
load("R-scripts/R2jags-objects/best-fitting-mods/a.alldata.mod.Rdata")

## Adult lifespan (lf)
load("R-scripts/R2jags-objects/best-fitting-mods/lf.arctic.mod.Rdata")

## Parasite development rate (PDR)
load("R-scripts/R2jags-objects/best-fitting-mods/PDR.arctic.mod.Rdata")

## Eggs per female per gonotrophic cycle (EFGC)
load("R-scripts/R2jags-objects/best-fitting-mods/EFGC.alldata.mod.Rdata")

## Egg viability (EV)
load("R-scripts/R2jags-objects/best-fitting-mods/EV.arctic.mod.Rdata")

## Larval-to-adult survival (pLA)
load("R-scripts/R2jags-objects/best-fitting-mods/pLA.arctic.mod.Rdata")

## Mosquito development rate (MDR)
load("R-scripts/R2jags-objects/best-fitting-mods/MDR.arctic.mod.Rdata")


#####  Pull out the derived/predicted values:
a.preds <- a.alldata.mod$BUGSoutput$sims.list$z.trait.mu.pred.pop ## Only get the global-level fit
bc.preds <- read_csv("data-processed/bc/bc.arctic.predictions.fullposts.csv")
bc.preds <- as.matrix(bc.preds)

lf.preds <- lf.arctic.mod$BUGSoutput$sims.list$z.trait.mu.pred.pop ## Only get the global-level fit
PDR.preds <- PDR.arctic.mod$BUGSoutput$sims.list$z.trait.mu.pred
EFGC.preds <- EFGC.alldata.mod$BUGSoutput$sims.list$z.trait.mu.pred.pop ## Only get the global-level fit
EV.preds <- EV.arctic.mod$BUGSoutput$sims.list$z.trait.mu.pred
pLA.preds <- pLA.arctic.mod$BUGSoutput$sims.list$z.trait.mu.pred
MDR.preds <- MDR.arctic.mod$BUGSoutput$sims.list$z.trait.mu.pred


## Pull out the full posterior distributions of TPC parameters
a.params.fullposts <- read.csv("data-processed/a/a.alldata.params.fullposts.csv")
bc.params.fullposts <- read.csv("data-processed/bc/bc.arctic.params.fullposts.csv")
lf.params.fullposts <- read.csv("data-processed/lf/lf.arctic.params.fullposts.csv")
PDR.params.fullposts <- read.csv("data-processed/PDR/PDR.arctic.params.fullposts.csv")
EFGC.params.fullposts <- read.csv("data-processed/EFGC/EFGC.alldata.params.fullposts.csv")
EV.params.fullposts <- read.csv("data-processed/EV/EV.arctic.params.fullposts.csv")
pLA.params.fullposts <- read.csv("data-processed/pLA/pLA.arctic.params.fullposts.csv")
MDR.params.fullposts <- read.csv("data-processed/MDR/MDR.arctic.params.fullposts.csv")



# 2.  Calculate S(T) -----------------------------------------------------------

## Columns = temp from 0 to 45ºC at a 0.1ºC interval, Rows = 15000 MCMC iterations
S.calc <- S(a.preds, bc.preds, lf.preds, PDR.preds, EFGC.preds, EV.preds, 
            pLA.preds, MDR.preds)

##### Temp sequence for derived quantity calculations
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)


# save all 15000 MCMC iterations of suitability calculation (451 row (temp), columns = temp and 15000 MCMC iterations)
S.preds.fullposts <- data.frame(Temp.xs, t(S.calc))
colnames(S.preds.fullposts) <- c("temp", paste0("iter", seq(1:nrow(S.calc))))

# Save output
write.csv(S.preds.fullposts, "data-processed/suitability/S.predictions.fullposts.csv")

# Get the Tmin, Topt and Tmax for each MCMC iteration
S.params.fullposts <- calcDerivedTPCParamPosteriors(S.calc, Temp.xs)
S.params.fullposts$iter <- seq(1:nrow(S.calc))
write.csv(S.params.fullposts, "data-processed/suitability/S.params.fullposts.csv")


# Get S mean, median, upper + lower CIs across temp gradient
S.out <- calcPostQuants(as.data.frame(S.calc), "S", Temp.xs)
S.out <- S.out %>% 
  dplyr::select(temperature, lowerCI, lowerQ, median, mean, upperQ, upperCI, trait)

head(S.out)

## Calculate relative S(T) 
# by scaling to the maximum median
S.out.median <- S.out %>% 
  mutate(scaled_lowerCI = S.out$lowerCI / max(S.out$median)) %>%
  mutate(scaled_lowerQ = S.out$lowerQ / max(S.out$median)) %>%
  mutate(scaled_median = S.out$median / max(S.out$median)) %>% 
  mutate(scaled_upperQ= S.out$upperQ / max(S.out$median)) %>% 
  mutate(scaled_upperCI = S.out$upperCI / max(S.out$median))
  
write.csv(S.out.median, "data-processed/suitability/S.predictions.summary.csv")




# Plot S
plot.S <- ggplot(data = S.out.median) +
  geom_ribbon(aes(x = temperature, ymin = scaled_lowerCI, ymax = scaled_upperCI),
              fill = "grey",
              alpha = 0.5) +
  geom_ribbon(aes(x = temperature, ymin = scaled_lowerQ, ymax = scaled_upperQ),
              fill = "grey",
              alpha = 0.7) +
  geom_line(aes(x = temperature, y = scaled_median), colour = "black", linewidth = 1) +

  # scale_x_continuous(limits = c(10, 35)) +
  scale_x_continuous(limits = c(10, 40)) +
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = "Suitability (S)") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))
  

plot.S

ggsave("figures/suitability.png", plot.S, width = 10.3, height = 5.6)



# 3. Calculate Tmin, Tmax and Topt (and CIs) for suitability --------------------------------------

S.params.summary <- extractDerivedTPC(as.data.frame(S.calc), "S", Temp.xs)

S.params.summary <- S.params.summary %>% 
  mutate(term = case_when(term == "cf.Tm" ~ "Tmax",
                          term == "cf.T0"~ "Tmin",
                          term == "Topt" ~ "Topt",
                          term == "Tbreadth" ~ "Tbreadth"))


S.params.summary$term <- factor(S.params.summary$term,
                         levels = c("Tmax", "Topt", "Tmin", "Tbreadth"))

# Save output
write.csv(S.params.summary, "data-processed/suitability/S.params.summary.csv")

# Plot
plot.S.params <- S.params.summary %>% 
  filter(term != "Tbreadth") %>% 
  ggplot() +
  geom_linerange(aes(xmin = lowerQ, xmax = upperQ, y = term, colour = term), 
                 linewidth = 1) +
  geom_linerange(aes(xmin = lowerCI, xmax = upperCI, y = term, colour = term), 
                 linewidth = 0.5) + #default size = 0.5
  geom_point(aes(x = median, y = term, colour = term)) +
  scale_x_continuous(limits = c(10, 40)) +
  labs(x = expression(paste("Temperature (", degree, "C)"))) +
  scale_y_discrete(labels=c("Tmin" = expression(paste("T"[min])), 
                            "Topt" = expression(paste("T"[opt])), 
                            "Tmax" = expression(paste("T"[max])))) +
  theme(axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

plot.S.params

plot.suitability <- plot_grid(plot.S, plot.S.params,
                              ncol = 1,
                              rel_heights = c(2,1),
                              align = "hv"
) +
  theme(panel.background = element_rect(fill = "white", color = NA))

plot.suitability


# 4. Sensitivity Analysis - partial derivatives --------------------------------

# Temperature levels and # MCMC steps
Temp.xs <- seq(0, 45, 0.1)
N.Temp.xs <-length(Temp.xs)
nMCMC <- nrow(S.calc)


##### Calculate trait means at each temperature
a.m <- colMeans(a.preds)
bc.m <- colMeans(bc.preds)
lf.m <- colMeans(lf.preds)
PDR.m <- colMeans(PDR.preds)
EFGC.m <- colMeans(EFGC.preds)
EV.m <- colMeans(EV.preds)
pLA.m <- colMeans(pLA.preds)
MDR.m <- colMeans(MDR.preds)


# Calculate sensitivity using partial derivatives
SA <- SensitivityAnalysis_pd(a.preds, bc.preds, lf.preds, PDR.preds, 
                             EFGC.preds, EV.preds, pLA.preds, MDR.preds,
                             a.params.fullposts, bc.params.fullposts, 
                             lf.params.fullposts, PDR.params.fullposts, 
                             EFGC.params.fullposts, EV.params.fullposts, 
                             pLA.params.fullposts, MDR.params.fullposts,
                             a.m, bc.m, lf.m, PDR.m, EFGC.m, EV.m, pLA.m, MDR.m)


# Get sensitivity posteriors for each term and summarize them
dS.da		<- calcPostQuants(as.data.frame(SA[[1]]), "a", Temp.xs)
dS.dbc		<- calcPostQuants(as.data.frame(SA[[2]]), "bc", Temp.xs)
dS.dlf		<- calcPostQuants(as.data.frame(SA[[3]]), "lf", Temp.xs)
dS.dPDR	<- calcPostQuants(as.data.frame(SA[[4]]), "PDR", Temp.xs)
dS.dEFGC		<- calcPostQuants(as.data.frame(SA[[5]]), "EFGC", Temp.xs)
dS.dEV 	<- calcPostQuants(as.data.frame(SA[[6]]),"EV",  Temp.xs)
dS.dpLA 	<- calcPostQuants(as.data.frame(SA[[7]]), "pLA", Temp.xs)
dS.dMDR		<- calcPostQuants(as.data.frame(SA[[8]]), "MDR", Temp.xs)
dS.dT		<- calcPostQuants(as.data.frame(SA[[9]]), "S", Temp.xs)



##### Plot results
plot.SA <- ggplot() +
  geom_line(data = dS.da, aes(x = temperature, y = median, colour = "a"), linewidth = 1) +
  geom_line(data = dS.dbc, aes(x = temperature, y = median, colour = "bc"), linewidth = 1) +
  geom_line(data = dS.dlf, aes(x = temperature, y = median, colour = "lf"), linewidth = 1) +
  geom_line(data = dS.dPDR, aes(x = temperature, y = median, colour = "PDR"), linewidth = 1) +
  geom_line(data = dS.dEFGC, aes(x = temperature, y = median, , colour = "EFGC"), linewidth = 1) +
  geom_line(data = dS.dEV, aes(x = temperature, y = median, , colour = "EV"), linewidth = 1) +
  geom_line(data = dS.dpLA, aes(x = temperature, y = median, colour = "pLA"), linewidth = 1) +
  geom_line(data = dS.dMDR, aes(x = temperature, y = median, colour = "MDR"), linewidth = 1) +
  geom_line(data = dS.dT, aes(x = temperature, y = median, , colour = "S"), linewidth = 1.2) +
  # scale_x_continuous(limits = c(10, 35)) +
  scale_x_continuous(limits = c(10, 40)) +
  labs(x = expression(paste("Temperature (", degree, "C)")),
       y = "Relative sensitivity") +
  
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
                                 "Pathogen development rate (PDR)", 
                                 "Eggs per female per\ngonotrophic cycle (EFGC)",  
                                 "Egg viability (EV)", 
                                 "Larval-to-adult survival (pLA)", 
                                 "Mosquito development rate (MDR)")) + 
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))


plot.SA

ggsave("figures/sensitivity.analysis.png", plot.SA, width = 10, height = 5)


plot.everything <- plot_grid(plot.S, plot.S.params, plot.SA,
                             ncol = 1,
                             labels = c(LETTERS[1:3]), 
                             rel_heights = c(2,1,2),
                             align = "hv"
                             ) +
  theme(panel.background = element_rect(fill = "white", color = NA))

plot.everything

ggsave("figures/Fig4-suitability.sensitivity.png", plot.everything, 
       width = 10, height = 10)




# 5. Packaging output for mapping ----------------------------------------------

head(S.out)

## Goal: create maps to compare the spatial distribution of days of thermal 
## suitability for transmission.

## We use two thermal suitability thresholds: S(T) > 0.001 and S(T) > 0.5, each 
## with a posterior probability greater than 0.975, 0.5, and 0.025.

## We will do that by scale the lower CI, median, and upper CI contours to 
## themselves get S(T) for mapping

S.out.scaled <- data.frame(temp = S.out$temperature,
                           scaled_lowerCI = S.out$lowerCI/max(S.out$lowerCI),
                           scaled_median = S.out$median/max(S.out$median),
                           scaled_upperCI = S.out$upperCI / max(S.out$upperCI)
                           )

head(S.out.scaled)

# Check output
plot(scaled_lowerCI ~ temp, data = S.out.scaled, type = "l")
plot(scaled_median ~ temp, data = S.out.scaled, type = "l")
plot(scaled_upperCI ~ temp, data = S.out.scaled, type = "l")

write.csv(S.out.scaled, "data-processed/suitability/S.predictions.scaled.csv")
