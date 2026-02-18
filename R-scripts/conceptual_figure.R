## Lilian Chan, University of Guelph
## Arctic vector-borne disease transmission suitability model
##
## Purpose: Create a conceptual figure for the study

##########
###### 0. Set-up workspace ----
##########

library(tidyverse)
library(ggsci)
library(ggpubr) # For ggarrange
library(grafify)


briere = function(T, T0, Tm, q){
  
  b <- c()
  
  for (i in 1:length(T)){
    if (T[i]>T0 && T[i]<Tm) # When trait value > 0
    {b[i] <- q * T[i] * (T[i]-T0) * (Tm-T[i])**0.5} # Briere function
    else {b[i] <- 0}
  }
  
  b # return output
  
}

quadratic = function(T, T0, Tm, q){
  
  b <- c()
  
  for (i in 1:length(T)){
    if (T[i]>T0 && T[i]<Tm) # When trait value > 0
    {b[i] <- -1 * q * (T[i]-T0) * (T[i] - Tm)} # Quadratic function
    else {b[i] <- 0}
  }
  
  b # return output
  
}


##########
###### 1. Figure 1 ----
##########

### Panel A: Estimate trait thermal responses ----

##### Generate the TPCs
Temp.xs <- seq(0, 45, 0.1)

bri <- briere(Temp.xs, T0 = 10, Tm = 30, q = 0.01)
quad <- quadratic(Temp.xs, T0 = 0, Tm = 25, q = 0.08)


##### Generate data points
set.seed(50) # Set a random seed for reproducibility of the simulation
temp1 <- c(2, 5, 7, 10, 13, 16, 18)
mean1 <- quadratic(temp1, T0 = 0, Tm = 25, q = 0.08)
sigma <- 0.8
data1 <- rnorm(length(temp1), mean = mean1, sd = sigma)

set.seed(50) # Set a random seed for reproducibility of the simulation
temp2 <- c(13, 16, 19, 21, 24, 26, 29)
mean2 <- briere(temp2, T0 = 10, Tm = 30, q = 0.01)
sigma <- 0.8
data2 <- rnorm(length(temp2), mean = mean2, sd = sigma)

lines <- data.frame(temp = Temp.xs, trait1 = quad, trait2 = bri)
lines <- pivot_longer(lines, cols = 2:3, names_to = "trait", values_to = "value")
quad_data <- data.frame(temp = temp1, data = data1)
bri_data <- data.frame(temp = temp2, data = data2)

panel_a <- lines %>% ggplot(aes(x = temp)) +
  geom_line(aes(y = value, colour = trait), size = 1) +
  geom_point(data = quad_data, aes(y = data1), colour = "#0072B2", size = 2) +
  geom_point(data = bri_data, aes(y = data2), colour = "#D55E00", size = 2) +
  xlim(0,35) +
  labs(x = "Temperature", y = "Traits") +
  scale_colour_manual(values = c("trait1" = "#0072B2", "trait2" = "#D55E00"),
                      name = element_blank(), # No legend title
                      labels = c("Trait 1 \n(e.g. egg viability) \n", "Trait 2 \n(e.g. biting rate)")) + 
  theme(legend.text = element_text(size = 8),
        legend.position = c(1.01, 1),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(0, 0, 0, 0),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", 
                                 size = 0.9, linetype = "solid"))

panel_a

ggsave("figures/conceptual_figure_panelA.png", panel_a, width = 3, height = 2.5)


### Panel B: Calculate thermal suitability ----

## We will use a normal distribution to illustrate the output of the thermal suitability model
y <- dnorm(seq(0, 35, 0.1), mean = 17.5, sd = 5)
suitability <- data.frame(temp = seq(0, 35, 0.1), suit = y/max(y)) # Scaled from 0 to 1

##### Generate data points for field data
set.seed(50) # Set a random seed for reproducibility of the simulation
temp3 <- c(5, 10, 12, 15, 21)
mean3 <- dnorm(temp3, mean = 17.5, sd = 5)
sigma <- 0.01
data3 <- rnorm(length(temp3), mean = mean3, sd = sigma)

field_data <- data.frame(temp = temp3, data = data3/max(y)) # scaled

panel_b <-suitability %>% ggplot(aes(x = temp)) +
  geom_line(aes(y = suit), colour = "black", size = 1) +
  geom_point(data = field_data, aes(y = data), colour = "azure4", size = 2.5) +
  xlim(0,35) +
  labs(x = "Temperature", y = "Suitability") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", 
                                 size = 0.9, linetype = "solid"))
panel_b
ggsave("figures/conceptual_figure_panelB.png", panel_b, width = 3, height = 2.5)
