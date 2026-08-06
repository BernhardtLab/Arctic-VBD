## Lilian Chan, University of Guelph
## Arctic vector-borne disease transmission suitability model
##
## Purpose: calculated the mean air temperature for June, July, and August 2025 
##          in Kugluktuk and Cambridge Bay. 
## 
## Table of content:
##    0. Set-up workspace
##    1. Load data
##  	2. 
##
## Inputs:
## data-raw/temperature-data/en_climate_daily_NU_2300903_2025_P1D.csv - 
##     Daily temperature at Kugluktuk
##
## data-raw/temperature-data/en_climate_daily_NU_2400601_2025_P1D.csv -
##     Daily temperature at Cambridge Bay
##
## Outputs: 
## figures/tableS10.docx - 
##      Supplementary table S10


# 0. Set-up workspace ----------------------------------------------------------

library(tidyverse)
library(janitor)
library(flextable)


#  1. Load data ----------------------------------------------------------------

## Look at the mean summer temperature of Kugluktuk and Cambridge Bay from 2025


## Kugluktuk data
KUG.data.2025 <- read_csv("data-raw/temperature-data/en_climate_daily_NU_2300903_2025_P1D.csv")

KUG.data.2025 <- KUG.data.2025 %>%
  clean_names() %>% 
  select(longitude_x, latitude_y, station_name, climate_id, date_time, year, 
         month, day, max_temp_c, min_temp_c, mean_temp_c) %>%
  mutate(station_name = "KUGLUKTUK")
  
## Cambridge Bay data
CBAY.data.2025 <- read_csv("data-raw/temperature-data/en_climate_daily_NU_2400601_2025_P1D.csv")

CBAY.data.2025 <- CBAY.data.2025 %>%
  clean_names() %>% 
  select(longitude_x, latitude_y, station_name, climate_id, date_time, year, 
         month, day, max_temp_c, min_temp_c, mean_temp_c) %>%
  mutate(station_name = "CAMBRIDGE BAY")


# Filter sampling period
KUG.data.2025.samp <- KUG.data.2025 %>%
  filter(date_time >= "2025-06-24" & date_time <= "2025-09-08")


CBAY.data.2025.samp <- CBAY.data.2025 %>%
  filter(date_time >= "2025-07-19" & date_time <= "2025-08-05")

temp.data.samp <- rbind(KUG.data.2025.samp, CBAY.data.2025.samp)

head(temp.data.samp)


#  2. Calculate mean temp ------------------------------------------------------

temp.data.samp.summary <- temp.data.samp %>% 
  filter(!is.na(mean_temp_c)) %>%
  group_by(station_name) %>% 
  summarise(mean_temp_sampling_period = mean(mean_temp_c))

temp.data.samp.summary

# Calculate Mean daily air temperatures for each month
temp.data <- rbind(KUG.data.2025, CBAY.data.2025)

## select sampling period
temp.data$year <- as.factor(temp.data$year)
temp.data$month <- as.factor(temp.data$month)

temp.data <- temp.data %>% 
  filter(!is.na(mean_temp_c)) %>% 
  filter(month %in% c("06", "07", "08", "09")) %>% 
  mutate(month2 = ifelse(month == "06", "Jun", 
                         ifelse(month == "07", "Jul", 
                                ifelse(month == "08", "Aug",
                                       "Sep")
                                )
                         )
         )
  
temp.data$month2 <- factor(temp.data$month2, levels = c("Jun", "Jul", "Aug", "Sep"))

temp.data.summary <- temp.data %>% 
  group_by(station_name, month2) %>% 
  summarise(mean_temp = mean(mean_temp_c))

temp.data.summary

# Change wide format
temp.data.summary.wide <- temp.data.summary %>% 
  pivot_wider(., names_from = month2,  values_from = mean_temp)

temp.data.summary.wide
temp.data.samp.summary

# Put the two table together
table <- inner_join(temp.data.samp.summary, temp.data.summary.wide, by = "station_name")
table[,-1] <- lapply(table[,-1], function(x) round(as.numeric(x), 1))

table <- flextable(table)

save_as_docx(table, path = "figures/tableS10.docx")

plot <- temp.data %>% ggplot(aes(x = mean_temp_c, fill = station_name)) +
  geom_histogram(binwidth = 1, color = "white", alpha = 0.8) +
  # Add a black dotted line showing the lower thermal limit for transmission
  geom_vline(aes(xintercept = 13.4), linetype = "dashed") +
  # mean daily temp for each month
  geom_vline(
    data = temp.data.summary,
    aes(xintercept = mean_temp),
    color = "red", linewidth = 1) +
  facet_grid(station_name ~ month2) +
  labs(
    x = "Temperature (°C)", y = "Count"
  ) +
  theme_bw() +
  theme(legend.position = "none")

plot



