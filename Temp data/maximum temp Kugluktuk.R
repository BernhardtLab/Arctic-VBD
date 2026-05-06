library(tidyverse)
library(janitor)

data.2014 <- read_csv("Temp data/en_climate_daily_NU_2300903_2014_P1D.csv")
data.2015 <- read_csv("Temp data/en_climate_daily_NU_2300903_2015_P1D.csv")
data.2016 <- read_csv("Temp data/en_climate_daily_NU_2300903_2016_P1D.csv")
data.2017 <- read_csv("Temp data/en_climate_daily_NU_2300903_2017_P1D.csv")
data.2018 <- read_csv("Temp data/en_climate_daily_NU_2300903_2018_P1D.csv")
data.2019 <- read_csv("Temp data/en_climate_daily_NU_2300903_2019_P1D.csv")
data.2020 <- read_csv("Temp data/en_climate_daily_NU_2300903_2020_P1D.csv")
data.2021 <- read_csv("Temp data/en_climate_daily_NU_2300903_2021_P1D.csv")
data.2022 <- read_csv("Temp data/en_climate_daily_NU_2300903_2022_P1D.csv")
data.2023 <- read_csv("Temp data/en_climate_daily_NU_2300903_2023_P1D.csv")
data.2024 <- read_csv("Temp data/en_climate_daily_NU_2300903_2024_P1D.csv")

data <- rbind(data.2015, data.2016, data.2017, data.2018,
              data.2019, data.2020, data.2021, data.2022, data.2023, data.2024)

data <- data %>% 
  clean_names() %>% 
  dplyr::select(latitude_y, longitude_x, station_name, date_time, year, month, 
                day, max_temp_c, min_temp_c, mean_temp_c) %>% 
  filter(!is.na(max_temp_c))

data$year <- as.factor(data$year)
data$month <- as.factor(data$month)

data <- data %>% 
  filter(month %in% c("06","07","08"))

# Find the highest maximum daily temp of each month each year
data.highest <- data %>% 
  group_by(year, month) %>% 
  summarise(max_temp = max(max_temp_c))

data.mean <- data %>% 
  group_by(year, month, station_name) %>% 
  summarise(mean_temp = max(mean_temp_c))

ggplot(data = data.mean, aes(x = year, y = mean_temp,group = station_name, color = station_name)) +
  geom_line() +
  facet_wrap(~month) +
  theme_bw()

