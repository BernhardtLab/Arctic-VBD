## Look at the mean summer temperature of Kugluktuk and Cambridge Bay from 2015-2025

library(tidyverse)
library(janitor)

## Kugluktuk data
KUG.data.2015 <- read_csv("Temp data/en_climate_daily_NU_2300903_2015_P1D.csv")
KUG.data.2016 <- read_csv("Temp data/en_climate_daily_NU_2300903_2016_P1D.csv")
KUG.data.2017 <- read_csv("Temp data/en_climate_daily_NU_2300903_2017_P1D.csv")
KUG.data.2018 <- read_csv("Temp data/en_climate_daily_NU_2300903_2018_P1D.csv")
KUG.data.2019 <- read_csv("Temp data/en_climate_daily_NU_2300903_2019_P1D.csv")
KUG.data.2020 <- read_csv("Temp data/en_climate_daily_NU_2300903_2020_P1D.csv")
KUG.data.2021 <- read_csv("Temp data/en_climate_daily_NU_2300903_2021_P1D.csv")
KUG.data.2022 <- read_csv("Temp data/en_climate_daily_NU_2300903_2022_P1D.csv")
KUG.data.2023 <- read_csv("Temp data/en_climate_daily_NU_2300903_2023_P1D.csv")
KUG.data.2024 <- read_csv("Temp data/en_climate_daily_NU_2300903_2024_P1D.csv")
KUG.data.2025 <- read_csv("Temp data/en_climate_daily_NU_2300903_2025_P1D.csv")


## Cambridge Bay data
CBAY.data.2015 <- read_csv("Temp data/en_climate_daily_NU_2400601_2015_P1D.csv")
CBAY.data.2016 <- read_csv("Temp data/en_climate_daily_NU_2400601_2016_P1D.csv")
CBAY.data.2017 <- read_csv("Temp data/en_climate_daily_NU_2400601_2017_P1D.csv")
CBAY.data.2018 <- read_csv("Temp data/en_climate_daily_NU_2400601_2018_P1D.csv")
CBAY.data.2019 <- read_csv("Temp data/en_climate_daily_NU_2400601_2019_P1D.csv")
CBAY.data.2020 <- read_csv("Temp data/en_climate_daily_NU_2400601_2020_P1D.csv")
CBAY.data.2021 <- read_csv("Temp data/en_climate_daily_NU_2400601_2021_P1D.csv")
CBAY.data.2022 <- read_csv("Temp data/en_climate_daily_NU_2400601_2022_P1D.csv")
CBAY.data.2023 <- read_csv("Temp data/en_climate_daily_NU_2400601_2023_P1D.csv")
CBAY.data.2024 <- read_csv("Temp data/en_climate_daily_NU_2400601_2024_P1D.csv")
CBAY.data.2025 <- read_csv("Temp data/en_climate_daily_NU_2400601_2025_P1D.csv")


temp.data <- rbind(KUG.data.2015, KUG.data.2016, KUG.data.2017, KUG.data.2018,
                   KUG.data.2019, KUG.data.2020, KUG.data.2021, KUG.data.2022,
                   KUG.data.2023, KUG.data.2024, KUG.data.2025,
                   CBAY.data.2015, CBAY.data.2016, CBAY.data.2017, 
                   CBAY.data.2018, CBAY.data.2019, CBAY.data.2020,
                   CBAY.data.2021, CBAY.data.2022, CBAY.data.2023, 
                   CBAY.data.2024, CBAY.data.2025)


temp.data <- temp.data %>% 
  clean_names() %>% 
  select(longitude_x, latitude_y, station_name, climate_id, date_time, year, 
         month, day, max_temp_c, min_temp_c, mean_temp_c)

write_csv(temp.data, "Temp data/station_data.csv")

## select summer months
temp.data$year <- as.factor(temp.data$year)
temp.data$month <- as.factor(temp.data$month)

temp.data <- temp.data %>% 
  filter(!is.na(mean_temp_c)) %>% 
  filter(month %in% c("05", "06", "07", "08", "09"))

temp.data.summary <- temp.data %>% 
  group_by(station_name, year, month) %>% 
  summarise(monthly_mean_temp = mean(mean_temp_c),
            num_observations = n())

temp.data.summary
temp.data.summary$month2 <- rep(c("May", "Jun", "Jul", "Aug", "Sep"), 22)


plot <- temp.data.summary %>% ggplot(aes(x = year, y = monthly_mean_temp)) +
  geom_line(aes(color = station_name, 
                group = interaction(station_name, month2)), 
            linewidth = 1.25) +
  scale_color_discrete(breaks = c("KUGLUKTUK A", "CAMBRIDGE BAY A"),
                       labels = c("Kugluktuk", "Cambridge Bay")) +
  labs(y = "Monthly mean temperature (ºC)") +
  facet_wrap(~factor(month2, c("May", "Jun", "Jul", "Aug", "Sep"))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))


plot

ggsave("Temp data/historical_temp.png", plot)



