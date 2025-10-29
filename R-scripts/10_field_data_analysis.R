## Lilian Chan, University of Guelph
## Arctic vector-borne disease transmission suitability model
##
## Purpose: To get an idea on how we collected the field data 
## 
## Table of content:
##    0. Set-up workspace
##
##    1. 


##########
###### 0. Set-up workspace ----
##########
library(tidyverse)
library(readxl)
library(janitor)
library(ggsci)
library(lubridate) # for time zone
library(RColorBrewer) # colour palette

# Load data
data <- read_csv("data/raw-data/field-data/guelph_data.csv") %>% 
  clean_names()


sites <- data %>% count(site) %>% arrange(desc(n))

sample.dates <- data %>% 
  group_by(site, sample_type_collection_method, date_collected) %>% 
  count()

coor <- data %>% 
  distinct(site, lat, long) %>% 
  arrange(site)

# write_csv(coor, "field-data/coordinates.csv")

sample <- data %>% 
  group_by(sample_replicate) %>% 
  count()


# HOBO ----

airport1 <- read_xlsx("data-raw/field-data/hobo/airport_road_1.xlsx", 
                      sheet = "Data") %>% 
  clean_names() %>% 
  mutate(site = "airport_road_1")


airport2 <- read_xlsx("data-raw/field-data/hobo/airport_road_2.xlsx", 
                      sheet = "Data") %>% 
  clean_names() %>% 
  mutate(site = "airport_road_2")


construction_shop <- read_xlsx("data-raw/field-data/hobo/construction_shop.xlsx", 
                      sheet = "Data") %>% 
  clean_names() %>% 
  mutate(site = "construction_shop")


dew_line_road <- read_xlsx("data-raw/field-data/hobo/dew_line_road.xlsx", 
                               sheet = "Data") %>% 
  clean_names() %>% 
  mutate(site = "dew_line_road")


end_of_pelly_road <- read_xlsx("data-raw/field-data/hobo/end_of_pelly_road.xlsx", 
                           sheet = "Data") %>% 
  clean_names() %>% 
  mutate(site = "end_of_pelly_road")


first_creek <- read_xlsx("data-raw/field-data/hobo/first_creek.xlsx", 
                           sheet = "Data") %>% 
  clean_names() %>% 
  mutate(site = "first_creek")



sewage_dump <- read_xlsx("data-raw/field-data/hobo/sewage_dump.xlsx", 
                         sheet = "Data") %>% 
  clean_names() %>% 
  mutate(site = "sewage_dump")


swimming_hole <- read_xlsx("data-raw/field-data/hobo/swimming_hole.xlsx", 
                         sheet = "Data") %>% 
  clean_names() %>% 
  mutate(site = "swimming_hole")


temp_data <- rbind(airport1, 
                   # airport2, # only one data point
                   construction_shop,
                   dew_line_road,
                   end_of_pelly_road,
                   first_creek,
                   sewage_dump,
                   swimming_hole)


## Convert the timezone to Mountain saving time (the one cambay uses)
temp_data$date_time_edt <- force_tz(temp_data$date_time_edt, "EST")
temp_data$date_time_edt <- as.POSIXct(temp_data$date_time_edt,
                                          tz = "MST")

# Rename the datetime column
colnames(temp_data)[2] <- "date_time_mst"

## The malaise trap at First creek was taken down on Aug 5 at around noon. 
## However, the HOBO continued to take measurements until Aug 6.
## Remove measurements after Aug 5 at around 12 pm
temp_data <- temp_data %>% filter(site != "first_creek" |
                       (site == "first_creek" & date_time_mst <= as.POSIXct("2025-08-05 12:00:00", tz = "MST"))
                     )


# write_csv(temp_data[,-1], "data/data-processed/temp_all.csv")
temp_data <- read_csv("data-processed/temp_all.csv")

## Since HOBOs were set up at different days, only use dates that all HOBOs have data on
## The first malaise trap we set up (and taken down) was sewage dump (Jul 19 - Aug 2);
## the last one was Dew Line Road (Jul 22- Aug 5)
## So we select data from Jul 23 00:00 to Aug 1 23:59

temp_summary_std <- temp_data %>% 
  mutate(site = as.factor(site)) %>% 
  mutate(date = as.Date(date_time_mst, tz = "MST")) %>% 
  filter(date >= "2025-07-23" & date < "2025-08-02")

temp_summary_std_max <- temp_summary_std %>% 
  group_by(site, date) %>% 
  summarise(daily_max = max(temperature_c)) %>% 
  group_by(site) %>% 
  summarise(mean_daily_max_std = mean(daily_max))


## Calculate the mean of 4-hourly maximum temperature
## Since measurements were taken every 5 mins, select the highest 12*5*4 = 48 measurements of each day
temp_summary_std_4hrsmax <- temp_summary_std %>% 
  group_by(site, date) %>% 
  slice_max(order_by = temperature_c, n = 48, with_ties = FALSE) %>% 
  group_by(site) %>% 
  summarise(mean_daily_4hrmax_std = mean(temperature_c))


temp_summary_std_mean_daily_mean <- temp_summary_std %>% 
  group_by(site, date) %>% 
  summarise(daily_mean = mean(temperature_c)) %>% 
  group_by(site) %>% 
  summarise(mean_daily_mean_std = mean(daily_mean))


## Repeat the steps above but use all dates
temp_summary <- temp_data %>% 
  mutate(site = as.factor(site)) %>% 
  mutate(date = as.Date(date_time_mst, tz = "MST"))

temp_summary_max <- temp_summary %>% 
  group_by(site, date) %>% 
  summarise(daily_max = max(temperature_c)) %>% 
  group_by(site) %>% 
  summarise(mean_daily_max = mean(daily_max))


## Calculate the mean of 4-hourly maximum temperature
## Since measurements were taken every 5 mins, select the highest 12*5*4 = 48 measurements of each day
temp_summary_4hrsmax <- temp_summary %>% 
  group_by(site, date) %>% 
  slice_max(order_by = temperature_c, n = 48, with_ties = FALSE) %>% 
  group_by(site) %>% 
  summarise(mean_daily_4hrmax = mean(temperature_c))


temp_summary_mean_daily_mean <- temp_summary %>% 
  group_by(site, date) %>% 
  summarise(daily_mean = mean(temperature_c)) %>% 
  group_by(site) %>% 
  summarise(mean_daily_mean = mean(daily_mean))

temp_summary <- left_join(temp_summary_std_mean_daily_mean, temp_summary_std_max,
                              by = "site") %>% 
  left_join(., temp_summary_std_4hrsmax, by = "site") %>% 
  left_join(., temp_summary_mean_daily_mean, by = "site") %>% 
  left_join(., temp_summary_max, by = "site") %>% 
  left_join(., temp_summary_4hrsmax, by = "site")
  
  


## Plotting
## Order the bars according to temperatures
neworder <- c("airport_road_1", "end_of_pelly_road", "sewage_dump", 
              "construction_shop", "first_creek", "swimming_hole", "dew_line_road")

## Time series
temporal_trend <- ggplot() +
  geom_line(data = temp_data, aes(x = date_time_mst, y = temperature_c, colour = factor(site, levels = neworder))) +
  scale_colour_brewer(palette = "Accent", 
                      name = "Site", 
                      labels = c("Airport Road 1", "End of Pelly Road",
                               "Sewage-dump", "Construction Shop", "First Creek",
                               "Swimming Hole", "Dew Line Road")) +
  labs(y = "Temperature ºC", x = "") +
  #scale_x_datetime(date_breaks = "3 days", date_minor_breaks = "1 day", date_labels = "%b %d") +
  theme_bw()

temporal_trend

# ggsave("figures/temperature/temp_time_series.png", temporal_trend, width = 10.9, height = 5.47)


## Create a labeller to change the facet_wrap label
summary_stat_names <- list(
  "mean_daily_max_std" = "Mean daily max",
  "mean_daily_4hrmax_std" = "Max daily 4-hr max",
  "mean_daily_mean_std" = "Mean daily mean"
)

summary_stat_labeller <- function(variable,value){
  return(summary_stat_names[value])
}

temp_summary_long <- pivot_longer(temp_summary, cols = 2:7, names_to = "summary_stat", values_to = "temperature_c") %>% 
  mutate(site = factor(site, levels = neworder))

temp_summary_std_plot <- temp_summary_long %>% 
  filter(summary_stat %in% c("mean_daily_max_std", "mean_daily_4hrmax_std", "mean_daily_mean_std")) %>% 
  ggplot() +
  geom_col(aes(x = site, y = temperature_c, 
               fill = site)) +
  facet_wrap(~factor(summary_stat, levels=c("mean_daily_max_std",
                                            "mean_daily_4hrmax_std",
                                            "mean_daily_mean_std")),
             labeller = summary_stat_labeller) +
  scale_fill_brewer(palette = "Accent", 
                    labels = c("Airport Road 1", "End of Pelly Road",
                               "Sewage-dump", "Construction Shop", "First Creek",
                               "Swimming Hole", "Dew Line Road")) +
  labs(y = "Temperature ºC") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

temp_summary_std_plot

# ggsave("figures/temperature/temp_summary_std_plot.png", temp_summary_std_plot, width = 10.9, height = 5.47)


temp_summary_plot <- temp_summary_long %>% 
  filter(summary_stat %in% c("mean_daily_max", "mean_daily_4hrmax", "mean_daily_mean")) %>% 
  ggplot() +
  geom_col(aes(x = site, y = temperature_c, 
               fill = site)) +
  facet_wrap(~factor(summary_stat, levels=c("mean_daily_max",
                                            "mean_daily_4hrmax",
                                            "mean_daily_mean")),
             labeller = summary_stat_labeller) +
  scale_fill_brewer(palette = "Accent", 
                    labels = c("Airport Road 1", "End of Pelly Road",
                               "Sewage-dump", "Construction Shop", "First Creek",
                               "Swimming Hole", "Dew Line Road")) +
  labs(y = "Temperature ºC") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

temp_summary_plot

# ggsave("figures/temperature/temp_summary_plot.png", temp_summary_plot, width = 10.9, height = 5.47)

