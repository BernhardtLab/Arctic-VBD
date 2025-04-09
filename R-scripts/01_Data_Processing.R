## Lilian Chan, University of Guelph
## Arctic vector-borne disease transmission suitability model
##
## Purpose: Process the raw data such that each trait has it own csv file
## 
## Table of content:
##    0. Set-up workspace
##    1. Biting rate (a)
##    2. Vector Competence (bc)
##    3. Parasite development rate (PDR)
##    4. Lifetime egg production (B)
##    5. Mosquito egg-to-adult survival (pEA)
##    6. Mosquito egg-to-adult development rate (MDR)
##    7. Adult mosquito lifespan (lf)

##########
###### 0. Set-up workspace ----
##########

library(tidyverse)
library(readxl)
library(janitor)

setwd("~/Documents/UofG/Arctic-VBD")

##########
###### 6. Mosquito egg-to-adult development rate (MDR) ----
##########

# Read data
## Aedes nigripes
MDR.nigripes <- read_excel("data/aedes_nigripes.Culler2015.xlsx", 
                           sheet = "Dev. Time")

MDR.nigripes <- MDR.nigripes %>% 
  clean_names() 

# Calculate development rate by 1/development time
MDR.nigripes <- MDR.nigripes %>% 
  mutate(trait = 1 / development_time_in_days) %>% 
  # Add new columns to provide more info
  mutate(trait_name = "MDR") %>% 
  mutate(genus = "aedes") %>% 
  mutate(species = "nigripes") %>% 
  mutate(citation = "Culler_2015_ProcRSocB") %>% 
  mutate(data_source = "raw_data") %>% 
  select(mean_temperature_during_development_c, trait, 
         trait_name, genus, species, citation, data_source, sex)

colnames(MDR.nigripes)[1] <- "temp"
colnames(MDR.nigripes)[8] <- "notes"

## Aedes sierrensis (for informative priors)
MDR.sierrensis <- read_csv("data/aedes_sierrensis.Couper2024.csv")

# Select relevant columns
MDR.sierrensis <- MDR.sierrensis %>% 
  clean_names() %>% 
  # Add new columns to provide more info
  mutate(trait_name = "MDR") %>% 
  mutate(genus = "aedes") %>% 
  mutate(species = "sierrensis") %>% 
  mutate(citation = "Couper_2024_ProcBiolSci") %>% 
  mutate(data_source = "raw_data") %>% # raw data from the paper
  # Combine info from population and sample_id into a new column called "notes"
  unite(population, sample_id, sep = "_", col = "notes") %>% 
  select(temp_treatment, juvenile_dev_rate, trait_name, genus, species, 
         citation, data_source, notes) %>% 
  filter(!is.na(juvenile_dev_rate))

# Rename columns
colnames(MDR.sierrensis)[1] <- "temp"
colnames(MDR.sierrensis)[2] <- "trait"

# Combine data from nigripes and sierrensis into a single dataframe
TraitData_MDR <- rbind(MDR.nigripes, MDR.sierrensis)

#write_csv(TraitData_MDR, "data/data-processed/TraitData_MDR.csv")

##########
###### 7. Adult mosquito lifespan (lf) ----
##########

## Aedes sierrensis (for informative priors)
MDR.sierrensis <- read_csv("data/aedes_sierrensis.Couper2024.csv")

# Select relevant columns
MDR.sierrensis <- MDR.sierrensis %>% 
  clean_names() %>% 
  # Add new columns to provide more info
  mutate(trait_name = "lf") %>% 
  mutate(genus = "aedes") %>% 
  mutate(species = "sierrensis") %>% 
  mutate(citation = "Couper_2024_ProcBiolSci") %>% 
  mutate(data_source = "raw_data") %>% # raw data from the paper
  # Combine info from population and sample_id into a new column called "notes"
  unite(population, sample_id, sep = "_", col = "notes") %>% 
  select(temp_treatment, AdultLifespan, trait_name, genus, species, 
         citation, data_source, notes) %>% 
  filter(!is.na(juvenile_dev_rate))

# Rename columns
colnames(MDR.sierrensis)[1] <- "temp"
colnames(MDR.sierrensis)[2] <- "trait"
