## Lilian Chan, University of Guelph
## Arctic vector-borne disease transmission suitability model
##
## Purpose: Process the raw data such that each trait has it own csv file
## 
## Table of content:
##    0. Set-up workspace
##    1. Biting rate (a)
##    2. Vector Competence (bc)
##       i) Infection efficiency (c)
##       ii) Transmission efficiency (b)
##
##    3. Parasite development rate (PDR)
##    4. Lifetime egg production (B)
##    5. Mosquito egg-to-adult survival (pEA)
##       i) Egg viability (EV)
##       ii) Larval survival (pLA)
##
##    6. Mosquito egg-to-adult development rate (MDR)
##    7. Adult mosquito lifespan (lf)

##########
###### 0. Set-up workspace ----
##########

library(tidyverse)
library(readxl)
library(janitor)



##########
###### 1. Biting rate (a) ----
##########

# Read data

## Aedes albopictus ----
a.albopictus.Delatte2009 <- read_csv("data-raw/aedes_albopictus.Delatte2009.csv") %>% 
  clean_names()
str(a.albopictus.Delatte2009)

a.albopictus.Delatte2009  <- a.albopictus.Delatte2009  %>% 
  filter(original_trait_def == "mean duration of gonotrophic cycle") %>% 
  # Add new columns to provide more info
  mutate(trait_name = "1/a") %>% 
  mutate(trait = original_trait_value) %>% 
  mutate(trait_def = "mean duration of gonotrophic cycle") %>% 
  mutate(citation = "Delatte_2009_JMedEntomol") %>% 
  select(trait_name, interactor1temp, trait, original_error_pos, 
         original_error_neg, original_error_unit, interactor1genus, 
         interactor1species, trait_def, citation, doi, figure_table, location)
  

colnames(a.albopictus.Delatte2009) <- c("trait_name", "temp", "trait", 
                                        "error_pos", "error_neg", "error_unit", 
                                        "genus", "species", "trait_def", 
                                        "citation", "doi", "data_source", "notes")



a.albopictus.Marini2020 <- read_csv("data-raw/aedes_albopictus.Marini2020.csv") %>% 
  clean_names() 

a.albopictus.Marini2020 <- a.albopictus.Marini2020 %>% 
  filter(original_trait_def == "mean duration of gonotrophic cycle") %>% 
  # Add new columns to provide more info
  mutate(trait_name = "1/a") %>% 
  mutate(trait = original_trait_value) %>% 
  mutate(trait_def = "mean duration of gonotrophic cycle") %>% 
  mutate(citation = "Marini_2020_Insects") %>% 
  select(trait_name, interactor1temp, trait, original_error_pos, 
         original_error_neg, original_error_unit, interactor1genus, 
         interactor1species, trait_def, citation, doi, figure_table, location)



colnames(a.albopictus.Marini2020) <- c("trait_name", "temp", "trait", 
                                       "error_pos", "error_neg", "error_unit", 
                                       "genus", "species", "trait_def", 
                                       "citation", "doi", "data_source", "notes")


## Aedes aegypti ----
a.aegypti.Goindin2015 <- read_csv("data-raw/aedes_aegypti.Goindin2015.csv") %>% 
  clean_names()


a.aegypti.Goindin2015  <- a.aegypti.Goindin2015  %>% 
  # Add new columns to provide more info
  mutate(trait_name = "1/a") %>% 
  mutate(trait = original_trait_value) %>% 
  mutate(trait_def = original_trait_def) %>% 
  mutate(citation = "Goindin_2015_PLoSOne") %>% 
  select(trait_name, interactor1temp, trait, interactor1genus, 
         interactor1species, trait_def, citation, doi, figure_table, location)


colnames(a.aegypti.Goindin2015) <- c("trait_name", "temp", "trait", 
                                        "genus", "species", "trait_def", 
                                        "citation","doi", "data_source", "notes")


a.aegypti <- read_csv("data-raw/aedes_aegypti.csv") %>% 
  clean_names() 

#a.aegypti$trait[a.aegypti$trait_name == "GCD"] <- 1/a.aegypti$trait[a.aegypti$trait_name == "GCD"]
a.aegypti$trait_name[a.aegypti$trait_name == "GCD"] <- "1/a"

TraitData_a <-  bind_rows(a.albopictus.Delatte2009, a.albopictus.Marini2020, 
                          a.aegypti.Goindin2015, a.aegypti)

# write_csv(TraitData_a, "data-processed/TraitData_a.csv")

## Plot raw data
plot.data.a <- TraitData_a %>%
  mutate(trait = ifelse(trait_name == "1/a", 1/trait, trait)) %>% 
  mutate(type = "non-Arctic") %>% 
  ggplot(aes(x = temp, y = trait)) +
  geom_point(aes(colour = species, shape = citation)) +
  labs(y = "Biting rate (1/days)", x = "Temperature ºC") +
  scale_colour_discrete(name = "Species", labels = c("Ae. aegypti",
                                                     "Ae. albopictus"
  )) +
  scale_shape_discrete(name = "Citation", labels = c("Delatte_2009)",
                                                     "Focks_1993",
                                                     "Focks_2006",
                                                     "Goindin_2015",
                                                     "Marini_2020",
                                                     "Morin_2015"
  )) +
  facet_grid(rows = vars(type)) +
  theme_bw()

plot.data.a

# ggsave("figures/raw_data/plot.data.a.png", plot.data.a, , width = 9.83, height = 6.17)



###### 2. Vector Competence (bc) ----
## Vector competence has two components: infection efficiency (c) and 
## transmission efficiency (b)



##########
###### 2i. Infection efficiency (c) ----
##########

## Ae. Trivittatus (transmitting Dirofilaria immitis); for informative priors) ----
c.trivittatus <- read_csv("data-raw/dirofilaria_immitis_x_aedes_trivittatus.Christensen1978.csv") %>% 
  clean_names() 

c.trivittatus <- c.trivittatus %>% 
  # Add new columns to provide more info
  mutate(trait_name = "c") %>% 
  mutate(trait_def = "infection rate %") %>% 
  # Convert development time to development rate
  mutate(trait = infection_rate_percent/100) %>% 
  mutate(genus = "aedes") %>% 
  mutate(species = "trivittatus ") %>% 
  mutate(paras.species = "Dirofilaria immitis") %>% 
  mutate(citation = "Christensen_1978_ProcHelmintholSocWash") %>% 
  mutate(doi = NA) %>% 
  mutate(data_source = "figure 2") %>% 
  select(trait_name, temp, trait, genus, species, paras.species, trait_def,
         citation, data_source, notes)

# write_csv(c.trivittatus, "data-processed/TraitData_c.csv")


## Plot raw data
plot.data.c <- c.trivittatus %>% 
mutate(type = "non-Arctic") %>% 
  ggplot() +
  geom_point(aes(x = temp, y = trait, colour = species)) +
  labs(y = "Infection probability", x = "Temperature ºC") +
  scale_color_discrete(name = "Species", label = "Ae. trivittatus") +
  facet_grid(rows = vars(type)) +
  theme_bw()

plot.data.c

# ggsave("figures/raw_data/plot.data.c.png", plot.data.c, , width = 9.83, height = 6.17)


##########
###### 2ii. Transmission efficiency (b) ----
##########


##########
###### 3. Parasite development rate (PDR) ----
##########

# Arctic species
## Varestrongylus eleguneniensis ----
PDR.eleguneniensis <- read_csv("data-raw/varestrongylus_eleguneniensis.Kafle2018.csv") %>% 
  clean_names()

PDR.eleguneniensis <- PDR.eleguneniensis %>% 
  # Add new columns to provide more info
  mutate(trait_name = "1/PDR") %>% 
  mutate(trait = dayi_l3) %>% 
  mutate(trait_def = "days first L3 observed") %>% 
  mutate(genus = "varestrongylus") %>% 
  mutate(species = "eleguneniensis") %>% 
  mutate(host.genus = "deroceras") %>% 
  mutate(host.species = "laeve") %>% 
  mutate(citation = "Kafle_2018_ParasitVectors") %>% 
  mutate(doi = "10.1186/s13071-018-2946-x") %>% 
  mutate(data_source = "table 1") %>% 
  mutate(notes = "Nematode infecting caribou and muskoxen in the Canadian Arctic; transmitted by gastropods") %>% 
  dplyr::select(trait_name, temp, trait, genus, species, host.genus, host.species, 
         trait_def, citation, doi, data_source, notes)



## At 8.5ºC, development did not occur after 101 days. Change development time to 1000 days
PDR.eleguneniensis[5,"trait"] <- "1000"
PDR.eleguneniensis$trait <- as.numeric(PDR.eleguneniensis$trait)
PDR.eleguneniensis[5,"notes"] <- "did not develop after 101 days"


## Setaria tundra ----
PDR.tundra <- read_csv("data-raw/setaria_tundra.Laaksonen2009.csv") %>% 
  clean_names()

PDR.tundra <- PDR.tundra %>% 
  # Add new columns to provide more info
  mutate(trait_name = "1/PDR") %>% 
  mutate(trait_def = "days first L3 observed") %>% 
  # Convert development time to development rate
  mutate(trait = parasite_development_time) %>% 
  mutate(genus = "setaria") %>% 
  mutate(species = "tundra") %>% 
  mutate(host.genus = "see notes") %>% 
  mutate(host.species = NA) %>% 
  mutate(citation = "Laaksonen_2009_ParasitVectors") %>% 
  mutate(doi = "10.1186/1756-3305-2-3") %>% 
  mutate(data_source = "result text") %>% 
  select(trait_name, temp, trait, genus, species, host.genus, host.species, 
         trait_def, citation, doi, data_source, notes)


## At 14.1ºC, development was not completed. Change development time to 1000 days
PDR.tundra[1, "trait"] <- 1000


## Dirofilaria immitis (in Ae. Trivittatus; for informing priors) ----
PDR.immitis <- read_csv("data-raw/dirofilaria_immitis_aedes_trivittatus.Christensen1978.csv") %>% 
  clean_names() 

PDR.immitis <- PDR.immitis %>% 
  # Add new columns to provide more info
  mutate(trait_name = "1/PDR") %>% 
  mutate(trait_def = "days first L3 observed") %>% 
  # Convert development time to development rate
  mutate(trait = parasite_development_time_days) %>% 
  mutate(genus = "dirofilaria") %>% 
  mutate(species = "immitis") %>% 
  mutate(host.genus = "aedes") %>% 
  mutate(host.species = "trivittatus") %>% 
  mutate(citation = "Christensen_1978_ProcHelmintholSocWash") %>% 
  mutate(doi = NA) %>% 
  mutate(data_source = "table 1") %>% 
  select(trait_name, temp, trait, genus, species, host.genus, host.species, 
         trait_def, citation, doi, data_source, notes)

## At 14.5ºC, development was not completed. Change development time to 1000 days
PDR.immitis[1,"trait"] <- 1000


## Dirofilaria immitis (in Ae. aegypti; for informative priors) ----
PDR.immitis.aegypti <- read_csv("data-raw/dirofilaria_immitis_aedes_aegypti.Ledesma2015.csv") %>% 
  clean_names() 

PDR.immitis.aegypti <- PDR.immitis.aegypti %>% 
  # Add new columns to provide more info
  mutate(trait_name = "1/PDR") %>% 
  mutate(trait_def = "days first L3 in Malpighian tubules") %>% 
  # Convert development time to development rate
  mutate(trait = days_post_infection_malpighian_tubules) %>% 
  mutate(genus = "dirofilaria") %>% 
  mutate(species = "immitis") %>% 
  mutate(host.genus = "aedes") %>% 
  mutate(host.species = "aegypti") %>% 
  mutate(citation = "Ledesma_2015_VetParasitol") %>% 
  mutate(doi = "10.1016/j.vetpar.2015.02.003") %>% 
  mutate(data_source = "table 1") %>% 
  mutate(notes = NA) %>% 
  select(trait_name, temp, trait, genus, species, host.genus, host.species, 
         trait_def, citation, doi, data_source, notes)


## lymphatic filarisis worms in Ae. polynesiensis
PDR.polynesiensis <- read_csv("data-raw/wuchereria_bancrofti_aedes_polynesiensis.Lardeux1997.csv") %>% 
  clean_names()

PDR.polynesiensis <- PDR.polynesiensis %>% 
  mutate(trait = trait) %>% 
  # Add new columns to provide more info
  mutate(trait_name = "1/PDR") %>% 
  mutate(trait_def = "days to appearance of L3 after experimental infection") %>% 
  mutate(error = sd) %>% 
  mutate(error_unit = "sd") %>% 
  mutate(genus = "wuchereria") %>% 
  mutate(species = "bancrofti") %>% 
  mutate(host.genus = "aedes") %>% 
  mutate(host.species = "polynesiensis") %>% 
  mutate(citation = "Ladeaux_1997_Parasitology") %>% 
  mutate(doi = "10.1017/s0031182096008359") %>% 
  mutate(data_source = "table 1") %>% 
  mutate(notes = NA) %>% 
  select(trait_name, temp, trait, error, error_unit, genus, species, host.genus, 
         host.species, trait_def, citation, doi, data_source, notes)

TraitData_PDR <- bind_rows(PDR.eleguneniensis, PDR.tundra, PDR.immitis, 
                           PDR.immitis.aegypti, PDR.polynesiensis)

## Reorder columns
TraitData_PDR <- TraitData_PDR %>% 
  relocate(error, .after = trait) %>% 
  relocate(error_unit, .after = error)


# write_csv(TraitData_PDR, "data-processed/TraitData_PDR.csv")



##########
###### 4. Lifetime egg production (B)----
##########

EFD.aedes <- read_csv("data/EFD.csv") %>% 
  clean_names()

# write_csv(TraitData_EFD, "data/data-processed/TraitData_EFD.csv")

## Plot raw data
plot.data.EFD <- EFD.aedes %>% 
  mutate(type = "non-Arctic") %>% 
  ggplot(aes(x = temp, y = trait)) +
  geom_point(aes(colour = species, shape = citation)) +
  labs(y = "Fecundity (as eggs per female per day)", x = "Temperature ºC") +
  scale_colour_discrete(name = "Species", labels = c("Ae. aegypti",
                                                     "Ae. albopictus"
  )) +
  # scale_shape_discrete(name = "Citation", labels = c("Beserra_2009)",
  #                                                    "Focks_1993",
  #                                                    "Focks_2006",
  #                                                    "Goindin_2015",
  #                                                    "Marini_2020",
  #                                                    "Morin_2015"
  # )) +
  facet_grid(rows = vars(type)) +
  theme_bw()

plot.data.EFD

# ggsave("figures/raw_data/plot.data.EFD.png", plot.data.EFD, , width = 9.83, height = 6.17)

###### 5. Mosquito egg-to-adult survival (pEA) ----
## pEA is broken down into two parts: egg viability (EV) and Larval survival (pLA)
## pEA = EV * pLA


##########
###### 5i. Egg viability (EV) ----
##########

# Aedes vexans ----
## Read data
EV.vexans <- read_csv("data-raw/aedes_vexans.McHaffey1972.csv") %>% 
  clean_names()

EV.vexans <- EV.vexans %>% 
  # Add new columns to provide more info
  mutate(trait_name = "EV") %>% 
  mutate(trait_def = "percent hatch") %>% 
  mutate(genus = "aedes") %>% 
  mutate(species = "vexans") %>% 
  select(!"host_code") %>% 
  relocate(genus, .after = "trait_2") %>% 
  relocate(species, .after = "genus") %>% 
  relocate(trait_def, .after = "species")

colnames(EV.vexans) <- c("trait_name", "temp", "trait", "error_pos", 
                         "error_neg", "trait2_name", "trait_2", "genus", 
                         "species", "trait_def", "citation", "doi", "data_source",
                         "notes")

## Plot raw data
plot.data.EV <- EV.vexans %>% 
  mutate(type = "Arctic") %>% 
  ggplot() +
  geom_point(aes(x = temp, y = trait, colour = species 
#                 colour = as.factor(trait_2)
                 )) +
  # geom_errorbar(aes(x = temp, ymin = trait - error_neg, ymax = trait + error_pos,
  #                   colour = as.factor(trait_2))) +
  labs(y = "Egg viability (%)", x = "Temperature ºC") +
#  scale_color_discrete(name = "Photoperiod") +
  facet_grid(rows = vars(type)) +
  theme_bw()

plot.data.EV

# ggsave("figures/raw_data/plot.data.EV.png", plot.data.EV, width = 9.83, height = 6.17)


##########
###### 5ii. Larval survival (pLA) ----
##########

# Aedes vexans ----
## Read data
pLA.vexans <- read_csv("data-raw/aedes_vexans.Brust1967.csv") %>% 
  clean_names()

pLA.vexans <- pLA.vexans %>% 
  # Add new columns to provide more info
  mutate(trait_def = "percent reaching the adult stage") %>% 
  mutate(genus = "aedes") %>% 
  mutate(species = "vexans") %>% 
  select(!"host_code") %>% 
  relocate(genus, .after = "trait_2") %>% 
  relocate(species, .after = "genus") %>% 
  relocate(trait_def, .after = "species")
  
colnames(pLA.vexans) <- c("trait_name", "temp", "trait", "error_pos", 
                          "error_neg", "trait2_name", "trait_2", "genus", 
                          "species", "trait_def", "citation", "doi",
                          "data_source", "notes")

pLA.vexans %>% 
  ggplot() +
  geom_point(aes(x = temp, y = trait)) +
  labs(y = "Larval survival (%)", x = "Temperature ºC") +
  theme_bw()


# Aedes spp. ----
pLA.aedes <- read_csv("data/pLA.csv") %>% 
  clean_names()

TraitData_pLA <- rbind(pLA.vexans, pLA.aedes)

# write_csv(TraitData_pLA, "data/data-processed/TraitData_pLA.csv")


##########
###### 6. Mosquito egg-to-adult development rate (MDR) ----
##########

# Read data
## Aedes nigripes ----
MDR.nigripes <- read_excel("data/aedes_nigripes.Culler2015.xlsx", 
                           sheet = "Dev. Time") %>% 
  clean_names() 

# Calculate development rate by 1/development time
MDR.nigripes <- MDR.nigripes %>% 
  mutate(trait = development_time_in_days) %>% 
  # Add new columns to provide more info
  mutate(trait_name = "1/MDR") %>% 
  mutate(trait_def = "days until emergence") %>% 
  mutate(genus = "aedes") %>% 
  mutate(species = "nigripes") %>% 
  mutate(citation = "Culler_2015_ProcRSocB") %>% 
  mutate(data_source = "raw data") %>% 
  select(trait_name, mean_temperature_during_development_c, trait, 
         genus, species, trait_def, citation, data_source, sex)

colnames(MDR.nigripes)[2] <- "temp"
colnames(MDR.nigripes)[9] <- "notes"


## Aedes sierrensis (for informative priors) ----
MDR.sierrensis <- read_csv("data/aedes_sierrensis.Couper2024.csv")

# Select relevant columns
MDR.sierrensis <- MDR.sierrensis %>% 
  clean_names() %>% 
  # Add new columns to provide more info
  mutate(trait_name = "MDR") %>% 
  mutate(trait_def = "1/days") %>% 
  mutate(genus = "aedes") %>% 
  mutate(species = "sierrensis") %>% 
  mutate(citation = "Couper_2024_ProcBiolSci") %>% 
  mutate(data_source = "raw data") %>% # raw data from the paper
  # Combine info from population and sample_id into a new column called "notes"
  unite(population, sample_id, sep = "_", col = "notes") %>% 
  select(trait_name, temp_treatment, juvenile_dev_rate, genus, species, 
         trait_def, citation, data_source, notes) %>% 
  filter(!is.na(juvenile_dev_rate))

# Rename columns
colnames(MDR.sierrensis)[2] <- "temp"
colnames(MDR.sierrensis)[3] <- "trait"


# Combine data from nigripes and sierrensis into a single dataframe
TraitData_MDR <- rbind(MDR.nigripes, MDR.sierrensis)

# write_csv(TraitData_MDR, "data/data-processed/TraitData_MDR.csv")


## Plot raw data
plot.data.MDR <- TraitData_MDR %>% 
  mutate(type = c(rep("Arctic", 75), rep("non-Arctic", 788))) %>% 
  ggplot(aes(x = round(temp,0), y = trait, colour = species)) +
  
  ## Since the Ae. sierrensis has many data, I will just plot the mean±SE
  #geom_point(data = ~filter(.x, type == "Arctic")) +
  stat_summary(fun = mean, geom = "point") +
  stat_summary(fun.data = "mean_se", geom = "errorbar") +
  
  labs(y = "Mosquito adult lifespan (days)", x = "Temperature ºC") +
  scale_colour_discrete(name = "Species", labels = c("Ae. vexans",
                                                     "Ae. sierrensis"
  )) +
  facet_grid(rows = vars(type)) +
  theme_bw()

plot.data.MDR

# ggsave("figures/raw_data/plot.data.MDR.png", plot.data.MDR, width = 9.83, height = 6.17)


##########
###### 7. Adult mosquito lifespan (lf) ----
##########

# Read data
## Aedes vexans ----
lf.vexans <- read_excel("data/aedes_vexans.Costello1971.xlsx") %>% 
  clean_names() 

lf.vexans <- lf.vexans %>% 
  # Choose the highest resource level (least restrictive)
  filter(adult_food == "honey+water" & relative_humidity == "80") %>% 
  # Only select female mosquitoes
  filter(sex == "F") %>% 
  # Add new columns to provide more info
  mutate(trait_name = "lf") %>% 
  mutate(trait_def = "days to 50% mortality") %>% 
  mutate(genus = "aedes") %>% 
  mutate(species = "vexans") %>% 
  mutate(citation = "Costello_1971_JEconEntomol") %>% 
  mutate(data_source = "table 1") %>% 
  mutate(notes = paste0("adult food: ", adult_food, 
                        "; RH: ", relative_humidity,
                        "; sex: ", sex)) %>% 
  select(trait_name, treatment_temperature, days_to_50_percent_mortality, 
         genus, species, trait_def, citation, data_source, notes)
  
  
colnames(lf.vexans)[2] <- "temp"
colnames(lf.vexans)[3] <- "trait"


## Aedes sierrensis (for informative priors) ----
lf.sierrensis <- read_csv("data/aedes_sierrensis.Couper2024.csv")

# Select relevant columns
lf.sierrensis <- lf.sierrensis %>% 
  clean_names() %>% 
  # Add new columns to provide more info
  mutate(trait_name = "lf") %>% 
  mutate(trait_def = "days") %>% 
  mutate(genus = "aedes") %>% 
  mutate(species = "sierrensis") %>% 
  mutate(citation = "Couper_2024_ProcBiolSci") %>% 
  mutate(data_source = "raw data") %>% # raw data from the paper
  # Combine info from population and sample_id into a new column called "notes"
  unite(population, sample_id, sep = "_", col = "notes") %>% 
  select(temp_treatment, adult_lifespan, trait_name, genus, species, trait_def,
         citation, data_source, notes) %>% 
  filter(!is.na(adult_lifespan))

# Rename columns
colnames(lf.sierrensis)[1] <- "temp"
colnames(lf.sierrensis)[2] <- "trait"


# Combine data from vexans and sierrensis into a single dataframe
TraitData_lf <- bind_rows(lf.vexans, lf.sierrensis)

# write_csv(TraitData_lf, "data/data-processed/TraitData_lf.csv")


## Plot raw data
plot.data.lf <- TraitData_lf %>% 
  mutate(type = c(rep("Arctic", 3), rep("non-Arctic", 787))) %>% 
  ggplot(aes(x = temp, y = trait, colour = citation)) +
  # geom_point(aes(colour = citation)) +
  
  ## Since the Ae. sierrensis has many data, I will just plot the mean±SE
  geom_point(data = ~filter(.x, type == "Arctic")) +
  geom_point(data = ~filter(.x, type == "Arctic")) +
  stat_summary(data = ~filter(.x, type == "non-Arctic"),
               fun = mean, geom = "point") +
  stat_summary(data = ~filter(.x, type == "non-Arctic"),
               fun.data = "mean_se", geom = "errorbar") +
  
  labs(y = "Mosquito adult lifespan (days)", x = "Temperature ºC") +
  scale_colour_discrete(name = "Species", labels = c("Ae. vexans",
                                                     "Ae. sierrensis (mean ± SE)"
                                                     )) +
  facet_grid(rows = vars(type)) +
  theme_bw()

plot.data.lf

# ggsave("figures/raw_data/plot.data.lf.png", plot.data.lf, width = 9.83, height = 6.17)
