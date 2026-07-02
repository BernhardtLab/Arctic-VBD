## Lilian Chan, University of Guelph
## Arctic vector-borne disease transmission suitability model
##
## Purpose: Clean and process raw trait data, then export each trait as a separate CSV file.
## 
## Table of content:
##    0. Set-up workspace
##    1. Biting rate (a)
##    2. Vector Competence (bc)
##    3. Pathogen development rate (PDR)
##    4. Egg viability (EV)
##    5. Larval survival (pLA)
##    6. Mosquito egg-to-adult development rate (MDR)
##    7. Adult mosquito lifespan (lf)
##    8. Eggs per female per gonotrophic cycle (EFGC)
##
##
## Inputs:
## csv and Excel files located in 'data-raw', containing published trait estimates
## across temperature gradients. Each file corresponds to a specific study or data source.
##
## Files include: a&EFGC&lf_aedes_spp.Sommerman1969.csv, a_Data_Mordecai2019.csv, 
## a_VecTrait.csv, bc&PDR_dirofilaria_immitis_aedes_vexans.Jankowski1976.csv,
## bc&PDR_dirofilaria_immitis_aedes_trivittatus.Christensen1978.csv, 
## EFGC_aedes_hexodontus.Barlow1955.csv, EV_aedes_triseriatus.Zimmerman2025.csv,
## EV_TraitData_Shocket2020.csv, EV_VecTrait.csv, Fecundity_Data_Mordecai2019.csv,
## lf_aedes_vexans.Costello1971.csv, lf_Data_Mordecai2019.csv, lf_VecTrait.csv,
## MDR&pLA_aedes_flavescens.Trpis1969.csv, MDR_aedes_nigripes.Culler2015.xlsx, 
## MDR_Data_Mordecai2019.csv, PDR_dirofilaria_immitis_aedes_triseriatus_vexans.Fortin1981.csv, 
## PDR_setaria_tundra.Laaksonen2009.csv, PDR_varestrongylus_eleguneniensis.Kafle2018.csv,
## PDR_wuchereria_bancrofti_aedes_polynesiensis.Lardeux1997.csv, pLA_Data_Mordecai2019.csv,
## pLA_VecTrait.csv
##
## Outputs:
## Processed trait datasets saved in 'data-processed/', organized by trait.
## Each output file follows the naming convention: TraitData_<trait>.csv
##
## Generated files:
## TraitData_a.csv, TraitData_bc.csv, TraitData_EFGC.csv, TraitData_EV.csv, 
## TraitData_lf.csv, TraitData_MDR.csv, TraitData_PDR.csv, TraitData_pLA.csv


# 0. Set-up workspace ----------------------------------------------------------
library(tidyverse)
library(readxl)
library(janitor)
library(sp) # For converting latitude and longitude from DMS to decimal


# 1. Biting rate (a) -----------------------------------------------------------

## Arctic species --------------------------------------------------------------
###### Ae. cinereus, Ae. communis, Ae. impiger, Ae. punctor in Alaska ######
a.aedes <- read_csv("data-raw/a&EFGC&lf_aedes_spp.Sommerman1969.csv") %>% 
  clean_names() %>% 
  filter(trait_name == "1/a") 

## Only include replicates with at least 8 individuals
a.aedes <- a.aedes %>% 
  # Separate the notes column into two columns: n and the number of replicates
  separate(notes, sep = " = ", into = c("n", "rep")) %>% 
  mutate(rep = as.numeric(rep)) %>% 
  # Only include replicates with at >=8 individuals
  filter(rep >= 8) %>% 
  # Combine n and rep back together
  unite(n, rep, col = "notes", sep = "=") %>% 
  mutate(type = "Arctic") %>% 
  relocate(type, .after = "longitude")


## Non-Arctic species (for informing priors) -------------------------------

# From VecTrait database
a.VecTrait <- read_csv("data-raw/a_VecTrait.csv") %>% 
  clean_names()

unique(a.VecTrait$interactor1) 
# Ae. aegypti is a tropical to subtropical species; Ae. albopictus is 
# adapted to cooler, temperate climate.


# From data compiled in Mordecaiet al. 2019
a.Mordecai2019 <- read_csv("data-raw/a_Data_Mordecai2019.csv") %>% 
  clean_names()

unique(a.Mordecai2019$host_code) 
# Ae. aegypti is a tropical to subtropical species, and Apse is Anopheles 
# pseudopunctipennis (not Aedes mosquitoes).
# Ae. albopictus is from Delatte et al. 2009, but the VecTrait database already 
# included this data.


###### Ae. albopictus ######
a.aalb <- a.VecTrait  %>% 
  filter(interactor1species == "albopictus") %>%
  filter(original_trait_def == "mean duration of gonotrophic cycle") %>% 
  # Add new columns to provide more info
  mutate(trait_name = "1/a") %>% 
  mutate(type = "non-Arctic") %>% 
  select(trait_name, interactor1temp, original_trait_value, original_error_pos, 
         original_error_neg, original_error_unit, original_trait_def, 
         interactor1genus, interactor1species, citation, doi, figure_table, 
         notes, type, latitude, longitude)


colnames(a.aalb) <- c("trait_name", "temp", "trait", "error_pos", 
                      "error_neg", "error_unit", "trait_def", "genus", 
                      "species", "citation", "doi", "data_source", "notes",
                      "type", "latitude", "longitude")
  


## Combine all data
TraitData_a <-  bind_rows(a.aedes,
                          a.aalb)

## Convert genotrophic cycle duration (1/a) to biting rate (a)
TraitData_a <- TraitData_a %>% 
  mutate(trait = ifelse(trait_name == "1/a", 1/trait, trait)) %>% 
  mutate(trait_name = ifelse(trait_name == "1/a", "a", trait_name)) %>% 
  mutate(notes = ifelse(is.na(notes), "converted from genotrophic cycle duration", 
                        paste0(notes, "; converted from genotrophic cycle duration")))

write_csv(TraitData_a, "data-processed/TraitData_a.csv")


## Plot raw data
plot.data.a <- TraitData_a %>%
  ggplot(aes(x = temp, y = trait)) +
  geom_point(aes(colour = species
                 #, shape = citation
                 ), size = 2) +
  labs(y = "Biting rate (1/days)", x = "Temperature ºC") +
  scale_colour_discrete(name = "Species", labels = c("Ae. albopictus",
                                                     "Ae. cinereus",
                                                     "Ae. communis",
                                                     "Ae. impiger",
                                                     "Ae. punctor"
  )) +
  facet_grid(rows = vars(type)) +
  theme_bw()

plot.data.a

ggsave("figures/raw_data/plot.data.a.png", plot.data.a, , width = 9.83, height = 6.17)




# 2. Vector Competence (bc) -----------------------------------------------
## Vector competence has two components: infection efficiency (c) and 
## transmission efficiency (b)

###### Ae. vexans and Ae. canadensis (transmitting Dirofilaria immitis) ######
bc.vexans <- read_csv("data-raw/bc&PDR_dirofilaria_immitis_aedes_vexans.Jankowski1976.csv") %>% 
  clean_names()

bc.vexans <- bc.vexans %>% 
  filter(trait_name == "bc") %>% 
  # Add new columns to provide more info
  mutate(type = "non-Arctic")


###### Ae. Trivittatus (transmitting Dirofilaria immitis) ######
# This data is infection efficiency (c)
bc.trivittatus <- read_csv("data-raw/bc&PDR_dirofilaria_immitis_aedes_trivittatus.Christensen1978.csv") %>% 
  clean_names()
  

bc.trivittatus <- bc.trivittatus %>% 
  filter(trait_name == "bc") %>% 
  # Add new columns to provide more info
  mutate(trait = trait/100) %>% # convert from % to proportion
  mutate(trait_def = "infection rate %") %>% 
  mutate(type = "non-Arctic")

TraitData_bc <- bind_rows(bc.vexans, bc.trivittatus)

write_csv(TraitData_bc, "data-processed/TraitData_bc.csv")


## Plot raw data
plot.data.bc <- TraitData_bc %>% 
  ggplot() +
  geom_point(aes(x = temp, y = trait, colour = host_species)) +
  labs(y = "vector competence", x = "Temperature ºC") +
  # scale_color_discrete(name = "Species", label = "Ae. trivittatus") +
  facet_grid(rows = vars(type)) +
  theme_bw()

plot.data.bc

ggsave("figures/raw_data/plot.data.bc.png", plot.data.bc, , width = 9.83, height = 6.17)




# 3. Parasite development rate (PDR) -------------------------------------------

## Arctic species --------------------------------------------------------------
###### Varestrongylus eleguneniensis ######
PDR.eleguneniensis <- read_csv("data-raw/PDR_varestrongylus_eleguneniensis.Kafle2018.csv") %>% 
  clean_names()

# This dataset is raw data downloaded from the paper
PDR.eleguneniensis <- PDR.eleguneniensis %>% 
  # Add new columns to provide more info
  mutate(trait_name = "PDR") %>% 
  mutate(trait = dev_rate) %>%
  mutate(error_pos = NA) %>% # this paper did not provide errors
  mutate(error_neg = NA) %>% 
  mutate(error_unit = NA) %>% 
  mutate(trait_def = "1/days first L3 observed") %>% 
  mutate(trait2_name = NA) %>% 
  mutate(trait2 = NA) %>% 
  mutate(paras_genus = "Varestrongylus") %>% 
  mutate(paras_species = "eleguneniensis") %>% 
  mutate(host_genus = "deroceras") %>% 
  mutate(host_species = "laeve") %>% 
  mutate(citation = "Kafle_2018_ParasitVectors") %>% 
  mutate(doi = "10.1186/s13071-018-2946-x") %>% 
  mutate(data_source = "table 1") %>% 
  mutate(notes = NA) %>% 
  mutate(type = "Arctic") %>% 
  mutate(latitude = NA) %>% 
  mutate(longitude = NA) %>% 
  select(trait_name, temp, trait, error_pos, error_neg, error_unit, trait_def,
         trait2_name, trait2, paras_genus, paras_species, host_genus, 
         host_species, citation, doi, data_source, notes, type, latitude, 
         longitude) 



## At 8.5ºC, development did not occur after 101 days. 
## Change development rate to 1/1000 days
PDR.eleguneniensis[5,"trait"] <- 1/1000
PDR.eleguneniensis[5,"notes"] <- "did not develop after 101 days"

# First-stage larvae of V. eleguneniensis were collected from northern Quebec (58.75°N, 68.55°W)
lat_quebec <- 58.75
lon_quebec <- -68.55
PDR.eleguneniensis$latitude <- lat_quebec
PDR.eleguneniensis$longitude <- lon_quebec


###### Setaria tundra ######
PDR.tundra <- read_csv("data-raw/PDR_setaria_tundra.Laaksonen2009.csv") %>% 
  clean_names()

PDR.tundra <- PDR.tundra %>% 
  # Add new columns to provide more info
  mutate(trait_name = "PDR") %>% 
  # Convert development time to development rate
  mutate(trait_def = "1/days first L3 observed") %>% 
  mutate(trait = 1/trait) %>% 
  mutate(type = "Arctic")


## At 14.1ºC, development was not completed. Change development time to 1000 days
PDR.tundra[1, "trait"] <- 1/1000


## Non-Arctic species (for informing priors) -------------------------------

###### Dirofilaria immitis (in Ae. Trivittatus) ######
PDR.immitis.trivittatus <- read_csv("data-raw/bc&PDR_dirofilaria_immitis_aedes_trivittatus.Christensen1978.csv") %>% 
  clean_names() 

PDR.immitis.trivittatus <- PDR.immitis.trivittatus %>% 
  filter(trait_name == "1/PDR") %>% 
  # Convert development time to development rate
  mutate(trait_name = "PDR") %>% 
  mutate(trait_def = "1/days first L3 observed") %>% 
  mutate(trait = 1/trait) %>% 
  mutate(type = "non-Arctic")
  
## At 14.5ºC, development was not completed. Change development rate to 1/1000 days
PDR.immitis.trivittatus[1,"trait"] <- 1/1000


###### Dirofilaria immitis (in Ae. canadensis and Ae. vexans) ###### 
PDR.immitis.vexans <- read_csv("data-raw/bc&PDR_dirofilaria_immitis_aedes_vexans.Jankowski1976.csv") %>% 
  clean_names() 


PDR.immitis.vexans <- PDR.immitis.vexans %>% 
  filter(trait_name == "1/PDR") %>% 
  # Convert development time to development rate
  mutate(trait_name = "PDR") %>% 
  mutate(trait_def = "1/days L3 larvae first observed") %>% 
  mutate(trait = 1/trait) %>% 
  mutate(type = "non-Arctic")
  


###### Dirofilaria immitis (in Ae. triseriatus ) ######
PDR.immitis.triseriatus <- read_csv("data-raw/PDR_dirofilaria_immitis_aedes_triseriatus_vexans.Fortin1981.csv") %>% 
  clean_names() 


PDR.immitis.triseriatus <- PDR.immitis.triseriatus %>% 
  # Convert development time to development rate
  mutate(trait_name = "PDR") %>% 
  mutate(trait_def = "1/days L3 larvae first observed in mouthparts") %>%
  mutate(trait = 1/trait) %>% 
  mutate(type = "non-Arctic")


## Development was not completed at 12, 14, and 16ºC. 
## Change development rate to 1/1000 days
PDR.immitis.triseriatus[1:3, "trait"] <- 1/1000



######  lymphatic filarisis worms (in Ae. polynesiensis) ###### 
PDR.bancrofti <- read_csv("data-raw/PDR_wuchereria_bancrofti_aedes_polynesiensis.Lardeux1997.csv") %>% 
  clean_names()

PDR.bancrofti <- PDR.bancrofti %>% 
  mutate(trait = trait) %>% 
  # Add new columns to provide more info
  # Convert development time to development rate
  mutate(error_pos = sd/trait^2) %>% 
  mutate(trait = 1/trait) %>% 
  mutate(error_unit = "sd") %>% 
  mutate(trait_def = "1/days to appearance of L3 after experimental infection") %>% 
  mutate(paras_genus = "Wuchereria") %>% 
  mutate(paras_species = "bancrofti") %>% 
  mutate(host_genus = "Aedes") %>% 
  mutate(host_species = "polynesiensis") %>% 
  mutate(citation = "Lardeux_1997_Parasitology") %>% 
  mutate(doi = "10.1017/s0031182096008359") %>% 
  mutate(data_source = "table 1") %>% 
  mutate(notes = NA) %>% 
  mutate(type = "non-Arctic") %>% 
  mutate(latitude = NA) %>% 
  mutate(longitude = NA) %>% 
  select(trait_name, temp, trait, error_pos, error_unit, paras_genus, 
         paras_species, host_genus, host_species, trait_def, citation, doi, 
         data_source, notes, type, latitude, longitude) 

# First-stage larvae of Wuchereria bancrofti were collected from French Polynesia
lat_frpolynesia <- -17.6
lon_frpolynesia <- -149.5
PDR.bancrofti$latitude <- lat_frpolynesia
PDR.bancrofti$longitude <- lon_frpolynesia


TraitData_PDR <- bind_rows(PDR.eleguneniensis, PDR.tundra, 
                           PDR.immitis.triseriatus, PDR.immitis.trivittatus,
                           PDR.immitis.vexans, PDR.bancrofti)


write_csv(TraitData_PDR, "data-processed/TraitData_PDR.csv")


# Plot the raw data
plot.data.PDR <- TraitData_PDR %>% 
  ggplot(aes(x = temp, y = trait)) +
  geom_point(aes(colour = host_species)) +
  labs(y = "Parasite development rate (1/days)", x = "Temperature ºC") +
  scale_colour_discrete(name = "Species", labels = c("D. immitis in Ae. canadensis",
                                                     "V. eleguneniensis",
                                                     "W. bancrofti",
                                                     "D. immitis in Ae. triseriatus",
                                                     "D. immitis in Ae. trivittatus",
                                                     "D. immitis in Ae. vexans",
                                                     "S. tundra"
  )) +
  facet_grid(rows = vars(type)) +
  theme_bw()

plot.data.PDR

ggsave("figures/raw_data/plot.data.PDR.png", plot.data.PDR, , width = 9.83, height = 6.17)




# 4. Egg viability (EV) ---------------------------------------------------

EV.Shocket2020 <- read_csv("data-raw/EV_TraitData_Shocket2020.csv") %>% 
  clean_names()

unique(EV.Shocket2020$host_code)

## Ae. vexans can be found in subarctic regions, while Ae. dorsalis and Ae. 
## nigromaculis are temperate species.


## Arctic species --------------------------------------------------------------
###### Ae. vexans ######

EV.vexans <- EV.Shocket2020 %>% 
  filter(host_code == "Avex") %>% 
  # Add new columns to provide more info
  mutate(trait_name = "EV") %>% 
  mutate(trait_def = "percent hatch") %>% 
  mutate(genus = "Aedes") %>% 
  mutate(species = "vexans") %>% 
  mutate(error_unit = "95% confidence interval") %>%  
  mutate(citation = "McHaffey_1972_JMedEntomol") %>% 
  mutate(doi = "10.1093/jmedent/9.6.564") %>% 
  mutate(type = "Arctic") %>% 
  mutate(latitude = NA) %>% 
  mutate(longitude = NA) %>% 
  select(trait_name, t, trait, error_pos_si, error_neg_si, error_unit, trait_def,
         trait2_name, trait_2, genus, species, citation, doi, figure, notes, 
         type, latitude, longitude) 



colnames(EV.vexans) <- c("trait_name", "temp", "trait", "error_pos", 
                         "error_neg", "error_unit", "trait_def", "trait2_name", 
                         "trait2", "genus", "species", "citation", "doi", 
                         "data_source", "notes", "type", "latitude", "longitude")

# Mosquitoes from McHaffey 1972 were collected from central Washington, USA
lat_washington <- 46
lon_washington <- -119
EV.vexans$latitude[EV.vexans$citation == "McHaffey_1972_JMedEntomol"] <- lat_washington
EV.vexans$longitude[EV.vexans$citation == "McHaffey_1972_JMedEntomol"] <- lon_washington


## Non-Arctic species (for informing priors) -----------------------------------

###### Ae. dorsalis ######

EV.dorsalis <- EV.Shocket2020 %>% 
  filter(host_code == "Ador") %>% 
  # Add new columns to provide more info
  mutate(trait_name = "EV") %>% 
  mutate(trait_def = "percent hatch") %>% 
  mutate(genus = "Aedes") %>% 
  mutate(species = "dorsalis") %>% 
  mutate(error_unit = "95% confidence interval") %>%  
  mutate(citation = "McHaffey_1970_JMedEntomol") %>% 
  mutate(doi = "10.1093/jmedent/7.6.631") %>% 
  mutate(type = "non-Arctic") %>% 
  mutate(latitude = NA) %>% 
  mutate(longitude = NA) %>% 
  select(trait_name, t, trait, error_pos_si, error_neg_si, error_unit, trait_def,
         trait2_name, trait_2, genus, species, citation, doi, figure, notes, 
         type, latitude, longitude) 

colnames(EV.dorsalis) <- c("trait_name", "temp", "trait", "error_pos", 
                           "error_neg", "error_unit", "trait_def", "trait2_name", 
                           "trait2", "genus", "species", "citation", "doi", 
                           "data_source", "notes", "type", "latitude", "longitude")


# Mosquitoes from McHaffey 1970 were collected in central Washington, USA
EV.dorsalis$latitude[EV.dorsalis$citation == "McHaffey_1970_JMedEntomol"] <- lat_washington
EV.dorsalis$longitude[EV.dorsalis$citation == "McHaffey_1970_JMedEntomol"] <- lon_washington


###### Ae. nigromaculis ######
EV.nigromaculis <- EV.Shocket2020 %>% 
  filter(host_code == "Anig") %>% 
  # Add new columns to provide more info
  mutate(trait_name = "EV") %>% 
  mutate(trait_def = "percent hatch") %>% 
  mutate(genus = "Aedes") %>% 
  mutate(species = "nigromaculis") %>% 
  mutate(error_unit = "95% confidence interval") %>%  
  mutate(citation = "McHaffey_1972_MosqNews") %>% 
  mutate(doi = "10.5281/zenodo.16126961") %>% 
  mutate(type = "non-Arctic") %>% 
  mutate(latitude = NA) %>% 
  mutate(longitude = NA) %>% 
  select(trait_name, t, trait, error_pos_si, error_neg_si, error_unit, trait_def,
         trait2_name, trait_2, genus, species, citation, doi, figure, notes, 
         type, latitude, longitude) 

colnames(EV.nigromaculis) <- c("trait_name", "temp", "trait", "error_pos", 
                               "error_neg", "error_unit", "trait_def", "trait2_name", 
                               "trait2", "genus", "species", "citation", "doi", 
                               "data_source", "notes", "type", "latitude", "longitude")


# Mosquitoes from McHaffey 1972 were collected in central Washington, USA
EV.nigromaculis$latitude[EV.nigromaculis$citation == "McHaffey_1972_MosqNews"] <- lat_washington
EV.nigromaculis$longitude[EV.nigromaculis$citation == "McHaffey_1972_MosqNews"] <- lon_washington


###### Ae. triseriatus ######
EV.triseriatus <- read_csv("data-raw/EV_aedes_triseriatus.Zimmerman2025.csv") %>% 
  clean_names()

EV.triseriatus <- EV.triseriatus %>% 
  # Add new columns to provide more info
  mutate(trait_name = "EV") %>% 
  mutate(error_unit = "IQR") %>%  
  mutate(trait_def = "total percentage of eggs hatching across all six flood rounds") %>% 
  mutate(citation = "Zimmerman_2025_JVectorEcol") %>% 
  mutate(doi = "10.52707/1081-1710-50.1-s1") %>% 
  mutate(data_source = "figure 1 insert") %>% 
  mutate(type = "non-Arctic") %>% 
  mutate(latitude = NA) %>% 
  mutate(longitude = NA) %>% 
  select(trait_name, temp, trait, error_pos, error_neg, error_unit, trait_def,
         genus, species, citation, doi, data_source, notes, type, latitude, 
         longitude) 


# Mosquitoes from Zimmerman 2025 were collected from Mercer County, NJ, USA
lat_mercer <- 40.3
lon_mercer <- -74.7
EV.triseriatus$latitude[EV.triseriatus$citation == "Zimmerman_2025_JVectorEcol"] <- lat_mercer
EV.triseriatus$longitude[EV.triseriatus$citation == "Zimmerman_2025_JVectorEcol"] <- lon_mercer


###### Ae. albopictus (from VecTrait database) ######
EV.VecTrait <- read_csv("data-raw/EV_VecTrait.csv") %>% 
  clean_names()

unique(EV.VecTrait$interactor1)

EV.VecTrait <- EV.VecTrait  %>% 
  # Only want egg stage survival
  filter(interactor1stage == "egg") %>% 
  # Add new columns to provide more info
  mutate(trait_name = "EV") %>% 
  mutate(trait2_name = ifelse(
    is.na(second_stressor_value), NA, # if no second_stressor, leave trait2_name empty
    paste0(second_stressor_def, " (", second_stressor_unit, ")"))) %>% # Combine second_stressor and the unit
  mutate(type = "non-Arctic") %>% 
  select(trait_name, interactor1temp, original_trait_value, original_error_pos, 
         original_error_neg, original_error_unit, original_trait_def, 
         trait2_name, second_stressor_value, interactor1genus, 
         interactor1species, citation, doi, figure_table, notes, type, 
         latitude, longitude)


colnames(EV.VecTrait) <- c("trait_name", "temp", "trait", "error_pos", 
                           "error_neg", "error_unit", "trait_def", "trait2_name", 
                           "trait2", "genus", "species", "citation", "doi", 
                           "data_source", "notes", "type", "latitude", "longitude")

## Convert percentage to proportion
EV.VecTrait <- EV.VecTrait %>% 
  mutate(trait = ifelse(
    trait_def == "percentage of individuals surviving life stage (emergence rate)",
    trait/100,
    trait
  )) %>% 
  ## Change trait_def to proportion as well
  mutate(trait_def = "proportion of individuals surviving life stage (emergence rate)")

## Change the trait2 column to character (so that it can combine with other dataset)
EV.VecTrait$trait2 <- as.character(EV.VecTrait$trait2) 


# Mosquitoes from Blagrove et al 2013 were from the Ascoli strain (Italy), where were was colonized from San Benedetto del Tronto, Italy (Blagrove et al. 2011)
lat_italy <- 42.9
lon_italy <- 13.9
EV.VecTrait$latitude[EV.VecTrait$doi == "10.1371/journal.pntd.0002152"] <- lat_italy
EV.VecTrait$longitude[EV.VecTrait$doi == "10.1371/journal.pntd.0002152"] <- lon_italy

TraitData_EV <- bind_rows(EV.vexans, EV.dorsalis, EV.nigromaculis, 
                          EV.triseriatus, EV.VecTrait)


write_csv(TraitData_EV, "data-processed/TraitData_EV.csv")


## Plot raw data
plot.data.EV <- TraitData_EV %>% 
  ggplot() +
  geom_point(aes(x = temp, y = trait, colour = species 
                 #                 colour = as.factor(trait_2)
  )) +
  # geom_errorbar(aes(x = temp, ymin = trait - error_neg, ymax = trait + error_pos,
  #                   colour = as.factor(trait_2))) +
  labs(y = "Egg viability (%)", x = "Temperature ºC") +
  scale_colour_discrete(name = "species", labels = c("Ae. albopictus",
                                                     "Ae. dorsalis",
                                                     "Ae. nigromaculis",
                                                     "Ae. triseriatus",
                                                     "Ae. vexans"
  )) +
  facet_grid(rows = vars(type)) +
  theme_bw()

plot.data.EV

ggsave("figures/raw_data/plot.data.EV.png", plot.data.EV, width = 9.83, height = 6.17)


# 5. Larval-to-adult survival (pLA) ---------------------------------------

## Arctic species ----------------------------------------------------------
pLA.Mordecai2019 <- read_csv("data-raw/pLA_Data_Mordecai2019.csv") %>% 
  clean_names()

unique(pLA.Mordecai2019$host_code)
## Ae. vexans is found in subarctic regions, while Ae. nigromaculis, Ae. 
## sollicitans, and Ae. triseriatus are temperate species.


###### Ae. vexans ######
pLA.vexans <- pLA.Mordecai2019 %>% 
  ## get data from Ae. vexans
  filter(host_code == "Avex") %>% 
  # select columns that we need
  select(trait_name, t, trait, error_pos_si, error_neg_si, trait2_name, 
         trait_2, citation, figure, notes) %>% 
  # Add new columns to provide more info
  mutate(error_unit = NA) %>% 
  relocate(error_unit, .after = "error_neg_si") %>% # rearrange columns
  mutate(trait_def = NA) %>% 
  relocate(trait_def, .after = "error_unit") %>% 
  mutate(genus = "Aedes") %>% 
  relocate(genus, .after = "trait_2") %>% 
  mutate(species = "vexans") %>% 
  relocate(species, .after = "genus") %>% 
  mutate(doi = NA) %>% 
  relocate(doi, .after = "citation") %>% 
  mutate(type = "Arctic") %>% 
  mutate(latitude = NA) %>% 
  mutate(longitude = NA)


## Rename columns
colnames(pLA.vexans) <- c("trait_name", "temp", "trait", "error_pos", "error_neg",
                         "error_unit", "trait_def", "trait2_name", "trait2",
                         "genus", "species", "citation","doi", "data_source", 
                         "notes", "type", "latitude", "longitude")

## Provide more info on the papers
### Brust 1967
pLA.vexans$doi[pLA.vexans$citation == "Brust_1967_TheCanadianEntomologist"] <-
  "10.4039/ent99986-9"

## Trpis and Shemanchuk 1970
pLA.vexans$doi[pLA.vexans$citation == "Trpis&Shemanchuk_1970_TheCanadianEntomologist"] <-
  "10.4039/ent1021048-8"


# Mosquitoes from Brust 1967 were collected near Winnipeg, Manitoba;
lat_winnipeg <- 49.9
lon_winnipeg <- -97.1
pLA.vexans$latitude[pLA.vexans$citation == "Brust_1967_TheCanadianEntomologist"] <- lat_winnipeg
pLA.vexans$longitude[pLA.vexans$citation == "Brust_1967_TheCanadianEntomologist"] <- lon_winnipeg


# Mosquitoes from Trpis and Shemanchuk 1970 were collected near Lethbridge, Alberta
lat_lethbridge <- 49.7
long_lethbridge <- -112.8
pLA.vexans$latitude[pLA.vexans$citation == "Trpis&Shemanchuk_1970_TheCanadianEntomologist"] <- lat_lethbridge
pLA.vexans$longitude[pLA.vexans$citation == "Trpis&Shemanchuk_1970_TheCanadianEntomologist"] <- long_lethbridge


###### Ae. flavescens ######
pLA.flavescens <- read_csv("data-raw/MDR&pLA_aedes_flavescens.Trpis1969.csv") %>% 
  clean_names()

pLA.flavescens <- pLA.flavescens %>% 
  filter(trait_name == "pLA") %>% 
  # Add new columns to provide more info
  mutate(type = "Arctic")


## Non-Arctic species (for informing priors) -------------------------------

###### Ae. nigromaculis ######
pLA.nigromaculis <- pLA.Mordecai2019 %>% 
  filter(host_code == "Anig") %>% 
  # select columns that we need
  select(trait_name, t, trait, error_pos_si, error_neg_si, trait2_name, 
         trait_2, citation, figure, notes) %>% 
  # Add new columns to provide more info
  mutate(genus = "Aedes") %>% 
  relocate(genus, .after = "trait_2") %>% 
  mutate(species = "nigromaculis") %>% 
  relocate(species, .after = "genus") %>% 
  mutate(doi = NA) %>% 
  relocate(doi, .after = "citation") %>% 
  mutate(type = "non-Arctic") %>% 
  mutate(latitude = NA) %>% 
  mutate(longitude = NA)


## Rename columns
colnames(pLA.nigromaculis) <- c("trait_name", "temp", "trait", "error_pos", 
                                "error_neg","trait2_name", "trait2", "genus", 
                                "species", "citation", "doi", "data_source", 
                                "notes", "type", "latitude", "longitude")


## Provide more info on the papers
### Brust 1967
pLA.nigromaculis$doi[pLA.nigromaculis$citation == "Brust_1967_TheCanadianEntomologist"] <-
  "10.4039/ent99986-9"

# Mosquitoes from Brust 1967 were collected near Winnipeg, Manitoba;
lat_winnipeg <- 49.9
lon_winnipeg <- -97.1
pLA.nigromaculis$latitude[pLA.nigromaculis$citation == "Brust_1967_TheCanadianEntomologist"] <- lat_winnipeg
pLA.nigromaculis$longitude[pLA.nigromaculis$citation == "Brust_1967_TheCanadianEntomologist"] <- lon_winnipeg



###### Ae. sollicitans ######
pLA.sollicitans <- pLA.Mordecai2019 %>% 
  filter(host_code == "Asol") %>% 
  # select columns that we need
  select(trait_name, t, trait, error_pos_si, error_neg_si, trait2_name, 
         trait_2, citation, figure, notes) %>% 
  # Add new columns to provide more info
  mutate(genus = "Aedes") %>% 
  relocate(genus, .after = "trait_2") %>% 
  mutate(species = "sollicitans") %>% 
  relocate(species, .after = "genus") %>% 
  mutate(doi = NA) %>% 
  relocate(doi, .after = "citation") %>% 
  mutate(type = "non-Arctic") %>% 
  mutate(latitude = NA) %>% 
  mutate(longitude = NA)


## Rename columns
colnames(pLA.sollicitans) <- c("trait_name", "temp", "trait", "error_pos", 
                                "error_neg","trait2_name", "trait2", "genus", 
                                "species", "citation", "doi", "data_source", 
                                "notes", "type", "latitude", "longitude")


# Shelton 1973 was conducted in Lake Charles, Louisiana
lat_lakecharles <- 30.2
lon_lakecharles <- -93.2
pLA.sollicitans$latitude[pLA.sollicitans$citation == "Shelton_1973_MosquitoNews"] <- lat_lakecharles
pLA.sollicitans$longitude[pLA.sollicitans$citation == "Shelton_1973_MosquitoNews"] <- lon_lakecharles


###### Ae. triseriatus ######
pLA.triseriatus <- pLA.Mordecai2019 %>% 
  filter(host_code == "Atri") %>% 
  # select columns that we need
  select(trait_name, t, trait, error_pos_si, error_neg_si, trait2_name, 
         trait_2, citation, figure, notes) %>% 
  # Add new columns to provide more info
  mutate(genus = "Aedes") %>% 
  relocate(genus, .after = "trait_2") %>% 
  mutate(species = "triseriatus") %>% 
  relocate(species, .after = "genus") %>% 
  mutate(doi = NA) %>% 
  relocate(doi, .after = "citation") %>% 
  mutate(type = "non-Arctic") %>% 
  mutate(latitude = NA) %>% 
  mutate(longitude = NA)


## Rename columns
colnames(pLA.triseriatus) <- c("trait_name", "temp", "trait", "error_pos", 
                               "error_neg", "trait2_name", "trait2", "genus", 
                               "species", "citation", "doi", "data_source", "notes", 
                               "type", "latitude", "longitude")


## Provide more info on the papers
### Teng and Apperson 2000
pLA.triseriatus$doi[pLA.triseriatus$citation == "TengApperson_2000_J.Med.Entomol."] <-
  "10.1603/0022-2585-37.1.40"

# Mosquitoes from Teng and Apperson 2000 were collected near the town of Cherokee, Jackson County, NC, USA
lat_cherokee <- 35.5
lon_cherokee <- -83.3
pLA.triseriatus$latitude[pLA.triseriatus$citation == "TengApperson_2000_J.Med.Entomol."] <- lat_cherokee
pLA.triseriatus$longitude[pLA.triseriatus$citation == "TengApperson_2000_J.Med.Entomol."] <- lon_cherokee


# Shelton 1973 was conducted in Lake Charles, Louisiana
pLA.triseriatus$latitude[pLA.triseriatus$citation == "Shelton_1973_MosquitoNews"] <- lat_lakecharles
pLA.triseriatus$longitude[pLA.triseriatus$citation == "Shelton_1973_MosquitoNews"] <- lon_lakecharles

# From VecTrait database
pLA.VecTrait <- read_csv("data-raw/pLA_VecTrait.csv") %>% 
  clean_names()

unique(pLA.VecTrait$interactor1) 


###### Ae. albopictus ######
pLA.albopictus <- pLA.VecTrait  %>% 
  filter(interactor1species == "albopictus") %>%
  filter(str_detect(interactor1stage, "juvenile")) %>% 
  # Add new columns to provide more info
  mutate(trait_name = "pLA") %>% 
  mutate(type = "non-Arctic") %>% 
  select(trait_name, interactor1temp, original_trait_value, original_error_pos, 
         original_error_neg, original_error_unit, original_trait_def, 
         interactor1genus, interactor1species, citation, doi, figure_table, 
         notes, type, latitude, longitude)

colnames(pLA.albopictus) <- c("trait_name", "temp", "trait", "error_pos", 
                              "error_neg", "error_unit", "trait_def", "genus", 
                              "species", "citation", "doi", "data_source", 
                              "notes", "type", "latitude", "longitude")

# Convert percentage to proportion
pLA.albopictus <- pLA.albopictus %>% 
  mutate(trait = ifelse(trait_def == "percent surviving life stage", trait/100, trait))
pLA.albopictus$trait_def[pLA.albopictus$trait_def == "percent surviving life stage"] = "proportion" 


TraitData_pLA <- bind_rows(pLA.vexans, pLA.flavescens,
                           pLA.nigromaculis, pLA.sollicitans, pLA.triseriatus, 
                           pLA.albopictus)


write_csv(TraitData_pLA, "data-processed/TraitData_pLA.csv")


# Plot the raw data
plot.data.pLA <- TraitData_pLA %>% 
  ggplot(aes(x = temp, y = trait)) +
  geom_point(aes(colour = species)) +
  labs(y = "Larval survival (%)", x = expression(paste("Temperature (", degree, "C)"))) +
  scale_colour_discrete(name = "Species", labels = c("Ae. albopictus",
                                                     "Ae. flavescens",
                                                     "Ae. nigromaculis",
                                                     "Ae. sollicitans",
                                                     "Ae. triseriatus",
                                                     "Ae. vexans"
  )) +
  facet_grid(rows = vars(type)) +
  theme_bw()

plot.data.pLA

ggsave("figures/raw_data/plot.data.pLA.png", plot.data.pLA, width = 9.83, height = 6.17)



# 6. Mosquito egg-to-adult development rate (MDR) -------------------------

## Arctic species ----------------------------------------------------------
###### Ae. nigripes ######
MDR.nigripes <- read_excel("data-raw/MDR_aedes_nigripes.Culler2015.xlsx", 
                           sheet = "Dev. Time") %>% 
  clean_names() 

# This dataset is raw data downloaded from the paper
MDR.nigripes <- MDR.nigripes %>% 
  mutate(trait = development_time_in_days) %>% 
  filter(sex == "female") %>% # Only female data is used
  # Add new columns to provide more info
  mutate(trait = 1/trait) %>% # convert from development time to development rate
  mutate(trait_name = "MDR") %>% 
  mutate(error_pos = NA) %>% 
  mutate(error_neg = NA) %>% 
  mutate(error_unit = NA) %>% 
  mutate(trait_def = "1/days until emergence") %>% 
  mutate(trait2_name = NA) %>% 
  mutate(trait2 = NA) %>% 
  mutate(genus = "Aedes") %>% 
  mutate(species = "nigripes") %>% 
  mutate(citation = "Culler_2015_ProcRSocB") %>% 
  mutate(doi = "10.1098/rspb.2015.1549") %>% 
  mutate(data_source = "raw data") %>% 
  mutate(notes = "converted from development time") %>% 
  mutate(type = "Arctic") %>% 
  select(trait_name, mean_temperature_during_development_c, trait, error_pos, 
         error_neg, error_unit, trait_def, trait2_name, trait2,
         genus, species, citation, doi, data_source, notes, type)


colnames(MDR.nigripes)[2] <- "temp"

read_excel("data-raw/MDR_aedes_nigripes.Culler2015.xlsx", 
                           sheet = "Pond Locations")
  
# Mosquitoes were collected in Black pond, 
lat_greenland <- as.numeric(char2dms("67d0'33.6\"N"))
lon_greenland <- as.numeric(char2dms("50d40'17.9\"W"))


MDR.nigripes <- MDR.nigripes %>% 
  mutate(latitude = lat_greenland) %>% 
  mutate(longitude = lon_greenland) 


###### Ae. flavescens ######
MDR.flavescens <- read_csv("data-raw/MDR&pLA_aedes_flavescens.Trpis1969.csv") %>% 
  clean_names()

MDR.flavescens <- MDR.flavescens %>% 
  filter(trait_name == "1/MDR") %>%
  # Convert from development time to development rate
  mutate(trait_name = "MDR") %>% 
  mutate(trait = 1/trait) %>% 
  mutate(trait_def = "1/median pre-adult development time (days)") %>% 
  # Add new columns to provide more info
  mutate(type = "Arctic")

## Development was not completed at 5 and 30ºC. 
## Change development rate to 1/1000 days-1
MDR.flavescens[c(1,6), "trait"] <- 1/1000


###### Ae. vexans ######
## Read Mordecai et al 2019 data
MDR.Mordecai2019 <- read_csv("data-raw/MDR_Data_Mordecai2019.csv") %>% 
  clean_names() %>% 
  # select columns that we need
  select(trait_name, t_c, trait, error_pos_si, error_neg_si, trait2_name, 
         trait_2, host_code, citation, figure, notes)
  
unique(MDR.Mordecai2019$host_code)
## Ae. vexans can be found in subarctic regions, while Ae. nigromaculis, Ae. 
## sollicitans, and Ae. triseriatus are temperate species.


MDR.vexans <- MDR.Mordecai2019 %>% 
  filter(host_code  == "Avex") %>% 
  # Convert from development time to development rate
  mutate(trait_name = "MDR") %>% 
  mutate(error_pos_si = error_pos_si/trait^2) %>% 
  mutate(error_neg_si = error_neg_si/trait^2) %>% 
  mutate(trait = 1/trait) %>% 
  # Add new columns to provide more info
  mutate(genus = "Aedes") %>% 
  relocate(genus, .after = "trait_2") %>% 
  mutate(species = "vexans") %>% 
  relocate(species, .after = "genus") %>% 
  mutate(type = "Arctic") %>% 
  select(!host_code) %>%   
  mutate(latitude = NA) %>% 
  mutate(longitude = NA)
  


## Rename columns
colnames(MDR.vexans) <- c("trait_name", "temp", "trait", "error_pos", 
                          "error_neg","trait2_name", "trait2", "genus", 
                          "species", "citation", "data_source", "notes", 
                          "type", "latitude", "longitude")


# Mosquitoes from Brust 1967 were collected near Winnipeg, Manitoba;
lat_winnipeg <- 49.9
lon_winnipeg <- -97.1
MDR.vexans$latitude[MDR.vexans$citation == "Brust_1967_CanEntomol"] <- lat_winnipeg
MDR.vexans$longitude[MDR.vexans$citation == "Brust_1967_CanEntomol"] <- lon_winnipeg

# Mosquitoes from Trpis and Shemanchuk 1970 were collected near Lethbridge, Alberta
lat_lethbridge <- 49.7
long_lethbridge <- -112.8
MDR.vexans$latitude[MDR.vexans$citation == "Trpis&Shemanchuk_1970_TheCanadianEntomologist"] <- lat_lethbridge
MDR.vexans$longitude[MDR.vexans$citation == "Trpis&Shemanchuk_1970_TheCanadianEntomologist"] <- long_lethbridge



## Non-Arctic species (for informing priors) -------------------------------

###### Ae. albopictus ######
MDR.albopictus <- MDR.Mordecai2019 %>% 
  filter(host_code  == "Aalb") %>% 
  # Add new columns to provide more info
  mutate(genus = "Aedes") %>% 
  relocate(genus, .after = "trait_2") %>% 
  mutate(species = "albopictus") %>% 
  relocate(species, .after = "genus") %>% 
  mutate(type = "non-Arctic") %>% 
  select(!host_code) %>% 
  mutate(latitude = NA) %>% 
  mutate(longitude = NA)


## Rename columns
colnames(MDR.albopictus) <- c("trait_name", "temp", "trait", "error_pos", 
                              "error_neg","trait2_name", "trait2", "genus", 
                              "species", "citation", "data_source", "notes", 
                              "type", "latitude", "longitude")


# Mosquitoes from Alto 2001 were collected from East St. Louis, IL
lat_stlouis <- 38.6
lon_stlouis <- -90.1
MDR.albopictus$latitude[MDR.albopictus$citation == "Alto_et_al_2001b_JMedEnto"] <- lat_stlouis
MDR.albopictus$longitude[MDR.albopictus$citation == "Alto_et_al_2001b_JMedEnto"] <- lon_stlouis

# Mosquitoes from Briegel 2001 were collected from Louisiana USA
lat_LA <- 30.7
lon_LA <- -92.1
MDR.albopictus$latitude[MDR.albopictus$citation == "Briegel_etal_2001_JMedEnto"] <- lat_LA
MDR.albopictus$longitude[MDR.albopictus$citation == "Briegel_etal_2001_JMedEnto"] <- lon_LA


# Mosquitoes from Calado 2002 were collected from Registro, São Paulo State, Brazil
lat_registro <- -24.5
lon_registro <- -47.9
MDR.albopictus$latitude[MDR.albopictus$citation == "Calado_Navarro-Silva_2002b_RevSaudaPublica"] <- lat_registro
MDR.albopictus$longitude[MDR.albopictus$citation == "Calado_Navarro-Silva_2002b_RevSaudaPublica"] <- lon_registro

# Mosquitoes from Delatte 2009 were collected from Saint Pierre, La Reunion Island, France
lat_stpierre <- -21.3
lon_stpierre <- 55.5
MDR.albopictus$latitude[MDR.albopictus$citation == "Delatte_etal_2009_JMedEnto"] <- lat_stpierre
MDR.albopictus$longitude[MDR.albopictus$citation == "Delatte_etal_2009_JMedEnto"] <- lon_stpierre

# Mosquitoes from Ezeakacha 2015 were collected in and around Hattiesburg, MS, USA
lat_Hattiesburg <- 31.3
lon_Hattiesburg <- -89.3
MDR.albopictus$latitude[MDR.albopictus$citation == "Ezeakacha_Dissertation_2015"] <- lat_Hattiesburg
MDR.albopictus$longitude[MDR.albopictus$citation == "Ezeakacha_Dissertation_2015"] <- lon_Hattiesburg


# Mosquitoes from Muturi 2011 were collected from Jacksonville, FL, USA
lat_jacksonville <- 30.4
lon_jacksonville <- -81.6
MDR.albopictus$latitude[MDR.albopictus$citation == "Muturi_et_al_2011_JMedEnto"] <- lat_jacksonville
MDR.albopictus$longitude[MDR.albopictus$citation == "Muturi_et_al_2011_JMedEnto"] <- lon_jacksonville


# Mosquitoes from Westbrook 2010 and thesis 2010 were collected from Palm Beach County, FL, USA
lat_palmbeach <- 26.7
lon_palmbeach <- -80.4
MDR.albopictus$latitude[MDR.albopictus$citation == "Westbrook_et_al_2010_VecBorn&ZooDis"] <- lat_palmbeach
MDR.albopictus$longitude[MDR.albopictus$citation == "Westbrook_et_al_2010_VecBorn&ZooDis"] <- lon_palmbeach
MDR.albopictus$latitude[MDR.albopictus$citation == "Westbrook_Thesis_2010"] <- lat_palmbeach
MDR.albopictus$longitude[MDR.albopictus$citation == "Westbrook_Thesis_2010"] <- lon_palmbeach


# Mosquitoes from Wiwatanaratanabutr 2006 were collected from Kanchanaburi  Province, western Thailand
lat_thai <- 14.8
lon_thai <- 99.0
MDR.albopictus$latitude[MDR.albopictus$citation == "Witwatana_Kittayaong_2006_Med&VetEnt"] <- lat_thai
MDR.albopictus$longitude[MDR.albopictus$citation == "Witwatana_Kittayaong_2006_Med&VetEnt"] <- lon_thai


# Mosquitoes from Yee 2017 (citation was wrong in the data) were collected from various locations
MDR.albopictus$citation[MDR.albopictus$citation == "Yee_et al_2016 JAE_in_review"] <- "Yee_etal_2017_Oikos"
MDR.albopictus$data_source[MDR.albopictus$citation == "Yee_etal_2017_Oikos"] <- "figure 3"



###### Ae. nigromaculis ######
MDR.nigromaculis <- MDR.Mordecai2019 %>% 
  filter(host_code  == "Anig") %>% 
  # convert from development time to development rate
  mutate(trait_name = "MDR") %>% 
  mutate(error_pos_si = error_pos_si/trait^2) %>% 
  mutate(error_neg_si = error_neg_si/trait^2) %>% 
  mutate(trait = 1/trait) %>% 
  # Add new columns to provide more info
  mutate(genus = "Aedes") %>% 
  relocate(genus, .after = "trait_2") %>% 
  mutate(species = "nigromaculis") %>% 
  relocate(species, .after = "genus") %>% 
  mutate(type = "non-Arctic") %>% 
  select(!host_code) %>% 
  mutate(latitude = NA) %>% 
  mutate(longitude = NA)


## Rename columns
colnames(MDR.nigromaculis) <- c("trait_name", "temp", "trait", "error_pos", 
                              "error_neg","trait2_name", "trait2", "genus", 
                              "species", "citation", "data_source", "notes", 
                              "type", "latitude", "longitude")

# Mosquitoes from Brust 1967 were collected from Winnipeg, Manitoba
MDR.nigromaculis$latitude[MDR.nigromaculis$citation == "Brust_1967_CanEntomol"] <- lat_winnipeg
MDR.nigromaculis$longitude[MDR.nigromaculis$citation == "Brust_1967_CanEntomol"] <- lon_winnipeg


###### Ae. sollicitans ######
MDR.sollicitans <- MDR.Mordecai2019 %>% 
  filter(host_code  == "Asol") %>% 
  # convert from development time to development rate
  mutate(trait_name = "MDR") %>% 
  mutate(error_pos_si = error_pos_si/trait^2) %>% 
  mutate(error_neg_si = error_neg_si/trait^2) %>% 
  mutate(trait = 1/trait) %>% 
  # Add new columns to provide more info
  mutate(genus = "Aedes") %>% 
  relocate(genus, .after = "trait_2") %>% 
  mutate(species = "sollicitans") %>% 
  relocate(species, .after = "genus") %>% 
  mutate(type = "non-Arctic") %>% 
  select(!host_code) %>% 
  mutate(latitude = NA) %>% 
  mutate(longitude = NA)


## Rename columns
colnames(MDR.sollicitans) <- c("trait_name", "temp", "trait", "error_pos", 
                                "error_neg","trait2_name", "trait2", "genus", 
                                "species", "citation", "data_source", "notes", 
                               "type", "latitude", "longitude")

# Shelton 1973 was conducted in Lake Charles, Louisiana
lat_lakecharles <- 30.2
lon_lakecharles <- -93.2
MDR.sollicitans$latitude[MDR.sollicitans$citation == "Shelton_1973_MosquitoNews"] <- lat_lakecharles
MDR.sollicitans$longitude[MDR.sollicitans$citation == "Shelton_1973_MosquitoNews"] <- lon_lakecharles


###### Ae. triseriatus ######
MDR.triseriatus <- MDR.Mordecai2019 %>% 
  filter(host_code  == "Atri") %>% 
  # Add new columns to provide more info
  mutate(genus = "Aedes") %>% 
  relocate(genus, .after = "trait_2") %>% 
  mutate(species = "triseriatus") %>% 
  relocate(species, .after = "genus") %>% 
  mutate(type = "non-Arctic") %>% 
  select(!host_code) %>% 
  mutate(latitude = NA) %>% 
  mutate(longitude = NA)


## Rename columns
colnames(MDR.triseriatus) <- c("trait_name", "temp", "trait", "error_pos", 
                               "error_neg","trait2_name", "trait2", "genus", 
                               "species", "citation", "data_source", "notes", 
                               "type", "latitude", "longitude")


# Convert from development time to development rate
MDR.triseriatus <- MDR.triseriatus %>% 
  mutate(trait = ifelse(trait_name == "1/MDR", 1/trait, trait)) %>% 
  mutate(trait_name = "MDR")

# Mosquitoes from Jalil 1972 were Alabama strain
lat_columbus <- 33.1
lon_columbus <- -86.6
MDR.triseriatus$latitude[MDR.triseriatus$citation == "Jalil 1972"] <- lat_columbus
MDR.triseriatus$longitude[MDR.triseriatus$citation == "Jalil 1972"] <- lon_columbus


MDR.triseriatus$latitude[MDR.triseriatus$citation == "Shelton_1973_MosquitoNews"] <- lat_lakecharles
MDR.triseriatus$longitude[MDR.triseriatus$citation == "Shelton_1973_MosquitoNews"] <- lon_lakecharles


# Combine data into a single dataframe
TraitData_MDR <- bind_rows(MDR.flavescens, MDR.nigripes, MDR.vexans, 
                           MDR.albopictus, MDR.nigromaculis,
                           MDR.sollicitans, MDR.triseriatus)


write_csv(TraitData_MDR, "data-processed/TraitData_MDR.csv")

TraitData_MDR %>% 
  filter(type == "Arctic") %>% 
  group_by(species, citation) %>% 
  mutate(unique_id = cur_group_id()) %>% 
  distinct(unique_id)

## Plot raw data
plot.data.MDR <- TraitData_MDR %>% 
  ggplot(aes(x = round(temp,0), y = trait, colour = species)) +
  geom_point() +
  labs(y = "Mosquito development rate (1/days)", x = "Temperature ºC") +
  scale_colour_discrete(name = "Species", labels = c("Ae. albopictus",
                                                     "Ae. flavescens",
                                                     "Ae. nigripes",
                                                     "Ae. nigromaculis",
                                                     "Ae. sollicitans",
                                                     "Ae. triseriatus",
                                                     "Ae. vexans"
  )) +
  facet_grid(rows = vars(type)) +
  theme_bw()

plot.data.MDR

ggsave("figures/raw_data/plot.data.MDR.png", plot.data.MDR, width = 9.83, height = 6.17)



# 7. Adult mosquito lifespan (lf) -----------------------------------------

## Arctic species ----------------------------------------------------------
###### Ae. vexans ######
lf.vexans <- read_csv("data-raw/lf_aedes_vexans.Costello1971.csv") %>% 
  clean_names() 

lf.vexans <- lf.vexans %>% 
  filter(sex == "F") %>% # Only select female mosquitoes
  filter(adult_food != "starved") %>% 
  # Add new columns to provide more info
  mutate(trait_name = "lf") %>% 
  mutate(error_pos = NA) %>% 
  mutate(error_neg = NA) %>% 
  mutate(error_unit = NA) %>% 
  mutate(trait_def = "days to 50% mortality") %>% 
  mutate(trait2_name = "adult food; RH") %>% 
  mutate(trait2 = paste0(adult_food, "; ", relative_humidity)) %>% 
  mutate(genus = "Aedes") %>% 
  mutate(species = "vexans") %>% 
  mutate(citation = "Costello_1971_JEconEntomol") %>% 
  mutate(doi = "10.1093/jee/64.1.324") %>% 
  mutate(data_source = "table 1") %>% 
  mutate(notes = NA) %>% 
  mutate(type = "Arctic") %>% 
  mutate(latitude = NA) %>% 
  mutate(longitude = NA) %>% 
  select(trait_name, treatment_temperature, days_to_50_percent_mortality, 
         error_pos, error_neg, error_unit, trait_def, trait2_name, trait2, 
         genus, species, citation, doi, data_source, notes, type, latitude, 
         longitude)
  
  
colnames(lf.vexans)[2] <- "temp"
colnames(lf.vexans)[3] <- "trait"

# Mosquitoes from Brust 1967 were collected near Winnipeg, Manitoba;
lat_winnipeg <- 49.9
lon_winnipeg <- -97.1
lf.vexans$latitude[lf.vexans$citation == "Costello_1971_JEconEntomol"] <- lat_winnipeg
lf.vexans$longitude[lf.vexans$citation == "Costello_1971_JEconEntomol"] <- lon_winnipeg


###### Ae. cinereus, Ae. communis, Ae. impiger, Ae. punctor in Alaska ######
lf.aedes <- read_csv("data-raw/a&EFGC&lf_aedes_spp.Sommerman1969.csv") %>%
  clean_names() %>%
  filter(trait_name == "lf") %>% 
  # Add new columns to provide more info
  mutate(trait_def = "average days alive") %>% 
  mutate(type = "Arctic")

# Convert trait2 to character
lf.aedes$trait2 <- as.character(lf.aedes$trait2)


## Non-Arctic species (for informing priors) -------------------------------

###### Ae. albopictus (from VecTrait) ######
lf.VecTrait <- read_csv("data-raw/lf_VecTrait.csv") %>% 
  clean_names()

unique(lf.VecTrait$interactor1)
# We will use data from Ae. albopictus

lf.VecTrait <- lf.VecTrait  %>% 
  filter(interactor1sex == "female") %>% # Only want female data
  filter(interactor1species == "albopictus") %>% 
  # Add new columns to provide more info
  mutate(trait_name = "lf") %>% 
  mutate(trait2_name = ifelse(
    is.na(second_stressor), NA, # if no second_stressor, leave trait2_name empty
    paste0(second_stressor, " (", second_stressor_unit, ")"))) %>% # Combine second_stressor and the unit
  mutate(type = "non-Arctic") %>% 
  select(trait_name, interactor1temp, original_trait_value, original_error_pos, 
         original_error_neg, original_error_unit, original_trait_def, 
         trait2_name, second_stressor_value, interactor1genus, 
         interactor1species, citation, doi, figure_table, notes, type, 
         latitude, longitude)


colnames(lf.VecTrait) <- c("trait_name", "temp", "trait", "error_pos", 
                          "error_neg", "error_unit", "trait_def", "trait2_name", 
                          "trait2", "genus", "species", "citation", "doi", 
                          "data_source", "notes", "type", "latitude", 
                          "longitude")

lf.VecTrait$trait2 <- as.character(lf.VecTrait$trait2) ## Change the trait2 column to character (so that it can combine with other dataset)

# The dataset also contains min/max duration of life stage. We only want the mean
lf.VecTrait <- lf.VecTrait %>% 
  filter(trait_def == "mean duration of life stage")



###### Ae. albopictus (from Mordecai et al. 2019) ######

lf.Mordecai2019 <- read_csv("data-raw/lf_Data_Mordecai2019.csv") %>% 
  clean_names() %>% 
  select(!series_id)
  
unique(lf.Mordecai2019$host_code) 
# Ae. taeniorhnchus and Ae. aegypti are tropical to subtropical spp.
# We will only use data from Ae. albopictus

lf.aalb <- lf.Mordecai2019 %>% 
  filter(host_code == "Aalb") %>% 
  # Add new columns to provide more info
  mutate(genus = "Aedes") %>% 
  mutate(species = "albopictus") %>% 
  mutate(type = "non-Arctic") %>% 
  mutate(latitude = NA) %>% 
  mutate(longitude = NA) %>%   
  select(trait_name, t, trait, error_pos_si, error_neg_si, trait_2, trait2_name,
         genus, species, citation, figure, notes, type, latitude, longitude)


## Rename columns
colnames(lf.aalb) <- c("trait_name", "temp", "trait", "error_pos", 
                       "error_neg","trait2_name", "trait2", "genus", # The values of trait2_name and trait2 are swapped for some reasons
                       "species", "citation", "data_source", "notes", 
                       "type", "latitude", "longitude")

## Change trait_name from 1/mu to lf (1/mu = lf), and from prop.dead to 1/lf
lf.aalb <- lf.aalb %>% 
  mutate(trait_name = ifelse(
    trait_name == "1/mu",
    "lf",
    ifelse(trait_name == "prop.dead", "1/lf", trait_name)))

## Remove mosquitoes that were starved
lf.aalb <- lf.aalb %>% 
  filter(is.na(trait2) | trait2 != "starved")

# Mosquitoes from Calado 2002 were collected from Registro, São Paulo State, Brazil
lat_registro <- -24.5
lon_registro <- -47.9
lf.aalb$latitude[lf.aalb$citation == "Calado_Navarro-Silva_2002_RevistaBrasileraDeEntomologia"] <- lat_registro
lf.aalb$longitude[lf.aalb$citation == "Calado_Navarro-Silva_2002_RevistaBrasileraDeEntomologia"] <- lon_registro

# Mosquitoes from Alto 2001 were collected from East St. Louis, IL, USA
lat_stlouis <- 38.6
lon_stlouis <- -90.1
lf.aalb$latitude[lf.aalb$citation == "Alto_et_al_2001b"] <- lat_stlouis
lf.aalb$longitude[lf.aalb$citation == "Alto_et_al_2001b"] <- lon_stlouis

# Mosquitoes from Ezeakacha 2015 were collected in and around Hattiesburg, MS, USA
lat_Hattiesburg <- 31.3
lon_Hattiesburg <- -89.3
lf.aalb$latitude[lf.aalb$citation == "Ezeakacha_Dissertation_2015"] <- lat_Hattiesburg
lf.aalb$longitude[lf.aalb$citation == "Ezeakacha_Dissertation_2015"] <- lon_Hattiesburg

## Since both VecTrait and Mordecai et al. 2019 database contains data from 
## Calado and Navarro-Silva 2002, we will remove the same data from VecTrait 
## database.
lf.VecTrait <- lf.VecTrait %>% 
  filter(citation != "Calado and Navarro-Silva. 2002. Influencia da temperatura sobre a longevidade fecundidade eatividade hematofagica de Aedes (Stegomyia) albopictus Skuse 1894 (Diptera: Culicidae) sob condicoes de laboratorio. Revista Brasileira de Entomologia 46: 93-98.")


# Combine data into a single dataframe
TraitData_lf <- bind_rows(lf.vexans, lf.aedes, 
                          lf.VecTrait, lf.aalb)


## Convert mortality rate (1/lf) to lifespan (lf)
TraitData_lf <- TraitData_lf %>% 
  mutate(trait = ifelse(trait_name == "1/lf", 1/trait, trait)) %>% 
  mutate(notes = ifelse(trait_name == "1/lf",
                        ifelse(is.na(notes),
                               "converted from mortality rate", 
                               paste0(notes, "; converted from mortality rate")
                               ),
                        notes)) %>%
  mutate(trait_name = "lf") 


write_csv(TraitData_lf, "data-processed/TraitData_lf.csv")


## Plot raw data
plot.data.lf <- TraitData_lf %>% 
  ggplot(aes(x = temp, y = trait, colour = species)) +
  geom_point(aes(colour = species)) +
  labs(y = "Mosquito adult lifespan (days)", x = "Temperature ºC") +
  scale_colour_discrete(name = "Species", labels = c("Ae. albopictus",
                                                     "Ae. cinereus", 
                                                     "Ae. communis",
                                                     "Ae. impiger",
                                                     "Ae. punctor",
                                                     "Ae. vexans"
  )) +
  facet_grid(rows = vars(type)) +
  theme_bw()

plot.data.lf

ggsave("figures/raw_data/plot.data.lf.png", plot.data.lf, width = 9.83, height = 6.17)




# 8. Eggs per female per gonotrophic cycle (EFGC) ------------------------------
## Arctic species --------------------------------------------------------------

##### Ae. hexodontus #####
EFGC.hexodontus <- read_csv("data-raw/EFGC_aedes_hexodontus.Barlow1955.csv") %>% 
  clean_names()

EFGC.hexodontus <- EFGC.hexodontus %>% 
  mutate(type = "Arctic")


##### Ae. cinereus, Ae. communis, Ae. impiger, Ae. punctor in Alaska #####
EFGC.aedes <- read_csv("data-raw/a&EFGC&lf_aedes_spp.Sommerman1969.csv") %>%
  clean_names() %>%
  filter(trait_name == "EFGC") 

EFGC.aedes <- EFGC.aedes %>% 
  ## Only include replicates with at least 8 individuals
  # Separate the notes column into two columns: n and the number of replicates
  separate(notes, sep = " = ", into = c("n", "rep")) %>% 
  mutate(rep = as.numeric(rep)) %>% 
  # Only include replicates with at >=8 individuals
  filter(rep >= 8) %>% 
  # Combine n and rep back together
  unite(n, rep, col = "notes", sep = "=") %>% 
  # Add new columns to provide more info
  mutate(type = "Arctic") %>% 
  relocate(type, .before = "latitude") 


## Non-Arctic species (for informing priors) -----------------------------------
EFGC.Mordecai2019 <- read_csv("data-raw/Fecundity_Data_Mordecai2019.csv") %>% 
  clean_names() %>% 
  # TFD was defined as eggs laid per female per gonotrophic cycle (number/female) in this dataset (see Mordecai 2017 Table A in S2 Text)
  filter(trait_name == "TFD") 

unique(EFGC.Mordecai2019$host_code) 
## Ae. albopictus is adapted to cooler, temperate climate.


###### Ae. albopictus ######
EFGC.albopictus <- EFGC.Mordecai2019 %>% 
  filter(host_code == "Aalb") %>% 
  # select columns that we need
  dplyr::select(trait_name, t, trait, error_pos_si, error_neg_si, trait2_name, 
                trait_2, citation, figure, notes) %>% 
  # Add new columns to provide more info
  mutate(trait_name = "EFGC") %>%  # Change from TFD to EFGC (eggs per female per gonotrophic cycle)
  mutate(error_unit = NA) %>% 
  relocate(error_unit, .after = "error_neg_si") %>% # rearrange columns
  mutate(trait_def = NA) %>% 
  relocate(trait_def, .after = "error_unit") %>% 
  mutate(genus = "Aedes") %>% 
  relocate(genus, .after = "trait_2") %>% 
  mutate(species = "albopictus") %>% 
  relocate(species, .after = "genus") %>% 
  mutate(doi = NA) %>% 
  relocate(doi, .after = "citation") %>% 
  mutate(type = "non-Arctic") %>% 
  mutate(latitude = NA) %>% 
  mutate(longitude = NA)


## Rename columns
colnames(EFGC.albopictus) <- c("trait_name", "temp", "trait", "error_pos", "error_neg",
                            "error_unit", "trait_def", "trait2_name", "trait2",
                            "genus", "species", "citation","doi", "data_source", 
                            "notes", "type", "latitude", "longitude")


## Provide more info on the papers
## Delatte et al. 2009
EFGC.albopictus$doi[EFGC.albopictus$citation == "Delatte_etal_2009_JMedEnto"] <-
  "10.1603/033.046.0105"

EFGC.albopictus$data_source[EFGC.albopictus$citation == "Delatte_etal_2009_JMedEnto"] <-
  "table 6"

# Mosquitoes from Delatte 2009 were collected from Saint Pierre, La Reunion Island, France
lat_stpierre <- -21.3
lon_stpierre <- 55.5
EFGC.albopictus$latitude[EFGC.albopictus$citation == "Delatte_etal_2009_JMedEnto"] <- lat_stpierre
EFGC.albopictus$longitude[EFGC.albopictus$citation == "Delatte_etal_2009_JMedEnto"] <- lon_stpierre


# Mosquitoes from Ezeakacha 2015 were collected in and around Hattiesburg, MS, USA
lat_Hattiesburg <- 31.3
lon_Hattiesburg <- -89.3
EFGC.albopictus$latitude[EFGC.albopictus$citation == "Ezeakacha_Dissertation_2015"] <- lat_Hattiesburg
EFGC.albopictus$longitude[EFGC.albopictus$citation == "Ezeakacha_Dissertation_2015"] <- lon_Hattiesburg

# Mosquitoes from Yee 2017 (citation was wrong in the data) were collected from various locations
EFGC.albopictus$citation[EFGC.albopictus$citation == "Yee_et al_2016 JAE_in_review"] <- "Yee_etal_2017_Oikos"


TraitData_EFGC <- bind_rows(EFGC.aedes, EFGC.hexodontus, EFGC.albopictus)

write_csv(TraitData_EFGC, "data-processed/TraitData_EFGC.csv")


## Plot raw data
plot.data.EFGC <- TraitData_EFGC %>% 
  ggplot() +
  geom_point(aes(x = temp, y = trait, colour = species)) +
  # xlim(c(10,35)) +
  labs(y = "Eggs per female per gonotrophic cycles", x = "Temperature ºC") +
  scale_colour_discrete(name = "Species", labels = c("Ae. albopictus",
                                                     "Ae. cinereus",
                                                     "Ae. communis",
                                                     "Ae. hexodontus",
                                                     "Ae. impiger",
                                                     "Ae. punctor"
  )) +
  facet_grid(rows = vars(type)) +
  theme_bw()


plot.data.EFGC

ggsave("figures/raw_data/plot.data.EFGC.png", plot.data.EFGC, width = 9.83, height = 6.17)
