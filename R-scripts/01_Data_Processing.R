## Lilian Chan, University of Guelph
## Arctic vector-borne disease transmission suitability model
##
## Purpose: Process the raw data and save each trait as it own csv file
## 
## Table of content:
##    0. Set-up workspace
##    1. Biting rate (a)
##    2. Vector Competence (bc)
##         i) Infection efficiency (c)
##        ii) Transmission efficiency (b)
##    3. Parasite development rate (PDR)
##    4. Egg viavility (EV)
##    5. Larval survival (pLA)
##    6. Mosquito egg-to-adult development rate (MDR)
##    7. Adult mosquito lifespan (lf)
##    8. Eggs per female per gonotrophic cycle (EFGC)



# 0. Set-up workspace -----------------------------------------------------
library(tidyverse)
library(readxl)
library(janitor)
library(ggpubr)



# 1. Biting rate (a) ------------------------------------------------------

## Arctic species ----------------------------------------------------------
###### Ae. cinereus, Ae. communis, Ae. impiger, Ae. punctor in Alaska ######
a.aedes <- read_csv("data-raw/aedes_spp.Sommerman1969.csv") %>% 
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
  unite(n, rep, col = "notes",) %>% 
  mutate(type = "Arctic")


## Non-Arctic species (for informing priors) -------------------------------
# From VecTrait database
a.VecTrait <- read_csv("data-raw/a_VecTrait.csv") %>% 
  clean_names()

###### Aedes albopictus ######
a.aalb <- a.VecTrait  %>% 
  filter(interactor1species == "albopictus") %>%
  filter(original_trait_def == "mean duration of gonotrophic cycle") %>% 
  # Add new columns to provide more info
  mutate(trait_name = "1/a") %>% 
  mutate(type = "non-Arctic") %>% 
  select(trait_name, interactor1temp, original_trait_value, original_error_pos, 
         original_error_neg, original_error_unit, original_trait_def, 
         interactor1genus, interactor1species, citation, doi, figure_table, 
         notes, type)


colnames(a.aalb) <- c("trait_name", "temp", "trait", "error_pos", 
                      "error_neg", "error_unit", "trait_def", "genus", 
                      "species", "citation", "doi", "data_source", "notes",
                      "type")
  

###### Aedes aegypti ######
a.aaeg <- a.VecTrait  %>% 
  filter(interactor1species == "aegypti") %>% 
  # Since this dataset provides individual-level data, we will calculate the mean at each temp
  group_by(interactor1temp, original_trait_def, interactor1genus, interactor1species, figure_table, location, doi) %>% 
  summarize(trait = mean(original_trait_value),
            error_pos = sd(original_trait_value),) %>% 
  # Add new columns to provide more info
  mutate(trait_name = "1/a") %>% 
  mutate(error_unit = "sd") %>% 
  mutate(citation = "Goindin_2015_PLoSOne") %>% 
  mutate(type = "non-Arctic") %>% 
  select(trait_name, interactor1temp, trait, error_pos, error_unit,
         original_trait_def, interactor1genus, interactor1species, citation, 
         doi, figure_table, location, type)


colnames(a.aaeg) <- c("trait_name", "temp", "trait", "error_pos",
                      "error_unit", "trait_def", "genus", "species", "citation",
                      "doi", "data_source", "notes", "type")

## Change trait_def from individual-level to mean
a.aaeg$trait_def <- ifelse(a.aaeg$trait_def == "individual-level duration of first gonotrophic cycle",
                           "mean duration of first gonotrophic cycle",
                           "mean duration of second gonotrophic cycle")


# From data compiled in Mordecaiet al. 2019
a.Mordecai2019 <- read_csv("data-raw/a_Data_Mordecai2019.csv") %>% 
  clean_names()

a.aaeg.Mordecai2019 <- a.Mordecai2019 %>% 
  filter(host_code == "Aaeg") %>% 
  # select columns that we need
  select(trait_name, t_c, trait, error_pos_si, error_neg_si, trait2_name, 
         trait_2, citation, figure, notes) %>% 
  # Add new columns to provide more info
  mutate(genus = "Aedes") %>% 
  relocate(genus, .after = "trait_2") %>% 
  mutate(species = "aegypti") %>% 
  relocate(species, .after = "genus") %>% 
  mutate(type = "non-Arctic")

  
## Change GCD (gonotrophic cycle duration) to 1/a
a.aaeg.Mordecai2019$trait_name[a.aaeg.Mordecai2019$trait_name == "GCD"] <- "1/a"


## Rename columns
colnames(a.aaeg.Mordecai2019) <- c("trait_name", "temp", "trait", "error_pos", 
                                   "error_neg", "trait2_name", "trait2", 
                                   "genus", "species", "citation","doi", 
                                   "data_source", "notes", "type")

## Provide more info on the papers
### Focks and Barrera 2006
a.aaeg.Mordecai2019$citation[a.aaeg.Mordecai2019$citation == "Focks_Barrera_2006_Research&TrainingTropicalDis_Geneva_Paper"] <- "Focks_2006_Research&TrainingTropicalDis_Geneva_Paper"
a.aaeg.Mordecai2019$data_source[a.aaeg.Mordecai2019$citation == "Focks_2006_Research&TrainingTropicalDis_Geneva_Paper"] <- "table 5"


### Focks et al. 1993
a.aaeg.Mordecai2019$citation[a.aaeg.Mordecai2019$citation == "Focks_et_al_1993a_JMedEntom"] <- "Focks_1993_JMedEntomol"
a.aaeg.Mordecai2019$data_source[a.aaeg.Mordecai2019$citation == "Focks_1993_JMedEntomol"] <- "figure 9"
a.aaeg.Mordecai2019$doi[a.aaeg.Mordecai2019$citation == "Focks_1993_JMedEntomol"] <- "10.1093/jmedent/30.6.1003"



## Combine all data
TraitData_a <-  bind_rows(a.aedes,
                          a.aalb)


write_csv(TraitData_a, "data-processed/TraitData_a.csv")


## Plot raw data
plot.data.a <- TraitData_a %>%
  mutate(trait = ifelse(trait_name == "1/a", 1/trait, trait)) %>% 
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
  # facet_grid(rows = vars(type)) +
  theme_bw()

plot.data.a

ggsave("figures/raw_data/plot.data.a.png", plot.data.a, , width = 9.83, height = 6.17)




# 2. Vector Competence (bc) -----------------------------------------------
## Vector competence has two components: infection efficiency (c) and 
## transmission efficiency (b)


## 2i. Infection efficiency (c) -----------------------------------------------

###### Ae. Trivittatus (transmitting Dirofilaria immitis) ######
c.trivittatus <- read_csv("data-raw/dirofilaria_immitis_aedes_trivittatus.Christensen1978.csv") %>% 
  clean_names() 

TraitData_c <- c.trivittatus %>% 
  # Add new columns to provide more info
  mutate(trait_name = "c") %>% 
  # Convert development time to development rate
  mutate(trait = infection_rate_percent/100) %>% 
  mutate(error_pos = NA) %>% # this paper did not provide errors
  mutate(error_neg = NA) %>% 
  mutate(error_unit = NA) %>% 
  mutate(trait_def = "infection rate %") %>% 
  mutate(trait2_name = NA) %>% 
  mutate(trait2 = NA) %>% 
  mutate(genus = "Aedes") %>% 
  mutate(species = "trivittatus") %>% 
  mutate(paras_genus = "dirofilaria") %>% 
  mutate(paras_species = "immitis") %>% 
  mutate(citation = "Christensen_1978_ProcHelmintholSocWash") %>% 
  mutate(doi = NA) %>% #No doi for this paper?
  mutate(data_source = "figure 2") %>% 
  mutate(type = "non-Arctic") %>% 
  select(trait_name, temp, trait, error_pos, error_neg, error_unit, trait_def, 
         trait2_name, trait2, genus, species, paras_genus, paras_species, 
         citation, doi, data_source, notes, type)

write_csv(TraitData_c, "data-processed/TraitData_c.csv")


## Plot raw data
plot.data.c <- TraitData_c %>% 
  ggplot() +
  geom_point(aes(x = temp, y = trait, colour = species)) +
  labs(y = "Infection probability", x = "Temperature ºC") +
  scale_color_discrete(name = "Species", label = "Ae. trivittatus") +
  facet_grid(rows = vars(type)) +
  theme_bw()

plot.data.c

ggsave("figures/raw_data/plot.data.c.png", plot.data.c, , width = 9.83, height = 6.17)


## 2ii. Transmission efficiency (b) -----------------------------------------------
# We couldn't find any data on this trait, so we excluded it from the model



# 3. Parasite development rate (PDR) --------------------------------------

## Arctic species ----------------------------------------------------------
###### Varestrongylus eleguneniensis ######
PDR.eleguneniensis <- read_csv("data-raw/varestrongylus_eleguneniensis.Kafle2018.csv") %>% 
  clean_names()


PDR.eleguneniensis <- PDR.eleguneniensis %>% 
  # Add new columns to provide more info
  mutate(trait_name = "1/PDR") %>% 
  mutate(trait = dayi_l3) %>% 
  mutate(error_pos = NA) %>% # this paper did not provide errors
  mutate(error_neg = NA) %>% 
  mutate(error_unit = NA) %>% 
  mutate(trait_def = "days first L3 observed") %>% 
  mutate(trait2_name = NA) %>% 
  mutate(trait2 = NA) %>% 
  mutate(genus = "varestrongylus") %>% 
  mutate(species = "eleguneniensis") %>% 
  mutate(host.genus = "deroceras") %>% 
  mutate(host.species = "laeve") %>% 
  mutate(citation = "Kafle_2018_ParasitVectors") %>% 
  mutate(doi = "10.1186/s13071-018-2946-x") %>% 
  mutate(data_source = "table 1") %>% 
  mutate(notes = "Nematode infecting caribou and muskoxen in the Canadian Arctic; transmitted by gastropods") %>% 
  mutate(type = "Arctic") %>% 
  select(trait_name, temp, trait, error_pos, error_neg, error_unit, trait_def,
         trait2_name, trait2, genus, species, host.genus, host.species, 
         citation, doi, data_source, notes, type)



## At 8.5ºC, development did not occur after 101 days. Change development time to 1000 days
PDR.eleguneniensis[5,"trait"] <- "1000"
PDR.eleguneniensis$trait <- as.numeric(PDR.eleguneniensis$trait)
PDR.eleguneniensis[5,"notes"] <- "did not develop after 101 days"


###### Setaria tundra ######
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
  mutate(type = "Arctic") %>% 
  select(trait_name, temp, trait, genus, species, host.genus, host.species, 
         trait_def, citation, doi, data_source, notes, type)


## At 14.1ºC, development was not completed. Change development time to 1000 days
PDR.tundra[1, "trait"] <- 1000



## Non-Arctic species (for informing priors) -------------------------------

###### Dirofilaria immitis (in Ae. Trivittatus) ######
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
  mutate(host.genus = "Aedes") %>% 
  mutate(host.species = "trivittatus") %>% 
  mutate(citation = "Christensen_1978_ProcHelmintholSocWash") %>% 
  mutate(doi = NA) %>% 
  mutate(data_source = "table 1") %>% 
  mutate(type = "non-Arctic") %>% 
  select(trait_name, temp, trait, genus, species, host.genus, host.species, 
         trait_def, citation, doi, data_source, notes, type)

## At 14.5ºC, development was not completed. Change development time to 1000 days
PDR.immitis[1,"trait"] <- 1000


###### Dirofilaria immitis (in Ae. aegypti) ###### 
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
  mutate(host.genus = "Aedes") %>% 
  mutate(host.species = "aegypti") %>% 
  mutate(citation = "Ledesma_2015_VetParasitol") %>% 
  mutate(doi = "10.1016/j.vetpar.2015.02.003") %>% 
  mutate(data_source = "table 1") %>% 
  mutate(notes = NA) %>% 
  mutate(type = "non-Arctic") %>% 
  select(trait_name, temp, trait, genus, species, host.genus, host.species, 
         trait_def, citation, doi, data_source, notes, type)


######  lymphatic filarisis worms (in Ae. polynesiensis) ###### 
PDR.polynesiensis <- read_csv("data-raw/wuchereria_bancrofti_aedes_polynesiensis.Lardeux1997.csv") %>% 
  clean_names()

PDR.polynesiensis <- PDR.polynesiensis %>% 
  mutate(trait = trait) %>% 
  # Add new columns to provide more info
  mutate(trait_name = "1/PDR") %>% 
  mutate(trait_def = "days to appearance of L3 after experimental infection") %>% 
  mutate(error_pos = sd) %>% 
  mutate(error_unit = "sd") %>% 
  mutate(genus = "wuchereria") %>% 
  mutate(species = "bancrofti") %>% 
  mutate(host.genus = "Aedes") %>% 
  mutate(host.species = "polynesiensis") %>% 
  mutate(citation = "Ladeaux_1997_Parasitology") %>% 
  mutate(doi = "10.1017/s0031182096008359") %>% 
  mutate(data_source = "table 1") %>% 
  mutate(notes = NA) %>% 
  mutate(type = "non-Arctic") %>% 
  select(trait_name, temp, trait, error_pos, error_unit, genus, species, host.genus, 
         host.species, trait_def, citation, doi, data_source, notes, type)

TraitData_PDR <- bind_rows(PDR.eleguneniensis, PDR.tundra, PDR.immitis, 
                           PDR.immitis.aegypti, PDR.polynesiensis)


write_csv(TraitData_PDR, "data-processed/TraitData_PDR.csv")


# Plot the raw data
plot.data.PDR <- TraitData_PDR %>% 
  ggplot(aes(x = temp, y = 1/trait)) +
  geom_point(aes(colour = species)) +
  labs(y = "Parasite development rate (1/days)", x = "Temperature ºC") +
  scale_colour_discrete(name = "Species", labels = c("W. bancrofti",
                                                     "V. eleguneniensis",
                                                     "D. immitis",
                                                     "S. tundra"
  )) +
  facet_grid(rows = vars(type)) +
  theme_bw()

plot.data.PDR

ggsave("figures/raw_data/plot.data.PDR.png", plot.data.PDR, , width = 9.83, height = 6.17)



# 4. Egg viability (EV) ---------------------------------------------------


## Arctic species ----------------------------------------------------------

###### Aedes vexans ######
EV.vexans <- read_csv("data-raw/TraitData_EV_Shocket2020.csv") %>% 
  clean_names()

EV.vexans <- EV.vexans %>% 
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
  select(trait_name, t, trait, error_pos_si, error_neg_si, error_unit, trait_def,
         trait2_name, trait_2, genus, species, citation, doi, figure, notes, type) 

colnames(EV.vexans) <- c("trait_name", "temp", "trait", "error_pos", 
                         "error_neg", "error_unit", "trait_def", "trait2_name", 
                         "trait2", "genus", "species", "citation", "doi", 
                         "data_source", "notes", "type")


## Non-Arctic species (for informing priors) -------------------------------
###### Aedes dorsalis ######
EV.dorsalis <- read_csv("data-raw/TraitData_EV_Shocket2020.csv") %>% 
  clean_names()

EV.dorsalis <- EV.dorsalis %>% 
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
  select(trait_name, t, trait, error_pos_si, error_neg_si, error_unit, trait_def,
         trait2_name, trait_2, genus, species, citation, doi, figure, notes, type) 

colnames(EV.dorsalis) <- c("trait_name", "temp", "trait", "error_pos", 
                         "error_neg", "error_unit", "trait_def", "trait2_name", 
                         "trait2", "genus", "species", "citation", "doi", 
                         "data_source", "notes", "type")





###### Aedes nigromaculis ######
EV.nigromaculis <- read_csv("data-raw/TraitData_EV_Shocket2020.csv") %>% 
  clean_names()
  
EV.nigromaculis <- EV.nigromaculis %>% 
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
  select(trait_name, t, trait, error_pos_si, error_neg_si, error_unit, trait_def,
         trait2_name, trait_2, genus, species, citation, doi, figure, notes, type) 

colnames(EV.nigromaculis) <- c("trait_name", "temp", "trait", "error_pos", 
                           "error_neg", "error_unit", "trait_def", "trait2_name", 
                           "trait2", "genus", "species", "citation", "doi", 
                           "data_source", "notes", "type")




###### Aedes triseriatus ######
EV.triseriatus <- read_csv("data-raw/aedes_triseriatus.Zimmerman2025.csv") %>% 
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
  select(trait_name, temp, trait, error_pos, error_neg, error_unit, trait_def,
         genus, species, citation, doi, data_source, notes, type) 


###### Aedes albopictus (from VecTrait database) ######
EV.VecTrait <- read_csv("data-raw/EV_VecTrait.csv") %>% 
  clean_names()


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
         interactor1species, citation, doi, figure_table, notes, type)


colnames(EV.VecTrait) <- c("trait_name", "temp", "trait", "error_pos", 
                           "error_neg", "error_unit", "trait_def", "trait2_name", 
                           "trait2", "genus", "species", "citation", "doi", 
                           "data_source", "notes", "type")

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


TraitData_EV <- bind_rows(EV.vexans, EV.dorsalis, EV.nigromaculis, EV.triseriatus, EV.VecTrait)


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
###### Ae. vexans ######
pLA.vexans <- read_csv("data-raw/pLA_Data_Mordecai2019.csv") %>% 
  clean_names()

pLA.vexans <- pLA.vexans %>% 
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
  mutate(type = "Arctic")


## Rename columns
colnames(pLA.vexans) <- c("trait_name", "temp", "trait", "error_pos", "error_neg",
                         "error_unit", "trait_def", "trait2_name", "trait2",
                         "genus", "species", "citation","doi", "data_source", 
                         "notes", "type")

## Provide more info on the papers
### Brust 1967
pLA.vexans$doi[pLA.vexans$citation == "Brust_1967_TheCanadianEntomologist"] <- "10.4039/ent99986-9"


## Non-Arctic species (for informing priors) -------------------------------
pLA.Mordecai2019 <- read_csv("data-raw/pLA_Data_Mordecai2019.csv") %>% 
  clean_names()


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
  mutate(type = "non-Arctic")


## Rename columns
colnames(pLA.nigromaculis) <- c("trait_name", "temp", "trait", "error_pos", 
                                "error_neg","trait2_name", "trait2", "genus", 
                                "species", "citation", "data_source", "notes",
                                "type")


## Provide more info on the papers
### Brust 1967
pLA.nigromaculis$doi[pLA.nigromaculis$citation == "Brust_1967_TheCanadianEntomologist"] <- "10.4039/ent99986-9"



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
  mutate(type = "non-Arctic")


## Rename columns
colnames(pLA.sollicitans) <- c("trait_name", "temp", "trait", "error_pos", 
                               "error_neg", "trait2_name", "trait2", "genus", 
                               "species", "citation", "data_source", "notes", 
                               "type")


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
  mutate(type = "non-Arctic")


## Rename columns
colnames(pLA.triseriatus) <- c("trait_name", "temp", "trait", "error_pos", 
                               "error_neg", "trait2_name", "trait2", "genus", 
                               "species", "citation", "data_source", "notes", 
                               "type")


TraitData_pLA <- bind_rows(pLA.vexans, pLA.nigromaculis, pLA.sollicitans, pLA.triseriatus)


write_csv(TraitData_pLA, "data-processed/TraitData_pLA.csv")


# Plot the raw data
plot.data.pLA <- TraitData_pLA %>% 
  ggplot(aes(x = temp, y = trait)) +
  geom_point(aes(colour = species, shape = citation)) +
  labs(y = "Larval survival (%)", x = expression(paste("Temperature (", degree, "C)"))) +
  scale_colour_discrete(name = "Species", labels = c("Ae. nigromaculis",
                                                     "Ae. sollicitans",
                                                     "Ae. triseriatus",
                                                     "Ae. vexans"
  )) +
  scale_shape_discrete(name = "Citation", labels = c("Brust 1967",
                                                     "Shelton 1973",
                                                     "Teng 2000",
                                                     "Trpis 1970")) +
  facet_grid(rows = vars(type)) +
  theme_bw()

plot.data.pLA



# 6. Mosquito egg-to-adult development rate (MDR) -------------------------

## Arctic species ----------------------------------------------------------
###### Aedes nigripes ######
MDR.nigripes <- read_excel("data-raw/aedes_nigripes.Culler2015.xlsx", 
                           sheet = "Dev. Time") %>% 
  clean_names() 


MDR.nigripes <- MDR.nigripes %>% 
  mutate(trait = development_time_in_days) %>% 
  filter(sex == "female") %>% # Only female data is used
  # Add new columns to provide more info
  mutate(trait_name = "1/MDR") %>% 
  mutate(error_pos = NA) %>% 
  mutate(error_neg = NA) %>% 
  mutate(error_unit = NA) %>% 
  mutate(trait_def = "days until emergence") %>% 
  mutate(trait2_name = NA) %>% 
  mutate(trait2 = NA) %>% 
  mutate(genus = "Aedes") %>% 
  mutate(species = "nigripes") %>% 
  mutate(citation = "Culler_2015_ProcRSocB") %>% 
  mutate(doi = "10.1098/rspb.2015.1549") %>% 
  mutate(data_source = "raw data") %>% 
  mutate(notes = NA) %>% 
  mutate(type = "Arctic") %>% 
  select(trait_name, mean_temperature_during_development_c, trait, error_pos, 
         error_neg, error_unit, trait_def, trait2_name, trait2,
         genus, species, citation, doi, data_source, notes, type)


colnames(MDR.nigripes)[2] <- "temp"


###### Aedes vexans ######
## Read Mordecai et al 2019 data
MDR.Mordecai2019 <- read_csv("data-raw/MDR_Data_Mordecai2019.csv") %>% 
  clean_names() %>% 
  # select columns that we need
  select(trait_name, t_c, trait, error_pos_si, error_neg_si, trait2_name, 
         trait_2, host_code, citation, figure, notes)
  

MDR.vexans <- MDR.Mordecai2019 %>% 
  filter(host_code  == "Avex") %>% 
  # Add new columns to provide more info
  mutate(genus = "Aedes") %>% 
  relocate(genus, .after = "trait_2") %>% 
  mutate(species = "vexans") %>% 
  relocate(species, .after = "genus") %>% 
  mutate(type = "Arctic") %>% 
  select(!host_code)


## Rename columns
colnames(MDR.vexans) <- c("trait_name", "temp", "trait", "error_pos", 
                          "error_neg","trait2_name", "trait2", "genus", 
                          "species", "citation", "data_source", "notes", "type")



## Non-Arctic species (for informing priors) -------------------------------

###### Aedes albopictus ######
MDR.albopictus <- MDR.Mordecai2019 %>% 
  filter(host_code  == "Aalb") %>% 
  # Add new columns to provide more info
  mutate(genus = "Aedes") %>% 
  relocate(genus, .after = "trait_2") %>% 
  mutate(species = "albopictus") %>% 
  relocate(species, .after = "genus") %>% 
  mutate(type = "non-Arctic") %>% 
  select(!host_code)


## Rename columns
colnames(MDR.albopictus) <- c("trait_name", "temp", "trait", "error_pos", 
                              "error_neg","trait2_name", "trait2", "genus", 
                              "species", "citation", "data_source", "notes", "type")


###### Aedes nigromaculis ######
MDR.nigromaculis <- MDR.Mordecai2019 %>% 
  filter(host_code  == "Anig") %>% 
  # Add new columns to provide more info
  mutate(genus = "Aedes") %>% 
  relocate(genus, .after = "trait_2") %>% 
  mutate(species = "nigromaculis") %>% 
  relocate(species, .after = "genus") %>% 
  mutate(type = "non-Arctic") %>% 
  select(!host_code)


## Rename columns
colnames(MDR.nigromaculis) <- c("trait_name", "temp", "trait", "error_pos", 
                              "error_neg","trait2_name", "trait2", "genus", 
                              "species", "citation", "data_source", "notes", "type")


###### Aedes sollicitans ######
MDR.sollicitans <- MDR.Mordecai2019 %>% 
  filter(host_code  == "Asol") %>% 
  # Add new columns to provide more info
  mutate(genus = "Aedes") %>% 
  relocate(genus, .after = "trait_2") %>% 
  mutate(species = "sollicitans") %>% 
  relocate(species, .after = "genus") %>% 
  mutate(type = "non-Arctic") %>% 
  select(!host_code)


## Rename columns
colnames(MDR.sollicitans) <- c("trait_name", "temp", "trait", "error_pos", 
                                "error_neg","trait2_name", "trait2", "genus", 
                                "species", "citation", "data_source", "notes", "type")


###### Aedes triseriatus ######
MDR.triseriatus <- MDR.Mordecai2019 %>% 
  filter(host_code  == "Asol") %>% 
  # Add new columns to provide more info
  mutate(genus = "Aedes") %>% 
  relocate(genus, .after = "trait_2") %>% 
  mutate(species = "triseriatus") %>% 
  relocate(species, .after = "genus") %>% 
  mutate(type = "non-Arctic") %>% 
  select(!host_code)


## Rename columns
colnames(MDR.triseriatus) <- c("trait_name", "temp", "trait", "error_pos", 
                               "error_neg","trait2_name", "trait2", "genus", 
                               "species", "citation", "data_source", "notes", "type")



###### Aedes sierrensis ######
MDR.sierrensis <- read_csv("data-raw/aedes_sierrensis.Couper2024.csv") %>% 
  clean_names() 

MDR.sierrensis <- MDR.sierrensis %>% 
  filter(!is.na(juvenile_dev_rate)) %>% 
  # Since this dataset provides individual-level data, we will calculate the mean at each temp
  group_by(temp_treatment) %>%
  summarize(trait = mean(juvenile_dev_rate),
            error_pos = sd(juvenile_dev_rate)) %>%
  mutate(trait_name = "MDR") %>% 
  mutate(error_unit = "sd") %>% 
  mutate(trait_def = "1/days") %>% 
  mutate(trait2_name = NA) %>% 
  mutate(trait2 = NA) %>% 
  mutate(genus = "Aedes") %>% 
  mutate(species = "sierrensis") %>% 
  mutate(citation = "Couper_2024_ProcBiolSci") %>% 
  mutate(doi = "10.1098/rspb.2023.2457") %>% 
  mutate(data_source = "calculated from raw data") %>% # raw data from the paper
  mutate(type = "non-Arctic") %>% 
  select(trait_name, temp_treatment, trait, 
                error_pos, error_unit, trait_def, trait2_name, trait2,
                genus, species, citation, doi, data_source, type)

# Rename columns
colnames(MDR.sierrensis)[2] <- "temp"


# Combine data from nigripes and sierrensis into a single dataframe
TraitData_MDR <- bind_rows(MDR.nigripes, MDR.vexans, MDR.albopictus,
                           MDR.nigromaculis, MDR.sierrensis, MDR.sollicitans, MDR.triseriatus)

write_csv(TraitData_MDR, "data-processed/TraitData_MDR.csv")


## Plot raw data
plot.data.MDR <- TraitData_MDR %>% 
  mutate(trait = ifelse(trait_name == "1/MDR", 1/trait, trait)) %>% 
  ggplot(aes(x = round(temp,0), y = trait, colour = species)) +
  geom_point() +
  labs(y = "Mosquito development rate (days)", x = "Temperature ºC") +
  scale_colour_discrete(name = "Species", labels = c("Ae. albopictus",
                                                     "Ae. nigripes",
                                                     "Ae. nigromaculis",
                                                     "Ae. sierrensis",
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
###### Aedes vexans ######
lf.vexans <- read_csv("data-raw/aedes_vexans.Costello1971.csv") %>% 
  clean_names() 

lf.vexans <- lf.vexans %>% 
  # Only select female mosquitoes
  filter(sex == "F") %>% 
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
  select(trait_name, treatment_temperature, days_to_50_percent_mortality, 
         error_pos, error_neg, error_unit, trait_def, trait2_name, trait2, 
         genus, species, citation, doi, data_source, notes, type)
  
  
colnames(lf.vexans)[2] <- "temp"
colnames(lf.vexans)[3] <- "trait"


###### Ae. cinereus, Ae. communis, Ae. impiger, Ae. punctor in Alaska ######
lf.aedes <- read_csv("data-raw/aedes_spp.sommerman1969.csv") %>%
  clean_names() %>%
  filter(trait_name == "lf") %>% 
  # Add new columns to provide more info
  mutate(trait_def = "average days alive") %>% 
  mutate(type = "Arctic")

# Convert trait2 to character
lf.aedes$trait2 <- as.character(lf.aedes$trait2)


## Non-Arctic species (for informing priors) -------------------------------
###### Aedes sierrensis ######
lf.sierrensis <- read_csv("data-raw/aedes_sierrensis.Couper2024.csv") %>% 
  clean_names()

# Select relevant columns
lf.sierrensis <- lf.sierrensis %>% 
  filter(!is.na(adult_lifespan)) %>% 
  # Since this dataset provides individual-level data, we will calculate the mean at each temp
  group_by(temp_treatment) %>%
  summarize(trait = mean(adult_lifespan),
            error_pos = sd(adult_lifespan)) %>%
  # Add new columns to provide more info
  mutate(trait_name = "lf") %>% 
  mutate(error_unit = "sd") %>% 
  mutate(trait_def = "days") %>% 
  mutate(genus = "Aedes") %>% 
  mutate(species = "sierrensis") %>% 
  mutate(citation = "Couper_2024_ProcBiolSci") %>% 
  mutate(doi = "10.1098/rspb.2023.2457") %>% 
  mutate(data_source = "calculated from raw data") %>% # raw data from the paper
  mutate(type = "non-Arctic") %>% 
  select(trait_name, temp_treatment, trait, error_pos, error_unit, trait_def, 
         genus, species, citation, doi, data_source, type)

# Rename columns
colnames(lf.sierrensis)[2] <- "temp"



###### Aedes albopictus (from VecTrait) ######
lf.VecTrait <- read_csv("data-raw/lf_VecTrait.csv") %>% 
  clean_names()


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
         interactor1species, citation, doi, figure_table, notes, type)


colnames(lf.VecTrait) <- c("trait_name", "temp", "trait", "error_pos", 
                          "error_neg", "error_unit", "trait_def", "trait2_name", 
                          "trait2", "genus", "species", "citation", "doi", 
                          "data_source", "notes", "type")

lf.VecTrait$trait2 <- as.character(lf.VecTrait$trait2) ## Change the trait2 column to character (so that it can combine with other dataset)


## Since Huxley et al. 2021 and Huxley et al. 2022 paper provide individual-level data, we will calculate the mean at each temp 
# lf.Huxley2021 <- lf.VecTrait %>% 
#   filter(citation == "Huxley et al. 2021. The effect of resource limitation on the temperature-dependance of mosquito fitness. Proc. R. Soc. B. 288: 20203217.") %>% 
#   group_by(trait_name, temp, trait2_name, trait2, genus, species, citation, doi) %>% 
#   summarize(trait = mean(trait),
#             error_pos = sd(trait)) %>% 
#   mutate(trait_def = "duration of life stage") %>% 
#   mutate(data_source = "calculated from supplementary material") %>% 
#   mutate(type = "non-Arctic")
# 
# 
# lf.Huxley2022 <- lf.VecTrait %>% 
#   filter(citation == "Huxley et al. 2022. Competition and resource depletion shape the thermal response of population fitness in Aedes aegypti. Commun. Biol. 5: 66.") %>% 
#   group_by(trait_name, temp, trait2_name, trait2, genus, species, citation, doi) %>% 
#   summarize(trait = mean(trait),
#             error_pos = sd(trait),) %>% 
#   mutate(trait_def = "duration of life stage") %>% 
#   mutate(data_source = "calculated from supplementary material") %>% 
#   mutate(type = "non-Arctic")
# 
# # Remove Huxley 2021 and Huxley 2022 data from lf.VecTrait
# lf.VecTrait <- lf.VecTrait %>% 
#   filter(citation != "Huxley et al. 2021. The effect of resource limitation on the temperature-dependance of mosquito fitness. Proc. R. Soc. B. 288: 20203217.") %>% 
#   filter(citation != "Huxley et al. 2022. Competition and resource depletion shape the thermal response of population fitness in Aedes aegypti. Commun. Biol. 5: 66.") 
  
###### Aedes albopictus (from Mordecai et al. 20190 ######

lf.Mordecai2019 <- read_csv("data-raw/AdultSurvival_Data_Mordecai2019.csv") %>% 
  clean_names() %>% 
  select(!series_id)
  

lf.aalb <- lf.Mordecai2019 %>% 
  filter(host_code == "Aalb") %>% 
  filter(trait_name == "1/mu") %>% 
  # Add new columns to provide more info
  mutate(genus = "Aedes") %>% 
  mutate(species = "albopictus") %>% 
  mutate(type = "non-Arctic") %>% 
  select(trait_name, t, trait, error_pos_si, error_neg_si, trait_2, trait2_name,
         genus, species, citation, figure, notes, type)


## Rename columns
colnames(lf.aalb) <- c("trait_name", "temp", "trait", "error_pos", 
                       "error_neg","trait2_name", "trait2", "genus", # The values of trait2_name and trait2 are swapped for some reasons
                       "species", "citation", "data_source", "notes", "type")



##### Aedes aegypti (from Mordecai et al. 2019) ######

# lf.aegypti <- lf.aegypti %>% 
#   ## get data from Ae. aegypti
#   filter(host_code == "Aaeg") %>% 
#   filter(trait_name == "1/mu") %>% 
#   # select columns that we need
#   select(trait_name, t, trait, error_pos_si, error_neg_si, trait2_name, 
#          trait_2, citation, figure, notes) %>% 
#   # Add new columns to provide more info
#   mutate(trait_name = "lf") %>% #change 1/mu to lf (lf = 1/mu)
#   mutate(genus = "Aedes") %>% 
#   relocate(genus, .after = "trait_2") %>% 
#   mutate(species = "aegypti") %>% 
#   relocate(species, .after = "genus") %>% 
#   mutate(type = "non-Arctic")
# 
# ## Rename columns
# colnames(lf.aegypti) <- c("trait_name", "temp", "trait", "error_pos", "error_neg",
#                           "trait2_name", "trait2", "genus", "species", 
#                           "citation", "data_source", "notes", "type")




# Combine data from vexans and sierrensis into a single dataframe
TraitData_lf <- bind_rows(lf.vexans, lf.aedes, lf.sierrensis, lf.VecTrait, 
                          lf.aalb)

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
                                                     "Ae. sierrensis",
                                                     "Ae. vexans"
  )) +
  facet_grid(rows = vars(type)) +
  theme_bw()

plot.data.lf

ggsave("figures/raw_data/plot.data.lf.png", plot.data.lf, width = 9.83, height = 6.17)




# 8. Eggs per female per gonotrophic cycle (EFGC) -------------------------

## Non-Arctic species (for informing priors) -------------------------------
###### Ae. albopictus ######
# EFGC is named as TFD in this dataset
EFGC.albopictus <- read_csv("data-raw/Fecundity_Data_Mordecai2019.csv") %>% 
  clean_names()


EFGC.albopictus <- EFGC.albopictus %>% 
  ## get data from Ae. albopictus
  filter(host_code == "Aalb") %>% 
  filter(trait_name == "TFD") %>% 
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
  mutate(type = "non-Arctic")


## Rename columns
colnames(EFGC.albopictus) <- c("trait_name", "temp", "trait", "error_pos", "error_neg",
                            "error_unit", "trait_def", "trait2_name", "trait2",
                            "genus", "species", "citation","doi", "data_source", 
                            "notes", "type")


TraitData_EFGC <- EFGC.albopictus

write_csv(TraitData_EFGC, "data-processed/TraitData_EFGC.csv")


## Plot raw data
plot.data.EFGC <- TraitData_EFGC %>% 
  ggplot() +
  geom_point(aes(x = temp, y = trait, colour = citation)) +
  # xlim(c(10,35)) +
  labs(y = "Eggs per female per gonotrophic cycles", x = "Temperature ºC") +
  scale_shape_discrete(name = "Citation", labels = c("Delatte 2009",
                                                    "Ezeakacha 2015",
                                                    "Yee 2016"
  )) +
  facet_grid(rows = vars(type), scales = "free_y") +
  theme_bw()



plot.data.EFGC

ggsave("figures/raw_data/plot.data.EFGC.png", plot.data.EFGC, width = 9.83, height = 6.17)
