## Lilian Chan, University of Guelph
## Arctic vector-borne disease transmission suitability model
##
## Purpose: Process the raw data such that each trait has it own csv file
## 
## Table of content:
##    0. Set-up workspace
##    1. Biting rate (a)
##    2. Vector Competence (bc)
##         i) Infection efficiency (c)
##        ii) Transmission efficiency (b)
##    3. Parasite development rate (PDR)
##    4. Lifetime egg production (B)
##    5. Mosquito egg-to-adult survival (pEA)
##         i) Egg viability (EV)
##        ii) Larval survival (pLA)
##    6. Mosquito egg-to-adult development rate (MDR)
##    7. Adult mosquito lifespan (lf)

##########
###### 0. Set-up workspace ----
##########

library(tidyverse)
library(readxl)
library(janitor)
library(ggpubr)



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
         original_error_neg, original_error_unit, trait_def, interactor1genus, 
         interactor1species, citation, doi, figure_table, location)
  

colnames(a.albopictus.Delatte2009) <- c("trait_name", "temp", "trait", 
                                        "error_pos", "error_neg", "error_unit", 
                                        "trait_def", "genus", "species", 
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
         original_error_neg, original_error_unit, trait_def, interactor1genus, 
         interactor1species, citation, doi, figure_table, location)



colnames(a.albopictus.Marini2020) <- c("trait_name", "temp", "trait", 
                                       "error_pos", "error_neg", "error_unit", 
                                       "trait_def", "genus", "species",
                                       "citation", "doi", "data_source", "notes")


## Aedes aegypti ----
a.aegypti.Goindin2015 <- read_csv("data-raw/aedes_aegypti.Goindin2015.csv") %>% 
  clean_names()


a.aegypti.Goindin2015  <- a.aegypti.Goindin2015  %>% 
  group_by(interactor1temp, original_trait_def, interactor1genus, interactor1species, figure_table, location, doi) %>% 
  summarize(trait = mean(original_trait_value)) %>% 
  # Add new columns to provide more info
  mutate(trait_name = "1/a") %>% 
  mutate(citation = "Goindin_2015_PLoSOne") %>% 
  select(trait_name, interactor1temp, trait, original_trait_def, interactor1genus, 
         interactor1species, citation, doi, figure_table, location)


colnames(a.aegypti.Goindin2015) <- c("trait_name", "temp", "trait", "trait_def", 
                                     "genus", "species", "citation","doi", 
                                     "data_source", "notes")


## Data from Marta
a.aegypti <- read_csv("data-raw/a_Data_fromShocket.csv") %>% 
  clean_names()


a.aegypti <- a.aegypti %>% 
  ## get data from Ae. aegypti
  filter(host_code == "Aaeg") %>% 
  # select columns that we need
  select(trait_name, t_c, trait, error_pos_si, error_neg_si, trait2_name, 
         trait_2, citation, figure, notes) %>% 
  # Add new columns to provide more info
  mutate(error_unit = NA) %>% 
  relocate(error_unit, .after = "error_neg_si") %>% # rearrange columns
  mutate(trait_def = NA) %>% 
  relocate(trait_def, .after = "error_unit") %>% 
  mutate(genus = "aedes") %>% 
  relocate(genus, .after = "trait_2") %>% 
  mutate(species = "aegypti") %>% 
  relocate(species, .after = "genus") %>% 
  mutate(doi = NA) %>% 
  relocate(doi, .after = "citation")

  
## Change GCD (gonotrophic cycle duration) to 1/a
a.aegypti$trait_name[a.aegypti$trait_name == "GCD"] <- "1/a"


## Rename columns
colnames(a.aegypti) <- c("trait_name", "temp", "trait", "error_pos", "error_neg",
                         "error_unit", "trait_def", "trait2_name", "trait2",
                         "genus", "species", "citation","doi", "data_source", 
                         "notes")

## Provide more info on the papers
### Focks and Barrera 2006
a.aegypti$citation[a.aegypti$citation == "Focks_Barrera_2006_Research&TrainingTropicalDis_Geneva_Paper"] <- "Focks_2006_Research&TrainingTropicalDis_Geneva_Paper"
a.aegypti$data_source[a.aegypti$citation == "Focks_2006_Research&TrainingTropicalDis_Geneva_Paper"] <- "table 5"


### Focks et al. 1993
a.aegypti$citation[a.aegypti$citation == "Focks_et_al_1993a_JMedEntom"] <- "Focks_1993_JMedEntomol"
a.aegypti$data_source[a.aegypti$citation == "Focks_1993_JMedEntomol"] <- "figure 9"
a.aegypti$doi[a.aegypti$citation == "Focks_1993_JMedEntomol"] <- "10.1093/jmedent/30.6.1003"


### Morin 2015
# Couldn't find the data from the paper. I will just skip the data from this paper
# a.aegypti <- a.aegypti %>% 
#   filter(citation != "Morin_et_al_2015")


## Aedes spp. in Alaska ----
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
  unite(n, rep, col = "notes",)
  

## Combine all data
TraitData_a <-  bind_rows(a.aedes,
                          a.albopictus.Delatte2009, 
                          a.albopictus.Marini2020, 
                          a.aegypti.Goindin2015, 
                          a.aegypti)


# write_csv(TraitData_a, "data-processed/TraitData_a.csv")

## Plot raw data
plot.data.a <- TraitData_a %>%
  mutate(trait = ifelse(trait_name == "1/a", 1/trait, trait)) %>% 
  mutate(type = c(rep("Arctic", 13), rep("non-Arctic", 30))) %>% 
  ggplot(aes(x = temp, y = trait)) +
  geom_point(aes(colour = species
                 #, shape = citation
                 ), size = 2) +
  labs(y = "Biting rate (1/days)", x = "Temperature ºC") +
  scale_colour_discrete(name = "Species", labels = c("Ae. aegypti",
                                                     "Ae. albopictus",
                                                     "Ae. cinereus",
                                                     "Ae. communis",
                                                     "Ae. impiger",
                                                     "Ae. punctor"
  )) +
  # scale_shape_discrete(name = "Citation", labels = c("Delatte_2009)",
  #                                                    "Focks_1993",
  #                                                    "Focks_2006",
  #                                                    "Goindin_2015",
  #                                                    "Marini_2020",
  #                                                    "Morin_2015"
  # )) +
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
c.trivittatus <- read_csv("data-raw/dirofilaria_immitis_aedes_trivittatus.Christensen1978.csv") %>% 
  clean_names() 

c.trivittatus <- c.trivittatus %>% 
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
  mutate(genus = "aedes") %>% 
  mutate(species = "trivittatus") %>% 
  mutate(paras_genus = "dirofilaria") %>% 
  mutate(paras_species = "immitis") %>% 
  mutate(citation = "Christensen_1978_ProcHelmintholSocWash") %>% 
  mutate(doi = NA) %>% #No doi for this paper?
  mutate(data_source = "figure 2") %>% 
  dplyr::select(trait_name, temp, trait, error_pos, error_neg, error_unit, trait_def, 
         trait2_name, trait2, genus, species, paras_genus, paras_species, 
         citation, doi, data_source, notes)

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



###### 4. Mosquito egg-to-adult survival (pEA) ----
## pEA is broken down into two parts: egg viability (EV) and Larval survival (pLA)
## pEA = EV * pLA


##########
###### 4i. Egg viability (EV) ----
##########

# Aedes vexans ----
## Read data
EV.vexans <- read_csv("data-raw/TraitData_EV_Shocket2020.csv") %>% 
  clean_names()

EV.vexans <- EV.vexans %>% 
  filter(host_code == "Avex") %>% 
  # Add new columns to provide more info
  mutate(trait_name = "EV") %>% 
  mutate(trait_def = "percent hatch") %>% 
  mutate(genus = "aedes") %>% 
  mutate(species = "vexans") %>% 
  mutate(error_unit = "95% confidence interval") %>%  
  mutate(citation = "McHaffey_1972_JMedEntomol") %>% 
  mutate(doi = "10.1093/jmedent/9.6.564") %>% 
  select(trait_name, t, trait, error_pos_si, error_neg_si, error_unit, trait_def,
         trait2_name, trait_2, genus, species, citation, doi, figure, notes) 

colnames(EV.vexans) <- c("trait_name", "temp", "trait", "error_pos", 
                         "error_neg", "error_unit", "trait_def", "trait2_name", 
                         "trait2", "genus", "species", "citation", "doi", 
                         "data_source", "notes")


# Aedes dorsalis ----
## Read data
EV.dorsalis <- read_csv("data-raw/TraitData_EV_Shocket2020.csv") %>% 
  clean_names()

EV.dorsalis <- EV.dorsalis %>% 
  filter(host_code == "Ador") %>% 
  # Add new columns to provide more info
  mutate(trait_name = "EV") %>% 
  mutate(trait_def = "percent hatch") %>% 
  mutate(genus = "aedes") %>% 
  mutate(species = "dorsalis") %>% 
  mutate(error_unit = "95% confidence interval") %>%  
  mutate(citation = "McHaffey_1970_JMedEntomol") %>% 
  mutate(doi = "10.1093/jmedent/7.6.631") %>% 
  select(trait_name, t, trait, error_pos_si, error_neg_si, error_unit, trait_def,
         trait2_name, trait_2, genus, species, citation, doi, figure, notes) 

colnames(EV.dorsalis) <- c("trait_name", "temp", "trait", "error_pos", 
                         "error_neg", "error_unit", "trait_def", "trait2_name", 
                         "trait2", "genus", "species", "citation", "doi", 
                         "data_source", "notes")





# Aedes nigromaculis ----
## Read data
EV.nigromaculis <- read_csv("data-raw/TraitData_EV_Shocket2020.csv") %>% 
  clean_names()
  
EV.nigromaculis <- EV.nigromaculis %>% 
  filter(host_code == "Anig") %>% 
  # Add new columns to provide more info
  mutate(trait_name = "EV") %>% 
  mutate(trait_def = "percent hatch") %>% 
  mutate(genus = "aedes") %>% 
  mutate(species = "nigromaculis") %>% 
  mutate(error_unit = "95% confidence interval") %>%  
  mutate(citation = "McHaffey_1972_MosqNews") %>% 
  mutate(doi = "10.5281/zenodo.16126961") %>% 
  select(trait_name, t, trait, error_pos_si, error_neg_si, error_unit, trait_def,
         trait2_name, trait_2, genus, species, citation, doi, figure, notes) 


colnames(EV.nigromaculis) <- c("trait_name", "temp", "trait", "error_pos", 
                           "error_neg", "error_unit", "trait_def", "trait2_name", 
                           "trait2", "genus", "species", "citation", "doi", 
                           "data_source", "notes")




# Aedes triseriatus ----
## Read data
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
  select(trait_name, temp, trait, error_pos, error_neg, error_unit, trait_def,
         genus, species, citation, doi, data_source, notes) 


## Aedes albopictus (from VecTraits) ----
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
  select(trait_name, interactor1temp, original_trait_value, original_error_pos, 
         original_error_neg, original_error_unit, original_trait_def, 
         trait2_name, second_stressor_value, interactor1genus, 
         interactor1species, citation, doi, figure_table, notes)


colnames(EV.VecTrait) <- c("trait_name", "temp", "trait", "error_pos", 
                           "error_neg", "error_unit", "trait_def", "trait2_name", 
                           "trait2", "genus", "species", "citation", "doi", 
                           "data_source", "notes")

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

# write_csv(TraitData_EV, "data-processed/TraitData_EV.csv")

 ## Plot raw data
plot.data.EV <- TraitData_EV %>% 
  mutate(type = c(rep("Arctic", 12), rep("non-Arctic", 29))) %>% 
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

# ggsave("figures/raw_data/plot.data.EV.png", plot.data.EV, width = 9.83, height = 6.17)

##########
###### 4ii. Larval survival (pLA) ----
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
###### 5. Mosquito egg-to-adult development rate (MDR) ----
##########

# Read data
## Aedes nigripes ----
MDR.nigripes <- read_excel("data-raw/aedes_nigripes.Culler2015.xlsx", 
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
MDR.sierrensis <- read_csv("data-raw/aedes_sierrensis.Couper2024.csv")

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
  
  labs(y = "Mosquito development rate (days)", x = "Temperature ºC") +
  scale_colour_discrete(name = "Species", labels = c("Ae. nigripes",
                                                     "Ae. sierrensis"
  )) +
  facet_grid(rows = vars(type)) +
  theme_bw()

plot.data.MDR

# ggsave("figures/raw_data/plot.data.MDR.png", plot.data.MDR, width = 9.83, height = 6.17)



##########
###### 6. Adult mosquito lifespan (lf) ----
##########

# Read data
## Aedes vexans ----
lf.vexans <- read_csv("data-raw/aedes_vexans.Costello1971.csv") %>% 
  clean_names() 

lf.vexans <- lf.vexans %>% 
  # Choose the highest resource level (least restrictive)
  #filter(adult_food == "honey+water" & relative_humidity == "80") %>% 
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
  mutate(genus = "aedes") %>% 
  mutate(species = "vexans") %>% 
  mutate(citation = "Costello_1971_JEconEntomol") %>% 
  mutate(doi = "10.1093/jee/64.1.324") %>% 
  mutate(data_source = "table 1") %>% 
  mutate(notes = NA) %>% 
  dplyr::select(trait_name, treatment_temperature, days_to_50_percent_mortality, 
         error_pos, error_neg, error_unit, trait_def, trait2_name, trait2, 
         genus, species, citation, doi, data_source, notes)
  
  
colnames(lf.vexans)[2] <- "temp"
colnames(lf.vexans)[3] <- "trait"


## Aedes spp in Alaska ----
lf.aedes <- read_csv("data-raw/aedes_spp.sommerman1969.csv") %>%
  clean_names() %>%
  filter(trait_name == "lf")

lf.aedes <- lf.aedes %>%
  # Add new columns to provide more info
  mutate(trait_def = "average days alive")

# Convert trait2 to character
lf.aedes$trait2 <- as.character(lf.aedes$trait2)



## Aedes sierrensis (for informative priors) ----
lf.sierrensis <- read_csv("data-raw/aedes_sierrensis.Couper2024.csv") %>% 
  clean_names()

# Select relevant columns
lf.sierrensis <- lf.sierrensis %>% 
  filter(!is.na(adult_lifespan)) %>% 
  group_by(temp_treatment) %>% 
  summarize(trait = mean(adult_lifespan)) %>% 
  # Add new columns to provide more info
  mutate(trait_name = "lf") %>% 
  mutate(trait_def = "days") %>% 
  mutate(genus = "aedes") %>% 
  mutate(species = "sierrensis") %>% 
  mutate(citation = "Couper_2024_ProcBiolSci") %>% 
  mutate(doi = "10.1098/rspb.2023.2457") %>% 
  mutate(data_source = "calculated from raw data") %>% # raw data from the paper
  # Combine info from population and sample_id into a new column called "notes"
  # unite(population, sample_id, sep = "_", col = "notes") %>% 
  dplyr::select(trait_name, temp_treatment, trait, trait_def, genus, species,
         citation, doi, data_source)

# Rename columns
colnames(lf.sierrensis)[2] <- "temp"
# colnames(lf.sierrensis)[3] <- "trait"



## Aedes aegypti (from VecTrait) ----
lf.VecTrait <- read_csv("data-raw/lf_VecTrait.csv") %>% 
  clean_names()


lf.VecTrait <- lf.VecTrait  %>% 
  filter(interactor1sex == "female") %>% 
  # Add new columns to provide more info
  mutate(trait_name = "lf") %>% 
  mutate(trait2_name = ifelse(
    is.na(second_stressor), NA, # if no second_stressor, leave trait2_name empty
    paste0(second_stressor, " (", second_stressor_unit, ")"))) %>% # Combine second_stressor and the unit
  select(trait_name, interactor1temp, original_trait_value, original_error_pos, 
         original_error_neg, original_error_unit, original_trait_def, 
         trait2_name, second_stressor_value, interactor1genus, 
         interactor1species, citation, doi, figure_table, notes)


colnames(lf.VecTrait) <- c("trait_name", "temp", "trait", "error_pos", 
                          "error_neg", "error_unit", "trait_def", "trait2_name", 
                          "trait2", "genus", "species", "citation", "doi", 
                          "data_source", "notes")

lf.VecTrait$trait2 <- as.character(lf.VecTrait$trait2) ## Change the trait2 column to character (so that it can combine with other dataset)


## Calculate mean lifespan for each temperature for Huxley et al. 2021 and Huxley et al. 2022 paper

lf.Huxley2021 <- lf.VecTrait %>% 
  filter(citation == "Huxley et al. 2021. The effect of resource limitation on the temperature-dependance of mosquito fitness. Proc. R. Soc. B. 288: 20203217.") %>% 
  group_by(trait_name, temp, trait2_name, trait2, genus, species, citation, doi) %>% 
  summarize(trait = mean(trait)) %>% 
  mutate(trait_def = "duration of life stage") %>% 
  mutate(data_source = "calculated from supplementary material")


lf.Huxley2022 <- lf.VecTrait %>% 
  filter(citation == "Huxley et al. 2022. Competition and resource depletion shape the thermal response of population fitness in Aedes aegypti. Commun. Biol. 5: 66.") %>% 
  group_by(trait_name, temp, trait2_name, trait2, genus, species, citation, doi) %>% 
  summarize(trait = mean(trait)) %>% 
  mutate(trait_def = "duration of life stage") %>% 
  mutate(data_source = "calculated from supplementary material")

lf.VecTrait <- lf.VecTrait %>% 
  filter(citation != "Huxley et al. 2021. The effect of resource limitation on the temperature-dependance of mosquito fitness. Proc. R. Soc. B. 288: 20203217.") %>% 
  filter(citation != "Huxley et al. 2022. Competition and resource depletion shape the thermal response of population fitness in Aedes aegypti. Commun. Biol. 5: 66.") 
  

## Aedes aegypti (from Shocket) ----
lf.aegypti <- read_csv("data-raw/AdultSurvival_Data_fromShocket.csv") %>% 
  clean_names() 

lf.aegypti <- lf.aegypti %>% 
  ## get data from Ae. aegypti
  filter(host_code == "Aaeg") %>% 
  filter(trait_name == "1/mu") %>% 
  # select columns that we need
  select(trait_name, t, trait, error_pos_si, error_neg_si, trait2_name, 
         trait_2, citation, figure, notes) %>% 
  # Add new columns to provide more info
  mutate(trait_name = "lf") %>% #change 1/mu to lf (lf = 1/mu)
  mutate(error_unit = NA) %>% 
  relocate(error_unit, .after = "error_neg_si") %>% # rearrange columns
  mutate(trait_def = NA) %>% 
  relocate(trait_def, .after = "error_unit") %>% 
  mutate(genus = "aedes") %>% 
  relocate(genus, .after = "trait_2") %>% 
  mutate(species = "aegypti") %>% 
  relocate(species, .after = "genus") %>% 
  mutate(doi = NA) %>% 
  relocate(doi, .after = "citation")

## Rename columns
colnames(lf.aegypti) <- c("trait_name", "temp", "trait", "error_pos", "error_neg",
                          "error_unit", "trait_def", "trait2_name", "trait2",
                          "genus", "species", "citation","doi", "data_source", 
                          "notes")


# Combine data from vexans and sierrensis into a single dataframe
TraitData_lf <- bind_rows(lf.vexans, lf.aedes, lf.sierrensis, lf.aegypti, lf.VecTrait, lf.Huxley2021, lf.Huxley2022)

# write_csv(TraitData_lf, "data-processed/TraitData_lf.csv")


## Plot raw data
plot.data.lf <- TraitData_lf %>% 
  mutate(type = c(rep("Arctic", 54), rep("non-Arctic", 97))) %>% 
  ggplot(aes(x = temp, y = trait, colour = species)) +
  geom_point(aes(colour = species), position = "jitter") +
  
  ## Since the Ae. sierrensis has many data, I will just plot the mean±SE
  # geom_point(data = ~filter(.x, type == "Arctic")) +
  # geom_point(data = ~filter(.x, type == "Arctic")) +
  # stat_summary(data = ~filter(.x, type == "non-Arctic"),
  #              fun = mean, geom = "point") +
  # stat_summary(data = ~filter(.x, type == "non-Arctic"),
  #              fun.data = "mean_se", geom = "errorbar") +
  labs(y = "Mosquito adult lifespan (days)", x = "Temperature ºC") +
  scale_colour_discrete(name = "Species", labels = c("Ae. aegypti",
                                                     "Ae. albopictus",
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

# ggsave("figures/raw_data/plot.data.lf.png", plot.data.lf, width = 9.83, height = 6.17)


##########
###### 7. Lifetime egg production (B)----
##########
# Read data
## Aedes hexodontus ----
B.hexodontus <- read_csv("data-raw/aedes_hexodontus.barlow1955.csv") %>% 
  clean_names()


## Aedes spp in Alaska ----
B.aedes.spp <- read_csv("data-raw/aedes_spp.sommerman1969.csv") %>% 
  clean_names() %>% 
  filter(trait_name == "B")

B.aedes.spp$trait2 <- as.character(B.aedes.spp$trait2)


## Aedes aegypti: EFD ----
EFD.aegypti <- read_csv("data-raw/Fecundity_Data_fromShocket.csv") %>% 
  clean_names()


EFD.aegypti <- EFD.aegypti %>% 
  ## get data from Ae. aegypti
  filter(host_code == "Aaeg") %>% 
  filter(trait_name == "EFD") %>% 
  # select columns that we need
  dplyr::select(trait_name, t, trait, error_pos_si, error_neg_si, trait2_name, 
         trait_2, citation, figure, notes) %>% 
  # Add new columns to provide more info
  mutate(error_unit = NA) %>% 
  relocate(error_unit, .after = "error_neg_si") %>% # rearrange columns
  mutate(trait_def = NA) %>% 
  relocate(trait_def, .after = "error_unit") %>% 
  mutate(genus = "aedes") %>% 
  relocate(genus, .after = "trait_2") %>% 
  mutate(species = "aegypti") %>% 
  relocate(species, .after = "genus") %>% 
  mutate(doi = NA) %>% 
  relocate(doi, .after = "citation")


## Rename columns
colnames(EFD.aegypti) <- c("trait_name", "temp", "trait", "error_pos", "error_neg",
                           "error_unit", "trait_def", "trait2_name", "trait2",
                           "genus", "species", "citation","doi", "data_source", 
                           "notes")

## Provide more info on the papers
### Beserra_2009
# EFD.aegypti$citation[EFD.aegypti$citation == "Beserra_2009"] <- "Focks_2006_Research&TrainingTropicalDis_Geneva_Paper"
# EFD.aegypti$data_source[EFD.aegypti$citation == "Focks_2006_Research&TrainingTropicalDis_Geneva_Paper"] <- "table 5"


## Aedes aegypti: lf ----
## For lf data for Ae. aegypti, we will filter from TraitData_lf.csv
lf.aegypti <- read_csv("data-processed/TraitData_lf.csv") %>% 
  clean_names() %>% 
  filter(species == "aegypti")

TraitData_B <- bind_rows(B.hexodontus, B.aedes.spp, EFD.aegypti, lf.aegypti)

# write_csv(TraitData_B, "data-processed/TraitData_B.csv")


## Plot raw data
plot.data.B <- TraitData_B %>% 
  mutate(type = c(rep("Arctic",6), rep("non-Arctic: EFD", 30), rep("non-Arctic: lf", 74))) %>% 
  ggplot() +
  geom_point(aes(x = temp, y = trait, colour = species)) +
  # xlim(c(10,35)) +
  labs(y = "lifetime egg production", x = "Temperature ºC") +
  scale_colour_discrete(name = "species", labels = c("Ae. aegypti",
                                                     "Ae. cinereus",
                                                     "Ae. communis",
                                                     "Ae. hexodontus",
                                                     "Ae. impiger",
                                                     "Ae. punctor"
  )) +
  facet_grid(rows = vars(type), scales = "free_y") +
  theme_bw()


plot.data.B


# ggsave("figures/raw_data/plot.data.fecundity.png", plot.data.B, width = 9.83, height = 6.17)


