## Code to preprocess microbiome sample data into easier to use format

library(tidyverse)
library(magrittr)
library(phyloseq)

#### Read in data ####
preprocess_metadata <- read_csv('../intermediate_files/preprocess_metadata.csv', 
                                show_col_types = FALSE)

microbiome_raw <- read_rds('../Data/ps_fl_tank.rds') %>%
  subset_taxa(Kingdom == "Bacteria" &
                Phylum != "Cyanobacteria" &
                !is.na(Phylum) &
                Family != "Mitochondria" &
                Class != "Chloroplast" & 
                Order != "Chloroplast")

#filter for at least 1000 reads
microbiome_raw <- prune_samples(sample_sums(microbiome_raw)>=1000, microbiome_raw)

#### Process Data ####
#Split into repeated measures portion and input portion
homogenate_samples <- sample_data(microbiome_raw) %>%
  as_tibble(rownames = 'sample_id') %>%
  filter(time %in% c('Acerv', 'T0')) %>%
  mutate(retain_sample = TRUE) %>%
  mutate(time = 'T0',
         exposure = case_when(is.na(exposure) ~ tank, 
                              exposure == 'REP1' ~ 'Field',
                              TRUE ~ exposure),
         final_disease_state = exposure,
         tank = case_when(str_detect(sample_id, 'REP1') ~ 'preTank',
                          str_detect(sample_id, 'Acerv_NA') ~ 'homogenate_fragment',
                          TRUE ~ tank))



exposure_samples <- sample_data(microbiome_raw) %>%
  as_tibble(rownames = 'sample_id') %>%
  select(sample_id) %>%
  inner_join(preprocess_metadata, by = 'sample_id') %>%
  group_by(tank, treatment, genotype) %>% 
  # mutate(retain_sample = n() == 2) %>% #only keep samples when they are replicated at T3 & T7
  mutate(retain_sample = TRUE) %>%
  ungroup %>%
  rename(exposure = treatment)

#Recombine and put data into sample_data slot of phyloseq object
processed_microbiome <- microbiome_raw
sample_data(processed_microbiome) <- bind_rows(homogenate_samples, exposure_samples) %>%
  group_by(genotype) %>%
  mutate(final_disease_state = ifelse(time == 'T0' & tank == 'preTank', 
                                      unique(final_disease_state[exposure == 'D']), final_disease_state)) %>%
  ungroup %>%
  column_to_rownames('sample_id')
processed_microbiome <- subset_samples(processed_microbiome, retain_sample)

#### Sample Summary Table ####
sample_data(processed_microbiome) %>%
  as_tibble %>%
  count(time, exposure, genotype, final_disease_state, tank) %>%
  filter(!tank %in% c('HOMO', 'homogenate_fragment')) %>%
  select(-exposure, -n) %>%
  pivot_wider(names_from = time, values_from = tank,
              values_fn = ~str_c(., collapse = ' + ')) %>%
  group_by(genotype) %>%
  mutate(T0 = if_else(is.na(T0), '', T0)) %>%
  summarise(across(everything(), ~str_c(., collapse = ' + '))) %>%
  mutate(T0 = str_remove(T0, ' \\+ ')) %>%
  relocate(final_disease_state, .after = everything()) %>%
  arrange(genotype) %>% View

#### Output preprocessed microbiome as rds ####
write_rds(processed_microbiome, '../intermediate_files/preprocess_microbiome.rds')


#### Simple Analysis ####
sample_data(processed_microbiome) %>%
  as_tibble %>%
  filter(time == 'T7') %>%
  count(exposure, final_disease_state) %>%
  pivot_wider(names_from = final_disease_state, 
              values_from = n) %>%
  column_to_rownames('exposure') %T>%
  print %>%
  chisq.test()
