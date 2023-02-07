## Code to preprocess metadata into easier to use format

#### Libraries ####
library(tidyverse)
library(magrittr)


#### Read in & Process Data ####
preprocess_metadata <- read_csv('../Data/Florida_Census_outcomes.csv', 
                                show_col_types = FALSE) %>%
  select(-Genotype_Code) %>%
  pivot_longer(cols = c(T3, T7),
               names_to = 'time',
               values_to = 'disease_state') %>%
  janitor::clean_names() %>%
  filter(!is.na(disease_state)) %>%
  filter(!tank %in% c('H2')) %>% # Remove tanks which crashed due to anoxia
  group_by(tank, treatment, genotype) %>%
  filter(n() == 2) %>%
  mutate(final_disease_state = disease_state[time == 'T7']) %>%
  ungroup %>%
  mutate(sample_id = str_c(time, treatment, tank, genotype, sep = '_')) %>%
  select(sample_id, everything())

#### Output preprocessed metadata ####
write_csv(preprocess_metadata, '../intermediate_files/preprocess_metadata.csv')

#### Simple analysis of outcome v treatment
preprocess_metadata %>%
  filter(time == 'T7') %>%
  count(treatment, final_disease_state) %>%
  pivot_wider(names_from = final_disease_state, 
              values_from = n) %>%
  column_to_rownames('treatment') %T>%
  print %>%
  chisq.test()
