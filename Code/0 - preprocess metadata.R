## Code to preprocess metadata into easier to use format

#### Libraries ####
library(tidyverse)

preprocess_metadata <- read_csv('../Data/Florida_Census_outcomes.csv', 
                                show_col_types = FALSE) %>%
  select(-Genotype_Code) %>%
  pivot_longer(cols = c(T3, T7),
               names_to = 'time',
               values_to = 'disease_state') %>%
  janitor::clean_names() %>%
  filter(!is.na(disease_state)) %>%
  group_by(tank, treatment, genotype) %>%
  filter(n() == 2) %>%
  mutate(final_disease_state = disease_state[time == 'T7']) %>%
  ungroup %>%
  mutate(sample_id = str_c(time, treatment, tank, genotype, sep = '_')) %>%
  select(sample_id, everything())

write_csv(preprocess_metadata, '../intermediate_files/preprocess_metadata.csv')

