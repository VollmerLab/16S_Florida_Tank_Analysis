## Code to preprocess metadata into easier to use format

#### Libraries ####
library(tidyverse)
library(magrittr)


#### Read in & Process Data ####
preprocess_metadata <- read_csv('../Data/acropora_master_data.csv', 
         show_col_types = FALSE) %>% 
  select(genotype, clone_group, starts_with('diseaseResistance'), 
         starts_with('fds'), starts_with('16s')) %>%
  filter(!if_all(starts_with('16s'), is.na)) %>%
  select(-diseaseResistance_H) %>% 
  rename(resistance = diseaseResistance_D,
         susceptability = diseaseResistance_categorical,
         `16s_T0_Field` = `16s_T0`) %>% 
  pivot_longer(cols = starts_with('16s'),
               names_to = c('time', 'tank'),
               names_pattern = '16s_(.*)_(.*)',
               values_to = 'reads') %>% 
  filter(!is.na(reads),
         !str_detect(genotype, 'homogenate')) %>%
  # mutate(tank2 = tank) %>%
  pivot_longer(cols = starts_with('fds'),
               names_to = 'tank2',
               values_to = 'final_disease_state',
               names_prefix = 'fds_') %>%
  filter((tank2 == tank | time == 'T0')) %>%
  
  # filter(time == 'T0') %>%
  mutate(final_disease_state = if_else(time == 'T0', NA_character_, final_disease_state)) %>%
  select(-tank2) %>%
  distinct %>% 
  mutate(treatment = str_extract(tank, '[DHF]'),
         sample_id = str_c(time, treatment, tank, genotype, sep = '_')) %>%
  select(sample_id, everything()) %>%
  rename(exposure = treatment)

# preprocess_metadata <- read_csv('../Data/Florida_Census_outcomes.csv', 
#                                 show_col_types = FALSE) %>%
#   select(-Genotype_Code) %>%
#   pivot_longer(cols = c(T3, T7),
#                names_to = 'time',
#                values_to = 'disease_state') %>%
#   janitor::clean_names() %>%
#   filter(!is.na(disease_state)) %>%
#   filter(!tank %in% c('H2')) %>% # Remove tanks which crashed due to anoxia
#   group_by(tank, treatment, genotype) %>%
#   filter(n() == 2) %>%
#   mutate(final_disease_state = disease_state[time == 'T7']) %>%
#   ungroup %>%
#   mutate(sample_id = str_c(time, treatment, tank, genotype, sep = '_')) %>%
#   select(sample_id, everything())

#### Get disease susceptibility for unlabeled samples ####
resistance_susceptable_data <- preprocess_metadata %>%
  select(genotype, resistance, susceptability) %>%
  distinct %>%
  mutate(is_resistant = susceptability == 'R')

resistance_model <- glm(is_resistant ~ resistance, 
    family = 'binomial',
    data = filter(resistance_susceptable_data, !is.na(susceptability)))


preprocess_metadata <- resistance_susceptable_data %>%
  mutate(pred_sus = predict(resistance_model, newdata = resistance_susceptable_data, type = 'response') %>%
           round) %>%
  mutate(pred_susceptability = c('S', 'R')[pred_sus + 1]) %>%
  select(genotype, pred_susceptability) %>%
  rename(susceptability = pred_susceptability) %>%
  full_join(select(preprocess_metadata, -susceptability), .,
            by = 'genotype')

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
