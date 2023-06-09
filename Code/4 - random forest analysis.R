## Goal - for disease exposed individuals predict final disease state based on T3 & T0 asv abundances
## Identify ASVs (and higher level groups most important to explaining final disease state

library(tidyverse)
library(tidymodels)
library(themis)
library(finetune)
library(tidytext)
library(vip)

#### Functions ####
variance_explained <- function(recipe, step_number){
  recipe$steps[[step_number]]$res$sdev^2 / sum(recipe$steps[[step_number]]$res$sdev^2)
}

pca_correlations_original_var <- function(pca){
  left_join(
    tidy(pca, 4, matrix = 'rotation') %>%
      select(terms, value, component),
    tidy(pca, 4, type = "variance") %>%
      filter(terms == 'variance') %>%
      select(component, value) %>%
      mutate(component = str_c('PC', component),
             sd = sqrt(value)) %>%
      select(-value),
    by = 'component'
  ) %>%
    mutate(r = value * sd) %>%
    mutate(t = r / sqrt((1 - r^2) / (15 - 2)),
           df = 15 - 2,
           p = 2 * pt(abs(t), df = df, lower.tail = FALSE))
}

#
#### Data ####
normalized_asv_counts <- read_csv('../intermediate_files/fully_preprocessed_samples.csv.gz', show_col_types = FALSE) %>%
  mutate(time = factor(time, ordered = TRUE)) %>%
  mutate(treatment = str_c(time, exposure, susceptability, sep = '_'),
         time_exposure = str_c(time, exposure, sep = '_'),
         timeC = str_extract(time, '[0-9]+') %>% as.numeric,
         across(Kingdom:Species, str_replace_na)) %>%
  mutate(asv_number = str_extract(asv_id, '[0-9]+') %>% as.integer)

#### Preprocess Metadata for Random Forest ####
field_samples <- normalized_asv_counts %>%
  select(sample_id, time, exposure, tank, genotype, final_disease_state) %>%
  distinct %>%
  filter(exposure == 'F') %>%
  select(genotype, time, sample_id) %>%
  pivot_wider(names_from = time,
              values_from = sample_id)

samples_through_time <- normalized_asv_counts %>%
  select(sample_id, time, exposure, tank, genotype, final_disease_state) %>%
  distinct %>%
  filter(exposure == 'D') %>%
  select(-exposure) %>%
  pivot_wider(names_from = time,
              values_from = sample_id) %>%
  filter(!is.na(T3), !is.na(T7)) %>%
  select(-tank) %>%
  left_join(field_samples, 
            by = 'genotype') %>%
  filter(!is.na(T0)) 

#### Preprocess ASVs for Random Forest ####
taxonomic_counts <- normalized_asv_counts %>%
  select(Kingdom:Genus, asv_id, sample_id, read_count, lib.size) %>%
  select(-Kingdom:-Class) %>%
  mutate(asv_id = str_c(Order, Family, Genus, asv_id, sep = ':')) %>%
  pivot_longer(cols = -c(sample_id, read_count, lib.size),
               names_to = 'taxonomic_level',
               values_to = 'taxon_id') %>%
  filter(taxon_id != 'NA') %>%
  group_by(sample_id, taxonomic_level, taxon_id, lib.size) %>%
  summarise(read_count = sum(read_count),
            .groups = 'drop') %>%
  mutate(cpm = read_count / lib.size * 1e6) %>%
  mutate(taxon = str_c(taxonomic_level, taxon_id, sep = '-'),
         .keep = 'unused') %>%
  select(-read_count) %>%
  pivot_wider(names_from = 'taxon',
              values_from = 'cpm') %>%
  rename(total_reads = lib.size) %>% #includes asvs not retained by filtering 
  
  select(sample_id, starts_with('asv_id')) %>%
  identity()

#### Recombine Data ####
prepared_dataset <- samples_through_time %>%
  pivot_longer(cols = starts_with('T'),
               names_to = 'timepoint',
               values_to = 'sample_id') %>%
  left_join(taxonomic_counts,
            by = 'sample_id') %>%
  select(-sample_id) %>%
  pivot_wider(names_from = timepoint,
              values_from = any_of(colnames(taxonomic_counts)),
              names_glue = "{timepoint}_{.value}") %>%
  mutate(final_disease_state = factor(final_disease_state))

#### Plot PCA Contributions ####
pca_results <- recipe(final_disease_state ~ ., data = prepared_dataset) %>%
  update_role(genotype, new_role = "ID") %>%
  step_nzv(all_predictors()) %>%
  step_log(all_predictors(), base = 2, offset = 0.5) %>%
  step_normalize(all_predictors()) %>%
  step_pca(all_predictors(), threshold = 0.75)

var_exp <- prep(pca_results) %>%
  variance_explained(4)

prep(pca_results) %>%
  tidy(4) %>%
  filter(component %in% str_c('PC', 1:6)) %>%
  group_by(component) %>%
  top_n(20, abs(value)) %>%
  ungroup() %>%
  mutate(terms = reorder_within(terms, abs(value), component)) %>%
  ggplot(aes(abs(value), terms, fill = value > 0)) +
  geom_col() +
  facet_wrap(~component, scales = "free_y") +
  scale_y_reordered() +
  labs(x = "Absolute value of contribution",
       y = NULL, fill = "Positive?")


significant_associations <- prep(pca_results) %>%
  pca_correlations_original_var() %>%
  group_by(component) %>%
  mutate(fdr = p.adjust(p, method = 'fdr')) %>%
  ungroup %>%
  filter(fdr < 0.05) %>%
  filter(component %in% str_c('PC', 1:2)) %>%
  pull(terms)

variable_loadings <- prep(pca_results) %>%
  tidy(4, matrix = "rotation") %>%
  filter(component %in% str_c('PC', 1:2)) %>%
  pivot_wider(names_from = "component", 
              values_from = "value") %>%
  filter(terms %in% significant_associations) %>%
  filter(str_detect(terms, 'asv')) %>%
  # select(terms) %>%
  mutate(time = str_extract(terms, '^T[037]'),
         asv_id = str_extract(terms, 'ASV_[0-9]+')) %>%
  select(-terms) %>%
  left_join(normalized_asv_counts %>%
              select(Order:Genus, asv_id) %>%
              distinct,
            by = 'asv_id') 

prep(pca_results) %>%
  juice() %>%
  ggplot(aes(x = PC1, y = PC2)) +
  
  geom_segment(data = variable_loadings %>% mutate(across(starts_with('PC'), ~. * 50)),
               aes(linetype = time, colour = Family),
               xend = 0, yend = 0, 
               arrow = arrow(angle = 20, ends = "first", 
                             type = "closed", 
                             length = grid::unit(8, "pt"))) +
  geom_text(data = variable_loadings %>% mutate(across(starts_with('PC'), ~. * 50)),
                  aes(colour = Family, label = asv_id), show.legend = FALSE) +
  
  geom_point(aes(shape = final_disease_state), alpha = 0.7, size = 2) +
  geom_text(aes(label = genotype), check_overlap = TRUE, hjust = "inward") +
  labs(color = NULL,
       x = str_c('PC1 (', scales::percent(var_exp[1]), ')'),
       y = str_c('PC2 (', scales::percent(var_exp[2]), ')'))


#### Train/Test Split ####
# not ideal because of sample size...but what can you do
# use 10 samples to build model and assess model on 5 (9 and 6 from 2/3 1/3)
# get parameters using 3 fold cv repeated 10 times
coral_split <- initial_split(prepared_dataset, prop = 2/3,
                             strata = final_disease_state)
coral_train <- training(coral_split)
coral_test <- testing(coral_split)

count(coral_train, final_disease_state)
count(coral_test, final_disease_state)

#### Random Forest of all data points ####

#### Random Forest Preprocessing ####
preprocess_recipie <- recipe(final_disease_state ~ ., data = coral_train) %>%
  update_role(genotype, new_role = "ID") %>%
  step_downsample(all_outcomes()) %>%
  step_nzv(all_predictors()) %>%
  step_log(all_predictors(), base = 2, offset = 0.5) %>%
  step_normalize(all_predictors()) %>%
  
  # step_pca(all_predictors(), threshold = 0.75) %>%
  identity()

prep(preprocess_recipie) %>%
  juice()


#### Build Random Forest Model ####
#not enough samples for any real parameter tuning 
forest_model <- rand_forest() %>%
  set_engine('ranger', importance = "impurity") %>%
  set_mode('classification') %>%
  set_args(mtry = tune(),
           trees = tune(),
           min_n = tune())

disease_wflow <- workflow() %>% 
  add_model(forest_model) %>% 
  add_recipe(preprocess_recipie)

#### Train Hyperparameters ####
coral_folds <- vfold_cv(coral_train, 
                        v = 3, repeats = 5, 
                        strata = final_disease_state)

tuning_params <- parameters(list(mtry = mtry(range = c(50, 500)), 
                                 trees = trees(range = c(1e3, 1e4)),
                                 min_n = min_n(range = c(1, 4))))


tune_aov <- tune_race_anova(disease_wflow,
                            resamples = coral_folds,
                            param_info = tuning_params,
                            metrics = metric_set(accuracy, roc_auc),
                            grid = 20,
                            control = control_race(verbose = TRUE,
                                                   verbose_elim = TRUE))
autoplot(tune_aov)


#### Finalize Hyper Parameters ####
final_rf <- finalize_workflow(disease_wflow,
  select_best(tune_aov, "accuracy"))

#### Assess Model Accuracy on Test Dataset ####
disease_fit <- final_rf %>%
  last_fit(coral_split)

disease_fit %>%
  collect_metrics()

disease_fit %>%
  collect_predictions()

#### Get Variable Importance ####
final_rf %>%
  fit(data = prepared_dataset) %>%
  extract_fit_parsnip() %>%
  vip(geom = "point",
      num_features = 10)

#### Plot Counts of top ASVs ####
final_rf %>%
  fit(data = prepared_dataset) %>%
  extract_fit_parsnip() %>% 
  vi() %>%
  arrange(-Importance) %>%
  slice(1:10) %>%
  mutate(asv_id = str_extract(Variable, 'ASV_[0-9]+')) %>%
  filter(!is.na(asv_id)) %>%
  left_join(normalized_asv_counts,
            by = 'asv_id') %>%
  mutate(asv_id = str_c(Family, Genus, asv_id, sep = ': ')) %>%
  
  ggplot(aes(x = time, y = log2_cpm, colour = exposure,
             shape = susceptability)) +
  stat_summary(fun.data = mean_se, position = position_dodge(0.5)) +
  facet_wrap(~asv_id)
 
