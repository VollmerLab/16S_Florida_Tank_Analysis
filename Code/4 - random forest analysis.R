## Goal - for disease exposed individuals predict final disease state based on T3 & T0 asv abundances
## Identify ASVs (and higher level groups most important to explaining final disease state

##TODO - add in dimensionality reduction options (e.g. PLS & UMAP)?
##TODO - add in other classifiers (e.g. naive bayes, discriminant analysis)?

library(tidyverse)
library(tidymodels)
library(themis)
library(finetune)
library(tidytext)
library(vip)
library(kernelshap)
library(shapviz)

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

#### NMDS Samples ####
library(vegan)

nmds <- metaMDS(select(prepared_dataset, -final_disease_state) %>%
          column_to_rownames('genotype'))

plot(nmds)

corr_asvs <- envfit(nmds, select(prepared_dataset, -final_disease_state) %>%
         column_to_rownames('genotype'))


asv_nmds <- scores(nmds, 'species') %>%
  as_tibble(rownames = 'asv_id') %>%
  filter(asv_id %in% names(corr_asvs$vectors$pvals[corr_asvs$vectors$pvals < 0.05])) %>%
  mutate(time = str_extract(asv_id, '^T[037]'),
         asv_id = str_extract(asv_id, 'ASV_[0-9]+')) %>%
  left_join(normalized_asv_counts %>%
              select(Order:Genus, asv_id) %>%
              distinct,
            by = 'asv_id') %>%
  mutate(asv_id = str_c(time, asv_id, sep = ': '))

scores(nmds, 'sites') %>%
  as_tibble(rownames = 'genotype') %>%
  left_join(select(samples_through_time, genotype, final_disease_state),
            by = 'genotype') %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  
  geom_point(data = asv_nmds,
             aes(colour = Family), size = 1) +
  geom_text(data = asv_nmds,
            aes(label = asv_id, colour = Family), 
            check_overlap = TRUE, hjust = "inward", 
            size = 3) +
  
  geom_point(aes(shape = final_disease_state), alpha = 1, size = 2) +
  geom_text(aes(label = genotype), check_overlap = TRUE, hjust = "inward") +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  labs(color = 'ASV Family',
       x = str_c('NMDS1'),
       y = str_c('NMDS2'),
       shape = 'Disease') +
  theme_classic() +
  theme(axis.text = element_text(colour = 'black'),
        axis.title = element_text(colour = 'black', size = 14),
        panel.background = element_rect(colour = 'black'))
ggsave('../Figures/nmds_diseaseExposed_genotype_asv.png', height = 10, width = 10)

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

#### Random Forest Preprocessing ####
preprocess_recipie <- recipe(final_disease_state ~ ., data = coral_train) %>%
  update_role(genotype, new_role = "ID") %>%
  step_downsample(all_outcomes()) %>%
  step_nzv(all_predictors()) %>%
  step_log(all_predictors(), base = 2, offset = 0.5) %>%
  step_normalize(all_predictors()) %>%
  identity()

#### Build Random Forest Model ####
forest_model <- rand_forest() %>%
  set_engine('ranger', importance = 'permutation', probability = TRUE) %>%
  set_mode('classification') %>%
  set_args(mtry = tune(),
           trees = tune(),
           min_n = tune())

disease_wflow <- workflow() %>% 
  add_model(forest_model) %>% 
  add_recipe(preprocess_recipie)

#### Train Hyperparameters ####
coral_folds <- vfold_cv(coral_train, 
                        v = 3, repeats = 10, 
                        strata = final_disease_state)

tuning_params <- parameters(list(mtry = mtry(range = c(10, 250)), 
                                 trees = trees(range = c(1e2, 1e4)),
                                 min_n = min_n(range = c(1, 4))))


# tune_aov <- tune_race_anova(disease_wflow,
#                             resamples = coral_folds,
#                             param_info = tuning_params,
#                             metrics = metric_set(accuracy),
#                             grid = 20,
#                             control = control_race(verbose = TRUE,
#                                                    verbose_elim = TRUE,
#                                                    alpha = 0.1,
#                                                    burn_in = 5))
# autoplot(tune_aov)

initial_grid_tune <- tune_grid(disease_wflow,
                               resamples = coral_folds,
                               param_info = tuning_params,
                               metrics = metric_set(accuracy),
                               grid = 10,
                               control = control_grid(verbose = TRUE))

bayes_tune <- tune_bayes(disease_wflow,
                         resamples = coral_folds,
                         param_info = tuning_params,
                         metrics = metric_set(accuracy),
                         iter = 100,
                         initial = initial_grid_tune,
                         objective = exp_improve(),
                         control = control_bayes(verbose = TRUE,
                                                 verbose_iter = TRUE,
                                                 no_improve = 10,
                                                 uncertain = 3))
autoplot(bayes_tune)

#### Finalize Hyper Parameters ####
final_rf <- finalize_workflow(disease_wflow,
  select_best(bayes_tune, "accuracy"))

#### Assess Model Accuracy on Test Dataset ####
disease_fit <- final_rf %>%
  last_fit(coral_split)

disease_fit %>%
  collect_metrics()

disease_fit %>%
  collect_predictions()

#### Fit Final Model for Variable Importance ####
fit_rf <- final_rf %>%
  fit(data = prepared_dataset)

#### Get Variable Importance ####
fit_rf %>%
  extract_fit_parsnip() %>% 
  vi() %>%
  mutate(asv_id = str_extract(Variable, 'ASV_[0-9]+'),
         importance_time = str_extract(Variable, '^T[037]')) %>%
  group_by(importance_time) %>%
  slice_max(n = 12, order_by = Importance) %>%
  ungroup %>%
  filter(Importance > 0) %>%
  mutate(Variable = reorder_within(Variable, Importance, 
                                   within = importance_time)) %>%
  ggplot(aes(x = Importance, y = Variable)) +
  geom_segment(xend = 0, aes(yend = Variable)) +
  geom_point(size = 2) +
  scale_y_reordered() +
  facet_wrap(~importance_time, scales = 'free_y',
             ncol = 1) +
  theme_classic() +
  labs(y = NULL) +
  theme(panel.background = element_rect(colour = 'black'),
        strip.background = element_blank())
ggsave('../Figures/random_forest_asvImportance.png', height = 15, width = 7)

#### SHAP Plots ####
shap_vals <- kernelshap(fit_rf, 
                        X = select(prepared_dataset, -final_disease_state), 
                        bg_X = prepared_dataset,
                        pred_fun = predict_function_rf,
                        verbose = TRUE)

shap_viz <- shapviz(shap_vals)

make_shap_force_plot <- function(shap_viz, max_asvs = 15){
  force_data <- shap_viz$S %>%
    as_tibble() %>%
    select(-genotype) %>%
    bind_cols(select(prepared_dataset, genotype, final_disease_state), .) %>%
    pivot_longer(cols = -c(genotype, final_disease_state),
                 names_to = 'asv_id',
                 values_to = 'shap_val') %>%
    group_by(genotype, final_disease_state) %>%
    arrange(-abs(shap_val), .by_group = TRUE) %>%
    mutate(a = row_number(-abs(shap_val))) %>%
    
    mutate(asv_id = case_when(a > max_asvs ~ "Other", TRUE ~ as.character(asv_id))) %>%
    group_by(genotype, final_disease_state, asv_id) %>%
    summarize(shap_val = sum(shap_val),
              .groups = 'drop') %>%
    mutate(direction = shap_val > 0) %>%
    mutate(direction = case_when(final_disease_state == 'D' ~ direction,
                                 TRUE ~ !direction)) %>%
    arrange(genotype, final_disease_state, direction, abs(shap_val)) %>%
    group_by(genotype, final_disease_state) %>%
    mutate(to = cumsum(shap_val),
           from = lag(to, default = 0)) %>%
    mutate(across(c(to, from), ~. + get_baseline(shap_viz))) %>%
    ungroup %>%
    mutate(asv_id = str_remove(asv_id, '_asv_id'))
  
  
  force_data %>%
    ggplot(aes(xmin = from, xmax = to, y = genotype, 
               # colour = factor(shap_val < 0, levels = c(FALSE, TRUE)),
               fill = asv_id)) +
    geom_vline(xintercept = shap_viz$baseline, linetype = 'dashed') +
    gggenes::geom_gene_arrow(show.legend = TRUE, aes(group = factor(shap_val < 0, levels = c(FALSE, TRUE))),
                             arrowhead_width = grid::unit(2, "mm"),
                             position = position_dodge(1)) +
    
    # geom_point(data = force_data %>%
    #                group_by(genotype, final_disease_state) %>%
    #                summarise(prob_d = sum(shap_val),
    #                          .groups = 'drop') %>%
    #                mutate(prob_d = prob_d + get_baseline(shap_viz)),
    #              inherit.aes = FALSE,
    #              aes(y = genotype, x = prob_d, colour = final_disease_state)) +
    geom_text(data = force_data %>%
                group_by(genotype, final_disease_state) %>%
                summarise(prob_d = sum(shap_val),
                          .groups = 'drop') %>%
                mutate(prob_d = prob_d + get_baseline(shap_viz)),
              inherit.aes = FALSE, hjust = 'outward',
              aes(y = genotype, x = prob_d, label = scales::percent(prob_d,
                                                                    accuracy = 0.1))) +
    
    facet_grid(final_disease_state ~ ., scales = 'free_y',
               labeller = labeller(final_disease_state = c('D' = 'Disease', 'H' = 'Healthy'))) +
    scale_x_continuous(labels = scales::percent_format(),
                       limits = c(0, 1)) +
    labs(y = NULL,
         x = 'Probability of Disease',
         fill = 'ASV ID') +
    guides(fill = guide_legend(title.position = 'top', title.hjust = 0.5, nrow = 5),
           colour = 'none') +
    theme_classic() +
    theme(legend.position = 'bottom',
          axis.text = element_text(colour = 'black', size = 12),
          panel.background = element_rect(colour = 'black'),
          strip.background = element_blank()) 
}

shap_plot_disease <- make_shap_force_plot(shap_viz, 15)
ggsave('../Figures/rf_shap_disease_expose.png',
       plot = shap_plot_disease, 
       height = 15,
       width = 10)

#### Plot Counts of top ASVs ####
fit_rf %>%
  extract_fit_parsnip() %>% 
  vi() %>%
  mutate(asv_id = str_extract(Variable, 'ASV_[0-9]+'),
         importance_time = str_extract(Variable, '^T[037]')) %>%
  slice_max(n = 12, order_by = Importance) %>%
  ungroup %>%
  filter(Importance > 0) %>%
  ungroup %>%
  filter(!is.na(asv_id)) %>%
  left_join(normalized_asv_counts,
            by = 'asv_id', relationship = 'many-to-many') %>%
  mutate(asv_id = str_c(Family, Genus, asv_id, sep = ': '),
         asv_id = fct_reorder(asv_id, -Importance)) %>%
  
  ggplot(aes(x = time, y = log2_cpm, colour = exposure,
             shape = susceptability)) +
  stat_summary(fun.data = mean_se, position = position_dodge(0.5)) +
  facet_wrap( ~ asv_id) +
  labs(y = 'log2(cpm)',
       x = NULL,
       colour = 'Exposure',
       shape = 'Susceptability') +
  theme_classic() +
  theme(panel.background = element_rect(colour = 'black'),
        strip.background = element_blank())
ggsave('../Figures/rf_important_asv_time.png', height = 15, width = 15)

#### SHAP Plot for Healthy exposed samples ####
healthy_through_time <- normalized_asv_counts %>%
  select(sample_id, time, exposure, tank, genotype, final_disease_state) %>%
  distinct %>%
  filter(exposure == 'H') %>%
  select(-exposure) %>%
  pivot_wider(names_from = time,
              values_from = sample_id) %>%
  filter(!is.na(T3), !is.na(T7)) %>%
  select(-tank) %>%
  left_join(field_samples, 
            by = 'genotype') %>%
  filter(!is.na(T0)) 

healthy_dataset <- healthy_through_time %>%
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

healthy_shap_vals <- kernelshap(fit_rf, 
                        X = select(prepared_dataset, -final_disease_state), 
                        bg_X = healthy_dataset,
                        pred_fun = predict_function_rf,
                        verbose = TRUE)

healthy_shap_viz <- shapviz(healthy_shap_vals)


shap_plot_healthy <- make_shap_force_plot(healthy_shap_viz, 15)
ggsave('../Figures/rf_shap_healthy_expose.png',
       plot = shap_plot_healthy, 
       height = 15,
       width = 10)
