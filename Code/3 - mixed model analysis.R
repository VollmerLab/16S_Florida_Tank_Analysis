##TODO - Set up/evaluate contrasts for relevant/interesting bacterial responses
#           -want to identify pathogens, opportunists, probiotics
##TODO - 

#### Libraries ####
library(tidyverse)
library(magrittr)
library(lmerTest)
library(emmeans)
library(multidplyr)
library(ggupset)
library(qvalue)
library(patchwork)

refit_models <- TRUE

#### Functions ####
cluster <- new_cluster(parallel::detectCores() - 1)
cluster_library(cluster, c('dplyr', 'lmerTest', 'emmeans', 'stringr', 'tidyr'))

fit_model <- function(formula, data, use_weights = TRUE){
  if(!use_weights){
    data$weights <- 1
  }
  
  full_model <- lmer(formula, 
                     weights = data$weight,
                     data = data, 
                     REML = TRUE,
                     control = variancePartition:::vpcontrol)
  
  
  main_formula <- as.character(formula)
  re_formula <- str_c(main_formula[2], main_formula[1], str_extract_all(main_formula[3], '\\(.*\\)')) %>%
    as.formula()
  
  re_model <- lmer(re_formula, 
                   weights = data$weight,
                   data = data, 
                   REML = TRUE,
                   control = variancePartition:::vpcontrol)
  tibble(model = list(full_model), re_model = list(re_model))
}

# model <- asv_models$model[[1]]; re_model <- asv_models$re_model[[1]]; random_anova <- asv_models$random_anova[[1]]
process_model <- function(model, re_model, random_anova){
  #create type 3 anova table with KR based p-values, marginal & conditional r2 and eta2 effect size
  #also output as single row all anova results sorted nicely
  #needs fit model and null model with only random effects 
  aov_tab <- anova(model, type = '3', ddf = 'Kenward-Roger')
  
  global_row <- anova(model, re_model) %>% 
    broom::tidy() %>%
    filter(term == 'model') %>%
    select(statistic, df, p.value) %>%
    rename(chisq = statistic,
           pvalue = p.value) %>%
    rename_with(~str_c(., '_global'))
  
  r2_row <- performance::r2(model) %>%
    as_tibble
  
  aov_row <- as_tibble(aov_tab, rownames = 'term') %>%
    mutate(term = str_replace(term, 'time_exposure', 'timeXexposure')) %>%
    rename(ss = 'Sum Sq',
           ms = 'Mean Sq',
           n.DF = NumDF,
           d.DF = DenDF,
           fvalue = 'F value',
           pvalue = 'Pr(>F)') %>%
    rowwise %>%
    mutate(eta2Partial = effectsize::F_to_eta2(fvalue, n.DF, d.DF, ci = NULL)$Eta2_partial) %>%
    tidyr::pivot_wider(names_from = term,
                       values_from = where(is.numeric),
                       names_vary = 'slowest') %>%
    rename_with(~str_replace_all(., ':', 'X'))
  
  varDecomp_row <- VarCorr(model) %>%
    as_tibble() %>%
    mutate(varComp = sdcor^2 / sum(sdcor^2)) %>%
    select(grp, varComp) %>%
    filter(grp != 'Residual') %>%
    left_join(random_anova %>%
                  as_tibble(rownames = 'term') %>%
                  mutate(term = str_extract(term, '\\| [0-9a-zA-Z]+') %>%
                           str_remove('\\| +')) %>%
                  filter(!is.na(term)) %>%
                  select(term, Df, LRT, `Pr(>Chisq)`) %>%
                  rename(df = Df,
                         chisq = LRT,
                         pvalue = `Pr(>Chisq)`),
                by = c('grp' = 'term')) %>%
    tidyr::pivot_wider(names_from = 'grp',
                       values_from = c('varComp', 'df', 'chisq', 'pvalue'),
                       names_vary = 'slowest') #%>%
    # rename_with(~str_replace_all(., '_', '.'))
  
  tibble::tibble(anova_table = list(aov_tab)) %>%
    bind_cols(global_row, r2_row, varDecomp_row, ., aov_row)
}

# posthoc <- asv_models$posthoc[[1]]
process_postHoc <- function(posthoc){
  post_row <- as_tibble(posthoc) %>%
    rename(tvalue = t.ratio,
           pvalue = p.value) %>%
    pivot_wider(names_from = c('contrast'),
                values_from = c('estimate', 'SE', 'df', 'tvalue', 'pvalue'),
                names_vary = 'slowest')
  post_row
}

safe_qvalue <- possibly(.f = ~qvalue(.)$qvalues, otherwise = NA_real_)

# df <- asv_models
reorder_columns <- function(df){
  p_cols <- str_subset(colnames(df), 'pvalue')
  fdr_cols <- str_replace(p_cols, 'pvalue', 'fdr')
  q_cols <- str_replace(p_cols, 'pvalue', 'qvalue')
  
  for(col_num in 1:length(p_cols)){
    df <- relocate(df, fdr_cols[col_num], q_cols[col_num], .after = p_cols[col_num])
  }
  df
}

p_adjust <- function(df, exclude_cols = NA_character_){
  exclude_cols <- if_else(is.na(exclude_cols), '@@@', exclude_cols)
  mutate(df, across(c(contains('pvalue'), -contains(exclude_cols)), ~p.adjust(., method = 'fdr'),
                .names = 'fdr_{.col}')) %>% 
    rename_with(~str_replace_all(., 'fdr_pvalue', 'fdr')) %>% 
    
    mutate(across(c(contains('pvalue'), -contains(exclude_cols)), safe_qvalue,
                  .names = 'qvalue_{.col}')) %>%
    rename_with(~str_replace_all(., 'qvalue_pvalue', 'qvalue')) %>%
    reorder_columns 
}

cluster_copy(cluster, c('fit_model', 'process_model', 'process_postHoc'))

#
#### Data ####
normalized_asv_counts <- read_csv('../intermediate_files/fully_preprocessed_samples.csv.gz', show_col_types = FALSE) %>%
  mutate(time = factor(time, ordered = TRUE)) %>%
  mutate(treatment = str_c(time, exposure, susceptability, sep = '_'),
         time_exposure = str_c(time, exposure, sep = '_'),
         timeC = str_extract(time, '[0-9]+') %>% as.numeric,
         across(Kingdom:Species, str_replace_na)) %>%
  mutate(asv_number = str_extract(asv_id, '[0-9]+') %>% as.integer) #%>%
  # filter(asv_number <= 200)

#
#### Plot raw data for an ASV ####
normalized_asv_counts %>%
  # filter(asv_id == 'ASV_3') %>%
  
  nest_by(across(c('asv_id', Kingdom:Species))) %>%
  ungroup %>%
  sample_n(1) %>%
  unnest(data) %>%
  
  ggplot(aes(x = time, y = log2_cpm, colour = exposure, shape = susceptability)) +
  geom_jitter() +
  stat_summary(fun.data = 'mean_se', position = position_dodge(1), size = 1.5)

#### Set-up Post-hoc contrasts ####
# model <- asv_models$model[[1]]
# emmeans(model, ~treatment)
#   -Time alone i.e. linear or quadratic
time_mainEffects_posthoc <- list('linear' = c(-1/2, -1/2, 0, 0, 0, 0, 1/4, 1/4, 1/4, 1/4),
                                 'quadratic' = c(1/2, 1/2, -2/4, -2/4, -2/4, -2/4, 1/4, 1/4, 1/4, 1/4),
                                 'early' = c(-1/2, -1/2, 1/4, 1/4, 1/4, 1/4, 0, 0, 0, 0),
                                 'late' = c(0, 0, -1/4, -1/4, -1/4, -1/4, 1/4, 1/4, 1/4, 1/4))

#   -Resistance alone i.e. difference between resistant & susceptabil across all time/exposure
resistance_mainEffects_posthoc <- list('R - S' = c(1/5, -1/5, 1/5, -1/5, 1/5, -1/5, 1/5, -1/5, 1/5, -1/5))

#   -Field vs Tank i.e. difference between field and average of healthy and diseased across all time/resistance
#   -Exposure alone i.e. difference between heatlhy and diseased across all resistance & T3/T7
exposure_mainEffects_posthoc <- list('Experimental - Field' = c(-1/2, -1/2, 1/8, 1/8, 1/8, 1/8, 1/8, 1/8, 1/8, 1/8),
                                     'D - H' = c(0, 0, 1/4, 1/4, -1/4, -1/4, 1/4, 1/4, -1/4, -1/4))
#   -Time x Resistance
#     -RvS (t3-t0), RvS (t7-t3) - name early/late resistance/susceptible
timeResistance_interaction_posthoc <- list('R(T3-T0) - S(T3-T0)' = c(-1, 0, 1/2, 0, 1/2, 0, 0, 0, 0, 0) - c(0, -1, 0, 1/2, 0, 1/2, 0, 0, 0, 0),
                                           'R(T7-T3) - S(T7-T3)' = c(0, 0, -1/2, 0, -1/2, 0, 1/2, 0, 1/2, 0) - c(0, 0, 0, -1/2, 0, -1/2, 0, 1/2, 0, 1/2))


#   -Time x Exposure
#     -HvD (t3-t0), HvD (t7-t3) - name early/late healthy/diseased
timeExposure_interaction_posthoc <- list('D(T3-T0)' = c(-1/2, -1/2, 1/2, 1/2, 0, 0, 0, 0, 0, 0),
                                         'H(T3-T0)' = c(-1/2, -1/2, 0, 0, 1/2, 1/2, 0, 0, 0, 0),
                                         'D(T3-T0) - H(T3-T0)' = c(-1/2, -1/2, 1/2, 1/2, 0, 0, 0, 0, 0, 0) - c(-1/2, -1/2, 0, 0, 1/2, 1/2, 0, 0, 0, 0),
                                         'D(T7-T3)' = c(0, 0, -1/2, -1/2, 0, 0, 1/2, 1/2, 0, 0),
                                         'H(T7-T3)' = c(0, 0, 0, 0, -1/2, -1/2, 0, 0, 1/2, 1/2),
                                         'D(T7-T3) - H(T7-T3)' = c(0, 0, -1/2, -1/2, 0, 0, 1/2, 1/2, 0, 0) - c(0, 0, 0, 0, -1/2, -1/2, 0, 0, 1/2, 1/2))

#   -Exposure x Resistance
#     -HSvDS, HSvDR, HRvDS, HRvDR
exposureResistance_interaction_posthoc <- list('HS - HR' = c(-1/3, 1/3, 0, 0, -1/3, 1/3, 0, 0, -1/3, 1/3),
                                               'HS - DS' = c(0, 0, 0, -1/2, 0, 1/2, 0, -1/2, 0, 1/2),
                                               'HS - DR' = c(-1/3, 1/3, -1/3, 0, 0, 1/3, -1/3, 0, 0, 1/3),
                                               'HR - DS' = c(1/3, -1/3, 0, -1/3, 1/3, 0, 0, -1/3, 1/3, 0),
                                               'HR - DR' = c(0, 0, -1/2, 0, 1/2, 0, -1/2, 0, 1/2, 0),
                                               'DS - DR' = c(-1/3, 1/3, -1/3, 1/3, 0, 0, -1/3, 1/3, 0, 0))

#   -Time x Resistance x Exposure
#       -DT7(S - R)
#       -ST7(D - H)

threeWay_interaction_posthoc <- list('DT7(S - R)' = c(0, 0, 0, 0, 0, 0, -1, 1, 0, 0),
                                     'ST7(D - H)' = c(0, 0, 0, 0, 0, 0, 0, 1, 0, -1))

posthoc_categories <- tibble(summarised_effect = c('time', 'resistance', 'exposure', 
                                                   'timeXresistance', 'timeXexposure', 'exposureXresistance',
                                                   'timeXresistanceXexposure'),
       contrasts = list(time_mainEffects_posthoc, resistance_mainEffects_posthoc, exposure_mainEffects_posthoc,
                     timeResistance_interaction_posthoc, timeExposure_interaction_posthoc, exposureResistance_interaction_posthoc,
                     threeWay_interaction_posthoc)) %>%
  unnest(contrasts) %>%
  mutate(contrast_name = names(contrasts)) 

#### Pathogen Pattern based Contrasts ####
# emmeans(asv_models$model[[1]], ~treatment)
posthoc_order <- c('T0.F.R', 'T0.F.S', 'T3.D.R', 'T3.D.S', 'T3.H.R', 'T3.H.S', 'T7.D.R', 'T7.D.S', 'T7.H.R', 'T7.H.S')
example_posthoc <- list('DS(T3-T0)' = c(0, -1, 0, 1, 0, 0, 0, 0, 0, 0),
                       'DS(T7-T3)' = c(0, 0, 0, -1, 0, 0, 0, 1, 0, 0),
                       'DS(T7-T0)' = c(0, -1, 0, 0, 0, 0, 0, 1, 0, 0),
                       'T0(DSvHR.DR.HS)' = c(-1, 1, 0, 0, 0, 0, 0, 0, 0, 0),
                       'T3(DSvHR.DR.HS)' = c(0, 0, -1/3, 1, -1/3, -1/3, 0, 0, 0, 0),
                       'T7(DSvHR.DR.HS)' = c(0, 0, 0, 0, 0, 0, -1/3, 1, -1/3, -1/3),
                       'T3(DSvDR)' = c(0, 0, -1, 1, 0, 0, 0, 0, 0, 0),
                       'T7(DSvDR)' = c(0, 0, 0, 0, 0, 0, -1, 1, 0, 0))

bacterial_growth <- list('(T3-T0)' = c(-1/2, -1/2, 1/4, 1/4, 1/4, 1/4, 0, 0, 0, 0),
                         '(T7-T3)' = c(0, 0, -1/4, -1/4, -1/4, -1/4, 1/4, 1/4, 1/4, 1/4),
                         '(T7-T0)' = c(-1/2, -1/2, 0, 0, 0, 0, 1/4, 1/4, 1/4, 1/4),
                         
                         'DS(T3-T0)' = c(0, -1, 0, 1, 0, 0, 0, 0, 0, 0),
                         'DS(T7-T3)' = c(0, 0, 0, -1, 0, 0, 0, 1, 0, 0),
                         'DS(T7-T0)' = c(0, -1, 0, 0, 0, 0, 0, 1, 0, 0),
                         
                         'DR(T3-T0)' = c(-1, 0, 1, 0, 0, 0, 0, 0, 0, 0),
                         'DR(T7-T3)' = c(0, 0, -1, 0, 0, 0, 1, 0, 0, 0),
                         'DR(T7-T0)' = c(-1, 0, 0, 0, 0, 0, 1, 0, 0, 0),
                         
                         'HS(T3-T0)' = c(0, -1, 0, 0, 0, 1, 0, 0, 0, 0),
                         'HS(T7-T3)' = c(0, 0, 0, 0, 0, -1, 0, 0, 0, 1),
                         'HS(T7-T0)' = c(0, -1, 0, 0, 0, 0, 0, 0, 0, 1),
                         
                         'HR(T3-T0)' = c(-1, 0, 0, 0, 1, 0, 0, 0, 0, 0),
                         'HR(T7-T3)' = c(0, 0, 0, 0, -1, 0, 0, 0, 1, 0),
                         'HR(T7-T0)' = c(-1, 0, 0, 0, 0, 0, 0, 0, 1, 0),
                         
                         'H(T3-T0)' = c(-1/2, -1/2, 0, 0, 1/2, 1/2, 0, 0, 0, 0),
                         'H(T7-T3)' = c(0, 0, 0, 0, -1/2, -1/2, 0, 0, 1/2, 1/2),
                         'H(T7-T0)' = c(-1/2, -1/2, 0, 0, 0, 0, 0, 0, 1/2, 1/2),
                         
                         'D(T3-T0)' = c(-1/2, -1/2, 1/2, 1/2, 0, 0, 0, 0, 0, 0),
                         'D(T7-T3)' = c(0, 0, -1/2, -1/2, 0, 0, 1/2, 1/2, 0, 0),
                         'D(T7-T0)' = c(-1/2, -1/2, 0, 0, 0, 0, 1/2, 1/2, 0, 0),
                         
                         'R(T3-T0)' = c(-1, 0, 1/2, 0, 1/2, 0, 0, 0, 0, 0),
                         'R(T7-T3)' = c(0, 0, -1/2, 0, -1/2, 0, 1/2, 0, 1/2, 0),
                         'R(T7-T0)' = c(-1, 0, 0, 0, -1/2, 0, 1/2, 0, 1/2, 0),
                         
                         'S(T3-T0)' = c(0, -1, 0, 1/2, 0, 1/2, 0, 0, 0, 0),
                         'S(T7-T3)' = c(0, 0, 0, -1/2, 0, -1/2, 0, 1/2, 0, 1/2),
                         'S(T7-T0)' = c(0, -1, 0, 0, 0, 0, 0, 1/2, 0, 1/2))

early_pathogen <- list('T3(DSvHR.DR.HS)' = c(0, 0, -1/3, 1, -1/3, -1/3, 0, 0, 0, 0),
                       'T3(DSvDR)' = c(0, 0, -1, 1, 0, 0, 0, 0, 0, 0))
continuous_pathogen <- list('T3(DSvHR.DR.HS)' = c(0, 0, -1/3, 1, -1/3, -1/3, 0, 0, 0, 0),
                            'T3(DSvDR)' = c(0, 0, -1, 1, 0, 0, 0, 0, 0, 0),
                            'T7(DSvHR.DR.HS)' = c(0, 0, 0, 0, 0, 0, -1/3, 1, -1/3, -1/3),
                            'T7(DSvDR)' = c(0, 0, 0, 0, 0, 0, -1, 1, 0, 0))
late_pathogen <- list('T7(DSvHR.DR.HS)' = c(0, 0, 0, 0, 0, 0, -1/3, 1, -1/3, -1/3),
                      'T7(DSvDR)' = c(0, 0, 0, 0, 0, 0, -1, 1, 0, 0))

early_opportunist <- list('T3(DS.DRvHS.HR)' = c(0, 0, 1/2, 1/2, -1/2, -1/2, 0, 0, 0, 0))
continuous_opportunist <- list('T3(DS.DRvHS.HR)' = c(0, 0, 1/2, 1/2, -1/2, -1/2, 0, 0, 0, 0),
                               'T7(DS.DRvHS.HR)' = c(0, 0, 0, 0, 0, 0, 1/2, 1/2, -1/2, -1/2))
late_opportunist <- list('T7(DS.DRvHS.HR)' = c(0, 0, 0, 0, 0, 0, 1/2, 1/2, -1/2, -1/2))


probiotic <- list('T0(SvR)' = c(-1, 1, 0, 0, 0, 0, 0, 0, 0, 0),
                  'T3(DSvDR)' = c(0, 0, -1, 1, 0, 0, 0, 0, 0, 0),
                  'T7(DSvDR)' = c(0, 0, 0, 0, 0, 0, -1, 1, 0, 0))

posthoc_categories <- tibble(microbial_signature = c('growth_comparisons',
                                                     'early_pathogen', 'continuous_pathogen', 'late_pathogen',
                                                     'early_opportunist', 'continuous_opportunist', 'late_opportunist',
                                                     'probiotic'),
                             contrasts = list(bacterial_growth,
                                              early_pathogen, continuous_pathogen, late_pathogen,
                                              early_opportunist, continuous_opportunist, late_opportunist,
                                              probiotic)) %>%
  unnest(contrasts) %>%
  mutate(contrast_name = names(contrasts)) %>%
  group_by(contrast_name, contrasts) %>%
  summarise(signatures = list(c(microbial_signature)),
            .groups = 'drop')

#### Model ASV counts ####
if(file.exists('../intermediate_files/mixed_model_results.rds.gz') & !refit_models){
  asv_models <- write_rds('../intermediate_files/mixed_model_results.rds.gz')
} else {
  cluster_copy(cluster, c('posthoc_categories'))
  
  asv_models <- normalized_asv_counts %>%
    nest_by(across(c('asv_id', Kingdom:Species))) %>%
    partition(cluster) %>%
    mutate(fit_model(log2_cpm ~ treatment + (1 | genotype) + (1 | tank),
                     data, 
                     use_weights = FALSE),
           random_anova = list(rand(model)),
           process_model(model, re_model, random_anova),
           posthoc = list(emmeans(model, ~treatment) %>%
                            contrast(method = posthoc_categories$contrasts, adjust = 'none'))) %>%
    collect() %>%
    select(-re_model, -ends_with('global')) %>%
    ungroup %>% 
    p_adjust %>%
    relocate(anova_table, .after = model) %>% 
    relocate(posthoc, .after = random_anova)
  
  write_rds(asv_models, '../intermediate_files/mixed_model_results.rds.gz')
  write_csv(select(asv_models, -where(is.list)), '../intermediate_files/mixed_model_results.csv.gz')
}

#### Main Effects Upset ####
asv_models %>%
  select(asv_id, starts_with('fdr')) %>% 
  mutate(across(starts_with('fdr'), ~. < 0.05)) %>%
  
  pivot_longer(cols = -asv_id,
               names_to = c('term'),
               values_to = 'significance',
               names_prefix = 'fdr_') %>%
  filter(significance) %>%
  group_by(asv_id) %>%
  summarise(terms = list(term),
            .groups = 'drop') %>%
  
  ggplot(aes(x = terms)) +
  geom_bar() +
  scale_x_upset() +
  theme_classic() +
  theme_combmatrix(combmatrix.label.make_space = TRUE)

#### Post-Hoc upset only on ASVs that have a significant treatment effect ####
significant_models <- asv_models %>%
  filter(fdr_treatment < 0.05) %>%
  # slice(26) %>%
  rowwise %>%
  mutate(process_postHoc(posthoc)) %>%
  ungroup() %>%
  p_adjust(exclude_cols = c('treatment', 'tank', 'genotype'))

bacterial_signature_asv <- significant_models %>%
  select(asv_id, starts_with('fdr')) %>% 
  select(-contains(c('treatment', 'tank', 'genotype'))) %>%
  mutate(across(starts_with('fdr'), ~. < 0.05)) %>%
  
  pivot_longer(cols = -asv_id,
               names_to = c('term'),
               values_to = 'significance') %>%
  mutate(term = str_remove(term, 'fdr_')) %>%
  
  #Select ASVs which fit all characteristics of any given bacterial signature
  left_join(posthoc_categories, by = c('term' = 'contrast_name')) %>%
  filter(signatures != 'growth_comparisons') %>%
  select(-contrasts) %>%
  unnest(signatures) %>%
  # mutate(signatures = if_else(signatures == 'growth_comparisons', term, signatures)) %>%
  # 
  # filter(signatures != 'growth_comparisons' | 
  #          (signatures == 'growth_comparisons' & significance)) %>%
  
  group_by(asv_id, signatures) %>%
  filter(all(significance)) %>%
  ungroup %>%
  select(-term) %>%
  distinct 

bacterial_signature_asv %>%
  group_by(asv_id) %>%
  summarise(terms = list(unique(signatures)),
            .groups = 'drop') %>%
  
  ggplot(aes(x = terms)) +
  geom_bar() +
  scale_x_upset() +
  theme_classic() +
  theme_combmatrix(combmatrix.label.make_space = TRUE)


# emily comp upset

testt <- bacterial_signature_asv %>% pivot_wider(names_from = signatures, values_from = significance)
testt[is.na(testt)] <- FALSE

testt <- testt %>% left_join(taxonomy_tibble, by = c('asv_id' = 'asv_names'))

upset(testt,
      colnames(select(testt, starts_with('late') | starts_with('early') | starts_with('cont'))), 
      
      base_annotations=list(
        'Intersection size'=intersection_size(
          mapping=aes(fill=Family)
        ) 
      ),
      queries=list(upset_query(set="early_opportunist", color="deepskyblue", fill = "deepskyblue"),
                   upset_query(set="late_opportunist", color="deepskyblue", fill = "deepskyblue"),
                   upset_query(set="continuous_opportunist", color="deepskyblue", fill = "deepskyblue"),
                   upset_query(set="early_pathogen", color="firebrick1", fill = "firebrick1"),
                   upset_query(set="late_pathogen", color="firebrick1", fill = "firebrick1")),
      
      name='asv_names', width_ratio=0.1, min_size = 0) +
  ggtitle("Disease-Associated Bacteria")




#### Plot Individual Groupings ####
the_plots <- bacterial_signature_asv %>%
  arrange(signatures) %>%
  group_by(asv_id) %>%
  summarise(signatures = str_c(signatures, collapse = ', ')) %>%
  inner_join(significant_models,
             by = 'asv_id') %>%
  rowwise %>%
  mutate(plot = list(emmeans(model, ~treatment) %>%
                       broom::tidy(conf.int = TRUE) %>%
                       separate(treatment, into = c('time', 'exposure', 'susceptability')) %>%
                       ggplot(aes(x = time, y = estimate, ymin = conf.low, ymax = conf.high, 
                                  colour = exposure,
                                  shape = susceptability)) +
                       geom_pointrange(position = position_dodge(0.5)) +
                       labs(title = str_c(Family, Genus, asv_id, sep = ': ')))) %>%
  group_by(signatures) %>%
  summarise(plots = list(wrap_plots(plot) + plot_layout(guides = 'collect') & plot_annotation(title = signatures)))

the_plots$plots[[1]]
the_plots$plots[[2]]
the_plots$plots[[3]]
the_plots$plots[[4]]
the_plots$plots[[5]]
the_plots$plots[[6]]

#### Work Zone ####
tmp <- filter(significant_models, 
              `fdr_DS(T7-T3)` < 0.05,
              `fdr_T7(DSvHR.DR.HS)` < 0.05,
              `fdr_T7(DSvDR)` < 0.05,
              `fdr_DS(T7-T0)` < 0.05) %>%
  filter(`estimate_DS(T7-T3)` > 0,
         `estimate_DS(T7-T0)` > 0,
         `estimate_T7(DSvDR)` > 0) %>%
  mutate(across(Kingdom:Species, str_replace_na)) 

tmp %>%
  rowwise %>%
  mutate(plot = list(emmeans(model, ~treatment) %>%
                       broom::tidy(conf.int = TRUE) %>%
                       separate(treatment, into = c('time', 'exposure', 'susceptability')) %>%
                       ggplot(aes(x = time, y = estimate, ymin = conf.low, ymax = conf.high, 
                                  colour = exposure,
                                  shape = susceptability)) +
                       geom_pointrange(position = position_dodge(0.5)) +
                       labs(title = str_c(Family, Genus, asv_id, sep = ': ')))) %>%
  pull(plot) %>%
  wrap_plots() +
  plot_layout(guides = 'collect')


filter(significant_models, 
       `fdr_DS(T7-T3)` < 0.05,
       `fdr_T7(DSvHR.DR.HS)` < 0.05,
       `fdr_T7(DSvDR)` < 0.05,
       # `fdr_DS(T7-T0)` < 0.05
       ) %>%
  filter(#`estimate_DS(T7-T3)` > 0,
         #`estimate_DS(T7-T0)` > 0,
         `estimate_T7(DSvDR)` < 0) %>%
  mutate(across(Kingdom:Species, str_replace_na)) %>%
  rowwise %>%
  mutate(plot = list(emmeans(model, ~treatment) %>%
                       broom::tidy(conf.int = TRUE) %>%
                       separate(treatment, into = c('time', 'exposure', 'susceptability')) %>%
                       ggplot(aes(x = time, y = estimate, ymin = conf.low, ymax = conf.high, 
                                  colour = exposure,
                                  shape = susceptability)) +
                       geom_pointrange(position = position_dodge(0.5)) +
                       labs(title = str_c(Family, Genus, asv_id, sep = ': ')))) %>%
  pull(plot) %>%
  wrap_plots() +
  plot_layout(guides = 'collect')

T3(DS.DRvHS.HR)
T7(DS.DRvHS.HR)
filter(significant_models, 
       # `fdr_DS(T7-T3)` < 0.05,
       `fdr_T3(DSvHR.DR.HS)` < 0.05,
       # `fdr_T3(DSvDR)` < 0.05,
       # `fdr_DS(T7-T0)` < 0.05,
       # `fdr_DS(T3-T0)` < 0.05,
       `fdr_T7(DSvHR.DR.HS)` < 0.05,
       `fdr_T7(DSvDR)` < 0.05,
) %>%
  filter(#`estimate_DS(T7-T3)` > 0,
    #`estimate_DS(T7-T0)` > 0,
    `estimate_T7(DSvDR)` > 0) %>%
  mutate(across(Kingdom:Species, str_replace_na)) %>%
  rowwise %>%
  mutate(plot = list(emmeans(model, ~treatment) %>%
                       broom::tidy(conf.int = TRUE) %>%
                       separate(treatment, into = c('time', 'exposure', 'susceptability')) %>%
                       ggplot(aes(x = time, y = estimate, ymin = conf.low, ymax = conf.high, 
                                  colour = exposure,
                                  shape = susceptability)) +
                       geom_pointrange(position = position_dodge(0.5)) +
                       labs(title = str_c(Family, Genus, asv_id, sep = ': ')))) %>%
  pull(plot) %>%
  wrap_plots() +
  plot_layout(guides = 'collect')




tmp <- significant_models %>%
  select(asv_id, starts_with('fdr')) %>% 
  select(-contains(c('treatment', 'tank', 'genotype'))) %>%
  mutate(across(starts_with('fdr'), ~. < 0.05)) %>%
  
  pivot_longer(cols = -asv_id,
               names_to = c('term'),
               values_to = 'significance') %>%
  mutate(term = str_remove(term, 'fdr_')) %>%
  left_join(posthoc_categories, by = c('term' = 'contrast_name')) %>%
  mutate(summarised_effect = if_else(is.na(summarised_effect), term, summarised_effect)) %>%
  group_by(asv_id, summarised_effect) %>%
  filter(any(significance)) %>%
  filter(summarised_effect == 'exposureXresistance') %>%
  ungroup %>%
  select(-contrasts) %>%
  group_by(asv_id, summarised_effect) %>%
  filter(significance) %>%
  summarise(sig_terms = str_c(term, collapse = ' & '),
            .groups = 'drop') %>%
  left_join(significant_models, by = 'asv_id') %>%
  mutate(across(Kingdom:Species, str_replace_na)) %>%
  rowwise %>%
  mutate(plot = list(emmeans(model, ~treatment) %>%
                       broom::tidy(conf.int = TRUE) %>%
                       separate(treatment, into = c('time', 'exposure', 'susceptability')) %>%
                       ggplot(aes(x = time, y = estimate, ymin = conf.low, ymax = conf.high, 
                                  colour = exposure,
                                  shape = susceptability)) +
                       geom_pointrange(position = position_dodge(0.5)) +
                       labs(title = str_c(Family, Genus, asv_id, sep = ': '),
                            subtitle = sig_terms)))



  
wrap_plots(tmp$plot) + plot_layout(guides = 'collect')



tmp <- significant_models %>%
  select(asv_id, starts_with('fdr')) %>% 
  select(-contains(c('treatment', 'tank', 'genotype'))) %>%
  mutate(across(starts_with('fdr'), ~. < 0.05)) %>%
  
  pivot_longer(cols = -asv_id,
               names_to = c('term'),
               values_to = 'significance') %>%
  mutate(term = str_remove(term, 'fdr_')) %>%
  left_join(posthoc_categories, by = c('term' = 'contrast_name')) %>%
  mutate(summarised_effect = if_else(is.na(summarised_effect), term, summarised_effect)) %>%
  group_by(asv_id, summarised_effect) %>%
  filter(any(significance)) %>%
  ungroup %>%
  group_by(asv_id) %>%
  filter(any(summarised_effect == 'timeXresistanceXexposure')) %>%
  ungroup %>%
  select(-contrasts) %>%
  group_by(asv_id) %>%
  filter(significance) %>%
  summarise(sig_terms = str_c(term, collapse = ' & '),
            .groups = 'drop') %>%
  left_join(significant_models, by = 'asv_id') %>%
  mutate(across(Kingdom:Species, str_replace_na)) %>%
  filter(`estimate_DT7(S - R)` > 0 & `estimate_ST7(D - H)` > 0) %>%
  rowwise %>%
  mutate(plot = list(emmeans(model, ~treatment) %>%
                       broom::tidy(conf.int = TRUE) %>%
                       separate(treatment, into = c('time', 'exposure', 'susceptability')) %>%
                       ggplot(aes(x = time, y = estimate, ymin = conf.low, ymax = conf.high, 
                                  colour = exposure,
                                  shape = susceptability)) +
                       geom_pointrange(position = position_dodge(0.5)) +
                       labs(title = str_c(Family, Genus, asv_id, sep = ': '),
                            subtitle = sig_terms)))

wrap_plots(tmp$plot) + 
  plot_layout(guides = 'collect')



tmp$model[[2]] %>%
  emmeans(~treatment) %>%
  contrast(method = c(exposure_mainEffects_posthoc))


reference_grid_2 <- tmp$model[[2]] %>%
  ref_grid(~treatment)

tmp1_grid <- add_grouping(reference_grid_2, newname = 'Experiment_Field', refname = 'treatment',
                          c('F', 'F', 'E', 'E', 'E', 'E', 'E', 'E', 'E', 'E'))
emmeans(tmp1_grid, ~Experiment_Field) %>%
  contrast('pairwise')

emmeans(tmp1_grid, ~Experiment_Field) %>%
  broom::tidy(conf.int = TRUE) %>%
  ggplot(aes(x = Experiment_Field, y = estimate, ymin = conf.low, ymax = conf.high)) +
  geom_pointrange()



tmp1_grid <- add_grouping(reference_grid_2, newname = 'H_D', refname = 'treatment',
                          c(NA, NA, 'D', 'D', 'H', 'H', 'D', 'D', 'H', 'H'))
emmeans(tmp1_grid, ~H_D) %>%
  contrast('pairwise')

emmeans(tmp1_grid, ~H_D) %>%
  broom::tidy(conf.int = TRUE) %>%
  ggplot(aes(x = H_D, y = estimate, ymin = conf.low, ymax = conf.high)) +
  geom_pointrange()



tmp <- significant_models %>%
  select(-contains(c('treatment', 'tank', 'genotype'))) %>%
  mutate(across(starts_with('fdr'), ~. < 0.05)) %>% 
  # select(asv_id, starts_with('fdr')) %>%
  filter(`fdr_D(T3-T0)` & `fdr_D(T3-T0) - H(T3-T0)`) 

reference_grid_2 <- tmp$model[[1]] %>%
  ref_grid(~treatment)
tmp1_grid <- add_grouping(reference_grid_2, newname = 'time_disease', refname = 'treatment',
                          c('T0_F', 'T0_F', 'T3_D', 'T3_D', 'T3_H', 
                            'T3_H', 'T7_D', 'T7_D', 'T7_H', 'T7_H'))
emmeans(tmp1_grid, ~time_disease) %>%
  contrast('pairwise')

emmeans(tmp1_grid, ~time_disease) %>%
  broom::tidy(conf.int = TRUE) %>%
  separate(time_disease, into = c('time', 'exposure')) %>%
  ggplot(aes(x = time, y = estimate, ymin = conf.low, ymax = conf.high,
             colour = exposure)) +
  geom_pointrange()
