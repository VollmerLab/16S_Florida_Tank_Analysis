##TODO - Set up/evaluate contrasts for relevant/interesting bacterial responses
#           -want to identify pathogens, opportunists, probiotics
##TODO - 

#### Libraries ####
library(magrittr)
library(lmerTest)
library(emmeans)
library(multidplyr)
library(ggupset)
library(qvalue)
library(patchwork)
library(relayer) #devtools::install_github("clauswilke/relayer")
library(ComplexUpset)
library(corrplot)
library(Hmisc)
library(broom.mixed)
library(tidyverse)

refit_models <- FALSE

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

run_posthoc <- function(model, contrast_list){
  em_out <- emmeans(model, ~treatment)
  
  contrast_list %>%
    rowwise(direction) %>%
    reframe(emmeans::contrast(em_out,
                              method = contrast$contrasts, 
                              adjust = 'none',
                              side = direction) %>%
              as_tibble)
}

# posthoc <- asv_models$posthoc[[1]]
process_postHoc <- function(posthoc){
  post_row <- as_tibble(posthoc) %>%
    dplyr::rename(tvalue = t.ratio,
           pvalue = p.value) %>%
    mutate(contrast = str_c(contrast, direction, sep = '_'), .keep = 'unused') %>%
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
normalized_asv_counts <- full_data %>% #read_csv('../intermediate_files/fully_preprocessed_samples.csv.gz', show_col_types = FALSE) %>%
  mutate(time = factor(time, ordered = TRUE)) %>%
  mutate(treatment = str_c(time, exposure, susceptability, sep = '_'),
         time_exposure = str_c(time, exposure, sep = '_'),
         timeC = str_extract(time, '[0-9]+') %>% as.numeric,
         across(Kingdom:Species, str_replace_na)) %>%
  mutate(asv_number = str_extract(asv_id, '[0-9]+') %>% as.integer) #%>%
  #filter(asv_id %in% otus_to_analyze) #otus in D, T3, T7
  # filter(asv_number <= 200)

taxonomy_tibble <- tax_table(microbiome_data) %>% 
  as.data.frame %>%
  as_tibble(rownames = "asv_names")

homogenate_data <- read_csv('../intermediate_files/homogenate_cpm.csv', show_col_types = FALSE)


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

#pathogens
early_pathogen <- list('T3(DSvHR.DR.HS)' = c(0, 0, -1/3, 1, -1/3, -1/3, 0, 0, 0, 0),
                       'T3(DSvDR)' = c(0, 0, -1, 1, 0, 0, 0, 0, 0, 0),
                       'DS(T3-T0)' = c(0, -1, 0, 1, 0, 0, 0, 0, 0, 0))
continuous_pathogen <- list('T3(DSvHR.DR.HS)' = c(0, 0, -1/3, 1, -1/3, -1/3, 0, 0, 0, 0),
                            'T3(DSvDR)' = c(0, 0, -1, 1, 0, 0, 0, 0, 0, 0),
                            'DS(T3-T0)' = c(0, -1, 0, 1, 0, 0, 0, 0, 0, 0),
                            'T7(DSvHR.DR.HS)' = c(0, 0, 0, 0, 0, 0, -1/3, 1, -1/3, -1/3),
                            'T7(DSvDR)' = c(0, 0, 0, 0, 0, 0, -1, 1, 0, 0),
                            'DS(T7-T0)' = c(0, -1, 0, 0, 0, 0, 0, 1, 0, 0))
late_pathogen <- list('T7(DSvHR.DR.HS)' = c(0, 0, 0, 0, 0, 0, -1/3, 1, -1/3, -1/3),
                      'T7(DSvDR)' = c(0, 0, 0, 0, 0, 0, -1, 1, 0, 0),
                      'DS(T7-T0)' = c(0, -1, 0, 0, 0, 0, 0, 1, 0, 0))


# early_pathogen <- list('T3(DSvHR.DR.HS)' = c(0, 0, -1/3, 1, -1/3, -1/3, 0, 0, 0, 0),
#                       'T3(DSvDR)' = c(0, 0, -1, 1, 0, 0, 0, 0, 0, 0),
#                       'linear' = c(-1/2, -1/2, 0, 0, 0, 0, 0, 1, 0, 0))
# continuous_pathogen <- list('T3(DSvHR.DR.HS)' = c(0, 0, -1/3, 1, -1/3, -1/3, 0, 0, 0, 0),
#                             'T3(DSvDR)' = c(0, 0, -1, 1, 0, 0, 0, 0, 0, 0),
#                             'T7(DSvHR.DR.HS)' = c(0, 0, 0, 0, 0, 0, -1/3, 1, -1/3, -1/3),
#                             'T7(DSvDR)' = c(0, 0, 0, 0, 0, 0, -1, 1, 0, 0),
#                             'linear' = c(-1/2, -1/2, 0, 0, 0, 0, 0, 1, 0, 0))
# late_pathogen <- list('T7(DSvHR.DR.HS)' = c(0, 0, 0, 0, 0, 0, -1/3, 1, -1/3, -1/3),
#                       'T7(DSvDR)' = c(0, 0, 0, 0, 0, 0, -1, 1, 0, 0),
#                       'linear' = c(-1/2, -1/2, 0, 0, 0, 0, 0, 1, 0, 0))

#opportunists
early_opportunist <- list('T3(DS.DRvHS.HR)' = c(0, 0, 1/2, 1/2, -1/2, -1/2, 0, 0, 0, 0),
                          'DS.DR(T3-T0)' = c(-1/2, -1/2, 1/2, 1/2, 0, 0, 0, 0, 0, 0))
continuous_opportunist <- list('T3(DS.DRvHS.HR)' = c(0, 0, 1/2, 1/2, -1/2, -1/2, 0, 0, 0, 0),
                               'T7(DS.DRvHS.HR)' = c(0, 0, 0, 0, 0, 0, 1/2, 1/2, -1/2, -1/2),
                               'DS.DR(T3-T0)' = c(-1/2, -1/2, 1/2, 1/2, 0, 0, 0, 0, 0, 0),
                               'DS.DR(T7-T0)' = c(-1/2, -1/2, 0, 0, 0, 0, 1/2, 1/2, 0, 0))
late_opportunist <- list('T7(DS.DRvHS.HR)' = c(0, 0, 0, 0, 0, 0, 1/2, 1/2, -1/2, -1/2),
                          'DS.DR(T7-T0)' = c(-1/2, -1/2, 0, 0, 0, 0, 1/2, 1/2, 0, 0))

#add in left tests
probiotic_t0 <- list('T0(SvR)' = c(-1, 1, 0, 0, 0, 0, 0, 0, 0, 0))
probiotic_t3 <- list('T3(DSvDR)' = c(0, 0, -1, 1, 0, 0, 0, 0, 0, 0))
probiotic_t7 <- list('T7(DSvDR)' = c(0, 0, 0, 0, 0, 0, -1, 1, 0, 0))

probiotic_t3_strict <- list('T3(DSvDR)' = c(0, 0, -1, 1, 0, 0, 0, 0, 0, 0),
                            'T3(DRvHR.DS.HS)' = c(0, 0, -1, 1/3, 1/3, 1/3, 0, 0, 0, 0),
                            'DR(T3-T0)' = c(0, 1, -1, 0, 0, 0, 0, 0, 0, 0))
probiotic_t7_strict <- list('T7(DSvDR)' = c(0, 0, 0, 0, 0, 0, -1, 1, 0, 0),
                            'T7(DRvHR.DS.HS)' = c(0, 0, 0, 0, 0, 0, -1, 1/3, 1/3, 1/3),
                            'DR(T7-T3)' = c(0, 0, 1, 0, 0, 0, -1, 0, 0, 0),
                            'DR(T7-T0)' = c(1, 0, 0, 0, 0, 0, -1, 0, 0, 0))

rev_probiotic_t0 <- list('T0(SvR)' = c(1, -1, 0, 0, 0, 0, 0, 0, 0, 0))

crasher_t3 <- list('DS(T3-T0)' = c(0, -1, 0, 1, 0, 0, 0, 0, 0, 0))
crasher_t7 <- list('DS(T7-T0)' = c(0, -1, 0, 0, 0, 0, 0, 1, 0, 0))

crasher_t3_strict <- list('DS(T3-T0)' = c(0, -1, 0, 1, 0, 0, 0, 0, 0, 0),
                          'T3(DSvHR.DR.HS)' = c(0, 0, -1/3, 1, -1/3, -1/3, 0, 0, 0, 0))
crasher_t7_strict <- list('DS(T7-T0)' = c(0, -1, 0, 0, 0, 0, 0, 1, 0, 0),
                          'T7(DSvHR.DR.HS)' = c(0, 0, 0, 0, 0, 0, -1/3, 1, -1/3, -1/3))

#bad commensalist def
commensalist <- list('T3vT7' = c(0, 0, 1/4, 1/4, 1/4, 1/4, -1/4, -1/4, -1/4, -1/4))

#Put tests for one/two sided and directionality into bins. Will combine after post-hoc into meaningful categorization
two_sided_tests <- tibble(microbial_signature = c('growth_comparisons'),
                          contrasts = list(bacterial_growth),
                          direction = '=') #=
right_tests <- tibble(microbial_signature = c('early_pathogen', 'continuous_pathogen', 'late_pathogen',
                                              'early_opportunist', 'continuous_opportunist', 'late_opportunist'),
                      contrasts = list(early_pathogen, continuous_pathogen, late_pathogen,
                                       early_opportunist, continuous_opportunist, late_opportunist),
                      direction = '>') #>
left_tests <- tibble(microbial_signature = c('crasher_t3', 'crasher_t7', 'probiotic_t3_strict', 'probiotic_t7_strict', 'rev_probiotic_t0'),
                     contrasts = list( crasher_t3, crasher_t7, probiotic_t3_strict, probiotic_t7_strict, rev_probiotic_t0),
                     direction = '<') #<

posthoc_categories <- bind_rows(two_sided_tests, right_tests, left_tests) %>%
  unnest(contrasts) %>%
  mutate(contrast_name = names(contrasts)) %>%
  group_by(contrast_name, contrasts, direction) %>%
  summarise(signatures = list(c(microbial_signature)),
            .groups = 'drop') %>%
  nest(contrast = -direction)


#### Model ASV counts ####

if(file.exists('../intermediate_files/mixed_model_results.rds.gz') & !refit_models){
  asv_models <- read_rds('../intermediate_files/mixed_model_results.rds.gz')
} else {
  cluster_copy(cluster, c('posthoc_categories', 'run_posthoc'))
  
  asv_models <- normalized_asv_counts %>%
    nest_by(across(c('asv_id', Kingdom:Species))) %>%
    partition(cluster) %>%
    mutate(fit_model(log2_cpm ~ treatment + (1 | genotype) + (1 | tank),
                     data, 
                     use_weights = FALSE),
           random_anova = list(rand(model)),
           process_model(model, re_model, random_anova),
           posthoc = list(run_posthoc(model, posthoc_categories))) %>%
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
  p_adjust(exclude_cols = c('treatment', 'tank', 'genotype')) #%>% #affects otu filter


bacterial_signature_asv <- significant_models %>%
  select(asv_id, starts_with('fdr')) %>% 
  select(-contains(c('treatment', 'tank', 'genotype'))) %>%
  mutate(across(starts_with('fdr'), ~. < 0.05)) %>%
  
  pivot_longer(cols = -asv_id,
               names_to = c('term'),
               values_to = 'significance') %>%
  mutate(term = str_remove(term, 'fdr_')) %>%
  
  #Select ASVs which fit all characteristics of any given bacterial signature
  mutate(direction = str_extract(term, '[><=]'),
         term = str_remove(term, '_[><=]')) %>%
  
  left_join(unnest(posthoc_categories, contrast), 
            by = c('term' = 'contrast_name', 'direction')) %>%
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


#### Complex Upset for Bacterial Strategies ####

#prep for comp upset
comp_upset_bac_strat <- bacterial_signature_asv %>% 
  pivot_wider(names_from = signatures, values_from = significance) %>%
  mutate(across(-asv_id, ~ifelse(is.na(.), FALSE, .))) %>% 
  left_join(taxonomy_tibble, by = c('asv_id' = 'asv_names'))


### Upset - bacterial strategies and ASV Presence Data

#  colnames(simple_comp_upset %>% select(-c(asv_id, colnames(taxonomy_tibble %>% select(-asv_names)))))

upset(
  comp_upset_bac_strat %>% left_join(venn_all_times_and_doses, by = c("asv_id" = "OTU")),
  colnames(comp_upset_bac_strat %>% left_join(venn_all_times_and_doses, by = c("asv_id" = "OTU")) %>%
             select(T0, D, H, probiotic_t7_strict, crasher_t7, crasher_t3, late_opportunist, early_opportunist, late_pathogen, early_pathogen)), 

  base_annotations=list(
    'Intersection size'=intersection_size(counts=T, text = aes(size = 6.5),
                                          bar_number_threshold = 25,
        mapping=aes(fill=Family, col = Family, label = parse_number(asv_id)), col = "gray10"
    ) +
      geom_text(size = 4, position = position_stack(vjust = 0.5), col = "gray10")
  ),
  matrix=(
    intersection_matrix(geom=geom_point(shape = "circle filled", size=3, stroke = 0.35, color = "gray20"))
    + scale_color_manual(
      values=c(
        #"T7" = "#650197",
        #"T3" = "#B21BFF",
        "T0" = "#B21BFF",
        "D" = "#AE0404",
        "H" = "#0FAB02",
        "late_pathogen" = "#FF1E1E",
        "early_pathogen" = "#FF9797",
        #"continuous_opportunist" = "royalblue4",
        "late_opportunist" = "deepskyblue3",
        "early_opportunist" = "lightskyblue",
        "crasher_t3" = "#FF9A00",
        "crasher_t7" = "#BD5E03",
        "probiotic_t7_strict" = "#49EACC"
      )
    )
  ),
  queries=list(
    #upset_query(set = "T7", fill = "#650197"),
    #upset_query(set = "T3", fill = "#B21BFF"), "#D98EFF"
    upset_query(set = "T0", fill = "#B21BFF"),
    upset_query(set = "D", fill = "#AE0404"),
    upset_query(set = "H", fill = "#0FAB02"),
    upset_query(set = "late_pathogen", fill = "#FF1E1E"),
    upset_query(set = "early_pathogen", fill = "#FF9797"),
    #upset_query(set = "continuous_opportunist", fill = "royalblue4"),
    upset_query(set = "late_opportunist", fill = "deepskyblue3"),
    upset_query(set = "early_opportunist", fill = "lightskyblue"),
    upset_query(set = "crasher_t3", fill = "#FF9A00"),
    upset_query(set = "crasher_t7", fill = "#BD5E03"),
    upset_query(set = "probiotic_t7_strict", fill = "#49EACC")
  ),
  name='ASVs', width_ratio=0.1, min_size = 0, sort_sets = FALSE,
  stripes = c("#F8EAFF", "#FFE4E4", "#E7FFE5", rep(c("gray91", "gray97"), 4))) +
  #stripes = c(rep(c("gray78", "gray87"), 2), "gray78", rep(c("gray97", "gray91"), 4))) +
  ggtitle("Bacterial Strategies/Presence Data") 


### Upset - just bacterial strategies

upset(
  comp_upset_bac_strat %>% left_join(venn_all_times_and_doses, by = c("asv_id" = "OTU")),
  colnames(comp_upset_bac_strat %>% left_join(venn_all_times_and_doses, by = c("asv_id" = "OTU")) %>%
             select(probiotic_t7_strict, crasher_t7, crasher_t3, late_opportunist, early_opportunist, late_pathogen, early_pathogen)), 
  
  base_annotations=list(
    'Intersection size'=intersection_size(counts=T, text = aes(size = 6.5),
                                          bar_number_threshold = 25,
                                          mapping=aes(fill=Family, col = Family, label = parse_number(asv_id)), col = "gray10"
    ) +
      geom_text(size = 4, position = position_stack(vjust = 0.5), col = "gray10")
  ),
  matrix=(
    intersection_matrix(geom=geom_point(shape = "circle filled", size=3, stroke = 0.35, color = "gray20"))
    + scale_color_manual(
      values=c(
        "late_pathogen" = "#FF1E1E",
        "early_pathogen" = "#FF9797",
        #"continuous_opportunist" = "royalblue4",
        "late_opportunist" = "deepskyblue3",
        "early_opportunist" = "lightskyblue",
        "crasher_t3" = "#FF9A00",
        "crasher_t7" = "#BD5E03",
        "probiotic_t7_strict" = "#49EACC"
      )
    )
  ),
  queries=list(
    upset_query(set = "late_pathogen", fill = "#FF1E1E"),
    upset_query(set = "early_pathogen", fill = "#FF9797"),
    #upset_query(set = "continuous_opportunist", fill = "royalblue4"),
    upset_query(set = "late_opportunist", fill = "deepskyblue3"),
    upset_query(set = "early_opportunist", fill = "lightskyblue"),
    upset_query(set = "crasher_t3", fill = "#FF9A00"),
    upset_query(set = "crasher_t7", fill = "#BD5E03"),
    upset_query(set = "probiotic_t7_strict", fill = "#49EACC")
  ),
  name='ASVs', width_ratio=0.1, min_size = 0, sort_sets = FALSE,
  stripes = c(rep(c("gray91", "gray97"), 4))) +
  #stripes = c(rep(c("gray78", "gray87"), 2), "gray78", rep(c("gray97", "gray91"), 4))) +
  ggtitle("Bacterial Strategies") 



#### Plot Bacterial Strategy Groupings ####

## Prep Data for Homogenate Dose Panel

homogenate_models <- homogenate_data %>%
  nest_by(asv_id) %>%
  mutate(em_model = list(emmeans(lm(log2_cpm ~ exposure, data = data),
                                 ~exposure))) %>%
  mutate(dose_pairwise = em_model %>%
           contrast('pairwise', adjust = 'fdr') %>%
           broom::tidy(conf.int = TRUE) %>% pull(p.value)) %>%
  mutate(homogenate_pred = list(em_model %>%
                                  broom::tidy(conf.int = TRUE) %>% mutate(time = "Dose", 
                                                                          graph_cat = "dose", c_time = ifelse(exposure == "D", -1.5, -1.2), susceptability = "na",
                                                                          facet_lab = "Doses") %>% 
                                  {. ->> set_one } %>% #save data to this var
                                  mutate(std.error = NA, df = NA, conf.low = NA, conf.high = NA, statistic = NA, p.value = dose_pairwise,
                                         c_time = c(-1.8, -0.9), susceptability = NA, exposure = NA) %>% #dummy set to change size of Dose Facet
                                  rbind(set_one) %>%
                                  {. ->> set_two } %>%
                                  arrange(desc(estimate)) %>%
                                  dplyr::slice(1) %>%
                                  mutate(exposure = "p_val", c_time = -1.35, estimate = estimate + 3) %>%
                                  rbind(set_two)
                                
  )) #recombine them

## Make Fancy Plots

the_plots <- bacterial_signature_asv %>%
  arrange(signatures) %>%
  group_by(asv_id) %>%
  summarise(signatures = str_c(signatures, collapse = ', ')) %>%
  mutate(signatures = case_when(signatures == "continuous_pathogen, early_pathogen, late_pathogen" ~ "continuous_pathogen",
                                signatures == "crasher_t3, crasher_t7" ~ "continuous_crasher",
                                signatures == "early_opportunist, early_pathogen" ~ "early_pathogen",
                                signatures == "late_opportunist, late_pathogen" ~ "late_pathogen",
                                signatures == "late_opportunist, probiotic_t7_strict" ~ "late_probiotic",
                                signatures == "crasher_t7" ~ "late_crasher",
                                signatures == "probiotic_t7_strict" ~ "late_probiotic",
                                signatures == "early_pathogen, late_opportunist" ~ "early_pathogen"
                                TRUE ~ signatures)) %>%
  inner_join(significant_models,
             by = 'asv_id') %>%
  rowwise %>%
  mutate(plot_info = list(emmeans(model, ~treatment) %>%
                            broom::tidy(conf.int = TRUE) %>%
                            separate(treatment, into = c('time', 'exposure', 'susceptability')) %>%
                            mutate(graph_cat = ifelse(time == "T0", NA, 
                                                      paste(exposure, susceptability, sep = "_"))) %>%
                            {. ->> intermed } %>%
                            mutate(graph_cat = ifelse(time == "T0", paste("H", susceptability, sep = "_"), 
                                                      graph_cat)) %>%
                            dplyr::slice(rep(1:2, 1)) %>%
                            rbind(intermed) %>%
                            mutate(graph_cat = ifelse(is.na(graph_cat), paste("D", susceptability, sep = "_"), 
                                                      graph_cat)) %>%
                            mutate(c_time = parse_number(time)) %>%
                            mutate(facet_lab = "Experimental") %>%
                            mutate(c_time = ifelse(time == "T0", ifelse(susceptability == "S", c_time - 0.1, c_time + 0.1),
                                                   case_when(graph_cat == "D_S" ~ c_time - 0.35,
                                                             graph_cat == "H_S" ~ c_time - 0.15,
                                                             graph_cat == "D_R" ~ c_time + 0.15,
                                                             graph_cat == "H_R" ~ c_time + 0.35))) %>%
                            mutate(graph_cat = factor(graph_cat, levels = c("D_S", "D_R", "H_S", "H_R"), labels = c("D_S", "D_R", "H_S", "H_R"))))) %>%
  left_join(homogenate_models, by = join_by(asv_id)) %>%
  mutate(plot_info = list(rbind(plot_info, homogenate_pred) %>% mutate(sig_p = ifelse(p.value < 0.05, "sig", "nonsig")))) %>%
  rowwise() %>%
  mutate(plot = list(
    ggplot(data = plot_info, aes(x = c_time, y = estimate, ymin = conf.low, ymax = conf.high)) +
      (geom_line(data = (plot_info %>% filter(graph_cat %in% c("D_S", "D_R"))), aes(colour1 = graph_cat, linetype = graph_cat)) %>%
         rename_geom_aes(new_aes = c("colour" = "colour1"))) + 
      (geom_line(data = (plot_info %>% filter(graph_cat %in% c("H_S", "H_R"))), aes(colour2 = graph_cat, linetype = graph_cat)) %>%
         rename_geom_aes(new_aes = c("colour" = "colour2"))) +
      (geom_errorbar(data = (plot_info %>% filter(graph_cat %in% c("D_S", "D_R"))), width = 0, aes(colour1 = graph_cat)) %>%
         rename_geom_aes(new_aes = c("colour" = "colour1"))) + 
      (geom_errorbar(data = (plot_info %>% filter(graph_cat %in% c("H_S", "H_R"))), width = 0, aes(colour2 = graph_cat)) %>%
         rename_geom_aes(new_aes = c("colour" = "colour2"))) +
      (geom_point(data = (plot_info %>% filter(graph_cat %in% c("D_S", "D_R"))), size = 3, aes(colour1 = graph_cat, pch = graph_cat)) %>%
         rename_geom_aes(new_aes = c("colour" = "colour1"))) + 
      (geom_point(data = (plot_info %>% filter(graph_cat %in% c("H_S", "H_R"))), size = 3, aes(colour2 = graph_cat, pch = graph_cat)) %>%
         rename_geom_aes(new_aes = c("colour" = "colour2"))) +
      
      geom_point(data = (plot_info %>% filter(exposure == "p_val")), size = 5, col = "gray20", pch = "*", aes(alpha = sig_p)) +
      
      (geom_point(data = (plot_info %>% filter(graph_cat == "dose" & susceptability == "na")), size = 3.7, aes(colour3 = exposure), shape = "diamond") %>%
         rename_geom_aes(new_aes = c("colour" = "colour3"))) +
      (geom_errorbar(data = (plot_info %>% filter(graph_cat == "dose" & susceptability == "na")), width = 0, aes(colour3 = exposure)) %>%
         rename_geom_aes(new_aes = c("colour" = "colour3"))) + 
      (geom_point(data = (plot_info %>% filter(graph_cat == "dose" & is.na(exposure))), size = 3, shape = 1, col = "black", alpha = 0)) +
      
      scale_color_manual(aesthetics = "colour1", values = c("#F75D5D", "#A70000"), guide = "legend", 
                         name = "Disease Exposed", labels = c("Susceptible", "Resistant")) +
      scale_color_manual(aesthetics = "colour3", values = c("#F10A0A", "#29A5B2"), guide = "legend", 
                         name = "Doses", labels = c("Diseased", "Healthy")) +
      scale_shape_manual(values = c(17, 16, 17, 16), guide = "none") +
      scale_color_manual(aesthetics = "colour2", values = c("#3DD8EA", "#048291"), guide = "legend", 
                         name = "Healthy Exposed", labels = c("Susceptible", "Resistant")) +
      guides(colour1 = guide_legend(
        override.aes=list(linetype = c(6, 1), shape = c(16, 17))),
        colour2 = guide_legend(
          override.aes=list(linetype = c(6, 1), shape = c(16, 17))),
        colour3 = guide_legend(
          override.aes=list(linetype = c(0, 0)))) +
      scale_x_continuous(breaks=c(0, 3, 7)) +
      scale_alpha_manual(values = c("sig" = 1, "nonsig" = 0), guide = "none") +
      scale_linetype_manual(values = c(1, 6, 1, 6), guide = "none") +
      theme_bw() +
      xlab("Time") +
      ylab(expression("Normalized log"[2]*" (cpm)")) +
      labs(title = str_c(Family, " ", Genus, " (", asv_id, ")", sep = "")) +
      facet_grid(cols = vars(facet_lab), space = "free", scales = "free_x")
  )) %>%
  group_by(signatures) %>%
  summarise(combo_plots = list(wrap_plots(plot) + plot_layout(guides = 'collect') & plot_annotation(title = signatures)))

#view the plots
the_plots$combo_plots[[1]]

#### Bac Strat NMDS and PCOA ####

clean_bacstrats <- bacterial_signature_asv %>% 
  group_by(asv_id) %>% 
  reframe(bacstrat = paste(list(signatures))) %>%
  mutate(bacstrat = case_when(bacstrat == "c(\"early_pathogen\", \"continuous_pathogen\", \"late_pathogen\")" ~ "continuous_pathogen",
                              bacstrat == 'c("crasher_t3", "crasher_t7")' ~ "continuous_crasher",
                              bacstrat == 'c("early_pathogen", "early_opportunist")' ~ "early_pathogen",
                              bacstrat == 'c("late_pathogen", "late_opportunist")' ~ "late_pathogen",
                              bacstrat == "c(\"probiotic_t7_strict\", \"late_opportunist\")" ~ "late_probiotic",
                              bacstrat == "crasher_t7" ~ "late_crasher",
                              bacstrat == "probiotic_t7_strict" ~ "late_probiotic",
                              bacstrat == 'c("early_pathogen", "late_opportunist")' ~ "early_pathogen",
                              TRUE ~ bacstrat))

asv_nmds <- full_data %>% 
  select(asv_id, sample_id, log2_cpm) %>%
  pivot_wider(names_from = sample_id, values_from = log2_cpm) %>%
  left_join(clean_bacstrats, by = join_by("asv_id")) %>%
  filter(!is.na(bacstrat)) %>%
  filter(!bacstrat %in% c("continuous_crasher", "late_crasher")) %>%
  select(-bacstrat) %>%
  column_to_rownames('asv_id') %>%
  t() %>%
  as.matrix()

test_nmds <- metaMDS(t(asv_nmds), distance = 'bray', k = 2, trymax = 100, autotransform = FALSE, verbose = TRUE)
sppscores(test_nmds) <- t(asv_nmds)

#shepard plot
plot(test_nmds$diss, test_nmds$dist)

suscep <- full_data %>% 
  select(asv_id, sample_id, resistance) %>%
  pivot_wider(names_from = sample_id, values_from = resistance) %>%
  left_join(clean_bacstrats, by = join_by("asv_id")) %>%
  filter(!is.na(bacstrat)) %>%
  select(-bacstrat) %>%
  column_to_rownames('asv_id') %>%
  t() %>%
  as.matrix()

#adonis2 fails but adonis works
perm.results <- vegan::adonis(asv_nmds ~ suscep, method="bray",perm=999)
perm.results$aov.tab #p value less than 0.05


#plot species or site alone
plot(test_nmds, "species")
orditorp(test_nmds, "species")


## NMDS with 95% CI around TIMEPOINTS

scores(test_nmds)$site %>%
  as_tibble(rownames = 'asv_names') %>%
  left_join(clean_bacstrats, by = c('asv_names' = 'asv_id')) %>%
  filter(!is.na(bacstrat)) %>%
  mutate(shape_determination = case_when(
    str_detect(bacstrat, "path")  ~ "pathogen",
    str_detect(bacstrat, "opp")  ~ "opportunist",
    str_detect(bacstrat, "crash")  ~ "crasher",
    str_detect(bacstrat, "pro")  ~ "probiotic"
  )) %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  stat_ellipse(data = as_tibble(scores(test_nmds)$species, rownames = aggregation_level) %>% 
                 mutate(exp_conditions = str_split(none, "_")) %>%
                 rowwise() %>%
                 mutate(time = exp_conditions[1], exposure = exp_conditions[2], resist = exp_conditions[3]), 
               geom = "polygon", alpha = 0.1, level = 0.95, aes(col = time, fill = time)) +
  geom_point(data = as_tibble(scores(test_nmds)$species, rownames = aggregation_level) %>% 
               mutate(exp_conditions = str_split(none, "_")) %>%
               rowwise() %>%
               mutate(time = exp_conditions[1], exposure = exp_conditions[2], resist = exp_conditions[3]),
             size = 1, aes(color = time, shape = exposure)) +
  geom_point(aes(col = bacstrat, shape = shape_determination), size = 4) + #alpha = shape_determination
  geom_text(aes(label = parse_number(asv_names)), size = 3) +
  scale_shape_manual(values = c(
    "opportunist" = 16,
    "pathogen" = 17,
    "crasher" = 15,
    "probiotic" = 8,
    "H" = 6,
    "D" = 4,
    "REP1" = 3
  )) +
  scale_fill_manual(values = c(
    "T0" = "gray80",
    "T3" = "gray40",
    "T7" = "gray10"
  ), guide = "none") +
  scale_color_manual(values = c(
    "continuous_crasher" = "sandybrown",
    "continuous_pathogen" = "firebrick4",
    "late_crasher" = "darkorange2",
    "late_probiotic" = "aquamarine",
    "early_pathogen" = "hotpink",
    "late_pathogen" = "firebrick1",
    "early_opportunist" = "deepskyblue3",
    "late_opportunist" = "royalblue2",
    "T0" = "gray80",
    "T3" = "gray40",
    "T7" = "gray10"
  )) +
  guides(color = guide_legend(override.aes = list(fill = NA))) +
  theme_bw()



facet_nmds_t3 <- scores(test_nmds)$site %>%
  as_tibble(rownames = 'asv_names') %>%
  left_join(clean_bacstrats, by = c('asv_names' = 'asv_id')) %>%
  filter(!is.na(bacstrat)) %>%
  mutate(shape_determination = case_when(
    str_detect(bacstrat, "path")  ~ "pathogen",
    str_detect(bacstrat, "opp")  ~ "opportunist",
    str_detect(bacstrat, "crash")  ~ "crasher",
    str_detect(bacstrat, "pro")  ~ "probiotic"
  )) %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  stat_ellipse(data = as_tibble(scores(test_nmds)$species, rownames = "sample_id") %>% 
                 left_join((full_data %>% select(sample_id, susceptability) %>% distinct()), by = join_by("sample_id")) %>%
                 mutate(exp_conditions = str_split(sample_id, "_")) %>%
                 rowwise() %>%
                 mutate(time = exp_conditions[1], exposure = exp_conditions[2], tank = exp_conditions[3]) %>%
                 mutate(exp_res = paste(exposure, susceptability, sep = "_")) %>%
                 filter(time == "T3"), 
               geom = "polygon", alpha = 0.1, level = 0.95, aes(col = exp_res, fill = exp_res, linetype = exp_res)) +
  geom_point(data = as_tibble(scores(test_nmds)$species, rownames = "sample_id") %>% 
               left_join((full_data %>% select(sample_id, susceptability) %>% distinct()), by = join_by("sample_id")) %>%
               mutate(exp_conditions = str_split(sample_id, "_")) %>%
               rowwise() %>%
               mutate(time = exp_conditions[1], exposure = exp_conditions[2], tank = exp_conditions[3]) %>%
               mutate(exp_res = paste(exposure, susceptability, sep = "_")),
             size = 1, aes(color = time, shape = exposure)) +
  geom_point(data = as_tibble(scores(test_nmds)$species, rownames = "sample_id") %>% 
               left_join((full_data %>% select(sample_id, susceptability) %>% distinct()), by = join_by("sample_id")) %>%
               mutate(exp_conditions = str_split(sample_id, "_")) %>%
               rowwise() %>%
               mutate(time = exp_conditions[1], exposure = exp_conditions[2], tank = exp_conditions[3]) %>%
               mutate(exp_res = paste(exposure, susceptability, sep = "_")) %>%
               filter(time == "T3"),
             size = 1, aes(color = exp_res, shape = exposure)) +
  geom_point(aes(col = bacstrat, shape = shape_determination), size = 4) + #alpha = shape_determination
  geom_text(aes(label = parse_number(asv_names)), size = 3) +
  scale_linetype_manual(values = c(
    "D_S" = "dashed",
    "H_S" = "dashed",
    "D_R" = "solid",
    "H_R" = "solid"
  ), guide = "none") +
  scale_shape_manual(values = c(
    "opportunist" = 16,
    "pathogen" = 17,
    "crasher" = 15,
    "probiotic" = 18,
    "H" = 6,
    "D" = 4,
    "REP1" = 7
  )) +
  scale_fill_manual(values = c(
    "T0" = "gray80",
    "T3" = "gray40",
    "T7" = "gray10",
    "D_S" = "#F75D5D",
    "H_S" = "#3DD8EA",
    "D_R" = "#A70000",
    "H_R" = "#048291"
  ), guide = "none") +
  scale_color_manual(values = c(
    "continuous_crasher" = "sandybrown",
    "continuous_pathogen" = "firebrick4",
    "late_crasher" = "darkorange2",
    "late_probiotic" = "aquamarine2",
    "early_pathogen" = "hotpink",
    "late_pathogen" = "firebrick1",
    "early_opportunist" = "deepskyblue3",
    "late_opportunist" = "royalblue2",
    "T0" = "gray80",
    "T3" = "gray55",
    "T7" = "gray55",
    "D_S" = "#F75D5D",
    "H_S" = "#3DD8EA",
    "D_R" = "#A70000",
    "H_R" = "#048291"
  )) +
  guides(color = guide_legend(override.aes = list(fill = NA))) +
  theme_bw() +
  ylim(-0.26, 0.245) +
  xlim(-0.25, 0.33) +
  ggtitle("T3")

facet_nmds_t7 <- scores(test_nmds)$site %>%
  as_tibble(rownames = 'asv_names') %>%
  left_join(clean_bacstrats, by = c('asv_names' = 'asv_id')) %>%
  filter(!is.na(bacstrat)) %>%
  mutate(shape_determination = case_when(
    str_detect(bacstrat, "path")  ~ "pathogen",
    str_detect(bacstrat, "opp")  ~ "opportunist",
    str_detect(bacstrat, "crash")  ~ "crasher",
    str_detect(bacstrat, "pro")  ~ "probiotic"
  )) %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  stat_ellipse(data = as_tibble(scores(test_nmds)$species, rownames = "sample_id") %>% 
                 left_join((full_data %>% select(sample_id, susceptability) %>% distinct()), by = join_by("sample_id")) %>%
                 mutate(exp_conditions = str_split(sample_id, "_")) %>%
                 rowwise() %>%
                 mutate(time = exp_conditions[1], exposure = exp_conditions[2], tank = exp_conditions[3]) %>%
                 mutate(exp_res = paste(exposure, susceptability, sep = "_")) %>%
                 filter(time == "T7"), 
               geom = "polygon", alpha = 0.1, level = 0.95, aes(col = exp_res, fill = exp_res, linetype = exp_res)) +
  geom_point(data = as_tibble(scores(test_nmds)$species, rownames = "sample_id") %>% 
               left_join((full_data %>% select(sample_id, susceptability) %>% distinct()), by = join_by("sample_id")) %>%
               mutate(exp_conditions = str_split(sample_id, "_")) %>%
               rowwise() %>%
               mutate(time = exp_conditions[1], exposure = exp_conditions[2], tank = exp_conditions[3]) %>%
               mutate(exp_res = paste(exposure, susceptability, sep = "_")),
             size = 1, aes(color = time, shape = exposure)) +
  geom_point(data = as_tibble(scores(test_nmds)$species, rownames = "sample_id") %>% 
               left_join((full_data %>% select(sample_id, susceptability) %>% distinct()), by = join_by("sample_id")) %>%
               mutate(exp_conditions = str_split(sample_id, "_")) %>%
               rowwise() %>%
               mutate(time = exp_conditions[1], exposure = exp_conditions[2], tank = exp_conditions[3]) %>%
               mutate(exp_res = paste(exposure, susceptability, sep = "_")) %>%
               filter(time == "T7"),
             size = 1, aes(color = exp_res, shape = exposure)) +
  geom_point(aes(col = bacstrat, shape = shape_determination), size = 4) + #alpha = shape_determination
  geom_text(aes(label = parse_number(asv_names)), size = 3) +
  scale_linetype_manual(values = c(
    "D_S" = "dashed",
    "H_S" = "dashed",
    "D_R" = "solid",
    "H_R" = "solid"
  ), guide = "none") +
  scale_shape_manual(values = c(
    "opportunist" = 16,
    "pathogen" = 17,
    "crasher" = 15,
    "probiotic" = 18,
    "H" = 6,
    "D" = 4,
    "REP1" = 7
  )) +
  scale_fill_manual(values = c(
    "T0" = "gray80",
    "T3" = "gray40",
    "T7" = "gray10",
    "D_S" = "#F75D5D",
    "H_S" = "#3DD8EA",
    "D_R" = "#A70000",
    "H_R" = "#048291"
  ), guide = "none") +
  scale_color_manual(values = c(
    "continuous_crasher" = "sandybrown",
    "continuous_pathogen" = "firebrick4",
    "late_crasher" = "darkorange2",
    "late_probiotic" = "aquamarine2",
    "early_pathogen" = "hotpink",
    "late_pathogen" = "firebrick1",
    "early_opportunist" = "deepskyblue3",
    "late_opportunist" = "royalblue2",
    "T0" = "gray80",
    "T3" = "gray55",
    "T7" = "gray55",
    "D_S" = "#F75D5D",
    "H_S" = "#3DD8EA",
    "D_R" = "#A70000",
    "H_R" = "#048291"
  )) +
  guides(color = guide_legend(override.aes = list(fill = NA))) +
  theme_bw() +
  ylim(-0.26, 0.245) +
  xlim(-0.25, 0.33) +
  ggtitle("T7")

wrap_plots(facet_nmds_t3, facet_nmds_t7) + plot_layout(guides = 'collect')

## PCOA with 95% CI around BACTERIAL STRATEGIES

bacstrat_dist_mat <- vegdist(t(asv_nmds))
bacstrat_pcoa <- cmdscale (bacstrat_dist_mat, eig = TRUE)
ordiplot (bacstrat_pcoa, display = 'sites', type = 'text')

bacstrat_pcoa$points %>%
  as_tibble(rownames = 'asv_names') %>%
  left_join(clean_bacstrats, by = c('asv_names' = 'asv_id')) %>%
  filter(!is.na(bacstrat)) %>%
  mutate(shape_determination = case_when(
    str_detect(bacstrat, "path")  ~ "pathogen",
    str_detect(bacstrat, "opp")  ~ "opportunist",
    str_detect(bacstrat, "crash")  ~ "crasher",
    str_detect(bacstrat, "pro")  ~ "probiotic"
  )) %>%
  ggplot(aes(x = V1, y = V2)) +
  stat_ellipse(geom = "polygon", alpha = 0.1, level = 0.95, aes(col = bacstrat, fill = bacstrat)) +
  geom_point(aes(col = bacstrat, shape = shape_determination), size = 4) + 
  geom_text(aes(label = parse_number(asv_names)), size = 3) +
  scale_shape_manual(values = c(
    "opportunist" = 16,
    "pathogen" = 17,
    "crasher" = 15,
    "probiotic" = 8,
    "H" = 6,
    "D" = 4,
    "REP1" = 3
  )) +
  scale_fill_manual(values = c(
    "continuous_crasher" = "sandybrown",
    "continuous_pathogen" = "firebrick4",
    "late_crasher" = "darkorange2",
    "late_probiotic" = "aquamarine",
    "early_pathogen" = "hotpink",
    "late_pathogen" = "firebrick1",
    "early_opportunist" = "deepskyblue3",
    "late_opportunist" = "royalblue2"
  ), guide = "none") +
  scale_color_manual(values = c(
    "continuous_crasher" = "sandybrown",
    "continuous_pathogen" = "firebrick4",
    "late_crasher" = "darkorange2",
    "late_probiotic" = "aquamarine",
    "early_pathogen" = "hotpink",
    "late_pathogen" = "firebrick1",
    "early_opportunist" = "deepskyblue3",
    "late_opportunist" = "royalblue2",
    "T0" = "gray80",
    "T3" = "gray40",
    "T7" = "gray10"
  )) +
  guides(color = guide_legend(override.aes = list(fill = NA))) +
  theme_bw() #+
  ylim(-0.3, 0.3) +
  xlim(-0.3, 0.25)


#### Emily Miscellaneous ####

    ## Access all ASVs in a given bacterial strategy:
  
asvs_by_signature <- clean_bacstrats %>%
  filter(bacstrat == "late_pathogen") %>%
  #filter(signatures %in% c("continuous_crasher", "late_crasher")) %>%
  pull(asv_id)
  
  
    ## heritability (random effect of genotype)

heritable_asvs <- asv_models %>%
  rowwise() %>%
  mutate(heritability = model %>% broom.mixed::tidy("ran_pars") 
         %>% filter(group == "genotype") %>% pull(estimate), .after = asv_id) %>%
  filter(heritability > 0)

#make formattable chart of heritability scores:
heritable_asvs %>% arrange(desc(heritability)) %>% filter(asv_id %in% asvs_by_signature) %>%
  select(Family, Genus, asv_id, heritability) %>% 
  mutate(heritability = round(heritability, 3)) %>%
  formattable(align = c("c", "c", "c", "c"), list(
    area(col = 4) ~ color_tile("lavender", "purple2")
  )) #%>%
  export_formattable("../Figures/put_pathogens_heritability.png")

  
    ##logfold change between D and H

#make formattable chart showing logfold change of D compared to H at Dose and T7
significant_models %>%
  filter(asv_id %in% asvs_by_signature) %>%
  select(asv_id, data) %>%
  unnest(data) %>%
  select(asv_id, sample_id, log2_cpm, exposure, time) %>%
  rbind(homogenate_data %>%
          mutate(time = "dose")) %>%
  filter(time %in% c("dose", "T7")) %>%
  group_by(asv_id, time, exposure) %>%
  reframe(ave_val = mean(log2_cpm)) %>%
  ungroup() %>%
  group_by(asv_id, time) %>%
  reframe(logfold_change = ave_val[exposure == "D"] - ave_val[exposure == "H"]) %>%
  mutate(logfold_change = round(logfold_change, 2)) %>%
  pivot_wider(names_from = time, values_from = logfold_change) %>%
  select(asv_id, dose, T7) %>%
  filter(!is.na(dose) & !is.na(T7)) %>%
  rename("asv_id" = "ASV", "dose" = "Dose Logfold Change", "T7" = "T7 Logfold Change") %>%
  formattable(align = c("l", "c", "c"), list(
    area(col = c(2,3)) ~ color_tile("white", "firebrick1")
  )) %>%
  export_formattable("../Figures/put_pathogens_logfold_table.png")


    ## lm models predicting w T0
  #TAKEAWAY: T0 abundance cannot predict disease resistance

t0_prediction_lms <- full_data %>%
  filter(time == "T0") %>%
  group_by(asv_id) %>%
  reframe(model = list(lm(log2_cpm ~ susceptability) %>% broom::tidy() %>% 
                         filter(term == "susceptability1"))) %>%
  unnest(model) %>%
  filter(p.value < 0.05)

bacterial_signature_asv %>% filter(asv_id %in% t0_prediction_lms$asv_id)


    ## Correlation Matrices

#corr matrix for 8 putative pathogen candidates only

asv_corr <- full_data %>% 
  select(asv_id, sample_id, Family, log2_cpm) %>%
  pivot_wider(names_from = sample_id, values_from = log2_cpm) %>%
  mutate(combo_name = paste(Family, asv_id, sep = "_"), .after = asv_id) %>%
  left_join(tgfff, by = join_by("asv_id")) %>%
  filter(!is.na(bacstrat)) %>%
  filter(bacstrat %in% c('c("late_pathogen", "late_opportunist")', "late_pathogen")) %>%
  select(-c(bacstrat, asv_id, Family)) %>%
  column_to_rownames('combo_name') %>%
  t() %>%
  as.matrix()

asv_corr_mat <- rcorr(asv_corr, type = c("pearson","spearman"))



corrplot(asv_corr_mat$r, type="upper", order="hclust", 
         p.mat = asv_corr_mat$P, sig.level = 0.05, insig = "blank",
         col.lim = c(0,1), col = rep(brewer.pal(11,"Spectral"), 2), tl.col = "gray20")

brewer.pal(11,"Spectral")
# corr matrix for All Colwelliaceae Thalassotaleas

colwell_corr <- full_data %>% 
  filter(Genus == "Thalassotalea") %>%
  select(asv_id, sample_id, log2_cpm) %>%
  pivot_wider(names_from = sample_id, values_from = log2_cpm) %>%
  column_to_rownames('asv_id') %>%
  t() %>%
  as.matrix()

colwell_test <- rcorr(colwell_corr, type = c("pearson","spearman"))

corrplot(colwell_test$r, type="upper", order="hclust", 
         p.mat = colwell_test$P, sig.level = 0.05, insig = "blank",
         #col.lim = c(-0.2,1), col = c(rep("gray30", 10), rainbow(20)), tl.col = "gray20")
         col.lim = c(-0.2,1), col = c(rep("gray30", 6), brewer.pal(11,"Spectral")), tl.col = "gray20")

#### alpha diversity ####

alpha_table <- microbiome::alpha(microbiome_data, index = "all") %>%
  as_tibble(rownames = 'sample_id') %>%
  inner_join(metadata, by = 'sample_id') %>%
  mutate(fragment_id = str_c(str_replace_na(exposure, 'NA'), tank, genotype, sep = '_'))

mod_alpha_tab <- alpha_table %>%
  filter(!tank %in% c("HOMO", "homogenate_fragment")) %>%
  mutate(treatment = str_c(time, exposure, susceptability, sep = '_')) %>%
  pivot_longer(cols = !c(colnames(metadata), "fragment_id", "treatment"),
               names_to = 'metric',
               values_to = 'alpha_div_value') %>%
  select(-c(retain_sample, final_disease_state, clone_group)) %>%
  nest_by(metric) %>%
  summarise(alpha_model = list(lmer(alpha_div_value ~ treatment + 
                                (1 | genotype) + (1 | tank), data = data))) %>% 
  rowwise() %>% 
  mutate(sig_terms = list(anova(alpha_model) %>% 
                            rownames_to_column(var = "sig_term") %>% 
                            as_tibble() %>% 
                            rename(`Pr(>F)` = "p_val") %>%
                            mutate(fdr_p_val = p.adjust(p_val, method = 'fdr')) %>%
                            filter(fdr_p_val < 0.05) %>%
                            filter(sig_term != "time") %>%
                            pull(sig_term))) %>%
  filter(length(sig_terms) > 0) %>%
  select(-sig_terms) %>%
  mutate(alpha_type = ifelse(metric %in% c("chao1", "observed"), "richness", str_extract(metric, "[^_]+")))
           
## figuring out the fig

alpha_graphs <- mod_alpha_tab %>%
  rowwise() %>%
  mutate(plot_info = list(emmeans(alpha_model, ~treatment) %>%
                          broom::tidy(conf.int = TRUE) %>%
                          separate(treatment, into = c('time', 'exposure', 'susceptability')) %>%
                          mutate(graph_cat = ifelse(time == "T0", NA, 
                                                    paste(exposure, susceptability, sep = "_"))) %>%
                          {. ->> intermed } %>%
                          mutate(graph_cat = ifelse(time == "T0", paste("H", susceptability, sep = "_"), 
                                                    graph_cat)) %>%
                          dplyr::slice(rep(1:2, 1)) %>%
                          rbind(intermed) %>%
                          mutate(graph_cat = ifelse(is.na(graph_cat), paste("D", susceptability, sep = "_"), 
                                                    graph_cat)) %>%
                          mutate(c_time = parse_number(time)) %>%
                          mutate(facet_lab = "Experimental") %>%
                          mutate(c_time = ifelse(time == "T0", ifelse(susceptability == "S", c_time - 0.1, c_time + 0.1),
                                                 case_when(graph_cat == "D_S" ~ c_time - 0.35,
                                                           graph_cat == "H_S" ~ c_time - 0.15,
                                                           graph_cat == "D_R" ~ c_time + 0.15,
                                                           graph_cat == "H_R" ~ c_time + 0.35))) %>%
                          mutate(graph_cat = factor(graph_cat, levels = c("D_S", "D_R", "H_S", "H_R"), labels = c("D_S", "D_R", "H_S", "H_R"))))) %>%
  rowwise() %>%
  mutate(plot = list(
    ggplot(data = plot_info, aes(x = c_time, y = estimate, ymin = conf.low, ymax = conf.high)) +
      (geom_line(data = (plot_info %>% filter(graph_cat %in% c("D_S", "D_R"))), aes(colour1 = graph_cat, linetype = graph_cat)) %>%
         rename_geom_aes(new_aes = c("colour" = "colour1"))) + 
      (geom_line(data = (plot_info %>% filter(graph_cat %in% c("H_S", "H_R"))), aes(colour2 = graph_cat, linetype = graph_cat)) %>%
         rename_geom_aes(new_aes = c("colour" = "colour2"))) +
      (geom_errorbar(data = (plot_info %>% filter(graph_cat %in% c("D_S", "D_R"))), width = 0, aes(colour1 = graph_cat)) %>%
         rename_geom_aes(new_aes = c("colour" = "colour1"))) + 
      (geom_errorbar(data = (plot_info %>% filter(graph_cat %in% c("H_S", "H_R"))), width = 0, aes(colour2 = graph_cat)) %>%
         rename_geom_aes(new_aes = c("colour" = "colour2"))) +
      (geom_point(data = (plot_info %>% filter(graph_cat %in% c("D_S", "D_R"))), size = 3, aes(colour1 = graph_cat, pch = graph_cat)) %>%
         rename_geom_aes(new_aes = c("colour" = "colour1"))) + 
      (geom_point(data = (plot_info %>% filter(graph_cat %in% c("H_S", "H_R"))), size = 3, aes(colour2 = graph_cat, pch = graph_cat)) %>%
         rename_geom_aes(new_aes = c("colour" = "colour2"))) +
      
      
      (geom_point(data = (plot_info %>% filter(graph_cat == "dose" & susceptability == "na")), size = 3.7, aes(colour3 = exposure), shape = "diamond") %>%
         rename_geom_aes(new_aes = c("colour" = "colour3"))) +
      (geom_errorbar(data = (plot_info %>% filter(graph_cat == "dose" & susceptability == "na")), width = 0, aes(colour3 = exposure)) %>%
         rename_geom_aes(new_aes = c("colour" = "colour3"))) + 
      (geom_point(data = (plot_info %>% filter(graph_cat == "dose" & is.na(exposure))), size = 3, shape = 1, col = "black", alpha = 0)) +
      
      scale_color_manual(aesthetics = "colour1", values = c("#F75D5D", "#A70000"), guide = "legend", 
                         name = "Disease Exposed", labels = c("Susceptible", "Resistant")) +
      scale_color_manual(aesthetics = "colour3", values = c("#F10A0A", "#29A5B2"), guide = "legend", 
                         name = "Doses", labels = c("Diseased", "Healthy")) +
      scale_shape_manual(values = c(17, 16, 17, 16), guide = "none") +
      scale_color_manual(aesthetics = "colour2", values = c("#3DD8EA", "#048291"), guide = "legend", 
                         name = "Healthy Exposed", labels = c("Susceptible", "Resistant")) +
      guides(colour1 = guide_legend(
        override.aes=list(linetype = c(6, 1), shape = c(16, 17))),
        colour2 = guide_legend(
          override.aes=list(linetype = c(6, 1), shape = c(16, 17))),
        colour3 = guide_legend(
          override.aes=list(linetype = c(0, 0)))) +
      scale_x_continuous(breaks=c(0, 3, 7)) +
      
      scale_linetype_manual(values = c(1, 6, 1, 6), guide = "none") +
      theme_bw() +
      xlab("Time") +
      ylab(expression("Normalized log"[2]*" (cpm)")) +
      labs(title = metric)
  )) %>%
  group_by(alpha_type) %>%
  summarise(combo_plots = list(wrap_plots(plot) + plot_layout(guides = 'collect') & plot_annotation(title = alpha_type)))


alpha_graphs$combo_plots[[5]]



#[WORK IN PROGRESS]

#### Simplified Complex Upset ####
simple_comp_upset <- clean_bacstrats %>% 
  mutate(sig = TRUE) %>%
  pivot_wider(names_from = bacstrat, values_from = sig) %>%
  mutate(across(-asv_id, ~ifelse(is.na(.), FALSE, .))) %>% 
  left_join(taxonomy_tibble, by = c('asv_id' = 'asv_names'))


### Upset - bacterial strategies

upset(
  simple_comp_upset,
  colnames(simple_comp_upset %>%
             select(starts_with("cont"), starts_with("late"), starts_with("early"))), 
  
  base_annotations=list(
    'Intersection size'=intersection_size(counts=T, text = aes(size = 6.5),
                                          bar_number_threshold = 25,
                                          mapping=aes(fill=Family, col = Family, label = parse_number(asv_id)), col = "gray10"
    ) +
      geom_text(size = 4, position = position_stack(vjust = 0.5), col = "gray10")
  ),
  matrix=(
    intersection_matrix(geom=geom_point(shape = "circle filled", size=3, stroke = 0.35, color = "gray20"))
    + scale_color_manual(
      values=c(
        "late_pathogen" = "#FF1E1E",
        "early_pathogen" = "#FF9797",
        "continuous_pathogen" = "firebrick4",
        "late_opportunist" = "deepskyblue3",
        "early_opportunist" = "lightskyblue",
        "continuous_crasher" = "#FF9A00",
        "late_crasher" = "#BD5E03",
        "late_probiotic" = "#49EACC"
      )
    )
  ),
  queries=list(
    upset_query(set = "late_pathogen", fill = "#FF1E1E"),
    upset_query(set = "early_pathogen", fill = "#FF9797"),
    upset_query(set = "continuous_pathogen", fill = "firebrick4"),
    upset_query(set = "late_opportunist", fill = "deepskyblue3"),
    upset_query(set = "early_opportunist", fill = "lightskyblue"),
    upset_query(set = "continuous_crasher", fill = "#FF9A00"),
    upset_query(set = "late_crasher", fill = "#BD5E03"),
    upset_query(set = "late_probiotic", fill = "#49EACC")
  ),
  name='ASVs', width_ratio=0.1, min_size = 0) + #, sort_sets = FALSE
  ggtitle("Bacterial Strategies") 

### Upset - bacterial strategies and Presence Data

upset(
  simple_comp_upset %>% left_join(venn_all_times_and_doses, by = c("asv_id" = "OTU")),
  colnames(simple_comp_upset %>% left_join(venn_all_times_and_doses, by = c("asv_id" = "OTU")) %>%
             select(T0, D, H, starts_with("cont"), starts_with("late"), starts_with("early"))),
  base_annotations=list(
    'Intersection size'=intersection_size(counts=T, text = aes(size = 6.5),
                                          bar_number_threshold = 25,
                                          mapping=aes(fill=Family, col = Family, label = parse_number(asv_id)), col = "gray10"
    ) +
      geom_text(size = 4, position = position_stack(vjust = 0.5), col = "gray10")
  ),
  matrix=(
    intersection_matrix(geom=geom_point(shape = "circle filled", size=3, stroke = 0.35, color = "gray20"))
    + scale_color_manual(
      values=c(
        #"T7" = "#650197",
        #"T3" = "#B21BFF",
        "T0" = "#B21BFF",
        "D" = "#AE0404",
        "H" = "#0FAB02",
        "late_pathogen" = "#FF1E1E",
        "early_pathogen" = "#FF9797",
        "continuous_pathogen" = "firebrick4",
        "late_opportunist" = "deepskyblue3",
        "early_opportunist" = "lightskyblue",
        "continuous_crasher" = "#FF9A00",
        "late_crasher" = "#BD5E03",
        "late_probiotic" = "#49EACC"
      )
    )
  ),
  queries=list(
    #upset_query(set = "T7", fill = "#650197"),
    #upset_query(set = "T3", fill = "#B21BFF"), "#D98EFF"
    upset_query(set = "T0", fill = "#B21BFF"),
    upset_query(set = "D", fill = "#AE0404"),
    upset_query(set = "H", fill = "#0FAB02"),
    upset_query(set = "late_pathogen", fill = "#FF1E1E"),
    upset_query(set = "early_pathogen", fill = "#FF9797"),
    upset_query(set = "continuous_pathogen", fill = "firebrick4"),
    upset_query(set = "late_opportunist", fill = "deepskyblue3"),
    upset_query(set = "early_opportunist", fill = "lightskyblue"),
    upset_query(set = "continuous_crasher", fill = "#FF9A00"),
    upset_query(set = "late_crasher", fill = "#BD5E03"),
    upset_query(set = "late_probiotic", fill = "#49EACC")
  ),
  name='ASVs', width_ratio=0.1, min_size = 0, sort_sets = FALSE,
  stripes = c("#F8EAFF", "#FFE4E4", "#E7FFE5", rep(c("gray91", "gray97"), 4))) +
  #stripes = c(rep(c("gray78", "gray87"), 2), "gray78", rep(c("gray97", "gray91"), 4))) +
  ggtitle("Bacterial Strategies/Presence Data") 


### Pathogens Only

pathogen_upset <- upset(
  simple_comp_upset %>% left_join(venn_all_times_and_doses, by = c("asv_id" = "OTU")) %>% 
    filter(continuous_pathogen | early_pathogen | late_pathogen),
  colnames(simple_comp_upset %>% left_join(venn_all_times_and_doses, by = c("asv_id" = "OTU")) %>%
             select(T0, D, H, continuous_pathogen, late_pathogen, early_pathogen)),
  base_annotations=list(
    'Intersection size'=intersection_size(counts=T, text = aes(size = 6.5),
                                          bar_number_threshold = 25,
                                          mapping=aes(fill=Family, col = Family, label = parse_number(asv_id)), col = "gray10"
    ) +
      geom_text(size = 4, position = position_stack(vjust = 0.5), col = "gray10") +
      ylim(0, 8) +
    scale_fill_manual(values = c(
      "Colwelliaceae" = "#CE2220",
      "Flavobacteriaceae" = "#E67F33",
      "Fokiniaceae" = "#7EB875",
      "Francisellaceae" = "#D0B541",
      "Oligoflexaceae" = "#57A2AC",
      "P13-46" = "#7C4942",
      "Hyphomonadaceae" = "#4E78C4",
      "Puniceicoccaceae" = "#B997C7", 
      "Rhodobacteraceae" = "#824D99",
      "Sphingomonadaceae" = "pink"
    ))
  ),
  matrix=(
    intersection_matrix(geom=geom_point(shape = "circle filled", size=3, stroke = 0.35, color = "gray20"))
    + scale_color_manual(
      values=c(
        #"T7" = "#650197",
        #"T3" = "#B21BFF",
        "T0" = "#B21BFF",
        "D" = "#AE0404",
        "H" = "#0FAB02",
        "late_pathogen" = "#FF1E1E",
        "early_pathogen" = "#FF9797",
        "continuous_pathogen" = "firebrick4"
      )
    )
  ),
  queries=list(
    #upset_query(set = "T7", fill = "#650197"),
    #upset_query(set = "T3", fill = "#B21BFF"), "#D98EFF"
    upset_query(set = "T0", fill = "#B21BFF"),
    upset_query(set = "D", fill = "#AE0404"),
    upset_query(set = "H", fill = "#0FAB02"),
    upset_query(set = "late_pathogen", fill = "#FF1E1E"),
    upset_query(set = "early_pathogen", fill = "#FF9797"),
    upset_query(set = "continuous_pathogen", fill = "firebrick4")
  ),
  name='ASVs', width_ratio=0.1, min_size = 0, sort_sets = FALSE,
  stripes = c("#F8EAFF", "#FFE4E4", "#E7FFE5", rep(c("gray91", "gray97"), 4))) +
  #stripes = c(rep(c("gray78", "gray87"), 2), "gray78", rep(c("gray97", "gray91"), 4))) +
  ggtitle("Putative Pathogen Candidates Only") 



opportunist_upset <- upset(
  simple_comp_upset %>% left_join(venn_all_times_and_doses, by = c("asv_id" = "OTU")) %>% mutate(Family = factor(Family)) %>%
    filter(early_opportunist | late_opportunist),
  colnames(simple_comp_upset %>% left_join(venn_all_times_and_doses, by = c("asv_id" = "OTU")) %>%
             select(T0, D, H, late_opportunist, early_opportunist)),
  base_annotations=list(
    'Intersection size'=intersection_size(counts=T, text = aes(size = 6.5),
                                          bar_number_threshold = 25,
                                          mapping=aes(fill=Family, col = Family, label = parse_number(asv_id)), col = "gray10"
    ) +
      geom_text(size = 4, position = position_stack(vjust = 0.5), col = "gray10") +
      scale_fill_manual(values = c(
        "Colwelliaceae" = "#CE2220",
        "Flavobacteriaceae" = "#E67F33",
        "Fokiniaceae" = "#7EB875",
        "Francisellaceae" = "#D0B541",
        "Oligoflexaceae" = "#57A2AC",
        "P13-46" = "#7C4942",
        "Hyphomonadaceae" = "#4E78C4",
        "Puniceicoccaceae" = "#B997C7", 
        "Rhodobacteraceae" = "#824D99",
        "Sphingomonadaceae" = "pink"
      ))
  ),
  
 
  
  
  matrix=(
    intersection_matrix(geom=geom_point(shape = "circle filled", size=3, stroke = 0.35, color = "gray20"))
    + scale_color_manual(
      values=c(
        #"T7" = "#650197",
        #"T3" = "#B21BFF",
        "T0" = "#B21BFF",
        "D" = "#AE0404",
        "H" = "#0FAB02",
        "late_opportunist" = "deepskyblue3",
        "early_opportunist" = "lightskyblue"
      )
    )
  ),
  queries=list(
    #upset_query(set = "T7", fill = "#650197"),
    #upset_query(set = "T3", fill = "#B21BFF"), "#D98EFF"
    upset_query(set = "T0", fill = "#B21BFF"),
    upset_query(set = "D", fill = "#AE0404"),
    upset_query(set = "H", fill = "#0FAB02"),
    upset_query(set = "late_opportunist", fill = "deepskyblue3"),
    upset_query(set = "early_opportunist", fill = "lightskyblue")
  ),
  name='ASVs', width_ratio=0.1, min_size = 0, sort_sets = FALSE,
  stripes = c("#F8EAFF", "#FFE4E4", "#E7FFE5", rep(c("gray91", "gray97"), 4))) +
  #stripes = c(rep(c("gray78", "gray87"), 2), "gray78", rep(c("gray97", "gray91"), 4))) +
  ggtitle("Putative Opportunist Candidates Only") 



wrap_plots(pathogen_upset, opportunist_upset) + 
  plot_layout(guides = 'collect')



#& plot_annotation(title = signatures)


#### JASON Work Zone ####
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
