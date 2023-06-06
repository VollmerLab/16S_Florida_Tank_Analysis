##TODO - version 1 - figure out post-hoc for polynomial changes through time etc.
##TODO - version 3 - figure out splitting time_exposure into time and exposure individually

##TODO - general - function to reorder columns so fdr goes after p-value

#### Libraries ####
library(tidyverse)
library(magrittr)
library(lmerTest)
library(emmeans)
library(multidplyr)
library(ggupset)
library(qvalue)
library(patchwork)

cluster <- new_cluster(parallel::detectCores() - 1)
cluster_library(cluster, c('dplyr', 'lmerTest', 'emmeans'))

#### Functions ####

#### Data ####
normalized_asv_counts <- read_csv('../intermediate_files/fully_preprocessed_samples.csv.gz', show_col_types = FALSE) %>%
  mutate(time = factor(time, ordered = TRUE)) %>%
  mutate(time_exposure = str_c(time, exposure, sep = '_'),
         timeC = str_extract(time, '[0-9]+') %>% as.numeric) %>%
  mutate(asv_number = str_extract(asv_id, '[0-9]+') %>% as.integer) %>%
  filter(asv_number <= 200)

#
#### Plot raw data for an ASV ####
normalized_asv_counts %>%
  # filter(asv_id == 'ASV_3') %>%
  
  nest_by(across(c('asv_id', Kingdom:Species))) %>%
  ungroup %>%
  slice(12) %>%
  unnest(data) %>%
  
  ggplot(aes(x = time, y = log2_cpm, colour = susceptability, shape = exposure)) +
  geom_jitter() +
  stat_summary(fun.data = 'mean_se', position = position_dodge(1), size = 1.5)

#### Method 1 - Duplicate Time 0 to have identical healthy & diseased rows then model as 3 fully crossed factors ####

method_1 <- normalized_asv_counts %>%
  # filter(asv_id == 'ASV_3') %>%
  nest_by(across(c('asv_id', 'exposure', 'susceptability', Kingdom:Species))) %>%
  pivot_wider(names_from = 'exposure', values_from = 'data') %>%
  rowwise %>%
  mutate(D = list(bind_rows(D, `F`)),
         H = list(bind_rows(H, `F`))) %>%
  select(-`F`) %>%
  pivot_longer(cols = c('D', 'H'),
               names_to = 'exposure',
               values_to = 'data') %>%
  unnest(data) %>%
  nest_by(across(c('asv_id', Kingdom:Species))) %>%
  partition(cluster) %>%
  mutate(full_model = list(lmer(log2_cpm ~ time * exposure * susceptability + (1 | genotype) + (1 | tank), 
                                # weights = data$weight,
                                data = data, 
                                REML = TRUE,
                                control = variancePartition:::vpcontrol)),
         anova_table = list(anova(full_model, type = '3', ddf = 'Kenward-Roger'))) %>%
  collect %>%
  rowwise %>%
  mutate(as_tibble(anova_table, rownames = 'param') %>%
           rename(ss = `Sum Sq`,
                  ms = `Mean Sq`,
                  nDF = NumDF,
                  dDF = DenDF,
                  fvalue = `F value`,
                  pvalue = `Pr(>F)`) %>%
           pivot_wider(names_from = param,
                       values_from = c(ss, ms, nDF, dDF, fvalue, pvalue),
                       names_vary = 'slowest') %>%
           rename_with(~str_replace_all(., ':', 'X'))) %>%
  ungroup %>%
  mutate(across(starts_with('pvalue'), ~p.adjust(., method = 'fdr'), 
                .names = 'fdr_{.col}')) %>% 
  rename_with(~str_remove(., 'pvalue_'), .cols = starts_with('fdr'))

anova(method_1$full_model[[1]], type = '3', ddf = 'Kenward-Roger')
#sig effect of time and susceptibility - no interaction = no difference in how S/R behaves across time

emmeans(tmp$full_model[[1]], ~time * exposure * susceptability) %>%
  broom::tidy(conf.int = TRUE) %>%
  ggplot(aes(x = time, y = estimate, ymin = conf.low, ymax = conf.high,
             colour = susceptability, shape = exposure)) +
  geom_pointrange(position = position_dodge(0.5))


method_1_upset <- method_1 %>%
  select(asv_id, starts_with('fdr')) %>%
  select(-fdr_time) %>%
  mutate(across(starts_with('fdr'), ~. < 0.05)) %>%
  
  pivot_longer(cols = -asv_id,
               names_to = c('term'),
               values_to = 'significance') %>%
  filter(significance) %>%
  mutate(term = str_remove(term, 'fdr_')) %>%
  group_by(asv_id) %>%
  summarise(terms = list(term),
            .groups = 'drop') %>%
  
  ggplot(aes(x = terms)) +
  geom_bar() +
  scale_x_upset() +
  theme_classic() +
  theme_combmatrix(combmatrix.label.make_space = TRUE) +
  labs(title = 'Method 1')

#
#### Method 2 - individual time curves for each exposure/resistance combination ####
#This method tests if there are changes through time for each exposure/resistance combination & if those changes significantly differe across exposure/resistance.
#It does not test if there are absolute differences in numbers between exposure/resistance combinations

## Step 1. Model change in time separately for each treatment
## Step 2. Extract time coefficients (with se) for each model 
## Step 3. Model if those coefficients differ by exposure & treatment

postprocess_time_model <- function(model){
  aov_tab <- anova(model, type = '3', ddf = 'Kenward-Roger')
  aov_row <- as_tibble(aov_tab, rownames = 'term') %>%
    rename(ss = 'Sum Sq',
           ms = 'Mean Sq',
           n.DF = NumDF,
           d.DF = DenDF,
           f.value = 'F value',
           pvalue = 'Pr(>F)') %>%
    rowwise %>%
    mutate(eta2Partial = effectsize::F_to_eta2(f.value, n.DF, d.DF, ci = NULL)$Eta2_partial) %>%
    tidyr::pivot_wider(names_from = term,
                       values_from = where(is.numeric),
                       names_vary = 'slowest') %>%
    rename_with(~str_replace_all(., ':', 'X'))
  
  coefficient_row <- bind_cols(bind_rows(fixef(model)) %>%
                                 rename_with(~str_c('mean_', .)) %>%
                                 rename_with(~str_remove_all(., '\\(|\\)')),
                               vcov(model) %>% diag %>% sqrt %>% bind_rows %>%
                                 rename_with(~str_c('se_', .)) %>%
                                 rename_with(~str_remove_all(., '\\(|\\)')))
  
  bind_cols(aov_row, coefficient_row)
}

tmp <- normalized_asv_counts %>%
  filter(asv_id == 'ASV_3') %>%
  nest_by(across(c('asv_id', 'exposure', 'susceptability', Kingdom:Species))) %>%
  pivot_wider(names_from = 'exposure', values_from = 'data') %>%
  rowwise %>%
  mutate(D = list(bind_rows(D, `F`)),
         H = list(bind_rows(H, `F`))) %>%
  select(-`F`) %>%
  pivot_longer(cols = c('D', 'H'),
               names_to = 'exposure',
               values_to = 'data') %>%
  rowwise %>%
  mutate(time_model = list(lmer(log2_cpm ~ time + (1 | genotype) + (1 | tank), 
                                weights = data$weight, #something kinda weird happening with weight inclusion here...possibly because weights assume exposure/resistance also in model
                                data = data, 
                                REML = TRUE,
                                control = variancePartition:::vpcontrol)),
         postprocess_time_model(time_model))

tmp %>%
  select(susceptability, exposure, data) %>%
  unnest(data) %>%
  
  ggplot(aes(x = time, y = log2_cpm)) +
  geom_jitter() +
  stat_summary(fun.data = mean_se, size = 2, colour = 'red') +
  facet_wrap(~susceptability + exposure)

tmp %>%
  select(susceptability, exposure, starts_with('mean'), starts_with('se')) %>%
  pivot_longer(cols = c(starts_with('mean'), starts_with('se')),
               names_to = c('.value', 'param'),
               names_pattern = '(.*)_(.*)') %>%
  ggplot(aes(x = susceptability, y = mean, ymin = mean - se, ymax = mean + se,
             colour = exposure)) +
  geom_pointrange(position = position_dodge(0.5)) +
  facet_wrap(~param)


library(emmeans)
emmeans(tmp$time_model[[1]], ~time) %>%
  contrast('poly')


tst_data <- tmp %>%
  select(susceptability, exposure, starts_with('mean'), starts_with('se')) 


tst_model <- brm(mean_Intercept | se(se_Intercept) ~ susceptability * exposure,
    data = tst_data,
    backend = 'cmdstanr')

emmeans(tst_model, ~ susceptability * exposure) %>%
  broom::tidy() %>%
  ggplot(aes(x = susceptability, y = estimate, ymin = lower.HPD, ymax = upper.HPD,
             colour = exposure)) +
  geom_pointrange(position = position_dodge(0.25)) +
  geom_pointrange(data = tst_data,
                  aes(y = mean_Intercept, ymin = mean_Intercept - se_Intercept, ymax = mean_Intercept + se_Intercept),
                  shape = 'square', linetype = 'dashed',
                  position = position_dodge(0.5))
emmeans(tst_model, ~ susceptability * exposure) %>%
  contrast('pairwise')


library(metafor)
tst_data
tst_model <- rma.mv(yi = mean_Intercept,
       V = se_Intercept,
       mods = ~ susceptability * exposure,
       data = tst_data,
       test = 'z')
anova(tst_model)
anova(tst_model, btt = 'susceptability')
summary(tst_model)
bind_cols(tst_data,
          as_tibble(predict(tst_model))) %>%
  
  ggplot(aes(x = susceptability, y = pred, ymin = ci.lb, ymax = ci.ub,
             colour = exposure)) +
  geom_pointrange(position = position_dodge(0.25)) +
  geom_pointrange(aes(y = mean_Intercept, ymin = mean_Intercept - se_Intercept, ymax = mean_Intercept + se_Intercept),
                  shape = 'square', linetype = 'dashed',
                  position = position_dodge(0.5))

#### Method 3 - concatenate time and exposure into single variable cross with resistance ####
posthoc_comparisons <- list('linear' = c(-1/2, 0, 0, 1/4, 1/4, -1/2, 0, 0, 1/4, 1/4),
                            'quadratic' = c(1/2, -2/4, -2/4, 1/4, 1/4, 1/2, -2/4, -2/4, 1/4, 1/4),
                            'T3 - T0' = c(-1/2, 1/4, 1/4, 0, 0, -1/2, 1/4, 1/4, 0, 0),
                            'T7 - T3' = c(0, -1/4, -1/4, 1/4, 1/4, 0, -1/4, -1/4, 1/4, 1/4),
                            
                            'R - S' = c(1/5, 1/5, 1/5, 1/5, 1/5, -1/5, -1/5, -1/5, -1/5, -1/5),
                            'linear R - S' = c(-1, 0, 0, 1/2, 1/2, 0, 0, 0, 0, 0) - c(0, 0, 0, 0, 0, -1, 0, 0, 1/2, 1/2),
                            'quadratic R - S' = c(1, -2/2, -2/2, 1/2, 1/2, 0, 0, 0, 0, 0) - c(0, 0, 0, 0, 0, 1, -2/2, -2/2, 1/2, 1/2),
                            'R(T3-T0) - S(T3-T0)' = c(-1, 1/2, 1/2, 0, 0, 0, 0, 0, 0, 0) - c(0, 0, 0, 0, 0, -1, 1/2, 1/2, 0, 0),
                            'R(T7-T3) - S(T7-T3)' = c(0, -1/2, -1/2, 1/2, 1/2, 0, 0, 0, 0, 0) - c(0, 0, 0, 0, 0, 0, -1/2, -1/2, 1/2, 1/2),
                            
                           
                            'D - H' = c(0, 1/4, -1/4, 1/4, -1/4, 0, 1/4, -1/4, 1/4, -1/4),
                            'linear D - H' = c(-1/2, 0, 0, 1/2, 0, -1/2, 0, 0, 1/2, 0) - c(-1/2, 0, 0, 0, 1/2, -1/2, 0, 0, 0, 1/2),
                            'quadratic D - H' = c(1/2, -2/2, 0, 1/2, 0, 1/2, -2/2, 0, 1/2, 0) - c(1/2, 0, -2/2, 0, 1/2, 1/2, 0, -2/2, 0, 1/2),
                            'D(T3-T0) - H(T3-T0)' = c(-1/2, 1/2, 0, 0, 0, -1/2, 1/2, 0, 0, 0) - c(-1/2, 0, 1/2, 0, 0, -1/2, 0, 1/2, 0, 0),
                            'D(T7-T3) - H(T7-T3)' = c(0, -1/2, 0, 1/2, 0, 0, -1/2, 0, 1/2, 0) - c(0, 0, -1/2, 0, 1/2, 0, 0, -1/2, 0, 1/2),
                            
                            'RH(T3-T0) - RD(T3-T0)' = c(-1, 0, 1, 0, 0, 0, 0, 0, 0, 0) - c(-1, 1, 0, 0, 0, 0, 0, 0, 0, 0),
                            'RH(T3-T0) - SH(T3-T0)' = c(-1, 0, 1, 0, 0, 0, 0, 0, 0, 0) - c(0, 0, 0, 0, 0, -1, 0, 1, 0, 0),
                            'RH(T3-T0) - SD(T3-T0)' = c(-1, 0, 1, 0, 0, 0, 0, 0, 0, 0) - c(0, 0, 0, 0, 0, -1, 1, 0, 0, 0),
                            'RD(T3-T0) - SH(T3-T0)' = c(-1, 1, 0, 0, 0, 0, 0, 0, 0, 0) - c(0, 0, 0, 0, 0, -1, 0, 1, 0, 0),
                            'RD(T3-T0) - SD(T3-T0)' = c(-1, 1, 0, 0, 0, 0, 0, 0, 0, 0) - c(0, 0, 0, 0, 0, -1, 1, 0, 0, 0),
                            'SH(T3-T0) - SD(T3-T0)' = c(0, 0, 0, 0, 0, -1, 0, 1, 0, 0) - c(0, 0, 0, 0, 0, -1, 1, 0, 0, 0),
                            
                            'RH(T7-T3) - RD(T7-T3)' = c(0, 0, -1, 0, 1, 0, 0, 0, 0, 0) - c(0, -1, 0, 1, 0, 0, 0, 0, 0, 0),
                            'RH(T7-T3) - SH(T7-T3)' = c(0, 0, -1, 0, 1, 0, 0, 0, 0, 0) - c(0, 0, 0, 0, 0, 0, 0, -1, 0, 1),
                            'RH(T7-T3) - SD(T7-T3)' = c(0, 0, -1, 0, 1, 0, 0, 0, 0, 0) - c(0, 0, 0, 0, 0, 0, -1, 0, 1, 0),
                            'RD(T7-T3) - SH(T7-T3)' = c(0, -1, 0, 1, 0, 0, 0, 0, 0, 0) - c(0, 0, 0, 0, 0, 0, 0, -1, 0, 1),
                            'RD(T7-T3) - SD(T7-T3)' = c(0, -1, 0, 1, 0, 0, 0, 0, 0, 0) - c(0, 0, 0, 0, 0, 0, -1, 0, 1, 0),
                            'SH(T7-T3) - SD(T7-T3)' = c(0, 0, 0, 0, 0, 0, 0, -1, 0, 1) - c(0, 0, 0, 0, 0, 0, -1, 0, 1, 0),
                            
                            'DT7(S - R)' = c(0, 0, 0, -1, 0, 0, 0, 0, 1, 0),
                            'ST7(D - H)' = c(0, 0, 0, 0, 0, 0, 0, 0, 1, -1),
                            
                            'R(T7-T3)' = c(0, -1/2, -1/2, 1/2, 1/2, 0, 0, 0, 0, 0),
                            'S(T7-T3)' = c(0, 0, 0, 0, 0, 0, -1/2, -1/2, 1/2, 1/2))

cluster_copy(cluster, c('posthoc_comparisons'))

# ANOVA <- method_3$anova_table[[1]]; posthoc <- method_3$posthoc[[1]]
post_process <- function(ANOVA, posthoc){
  aov_row <- as_tibble(ANOVA, rownames = 'param') %>%
    mutate(param = str_replace(param, 'time_exposure', 'timeXexposure')) %>%
    rename(ss = `Sum Sq`,
           ms = `Mean Sq`,
           nDF = NumDF,
           dDF = DenDF,
           fvalue = `F value`,
           pvalue = `Pr(>F)`) %>%
    pivot_wider(names_from = param,
                values_from = c(ss, ms, nDF, dDF, fvalue, pvalue),
                names_vary = 'slowest') %>%
    rename_with(~str_replace_all(., ':', 'X'))
  
  post_row <- as_tibble(posthoc) %>%
    # filter(contrast %in% c('linear', 'quadratic', 'D - H')) %>%
    # mutate(contrast = str_replace_all(contrast, c('linear' = 'timeL',
    #                                               'quadratic' = 'timeQ',
    #                                               'D - H' = 'exposure'))) %>%
    # select(-estimate, -SE) %>%
    rename(tvalue = t.ratio,
           pvalue = p.value) %>%
    pivot_wider(names_from = c('contrast'),
                values_from = c('estimate', 'SE', 'df', 'tvalue', 'pvalue'))
  bind_cols(aov_row, post_row)
}

method_3 <- normalized_asv_counts %>%
  # filter(asv_id == 'ASV_3') %>%
  
  nest_by(across(c('asv_id', Kingdom:Species))) %>%
  partition(cluster) %>%
  mutate(full_model = list(lmer(log2_cpm ~ time_exposure * susceptability + (1 | genotype) + (1 | tank), 
                                # weights = data$weight,
                                data = data, 
                                REML = TRUE,
                                control = variancePartition:::vpcontrol)),
         anova_table = list(anova(full_model, type = '3', ddf = 'Kenward-Roger')),
         posthoc = list(emmeans(full_model, ~time_exposure * susceptability) %>%
                          contrast(method = posthoc_comparisons, adjust = 'none'))) %>%
  collect %>%
  rowwise %>%
  mutate(post_process(anova_table, posthoc)) %>%
  ungroup %>%
  mutate(across(starts_with('pvalue'), ~p.adjust(., method = 'fdr'), 
                .names = 'fdr_{.col}')) %>% 
  rename_with(~str_remove(., 'pvalue_'), .cols = starts_with('fdr'))

method_3_upset <- method_3 %>%
  select(asv_id, starts_with('fdr')) %>% 
  select(-fdr_linear, -fdr_quadratic, -'fdr_T3 - T0', -'fdr_T7 - T3') %>%
  select(-contains('linear'), -contains('quadratic')) %>%
  # # select(asv_id, fdr_timeXexposure, fdr_susceptability, fdr_timeXexposureXsusceptability,
  # #        fdr_linear, fdr_quadratic, `fdr_D - H`) %>%
  # # slice(12) %>%
  # select(-contains('time')) %>%
  mutate(across(starts_with('fdr'), ~. < 0.05)) %>%
  
  pivot_longer(cols = -asv_id,
               names_to = c('term'),
               values_to = 'significance') %>%
  filter(significance) %>%
  mutate(term = str_remove(term, 'fdr_')) %>%
  filter(!term %in% c('D - H', 'R - S', 'susceptability')) %>%
  mutate(term = str_replace_all(term, c('D\\(T3-T0\\) - H\\(T3-T0\\)' = 'Early Exposure',
                                        'D\\(T7-T3\\) - H\\(T7-T3\\)' = 'Late Exposure',
                                        'R\\(T3-T0\\) - S\\(T3-T0\\)' = 'Early Resistance',
                                        'R\\(T7-T3\\) - S\\(T7-T3\\)' = 'Late Resistance'))) %>%
  group_by(asv_id) %>%
  summarise(terms = list(term),
            .groups = 'drop') %>%
  
  ggplot(aes(x = terms)) +
  geom_bar() +
  scale_x_upset() +
  theme_classic() +
  theme_combmatrix(combmatrix.label.make_space = TRUE) +
  labs(title = 'Method 3')
method_3_upset

method_3 %>%
  select(asv_id, starts_with('fdr')) %>% 
  select(-fdr_linear, -fdr_quadratic, -'fdr_T3 - T0', -'fdr_T7 - T3') %>%
  select(-contains('linear'), -contains('quadratic')) %>%
  # select(asv_id, fdr_timeXexposure, fdr_susceptability, fdr_timeXexposureXsusceptability,
  #        fdr_linear, fdr_quadratic, `fdr_D - H`) %>%
  # slice(12) %>%
  select(-contains('time')) %>%
  mutate(across(starts_with('fdr'), ~. < 0.05)) %>%
  pivot_longer(cols = -asv_id,
               names_to = c('term'),
               values_to = 'significance') %>%
  filter(significance)  %>%
  left_join(select(normalized_asv_counts, asv_id, Kingdom:Species) %>%
              distinct,
            by = 'asv_id') %>%
  group_by(across(c(asv_id, Kingdom:Species))) %>%
  summarise(term = str_c(term, collapse = ';;'),
            .groups = 'drop') %>% 
  count(term, Family, Genus) %>%
  filter(term != 'fdr_D(T7-T3) - H(T7-T3)') %>%
  arrange(-n)


tmp <- method_3 %>% 
  # filter(`fdr_R(T7-T3) - S(T7-T3)` < 0.05 | `fdr_R(T3-T0) - S(T3-T0)` < 0.05) %>%
  select(asv_id, starts_with('fdr')) %>%
  mutate(across(starts_with('fdr'), ~. < 0.05)) %>%
  pivot_longer(cols = starts_with('fdr'),
               names_to = c('term'),
               values_to = 'significance') %>% 
  filter(significance) %>%
  filter(term %in% c('fdr_R(T7-T3) - S(T7-T3)', 'fdr_R(T3-T0) - S(T3-T0)')) %>%
  group_by(asv_id) %>%
  summarise(term = str_c(term, collapse = '; ')) %>%
  left_join(method_3, by = 'asv_id') %>%
  mutate(across(asv_id:Species, ~str_replace_na(., 'NA'))) %>%
  rowwise %>%
  mutate(plot = list(full_model %>%
                       emmeans(~time_exposure * susceptability) %>%
                       broom::tidy(conf.int = TRUE) %>%
                       separate(time_exposure, into = c('time', 'exposure')) %>%
                       ggplot(aes(x = time, y = estimate, ymin = conf.low, ymax = conf.high,
                                  colour = susceptability, shape = exposure)) +
                       geom_pointrange(position = position_dodge(0.5)) +
                       labs(title = str_c(Family, Genus, asv_id, sep = ':'),
                            subtitle = term)))

tmp$plot[[7]]
select(tmp, asv_id, term, Genus, Family, 
       'estimate_R(T7-T3) - S(T7-T3)',
       'estimate_R(T3-T0) - S(T3-T0)')
colnames(tmp)



tmp <- method_3 %>% 
  # filter(`fdr_R(T7-T3) - S(T7-T3)` < 0.05 | `fdr_R(T3-T0) - S(T3-T0)` < 0.05) %>%
  select(asv_id, starts_with('fdr')) %>%
  mutate(across(starts_with('fdr'), ~. < 0.05)) %>%
  select(-contains('linear'), -contains('quadratic')) %>%
  select(asv_id, contains('H'), contains('R')) %>% 
  select(-`fdr_D - H`, -fdr_timeXexposure, -fdr_susceptability, -fdr_timeXexposureXsusceptability,
         -`fdr_T3 - T0`, -`fdr_T7 - T3`, -`fdr_R - S`) %>%
  pivot_longer(cols = starts_with('fdr'),
               names_to = c('term'),
               values_to = 'significance') %>% 
  filter(significance) %>%
  group_by(asv_id) %>%
  summarise(term = str_c(term, collapse = '; '),
            n_term = n(),
            .groups = 'drop') %>%
  filter(n_term > 1) %>%
  
  left_join(method_3, by = 'asv_id') %>%
  mutate(across(asv_id:Species, ~str_replace_na(., 'NA'))) %>%
  rowwise %>%
  mutate(plot = list(full_model %>%
                       emmeans(~time_exposure * susceptability) %>%
                       broom::tidy(conf.int = TRUE) %>%
                       separate(time_exposure, into = c('time', 'exposure')) %>%
                       ggplot(aes(x = time, y = estimate, ymin = conf.low, ymax = conf.high,
                                  colour = susceptability, shape = exposure)) +
                       geom_pointrange(position = position_dodge(0.5)) +
                       labs(title = str_c(Family, Genus, asv_id, sep = ':'),
                            subtitle = term)))
select(tmp, asv_id, term, Genus, Family)


plotting_emmeans <- function(model, sig_terms){
  sig_terms <- c('exposure', 'susceptability')
  if(length(sig_terms) == 1){
    if(sig_terms == 'exposure'){
      em_form <- ~time_exposure
      
      out <- emmeans(model, em_form) %>%
        broom::tidy(conf.int = TRUE) %>%
        separate(time_exposure, into = c('time', 'exposure'))
      
    } else if(sig_terms == 'susceptability'){
      em_form <- ~time_exposure * susceptability
      
      out <- emmeans(model, em_form) %>%
        contrast(method = list('T0_R' = c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                               'T0_S' = c(0, 0, 0, 0, 0, 1, 0, 0, 0, 0),
                               'T3_R' = c(0, 1/2, 1/2, 0, 0, 0, 0, 0, 0, 0),
                               'T3_S' = c(0, 0, 0, 0, 0, 0, 1/2, 1/2, 0, 0),
                               'T7_R' = c(0, 0, 0, 1/2, 1/2, 0, 0, 0, 0, 0),
                               'T7_S' = c(0, 0, 0, 0, 0, 0, 0, 0, 1/2, 1/2)),
                 adjust = 'none') %>%
        broom::tidy(conf.int = TRUE) %>%
        select(-term) %>%
        separate(contrast, into = c('time', 'susceptability'))
    }
  } else {
    em_form <- ~time_exposure * susceptability
    out <- emmeans(model, em_form) %>%
      broom::tidy(conf.int = TRUE) %>%
      separate(time_exposure, into = c('time', 'exposure'))
  }
  out
}

make_plot <- function(plot_data, sig_terms, TITLE, SUBTITLE){
  sig_terms <- c('exposure', 'susceptability')
  if(length(sig_terms) == 1){
    if(sig_terms == 'exposure'){
      out <- ggplot(plot_data, 
                    aes(x = time, y = estimate,
                        ymin = conf.low, ymax = conf.high,
                        colour = exposure)) +
        scale_colour_manual(values = c('F' = 'black', 'D' = 'red', 'H' = 'blue'))
      
    } else if(sig_terms == 'susceptability'){
      out <- ggplot(plot_data, 
                    aes(x = time, y = estimate,
                        ymin = conf.low, ymax = conf.high,
                        shape = susceptability))
    }
  } else {
    out <- ggplot(plot_data, 
                  aes(x = time, y = estimate,
                      ymin = conf.low, ymax = conf.high,
                      colour = exposure, shape = susceptability)) +
      scale_colour_manual(values = c('F' = 'black', 'D' = 'red', 'H' = 'blue'))
  }
  out +
    geom_pointrange(position = position_dodge(0.5)) +
    labs(x = NULL,
         y = 'log2_cpm',
         title = TITLE,
         subtitle = SUBTITLE) +
    theme_classic() +
    theme(axis.text = element_text(colour = 'black', size = 12),
          axis.title = element_text(colour = 'black', size = 16),
          panel.background = element_rect(colour = 'black'))
}

tmp <- method_3 %>%
  # select(asv_id, contains('(T7-T3)'), contains('(T3-T0)')) %>%
  select(asv_id, starts_with('fdr')) %>%
  pivot_longer(cols = starts_with('fdr'),
               names_to = 'term',
               values_to = 'fdr',
               names_prefix = 'fdr_') %>%
  filter(fdr < 0.05) %>%
  group_by(asv_id) %>%
  summarise(term = str_c(term, collapse = '; '),
            n_term = n(),
            .groups = 'drop') %>%
  mutate(em_form = case_when(str_detect(term, 'D\\(') & str_detect(term, 'R\\(') ~ list(c('exposure', 'susceptability')),
                             str_detect(term, 'D\\(') ~ list('exposure'),
                             str_detect(term, 'R\\(') ~ list('susceptability'))) %>% 
  left_join(select(method_3, asv_id:Species, data, full_model, starts_with('estimate')),
            by = 'asv_id') %>%
  rowwise %>%
  mutate(plot_data = list(plotting_emmeans(full_model, em_form))) %>%
  ungroup %>%
  mutate(across(Kingdom:Species, ~str_replace_na(.))) %>%
  mutate(term = str_replace_all(term, c('D\\(T3-T0\\) - H\\(T3-T0\\)' = 'Early Exposure',
                                        'D\\(T7-T3\\) - H\\(T7-T3\\)' = 'Late Exposure',
                                        'R\\(T3-T0\\) - S\\(T3-T0\\)' = 'Early Resistance',
                                        'R\\(T7-T3\\) - S\\(T7-T3\\)' = 'Late Resistance'))) %>%
  rowwise %>%
  mutate(plot = list(make_plot(plot_data, em_form, 
                               str_c(Family, Genus, asv_id, sep = ': '),
                               term))) %>%
  ungroup
  
method_3 %>%
  filter(asv_id == 'ASV_85') %>%
  select(starts_with('fdr')) %>%
  pivot_longer(cols = everything()) %>%
  filter(value < 0.05)



tmp %>%
  # filter(str_detect(term, 'RD\\(T7-T3\\) - SD\\(T7-T3\\)')) %>%
  filter(str_detect(term, 'DT7\\(S - R\\)')) %>%
  filter(`estimate_DT7(S - R)` > 0) %>% 
  filter(str_detect(term, 'timeXexposure')) %>%
  filter(str_detect(term, 'ST7\\(D - H\\)'))
  # filter(str_detect(term, 'T7\\(DS - H\\)')) %>%
  # filter(str_detect(term, 'T7 - T3')) %>%
  # filter(`estimate_T7 - T3` > 0) %>%
  # filter(str_detect(term, 'S\\(T7-T3\\)')) %>%
  # filter(`estimate_S(T7-T3)` > 0) %>%
  pull(plot) %>%
  wrap_plots() +
  plot_layout(guides = 'collect')

tmp %>%
  filter(str_detect(term, 'timeXexposureXsusceptability')) %>% 
  # filter(Family == 'Colwelliaceae') %>%
  # pull(term)
  rowwise %>%
  mutate(keep = is.null(em_form),
         keep = keep | length(em_form) == 2) %>%
  filter(keep) %>%
  pull(plot) %>%
  wrap_plots()


tmp %>%
  mutate(across(c(Family, Genus), str_replace_na)) %>%
  count(Family, Genus, term) %>%
  pivot_wider(names_from = 'term',
              values_from = 'n',
              values_fill = 0L) %>% 
  rowwise %>%
  mutate(total_sig = sum(c_across(where(is.integer)))) %>%
  ungroup %>%
  left_join(method_3 %>%
              mutate(across(c(Family, Genus), str_replace_na)) %>%
              count(Family, Genus) %>%
              rename(total = n),
            by = c('Family', 'Genus')) %>%
  mutate(pct_sig = total_sig / total) %>%
  arrange(-pct_sig) %>% View

blah <- tmp %>%
  filter(Genus == 'Thalassotalea') %>%
  summarise(plots = list(wrap_plots(plot) +
                           plot_layout(guides = 'collect') &
                           plot_annotation(title = Genus)),
            n = n())
blah$plots[[1]]

tmp %>%
  filter(term == 'Late Exposure') %>%
  ggplot(aes(x = `estimate_D(T7-T3) - H(T7-T3)`)) +
  geom_histogram(bins = 30) 

tmp %>%
  filter(term == 'Late Exposure') %>%
  filter(abs(`estimate_D(T7-T3) - H(T7-T3)`) > 3) %>%
  summarise(plots = list(wrap_plots(plot) +
                           plot_layout(guides = 'collect')),
            n = n()) %>%
  pull(plots) %>%
  pluck(1)

tmp %>%
  group_by(term) %>%
  summarise(n = n())

tst <- tmp %>%
  filter(term != 'Late Exposure') %>%
  group_by(term) %>%
  summarise(plots = list(wrap_plots(plot) +
                           plot_layout(guides = 'collect') &
                           plot_annotation(title = term)),
            n = n())
tst$plots[[6]]

method_3 %>%
  select(asv_id, starts_with('fdr')) %>% 
  select(-fdr_linear, -fdr_quadratic, -'fdr_T3 - T0', -'fdr_T7 - T3') %>%
  select(-contains('linear'), -contains('quadratic')) %>%
  # select(asv_id, fdr_timeXexposure, fdr_susceptability, fdr_timeXexposureXsusceptability,
  #        fdr_linear, fdr_quadratic, `fdr_D - H`) %>%
  # slice(12) %>%
  select(-contains('time')) %>%
  mutate(across(starts_with('fdr'), ~. < 0.05)) %>%
  
  pivot_longer(cols = -asv_id,
               names_to = c('term'),
               values_to = 'significance') %>%
  filter(significance) %>%
  group_by(asv_id) %>%
  summarise(terms = list(term),
            .groups = 'drop') %>%
  
  ggplot(aes(x = terms)) +
  geom_bar() +
  scale_x_mergelist(sep = '@') +
  theme_classic() +
  axis_combmatrix(sep = "@")
  theme_combmatrix(combmatrix.label.make_space = TRUE) +
  labs(title = 'Method 3')

emmeans(method_3$full_model[[1]], ~time_exposure * susceptability) %>%
  broom::tidy(conf.int = TRUE) %>%
  separate(time_exposure, into = c('time', 'exposure'), remove = FALSE) %>% 
  mutate(time_exposure = factor(time_exposure, levels = c('T7_D', 'T3_D', 'T0_F', 'T3_H', 'T7_H'))) %>%
  arrange(time_exposure) %>%
  ggplot(aes(x = time, y = estimate, ymin = conf.low, ymax = conf.high,
             colour = susceptability)) +
  # geom_line(data = . %>% filter(exposure != 'H'),
  #           aes(group = interaction(susceptability)),
  #           position = position_dodge(0.5)) +
  # geom_line(data = . %>% filter(exposure != 'D'),
  #           aes(group = interaction(susceptability)),
  #           position = position_dodge(0.5)) +
  geom_path(aes(group = interaction(susceptability, exposure)),
            position = position_dodge(0.5)) +
  geom_pointrange(aes(shape = exposure), position = position_dodge(0.5))

tmp_ref <- ref_grid(tmp$full_model[[1]])

#time - polynomial T0, T3, T7 average across time
tst_ref <- add_grouping(tmp_ref, newname = 'time', refname = 'time_exposure', newlevs = c('T0', 'T3', 'T3', 'T7', 'T7'))
emmeans(tst_ref, ~time) %>% contrast('poly')
emmeans(tst_ref, ~time) %>% contrast('consec')
emmeans(tst_ref, ~susceptability) %>% contrast('pairwise')

emmeans(tst_ref, ~time | susceptability) %>% 
  contrast('poly')
emmeans(tst_ref, ~time | susceptability) %>%
  contrast('consec')

emmeans:::poly.emmc(1:3)
emmeans:::pairwise.emmc(1:2)
emmeans(tst_ref, ~time * susceptability) %>%
  #first 3 numbers are coding for polynomial in R second 3 numbers are coding for polynomial in S
  contrast(method = list('linear R-S' = c(-1, 0, 1, 0, 0, 0) - c(0, 0, 0, -1, 0, 1),
                         'quadratic R-S' = c(1, -2, 1, 0, 0, 0) - c(0, 0, 0, 1, -2, 1)))


emmeans:::consec.emmc(1:3)
emmeans(tst_ref, ~time * susceptability) %>%
  #first 3 numbers are coding for polynomial in R second 3 numbers are coding for polynomial in S
  contrast(method = list('T3-T0 R-S' = c(-1, 1, 0, 0, 0, 0) - c(0, 0, 0, -1, 1, 0),
                         'T7-T3 R-S' = c(0, -1, 1, 0, 0, 0) - c(0, 0, 0, 0, -1, 1)))


compare_coefs <- function(model){
  #https://stats.stackexchange.com/questions/93540/testing-equality-of-coefficients-from-two-different-regressions
  blah <- emmeans(tst_ref, ~time | susceptability) %>% 
    contrast('poly') %>%
    broom::tidy() %>%
    group_split(susceptability)
  
  full_join(blah[[1]] %>%
              select(contrast, estimate, std.error) %>%
              rename_with(.cols = where(is.numeric), ~str_c(., '_R')),
            
            blah[[2]] %>%
              select(contrast, estimate, std.error) %>%
              rename_with(.cols = where(is.numeric), ~str_c(., '_S')),
            by = 'contrast') %>%
    mutate(diff = estimate_R - estimate_S,
           Z = (estimate_R - estimate_S) / sqrt(std.error_R^2 + std.error_S^2),
           p = 2 * pnorm(-abs(Z)))
}
compare_coefs(tst_ref)
         
         
#exposure - H/D average across T3 & T7
tst_ref2 <- add_grouping(tmp_ref, newname = 'exposure', refname = 'time_exposure', newlevs = c(NA, 'D', 'H', 'D', 'H'))
emmeans(tst_ref2, ~exposure) %>% contrast('pairwise')

emmeans:::poly.emmc(1:3)


#All contrasts 
emmeans(method_3$full_model[[1]], ~time_exposure * susceptability) %>%
  contrast(method = posthoc_comparisons, 
           adjust = 'none')



#### Method 4 - Time as a Random Effect ####
tmp <- normalized_asv_counts %>%
  filter(asv_id == 'ASV_3') %>%
  
  nest_by(across(c('asv_id', Kingdom:Species))) %>%
  mutate(full_model = list(lmer(log2_cpm ~ exposure * susceptability + (1 | time) + (1 | genotype) + (1 | tank), 
                                weights = data$weight,
                                data = data, 
                                REML = TRUE,
                                control = variancePartition:::vpcontrol)))
anova(tmp$full_model[[1]], type = '3', ddf = 'Kenward-Roger')
summary(tmp$full_model[[1]])
rand(tmp$full_model[[1]]) #works inside function

library(bootpredictlme4)

predict.merMod
pred_dat <- expand_grid(time = c('T0', 'T3', 'T7'),
            susceptability = c('R', 'S'),
            exposure = c('F', 'H', 'D')) %>%
  filter(!(time != 'T0' & exposure == 'F')) %>%
  filter(!(time == 'T0' & exposure != 'F'))

predict_lmer <- function(model, SIM){
  preds <- predict(model, nsim = SIM, 
                   newdata = pred_dat, se.fit = TRUE,
                   re.form = ~(1 | time))
  
  
  bind_cols(pred_dat, preds[1:3] %>% bind_cols) %>%
    bind_cols(t(preds[[4]])) %>%
    rename(lower = '2.5%',
           upper = '97.5%')
}

predict_lmer(tmp$full_model[[1]], 1000) %>%
  
  ggplot(aes(x = time, y = fit, ymin = lower, ymax = upper,
             colour = susceptability, shape = exposure)) +
  geom_pointrange(position = position_dodge(0.5))

emmeans(tmp$full_model[[1]], ~susceptability)

#Method 5 - Model Differences between T0 & T3 and T3 & T7 as mini time-series 
normalized_asv_counts %>%
  filter(asv_id == 'ASV_3') %>%
  nest_by(across(c('asv_id', 'exposure', 'susceptability', Kingdom:Species))) %>%
  pivot_wider(names_from = 'exposure', values_from = 'data') %>%
  rowwise %>%
  mutate(D = list(bind_rows(D, `F`)),
         H = list(bind_rows(H, `F`))) %>%
  select(-`F`) %>%
  pivot_longer(cols = c('D', 'H'),
               names_to = 'exposure',
               values_to = 'data') %>%
  unnest(data) %>%
  select(-reads, -time_c, -time_exposure, -sample_id, -weight) %>%
  group_by(genotype) %>%
  
  pivot_wider(names_from = time,
              values_from = log2_cpm) %>% 
  filter(any(!is.na(T0))) %>%
  filter(genotype == 'M5')
