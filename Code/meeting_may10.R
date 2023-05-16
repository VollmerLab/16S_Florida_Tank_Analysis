library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(qvalue)
library(afex)
library(multidplyr)
library(patchwork)

#### Functions ####
compare_coefs <- function(healthy_regression, disease_regression){
  #https://stats.stackexchange.com/questions/93540/testing-equality-of-coefficients-from-two-different-regressions
  healthy_poly <- emmeans(healthy_regression, ~time) %>%
    contrast('poly')
  
  disease_poly <- emmeans(disease_regression, ~time) %>%
    contrast('poly')
  
  full_join(broom::tidy(healthy_poly) %>%
              select(contrast, estimate, std.error) %>%
              rename_with(.cols = where(is.numeric), ~str_c(., '_healthy')),
            
            broom::tidy(disease_poly) %>%
              select(contrast, estimate, std.error) %>%
              rename_with(.cols = where(is.numeric), ~str_c(., '_disease')),
            by = 'contrast') %>%
    mutate(Z = (estimate_healthy - estimate_disease) / sqrt(std.error_healthy^2 + std.error_disease^2),
           p = 2 * pnorm(-abs(Z)))
}

make_ts_plot_data <- function(healthy_regression, disease_regression){
  bind_rows(healthy = broom::tidy(emmeans(healthy_regression, ~ time), conf.int = TRUE),
            disease = broom::tidy(emmeans(disease_regression, ~ time), conf.int = TRUE),
            .id = 'exposure') 
}

clean_afex <- function(model){
  #model is an afex::mixed model
  model$anova_table %>%
    as_tibble(rownames = 'param') %>%
    janitor::clean_names() %>%
    rename_with(~str_replace_all(., '_df', 'DF')) %>%
    rename(pvalue = pr_f) %>%
    mutate(param = str_replace(param, ':', 'X'),
           param = str_replace(param, 'final_disease_state', 'finalDisease')) %>%
    pivot_wider(names_from = 'param',
                values_from = where(is.numeric),
                names_vary = 'slowest')
}


#### Data ####
microbe_data <- read_csv('C:/Users/jdsel/Documents/Google Drive/Research/Vollmer Lab PostDoc/Emily_Microbiome/intermediate_files/pre_model_data.csv',
                         show_col_types = FALSE) %>%
  filter(!tank %in% c('HOMO', 'homogenate_fragment')) %>%
  mutate(time = factor(time, levels = c('T0', 'T3', 'T7')))

#for emily
microbe_data <- read_csv('../intermediate_files/pre_model_data.csv',
                         show_col_types = FALSE) %>%
  filter(!tank %in% c('HOMO', 'homogenate_fragment')) %>%
  mutate(time = factor(time, levels = c('T0', 'T3', 'T7')))

count(microbe_data, exposure)

microbe_data %>%
  select(sample_id:final_disease_state) %>%
  distinct() %>%
  count(time, exposure, final_disease_state)

microbe_data %>%
  select(sample_id:final_disease_state) %>%
  distinct() %>%
  filter(exposure == 'H', final_disease_state == 'D')

healthy_regression <- lmer(value ~ time + (1 | genotype) + (1 | tank),
                           data = filter(microbe_data, asv_names == 'ASV_1', exposure != 'D'))

disease_regression <- lmer(value ~ time + (1 | genotype) + (1 | tank),
                           data = filter(microbe_data, asv_names == 'ASV_1', exposure != 'H'))

compare_coefs(healthy_regression, disease_regression)


microbe_data %>%
  filter(tank != 'homogenate_fragment')

tst <- microbe_data %>%
  nest(data = -c(asv_names)) 


cluster <- new_cluster(parallel::detectCores() - 1)
cluster_library(cluster, c('dplyr', 'stringr', 'lme4', 'emmeans'))
cluster_copy(cluster, c('compare_coefs'))


tst <- microbe_data %>%
  nest(data = -c(asv_names)) %>%
  # sample_n(10) %>%
  rowwise %>%
  #partition(cluster) %>%
  mutate(healthy_regression = list(lmer(value ~ time + (1 | genotype) + (1 | tank),
                                        data = filter(data, exposure != 'D'),
                                        control = variancePartition:::vpcontrol)),
         disease_regression = list(lmer(value ~ time + (1 | genotype) + (1 | tank),
                                        data = filter(data, exposure != 'H'),
                                        control = variancePartition:::vpcontrol))) %>%
  mutate(coef_comp = list(compare_coefs(healthy_regression, disease_regression))) %>%
  #collect %>%
  ungroup 


tst_p <- tst %>%
  select(asv_names, coef_comp) %>%
  unnest(coef_comp) %>%
  group_by(contrast) %>%
  mutate(fdr = p.adjust(p, 'fdr'),
         qvalue = qvalue(p)$qvalues) %>%
  ungroup

tst_pre_plot <- tst_p %>%
  group_by(asv_names) %>%
  filter(any(fdr < 0.05)) %>%
  ungroup %>%
  mutate(is_significant = fdr < 0.05) %>%
  select(asv_names, contrast, is_significant) %>%
  pivot_wider(names_from = contrast,
              values_from = is_significant) %>%
  left_join(tst,
            by = 'asv_names') %>%
  rowwise(asv_names, linear, quadratic) %>%
  summarise(make_ts_plot_data(healthy_regression, disease_regression),
            .groups = 'drop') %>%
  nest(data = -c(linear, quadratic))
  

first_iter_model <- tst_p %>%
  group_by(asv_names) %>%
  filter(any(fdr < 0.05)) %>%
  ungroup %>%
  mutate(is_significant = fdr < 0.05) %>%
  select(asv_names, contrast, is_significant) %>%
  pivot_wider(names_from = contrast,
              values_from = is_significant) %>%
  left_join(select(microbe_data, asv_names, Order:Species) %>%
              distinct,
            by = 'asv_names') #%>% View



tst_pre_plot$data[[1]] %>%
  ggplot(aes(x = time, y = estimate, ymin = conf.low, ymax = conf.high,
             colour = exposure)) +
  geom_pointrange(position = position_dodge(0.5)) +
  facet_wrap(~asv_names)

blah <- tst_pre_plot %>%
  rowwise %>%
  mutate(data = list(data %>% 
                       left_join(select(microbe_data, asv_names, Order:Species) %>%
                                   distinct,
                                 by = 'asv_names')),
         plot = list(ggplot(data, aes(x = time, y = estimate, ymin = conf.low, ymax = conf.high,
                                      colour = exposure)) +
                       geom_pointrange(position = position_dodge(0.5)) +
                       facet_wrap(Order + Family + Genus ~ asv_names, scales = 'free_y') +
                       labs(title = str_c('linear: ', linear, '; quadratic: ', quadratic)))) %>%
  ungroup 

blah$plot[[3]]

#%>%
  summarise(plot = list(wrap_plots(plot))) %>%
  pull(plot) %>%
  pluck(1)
  
tst_p2$data[[1]] %>%
  filter(genotype == 'K1')

tst_p2 <- tst_p %>%
  group_by(asv_names) %>%
  #filter(any(fdr < 0.05)) %>%
  select(asv_names) %>%
  ungroup %>%
  distinct %>%
  left_join(tst,
            by = 'asv_names') %>%
  rowwise %>%
  mutate(final_outcome_model = list(mixed(value ~ time * final_disease_state + (1 | genotype) + 
                                            (1 | tank),
                                          data = filter(data, exposure == 'D'),
                                          method = 'KR',
                                          control = variancePartition:::vpcontrol))) 

second_iter_model <- tst_p2 %>%
  rowwise(asv_names, data, ends_with('regression'), ends_with('model')) %>%
  summarise(clean_afex(final_outcome_model),
            .groups = 'drop') %>%
  mutate(across(starts_with('pvalue'), p.adjust, method = 'fdr')) %>%
  mutate(across(starts_with('pvalue'), ~. < 0.05)) %>%
  select(-ends_with('time')) %>%
  filter(if_any(.cols = starts_with('pvalue'))) %>%
  rowwise(asv_names) %>%
  summarise(emmeans(final_outcome_model, ~time * final_disease_state) %>%
              broom::tidy(conf.int = TRUE)) %>%
  left_join(select(microbe_data, asv_names, Order:Species) %>%
              distinct,
            by = 'asv_names') 

second_iter_model %>%
  ggplot(aes(x = time, y = estimate, ymin = conf.low, ymax = conf.high,
             colour = final_disease_state)) +
  geom_pointrange(position = position_dodge(0.5)) +
  facet_wrap(Order + Family + Genus ~ asv_names, scales = 'free_y')


#### Significant ASVs ####

plot_method_asvs <- blah %>% select(-plot) %>%
  unnest(cols = c(data)) %>%
  mutate(terms = paste(ifelse(linear, "L", "n"), ifelse(quadratic, "Q", "n"), sep = ""), .after = asv_names) %>%
  #select(asv_names, terms) %>%
  select(asv_names, linear, quadratic, terms) %>%
  distinct()

write_rds(plot_method_asvs, "../intermediate_files/plot_method_asvs.rds")



first_iter_model

second_iter_model

write_rds(list(first_iter_model, second_iter_model, tst_p2), "../intermediate_files/both_jasons_models.rds")

#### ####


  