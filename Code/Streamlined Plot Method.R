#Streamlined version of Jason's Plot-Based Significance Model

setwd("~/Desktop/Screenshots/Career/Vollmer Lab/GitHub/16S_Florida_Tank_Analysis/Code")
rerun_model <- TRUE

#### Packages ####
library(lme4)
library(lmerTest)
library(emmeans)
library(qvalue)
library(afex)
library(multidplyr)
library(patchwork)
library(tidyverse)

select <- dplyr::select

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

#### Read in Data ####

microbe_data <- read_csv('../intermediate_files/pre_model_data.csv',
                         show_col_types = FALSE) %>%
  filter(!tank %in% c('HOMO', 'homogenate_fragment')) %>%
  mutate(time = factor(time, levels = c('T0', 'T3', 'T7')))

#### First Iteration - Exposure ####
control_options <- variancePartition:::vpcontrol
# control_options$optimizer <- 'bobyqa'

microbe_data %>%
  mutate(time = factor(time, ordered = TRUE)) %>%
  group_by(asv_names, exposure) %>%
  mutate(total = n(), .after = time) %>%
  filter(value < 7.13) %>%
  reframe(zer = n(), total = total) %>%
  mutate(prop = zer/total) %>%
  distinct() %>%
  ungroup() %>%
  group_by(asv_names) %>%
  filter(prop > 0.90, exposure != "Field") %>%
    pull(asv_names) -> many_zeros

if(file.exists("../intermediate_files/plot_model_iter1.rds") & !rerun_model){
  iter1 <- read_rds("../intermediate_files/plot_model_iter1.rds")
} else {
  iter1 <- microbe_data %>%
    # group_by(asv_names) %>%
    # filter(n_distinct(sample_id[value > 7.13]) > 24) %>%
    mutate(time = factor(time, ordered = TRUE)) %>%
    
    #add tiny random number to deal with constant values in a sample
    #error was in 'ASV_1634': emmeans -> emm_basis -> pbkrtest::vcovAdj.lmerMod -> 
    #pbkrtest:::vcovAdj_internal -> forceSymmetric(2 * solve(IE2))
    
    mutate(value = value + rnorm(nrow(.), 0.0001, 0.001)) %>%
    nest(data = -c(asv_names)) %>%
    rowwise %>%
    #partition(cluster) %>%
    mutate(healthy_regression = list(lmer(value ~ time + (1 | genotype) + (0 + dummy(time, 'T0') | tank),
                                          data = filter(data, exposure != 'D'),
                                          control = control_options)),
           disease_regression = list(lmer(value ~ time + (1 | genotype) + (0 + dummy(time, 'T0') | tank),
                                          data = filter(data, exposure != 'H'),
                                          control = control_options))) %>%
    mutate(coef_comp = list(compare_coefs(healthy_regression, disease_regression))) %>%
    #collect %>%
    ungroup
  write_rds(iter1, "../intermediate_files/plot_model_iter1.rds")
}


iter1_pvals <- iter1 %>%
  select(asv_names, coef_comp) %>%
  unnest(coef_comp) %>%
  group_by(contrast) %>%
  mutate(fdr = p.adjust(p, 'fdr'),
         qvalue = qvalue(p)$qvalues) %>%
  ungroup

iter1_pre_plot <- iter1_pvals %>%
  group_by(asv_names) %>%
  filter(any(fdr < 0.05)) %>%
  ungroup %>%
  mutate(is_significant = fdr < 0.05) %>%
  select(asv_names, contrast, is_significant) %>%
  pivot_wider(names_from = contrast,
              values_from = is_significant) %>%
  left_join(iter1,
            by = 'asv_names') %>%
  rowwise(asv_names, linear, quadratic) %>%
  summarise(make_ts_plot_data(healthy_regression, disease_regression),
            .groups = 'drop') %>%
  nest(data = -c(linear, quadratic))


iter1_model <- iter1_pvals %>%
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

iter1_plots <- iter1_pre_plot %>%
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
iter1_plots$plot[[3]]
#### Second Iteration - Final Disease State ####

iter2 <- iter1_pvals %>%
  group_by(asv_names) %>%
  filter(any(fdr < 0.05)) %>%
  select(asv_names) %>%
  ungroup %>%
  distinct %>%
  left_join(iter1,
            by = 'asv_names') %>%
  rowwise %>%
  mutate(final_outcome_model = list(mixed(value ~ time * final_disease_state + 
                                            (0 + dummy(time, 'T0') | fragment_id),
                                          data = filter(data, exposure != 'H', !is.na(final_disease_state)),
                                          method = 'KR',
                                          control = variancePartition:::vpcontrol))) 

write_rds(iter2, "../intermediate_files/plot_model_iter2.rds")

iter2_model <- iter2 %>%
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

iter2_model %>%
  ggplot(aes(x = time, y = estimate, ymin = conf.low, ymax = conf.high,
             colour = final_disease_state)) +
  geom_pointrange(position = position_dodge(0.5)) +
  facet_wrap(Order + Family + Genus ~ asv_names, scales = 'free_y')

#### Plot Individual ASVs ####

#Iteration 1
iter1_indiv_plots <- iter1_pre_plot %>%
  #mutate(sig_aspects = paste(ifelse(linear, "L", "n"), ifelse(quadratic, "Q", "n"), sep = "")) %>%
  unnest(data) %>%
  ungroup() %>%
  nest_by(asv_names) %>%
  rowwise() %>%
  mutate(first_plot = list(ggplot(data = data) +
                             geom_pointrange(aes(x = time, y = estimate, ymin = conf.low, ymax = conf.high,
                                                 colour = exposure), position = position_dodge(0.5)) +
                             labs(title = paste("First - ", asv_names, " (linear: ", data$linear, ", quadratic: ", data$quadratic , ")", sep = ""))))
#Iteration 2
iter2_indiv_plots <- iter2_model %>%
  ungroup() %>%
  nest_by(asv_names) %>%
  rowwise() %>%
  mutate(second_plot = list(ggplot(data = data, aes(x = time, y = estimate, ymin = conf.low, ymax = conf.high,
                                                    colour = final_disease_state)) +
                              geom_pointrange(position = position_dodge(0.5)) +
                              labs(title = paste("Second - ", asv_names, sep = ""))))



jason_model_plots <- iter1_indiv_plots %>%
  select(asv_names, first_plot) %>%
  full_join(iter2_indiv_plots %>% select(asv_names, second_plot))

write_rds(jason_model_plots, "../intermediate_files/comp_models_plots_1_2.rds")

#### Significant ASVs ####

plot_method_asvs <- iter1_plots %>% 
  select(-plot) %>%
  unnest(cols = c(data)) %>%
  mutate(terms = paste(ifelse(linear, "L", "n"), ifelse(quadratic, "Q", "n"), sep = ""), .after = asv_names) %>%
  #select(asv_names, terms) %>%
  select(asv_names, linear, quadratic, terms) %>%
  distinct()

write_rds(plot_method_asvs, "../intermediate_files/plot_method_asvs.rds")

write_rds(list(iter1_model, iter2_model, iter2), "../intermediate_files/both_jasons_models.rds")

#### ####

