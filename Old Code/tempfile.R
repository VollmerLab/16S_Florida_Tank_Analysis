tst <- normalized_asv_counts %>%
  nest_by(across(c('asv_id', Kingdom:Species))) %>%
  filter(asv_id == 'ASV_202') %>%
  mutate(fit_model(log2_cpm ~ treatment + (1 | genotype) + (1 | tank),
                   data, 
                   use_weights = FALSE),
         random_anova = list(rand(model)),
         process_model(model, re_model, random_anova),
         posthoc = list(emmeans(model, ~treatment) %>%
                          contrast(method = posthoc_categories$contrasts, adjust = 'none'))) 


emmeans(tst$model[[1]], ~treatment) %>%
  contrast(method = list('T7(DSvDR)' = c(0, 0, 0, 0, 0, 0, -1, 1, 0, 0),
                         'T7(DRvDS)' = c(0, 0, 0, 0, 0, 0, 1, -1, 0, 0)), 
           adjust = 'none', side = '>')



two_sided_tests <- tibble(microbial_signature = c('growth_comparisons'),
                          contrasts = list(bacterial_growth),
                          direction = '=') #=
right_tests <- tibble(microbial_signature = c('early_pathogen', 'continuous_pathogen', 'late_pathogen'),
                      contrasts = list(early_pathogen, continuous_pathogen, late_pathogen),
                      direction = '>') #>
left_tests <- tibble(microbial_signature = c('probiotic'),
                    contrasts = list(probiotic),
                    direction = '<') #<

posthoc_categories <- bind_rows(two_sided_tests, right_tests, left_tests) %>%
  unnest(contrasts) %>%
  mutate(contrast_name = names(contrasts)) %>%
  group_by(contrast_name, contrasts, direction) %>%
  summarise(signatures = list(c(microbial_signature)),
            .groups = 'drop') %>%
  nest(contrast = -direction)


asv_models <- normalized_asv_counts %>%
  nest_by(across(c('asv_id', Kingdom:Species))) %>%
  
  ungroup %>%
  sample_n(5)


blat <- asv_models %>%
  rowwise %>%
  
  partition(cluster) %>%
  mutate(fit_model(log2_cpm ~ treatment + (1 | genotype) + (1 | tank),
                   data, 
                   use_weights = FALSE),
         random_anova = list(rand(model)),
         process_model(model, re_model, random_anova),
         posthoc = list(run_posthoc(model, posthoc_categories))) %>%
  collect

blat$posthoc[[1]]

contrast_list <- posthoc_categories; model <- blat$model[[1]]




blat %>%
  select(-re_model:-eta2Partial_treatment) %>%
  expand_grid(posthoc_categories) %>%
  rowwise %>%
  partition(cluster) %>%
  mutate(posthoc = list(emmeans(model, ~treatment) %>%
                          contrast(method = contrast$contrasts, 
                                   adjust = 'none',
                                   side = direction))) %>%
  collect()

,
         posthoc = list(emmeans(model, ~treatment) %>%
                          contrast(method = posthoc_categories$contrasts, adjust = 'none'))) %>%
  collect() 
