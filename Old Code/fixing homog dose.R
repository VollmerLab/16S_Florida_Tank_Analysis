homogenate_models_test <- homogenate_data %>%
  filter(asv_id %in% c("ASV_111", "ASV_117", "ASV_17", "ASV_175", "ASV_65")) %>%
  nest_by(asv_id) %>%
  mutate(model = list(lm(log2_cpm ~ exposure, data = data))) %>%
  mutate(homogenate_pred = list(model %>%
            broom::tidy(conf.int = TRUE) %>%
            mutate(term = case_when(term == "(Intercept)" ~ "D",
                                    term == "exposureH" ~ "H")) %>%
            mutate(estimate = ifelse(term == "H", estimate[term == "H"] + estimate[term == "D"], estimate),
                   std.error = NA,
                   statistic = ifelse(term == "H", statistic[term == "H"] + statistic[term == "D"], statistic),
                   conf.low = ifelse(term == "H", conf.low[term == "H"] + conf.low[term == "D"], conf.low),
                   conf.high = ifelse(term == "H", conf.high[term == "H"] + conf.high[term == "D"], conf.high),
                   p.value = ifelse(term == "D", p.value[term == "H"], p.value),
                   df = NA) %>%
            rename("exposure" = "term") %>%
            mutate(time = "Dose", 
                   graph_cat = "dose", 
                   c_time = ifelse(exposure == "D", -1.2, -1.5), 
                   final_disease_state = "na",
                   facet_lab = "Doses") %>%
            {. ->> set_one } %>% #save data to this var
            mutate(std.error = NA, df = NA, conf.low = NA, conf.high = NA, statistic = NA,
                   c_time = c(-1.8, -0.9), final_disease_state = NA, exposure = NA) %>% #dummy set to change size of Dose Facet
            rbind(set_one) %>%
            {. ->> set_two } %>%
            arrange(desc(estimate)) %>%
            dplyr::slice(1) %>%
            mutate(exposure = "p_val", c_time = -1.35, estimate = estimate + 3) %>%
            rbind(set_two)
            )) %>%
  mutate(homog_p_val = homogenate_pred %>% filter(exposure == "H") %>% pull(p.value)) %>%
  ungroup %>%
  mutate(adj_homog_p_val = p.adjust(homog_p_val, method = 'fdr')) %>%
  unnest(homogenate_pred) %>%
  mutate(p.value = adj_homog_p_val) %>%
  nest(homogenate_pred = -c(asv_id, data, model, homog_p_val, adj_homog_p_val))

homogenate_models_test$homogenate_pred[[1]]


#%>%
  mutate(p.valueee = dose_pairwise %>% filter(term == "exposureH") %>% pull(p.value))
  ungroup %>%
  mutate(dose_pairwise_adj = p.adjust(dose_pairwise, method = 'fdr')) %>%
  rowwise(asv_id) %>%
  mutate(homogenate_pred = list(model %>%
                                  broom::tidy(conf.int = TRUE) %>% mutate(time = "Dose", 
                                                                          graph_cat = "dose", c_time = ifelse(exposure == "D", -1.2, -1.5), final_disease_state = "na",
                                                                          facet_lab = "Doses") %>% 
                                  {. ->> set_one } %>% #save data to this var
                                  mutate(std.error = NA, df = NA, conf.low = NA, conf.high = NA, statistic = NA, p.value = dose_pairwise,
                                         c_time = c(-1.8, -0.9), final_disease_state = NA, exposure = NA) %>% #dummy set to change size of Dose Facet
                                  rbind(set_one) %>%
                                  {. ->> set_two } %>%
                                  arrange(desc(estimate)) %>%
                                  dplyr::slice(1) %>%
                                  mutate(exposure = "p_val", c_time = -1.35, estimate = estimate + 3) %>%
                                  rbind(set_two)
  ))



homogenate_models <- homogenate_data %>%
  filter(asv_id %in% c("ASV_111", "ASV_117", "ASV_17", "ASV_175")) %>%
  nest_by(asv_id) %>%
  mutate(em_model = list(emmeans(lm(log2_cpm ~ exposure, data = data),
                                 ~exposure))) %>%
  mutate(dose_pairwise = em_model %>%
           contrast('pairwise', adjust = 'none') %>%
           broom::tidy(conf.int = TRUE) %>% pull(p.value)) %>%
  ungroup %>%
  mutate(dose_pairwise_adj = p.adjust(dose_pairwise, method = 'fdr')) %>%
  rowwise(asv_id) %>%
  mutate(homogenate_pred = list(em_model %>%
                                  broom::tidy(conf.int = TRUE) %>% mutate(time = "Dose", 
                                                                          graph_cat = "dose", c_time = ifelse(exposure == "D", -1.2, -1.5), final_disease_state = "na",
                                                                          facet_lab = "Doses") %>% 
                                  {. ->> set_one } %>% #save data to this var
                                  mutate(std.error = NA, df = NA, conf.low = NA, conf.high = NA, statistic = NA, p.value = dose_pairwise,
                                         c_time = c(-1.8, -0.9), final_disease_state = NA, exposure = NA) %>% #dummy set to change size of Dose Facet
                                  rbind(set_one) %>%
                                  {. ->> set_two } %>%
                                  arrange(desc(estimate)) %>%
                                  dplyr::slice(1) %>%
                                  mutate(exposure = "p_val", c_time = -1.35, estimate = estimate + 3) %>%
                                  rbind(set_two)
                                
  )) #recombine them
