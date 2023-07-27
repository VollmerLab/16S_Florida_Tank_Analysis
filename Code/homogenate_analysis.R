homogenate_data <- read_csv('../intermediate_files/homogenate_cpm.csv', show_col_types = FALSE)

library(patchwork)
library(emmeans)


  mutate(plot = list(ggplot(data = homogenate_pred,
                            aes(x = exposure, y = estimate, colour = exposure,
                                ymin = conf.low, ymax = conf.high)) +
                       geom_pointrange(position = position_dodge(0.5)) +
                       geom_jitter(data = data, inherit.aes = FALSE,
                                  aes(x = exposure, y = log2_cpm, colour = exposure), alpha = 0.6) +
                       theme_bw() +
                       scale_color_manual(values = c("#F75D5D", "#3DD8EA"), guide = "none") +
                       labs(title = asv_id))) #%>%
  
  
  pull(plot) %>%
  wrap_plots()
  
  
  homogenate_models <- homogenate_data %>%
    filter(asv_id %in% c('ASV_142', 'ASV_202', 'ASV_84', 'ASV_96', 'ASV_165', 'ASV_85',"ASV_132", "ASV_439")) %>% 
    nest_by(asv_id) %>%
    mutate(dose_pairwise = (emmeans(lm(log2_cpm ~ exposure, data = data),
                                          ~exposure) %>%
                                    contrast('pairwise', adjust = 'fdr') %>%
                                    broom::tidy(conf.int = TRUE) %>% pull(p.value))) %>%
    mutate(homogenate_pred = list(emmeans(lm(log2_cpm ~ exposure, data = data),
                                          ~exposure) %>%
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
  
  
  
  homogenate_models <- homogenate_data %>%
    #filter(asv_id %in% c('ASV_142', 'ASV_202', 'ASV_84', 'ASV_96', 'ASV_165', 'ASV_85',"ASV_132", "ASV_439")) %>% 
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
x  
  
