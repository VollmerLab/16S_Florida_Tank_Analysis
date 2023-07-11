homogenate_data <- read_csv('../intermediate_files/homogenate_cpm.csv', show_col_types = FALSE)

library(patchwork)
library(emmeans)
homogenate_models <- homogenate_data %>%
  filter(asv_id %in% c('ASV_3', 'ASV_202', 'ASV_84', 'ASV_96', 'ASV_165', 'ASV_85')) %>% 
  nest_by(asv_id) %>%
  mutate(homogenate_pred = list(emmeans(lm(log2_cpm ~ exposure, data = data),
                                        ~exposure) %>%
                                  broom::tidy(conf.int = TRUE) %>% mutate(time = "Dose", 
                                  graph_cat = "dose", c_time = ifelse(exposure == "D", -1.5, -1.3), susceptability = "na",
                                  facet_lab = "Doses"))) %>%
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
  
  
  homogenate_models <- homogenate_models %>% rename(plot = "dose_plot")
lim_homogenate_models <- homogenate_models %>% select(asv_id, homogenate_pred) 
