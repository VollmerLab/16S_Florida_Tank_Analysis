the_plots12 <- bacterial_signature_asv %>%
  #filter(asv_id %in% c("ASV_202", "ASV_85", "ASV_142", "ASV_132", "ASV_439")) %>%
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
  #mutate(plot_info = list(rbind(plot_info, buffer_space))) %>%
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

the_plots12$plot[[1]]
the_plots12$combo_plots[[8]]
