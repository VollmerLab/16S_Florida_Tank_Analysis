
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

rick_plots <- asv_models %>%
  filter(Order == "Rickettsiales") %>%
  group_by(Genus) %>%
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
  group_by(Genus) %>%
  summarise(combo_plots = list(wrap_plots(plot) + plot_layout(guides = 'collect') & plot_annotation(title = Genus)))

#view the plots
rick_plots$combo_plots[[2]]



#any sig Rickettsias?

bacterial_signature_asv %>% filter(asv_id %in% (significant_models %>% 
                                   filter(Order == "Rickettsiales") %>% pull(asv_id)))
#one crasher
#one early pathogen






