asvs_for_logfold <- bacterial_signature_asv %>%
  #filter(asv_id == "ASV_65") %>%
  group_by(asv_id) %>%
  reframe(signatures = str_c(signatures, collapse = ', ')) %>%
  mutate(grouped_signatures = case_when(str_detect(signatures, "Aquaria") ~ "Tank Effect",
                                        str_detect(signatures, "Field") ~ "Tank Effect",
                                        str_detect(signatures, "DiseaseOutcome") & str_detect(signatures, "DD_early") & str_detect(signatures, "DD_vs_DH_early") & asv_id %in% sig_homog_asv_list ~ "Putative Early Pathogens",
                                        str_detect(signatures, "DiseaseOutcome") & str_detect(signatures, "DD_late") & str_detect(signatures, "DD_vs_DH_late") & asv_id %in% sig_homog_asv_list ~ "Pathogens",
                                        str_detect(signatures, "DiseaseOutcome") & str_detect(signatures, "DD_late") & str_detect(signatures, "DD_vs_DH_late") & !asv_id %in% sig_homog_asv_list~ "White Band Disease-Associated Opportunists",
                                        str_detect(signatures, "DiseaseOutcome") ~ "Unlikely Pathogens",
                                        str_detect(signatures, "HealthyOutcome") ~ "Healthy-Associated",
                                        TRUE ~ signatures)) %>%
  #don't bother plotting tank effect, etc...
  filter(grouped_signatures %in% c("Putative Early Pathogens", "Pathogens", "White Band Disease-Associated Opportunists", "Healthy-Associated")) %>%
  group_by(asv_id, grouped_signatures) %>%
  distinct() %>%
  inner_join(significant_models,
             by = 'asv_id') %>%
  mutate(across(Family:Species, ~ifelse(str_detect(., "NA"), NA, .))) %>%
  mutate(taxonomy_for_graph = case_when(is.na(Order) ~ Class,
                                        is.na(Family) ~ Order,
                                        is.na(Genus) ~ Family,
                                        is.na(Species) ~ str_c(Family, Genus, sep = " "),
                                        TRUE ~ str_c(Family, Species, sep = " "))) %>%
  mutate(formatted_ASV = str_c("(ASV ", parse_number(asv_id), ")", sep = "")) %>%
  mutate(prep = ifelse(is.na(Family),
                       list(substitute(taxonomy_for_graph~italic(species_abbrev)~formatted_ASV, list(taxonomy_for_graph = taxonomy_for_graph, species_abbrev = "sp.", formatted_ASV = formatted_ASV))),
                       ifelse(is.na(Species),
                              list(substitute(italic(taxonomy_for_graph)~italic(species_abbrev)~formatted_ASV, list(taxonomy_for_graph = taxonomy_for_graph, species_abbrev = "sp.", formatted_ASV = formatted_ASV))),
                              list(substitute(italic(taxonomy_for_graph)~formatted_ASV, list(taxonomy_for_graph = taxonomy_for_graph, formatted_ASV = formatted_ASV))))
                       
  ))

#### BAD, T7-based t-tests ####
  
significant_models %>%
  unnest(data)
  select(data, asv_id, starts_with('fdr')) %>% 
  select(-contains(c('treatment', 'tank', 'genotype'))) %>%
  mutate(across(starts_with('fdr'), ~. < 0.05)) %>%
  
  pivot_longer(cols = -c(asv_id, data),
               names_to = c('term'),
               values_to = 'significance') %>%
  mutate(term = str_remove(term, 'fdr_')) %>%
  
  #Select ASVs which fit all characteristics of any given bacterial signature
  mutate(direction = str_extract(term, '[><=]'),
         term = str_remove(term, '_[><=]')) %>%
  filter(asv_id %in% asvs_for_logfold$asv_id, term %in% c("DD_vs_DH_late", "outcome")) #%>%
  
  
  look <- significant_models %>%
    filter(asv_id %in% asvs_for_logfold$asv_id) %>%
    unnest(data) %>%
    filter(time == "T7") %>%
    group_by(asv_id) %>%
    summarise_at(vars(log2_cpm), list(outcome_t_test = ~ list(t.test(. ~ final_disease_state))))
  
# dose
  
logfold_dose_1 <- homogenate_data %>%
  filter(asv_id %in% asvs_for_logfold$asv_id) %>%
  group_by(asv_id) %>%
  filter(!asv_id %in% c("ASV_301", "ASV_372", "ASV_458")) %>% #data are essentially constant, can't run t-test
  summarise_at(vars(log2_cpm), list(dose_t_test = ~ list(t.test(. ~ exposure)))) %>%
  rowwise() %>%
  mutate(t_stderr = dose_t_test$stderr) %>%
  mutate(dose_t_test = broom::tidy(dose_t_test)) %>%
  unnest(dose_t_test) %>%
  rename("D" = "estimate1", "H" = "estimate2", "dose_D_minus_H" = "estimate") %>%
  select(-c(parameter, conf.low, conf.high, method, alternative)) %>%
  mutate(fdr_dose = p.adjust(p.value)) %>%
  select(asv_id, dose_D_minus_H, fdr_dose, t_stderr) %>%
  pivot_longer(dose_D_minus_H, names_to = "comparison", values_to = "estimate") %>%
  rename("fdr_pval" = "fdr_dose")

logfold_dose <- logfold_dose_1 %>% ungroup() %>% 
  dplyr::slice(1:3) %>% mutate(across(everything(), ~NA)) %>% 
  mutate(asv_id = c("ASV_301", "ASV_372", "ASV_458"), fdr_pval = c(1, 1, 1), 
         comparison = c("dose_D_minus_H", "dose_D_minus_H", "dose_D_minus_H"), estimate = c(0, 0, 0), t_stderr = c(0, 0, 0)) %>%
  rbind(logfold_dose_1)
  
#outcome

logfold_outcome <- significant_models %>%
  filter(asv_id %in% asvs_for_logfold$asv_id) %>%
  unnest(data) %>%
  filter(time == "T7") %>%
  group_by(asv_id) %>%
  summarise_at(vars(log2_cpm), list(outcome_t_test = ~ list(t.test(. ~ final_disease_state)))) %>%
  rowwise() %>%
  mutate(t_stderr = outcome_t_test$stderr) %>%
  mutate(outcome_t_test = broom::tidy(outcome_t_test)) %>%
  unnest(outcome_t_test) %>%
  rename("D" = "estimate1", "H" = "estimate2", "outcome_D_minus_H" = "estimate") %>%
  select(-c(parameter, conf.low, conf.high, method, alternative)) %>%
  mutate(fdr_outcome = p.adjust(p.value)) %>%
  select(asv_id, outcome_D_minus_H, fdr_outcome, t_stderr) %>%
  pivot_longer(outcome_D_minus_H, names_to = "comparison", values_to = "estimate") %>%
  rename("fdr_pval" = "fdr_outcome")
  
#exposure

logfold_exposure <- significant_models %>%
  filter(asv_id %in% asvs_for_logfold$asv_id) %>%
  unnest(data) %>%
  filter(time == "T7") %>%
  group_by(asv_id) %>%
  summarise_at(vars(log2_cpm), list(exposure_t_test = ~ list(t.test(. ~ exposure)))) %>%
  rowwise() %>%
  mutate(t_stderr = exposure_t_test$stderr) %>%
  mutate(exposure_t_test = broom::tidy(exposure_t_test)) %>%
  unnest(exposure_t_test) %>%
  rename("D" = "estimate1", "H" = "estimate2", "exposure_D_minus_H" = "estimate") %>%
  select(-c(parameter, conf.low, conf.high, method, alternative)) %>%
  mutate(fdr_exposure = p.adjust(p.value)) %>%
  select(asv_id, exposure_D_minus_H, fdr_exposure, t_stderr) %>%
  pivot_longer(exposure_D_minus_H, names_to = "comparison", values_to = "estimate") %>%
  rename("fdr_pval" = "fdr_exposure")

#exposure on outcome

logfold_DD_DH <- significant_models %>%
  filter(asv_id %in% asvs_for_logfold$asv_id) %>%
  unnest(data) %>%
  filter(treatment %in% c("T7_D_D", "T7_D_H")) %>%
  group_by(asv_id) %>%
  summarise_at(vars(log2_cpm), list(DD_DH_t_test = ~ list(t.test(. ~ treatment)))) %>%
  rowwise() %>%
  mutate(t_stderr = DD_DH_t_test$stderr) %>%
  mutate(DD_DH_t_test = broom::tidy(DD_DH_t_test)) %>%
  unnest(DD_DH_t_test) %>%
  rename("D_D" = "estimate1", "D_H" = "estimate2", "DD_minus_DH" = "estimate") %>%
  select(-c(parameter, conf.low, conf.high, method, alternative)) %>%
  mutate(fdr_DD_DH = p.adjust(p.value)) %>%
  select(asv_id, DD_minus_DH, fdr_DD_DH, t_stderr) %>%
  pivot_longer(DD_minus_DH, names_to = "comparison", values_to = "estimate") %>%
  rename("fdr_pval" = "fdr_DD_DH")
   
#combine
for_plot <- logfold_outcome %>%
  rbind(logfold_DD_DH, logfold_dose) %>% #logfold_exposure
  mutate(graph_pval = str_c(comparison, ifelse(fdr_pval < 0.05, "sig", "nonsig"), sep = "_")) %>%
  left_join(asvs_for_logfold %>% select(asv_id, grouped_signatures), by = join_by(asv_id)) %>%
  left_join(the_plots1 %>% ungroup() %>% select(asv_id, prep), by = join_by(asv_id)) %>%
  mutate(grouped_signatures = factor(grouped_signatures, levels = c("Putative Late Pathogens",
                                            "Specialized Opportunists", "Healthy Associated"))) %>%
  arrange(grouped_signatures) %>%
  mutate(asv_id = factor(asv_id), ordered = TRUE)
  #left_join(significant_models %>% select(asv_id, Domain:Species), by = join_by(asv_id)) %>%



logfold_change_plot <- ggplot(for_plot) +
  geom_hline(yintercept = 0) +
  geom_pointrange(aes(x = asv_id, y = estimate, fill = graph_pval, col = comparison, pch = comparison, ymin = estimate - t_stderr, ymax = estimate + t_stderr), 
             size = 1, stroke = 1, position = position_dodge(width = 0.7)) + #position = position_dodge(width = 1)
  scale_shape_manual(values = c(24, 22, 21, 23)) +
  scale_color_manual(values = c("outcome_D_minus_H" = "#A70000",
                                "exposure_D_minus_H" = "#4722B5",
                                "dose_D_minus_H" = "#DF369B", 
                                "DD_minus_DH" = "#EC8E00")) +
  scale_fill_manual(values = c("outcome_D_minus_H_nonsig" = "white", 
                               "outcome_D_minus_H_sig" = "#A70000",
                               "exposure_D_minus_H_nonsig" = "white", 
                               "exposure_D_minus_H_sig" = "#4722B5",
                               "dose_D_minus_H_nonsig" = "white", 
                               "dose_D_minus_H_sig" = "#DF369B",
                               "DD_minus_DH_nonsig" = "white", 
                               "DD_minus_DH_sig" = "#EC8E00"), guide = "none") +
  coord_flip() +
  theme_bw() +
    theme(strip.background = element_blank(),
          strip.text.y = element_blank(),
          legend.position = "none") +
  scale_x_discrete(breaks = for_plot$asv_id, labels = for_plot$prep) +
  facet_grid(rows = vars(grouped_signatures), space = "free", scales = "free") +
  xlab(NULL) +
  ylab("Logfold Difference (D-H)")


  

  
relayer_logfold_change <- ggplot(for_plot, aes(x = asv_id, y = estimate, ymin = estimate - t_stderr, ymax = estimate + t_stderr)) +
    geom_hline(yintercept = 0) +
    #outcome main effect
    (geom_point(data = (for_plot %>% filter(comparison == "outcome_D_minus_H")), aes(colour1 = graph_pval, fill1 = graph_pval), pch = 23, size = 3, stroke = 1) %>%
       rename_geom_aes(new_aes = c("colour" = "colour1", "fill" = "fill1"))) + 
    scale_color_manual(aesthetics = "colour1", values = c("#A70000", "#A70000"), guide = "legend", 
                       name = "Outcome Main Effect", breaks = c("outcome_D_minus_H_sig", "outcome_D_minus_H_nonsig"), labels = c("Significant", "Non-Significant")) +
    scale_fill_manual(aesthetics = "fill1", values = c("#A70000", "transparent"), guide = "legend", 
                       name = "Outcome Main Effect", breaks = c("outcome_D_minus_H_sig", "outcome_D_minus_H_nonsig"), labels = c("Significant", "Non-Significant")) +
    #exposure main effect
    (geom_point(data = (for_plot %>% filter(comparison == "exposure_D_minus_H")), aes(colour2 = graph_pval, fill2 = graph_pval), pch = 21, size = 3, stroke = 1) %>%
       rename_geom_aes(new_aes = c("colour" = "colour2", "fill" = "fill2"))) + 
    scale_color_manual(aesthetics = "colour2", values = c("#4722B5", "#4722B5"), guide = "legend", 
                       name = "Exposure Main Effect", breaks = c("exposure_D_minus_H_sig", "exposure_D_minus_H_nonsig"), labels = c("Significant", "Non-Significant")) +
    scale_fill_manual(aesthetics = "fill2", values = c("#4722B5", "transparent"), guide = "legend", 
                      name = "Exposure Main Effect", breaks = c("exposure_D_minus_H_sig", "exposure_D_minus_H_nonsig"), labels = c("Significant", "Non-Significant")) +
    #Dose main effect
    (geom_point(data = (for_plot %>% filter(comparison == "dose_D_minus_H")), aes(colour3 = graph_pval, fill3 = graph_pval), pch = 22, size = 3, stroke = 1) %>%
      rename_geom_aes(new_aes = c("colour" = "colour3", "fill" = "fill3"))) + 
    scale_color_manual(aesthetics = "colour3", values = c("#DF369B", "#DF369B"), guide = "legend", 
                       name = "Homogenate Doses (D vs. H)", breaks = c("dose_D_minus_H_sig", "dose_D_minus_H_nonsig"), labels = c("Significant", "Non-Significant")) +
    scale_fill_manual(aesthetics = "fill3", values = c("#DF369B", "transparent"), guide = "legend", 
                      name = "Homogenate Doses (D vs. H)", breaks = c("dose_D_minus_H_sig", "dose_D_minus_H_nonsig"), labels = c("Significant", "Non-Significant")) +
    #exposure main effect
    (geom_point(data = (for_plot %>% filter(comparison == "DD_minus_DH")), aes(colour4 = graph_pval, fill4 = graph_pval), pch = 24, size = 3, stroke = 1) %>%
      rename_geom_aes(new_aes = c("colour" = "colour4", "fill" = "fill4"))) + 
    scale_color_manual(aesthetics = "colour4", values = c("#EC8E00", "#EC8E00"), guide = "legend", 
                       name = "DD_minus_DH", breaks = c("DD_minus_DH_sig", "DD_minus_DH_nonsig"), labels = c("Significant", "Non-Significant")) +
    scale_fill_manual(aesthetics = "fill4", values = c("#EC8E00", "transparent"), guide = "legend", 
                      name = "DD_minus_DH", breaks = c("DD_minus_DH_sig", "DD_minus_DH_nonsig"), labels = c("Significant", "Non-Significant")) +
    coord_flip() +
    theme_bw() +
    theme(strip.background = element_blank(),
          strip.text.y = element_blank()) +
    scale_x_discrete(breaks = for_plot$asv_id, labels = for_plot$prep) +
    facet_grid(rows = vars(grouped_signatures), space = "free", scales = "free") +
    xlab(NULL) +
    ylab("Logfold Difference (D-H)") +
    guides(#outcome
           color1 = guide_legend(order = 2),
           fill1 = guide_legend(order = 2),
           #exposure
           color2 = guide_legend(order = 1),
           fill2 = guide_legend(order = 1),
           #dose
           color3 = guide_legend(order = 4),
           fill3 = guide_legend(order = 4),
           #DD_DH
           color4 = guide_legend(order = 3),
           fill4 = guide_legend(order = 3))


relayer_logfold_legend <- cowplot::get_legend(relayer_logfold_change)    

((logfold_change_plot + ggtitle("T7 Only")) | relayer_logfold_legend) + plot_layout(widths = c(5, 1))






##### FROM MODEL ####

logfold_model_metrics <- significant_models %>% select(asv_id, (contains("estimate") | 
                                                                  contains("SE")) & (contains("DD_vs_DH_late") | contains("exposure") | 
                                                                                       contains("outcome"))) %>%
  pivot_longer(cols = contains("estimate"), names_prefix = "estimate_", names_to = "comparison", values_to = "estimate") %>%
  mutate(comparison = str_remove(comparison, '_[><=]')) %>%
  pivot_longer(cols = contains("SE"), names_prefix = "SE_", names_to = "comparison2", values_to = "t_stderr") %>%
  mutate(comparison2 = str_remove(comparison2, '_[><=]')) %>%
  filter(comparison == comparison2) %>%
  select(-comparison2)

model_for_plot <- significant_models %>%
  select(asv_id, starts_with('fdr')) %>% 
  select(-contains(c('treatment', 'tank', 'genotype'))) %>%
  mutate(across(starts_with('fdr'), ~ifelse(. < 0.05, "sig", "nonsig"))) %>%
  
  pivot_longer(cols = -asv_id,
               names_to = c('comparison'),
               values_to = 'significance') %>%
  mutate(comparison = str_remove(comparison, 'fdr_')) %>%
  mutate(comparison = str_remove(comparison, '_[><=]')) %>%
  filter(asv_id %in% asvs_for_logfold$asv_id) %>%
  filter(comparison %in% c("DD_vs_DH_late", "exposure", "outcome")) %>%
  left_join(logfold_model_metrics, by = join_by(asv_id, comparison)) %>%
  rbind(logfold_dose %>%
               mutate(fdr_pval = ifelse(fdr_pval < 0.05, "sig", "nonsig"),
                      comparison = "dose") %>%
          rename("significance" = "fdr_pval")) %>%
  mutate(graph_pval = str_c(comparison, significance, sep = "_")) %>%
  filter(comparison != "exposure") %>%
  left_join(asvs_for_logfold %>% select(asv_id, grouped_signatures), by = join_by(asv_id)) %>%
  left_join(asvs_for_logfold %>% ungroup() %>% select(asv_id, prep), by = join_by(asv_id)) %>%
  mutate(grouped_signatures = factor(grouped_signatures, 
                                     levels = c("Pathogens", "White Band Disease-Associated Opportunists", 
                                                "Healthy-Associated"))) %>%
  arrange(grouped_signatures) %>%
  mutate(asv_id = factor(asv_id), ordered = TRUE)

altered_data_make_legend <- for_plot %>%
  filter(asv_id %in% c("ASV_65", "ASV_134")) %>%
  mutate(graph_pval = ifelse((asv_id == "ASV_134" & comparison != "dose_D_minus_H"), str_replace(graph_pval, "sig", "nonsig"), graph_pval))

relayer_logfold_change_for_legend <- ggplot(altered_data_make_legend, aes(x = asv_id, y = estimate, ymin = estimate - t_stderr, ymax = estimate + t_stderr)) +
  geom_hline(yintercept = 0) +
  #outcome main effect
  (geom_point(data = (altered_data_make_legend %>% filter(comparison == "outcome_D_minus_H")), aes(colour1 = graph_pval, fill1 = graph_pval), pch = 23, size = 3, stroke = 1) %>%
     rename_geom_aes(new_aes = c("colour" = "colour1", "fill" = "fill1"))) + 
  scale_color_manual(aesthetics = "colour1", values = c("#A70000", "#A70000"), guide = "legend", 
                     name = "Outcome Main Effect", breaks = c("outcome_D_minus_H_sig", "outcome_D_minus_H_nonsig"), labels = c("Significant", "Non-Significant")) +
  scale_fill_manual(aesthetics = "fill1", values = c("#A70000", "transparent"), guide = "legend", 
                    name = "Outcome Main Effect", breaks = c("outcome_D_minus_H_sig", "outcome_D_minus_H_nonsig"), labels = c("Significant", "Non-Significant")) +
  #Dose main effect
  (geom_point(data = (altered_data_make_legend %>% filter(comparison == "dose_D_minus_H")), aes(colour3 = graph_pval, fill3 = graph_pval), pch = 22, size = 4, stroke = 1) %>%
     rename_geom_aes(new_aes = c("colour" = "colour3", "fill" = "fill3"))) + 
  scale_color_manual(aesthetics = "colour3", values = c("#732dcf", "#732dcf"), guide = "legend", 
                     name = "Homogenate Doses", breaks = c("dose_D_minus_H_sig", "dose_D_minus_H_nonsig"), labels = c("Significant", "Non-Significant")) +
  scale_fill_manual(aesthetics = "fill3", values = c("#732dcf", "transparent"), guide = "legend", 
                    name = "Homogenate Doses", breaks = c("dose_D_minus_H_sig", "dose_D_minus_H_nonsig"), labels = c("Significant", "Non-Significant")) +
  #exposure main effect
  (geom_point(data = (altered_data_make_legend %>% filter(comparison == "DD_minus_DH")), aes(colour4 = graph_pval, fill4 = graph_pval), pch = 24, size = 3, stroke = 1) %>%
     rename_geom_aes(new_aes = c("colour" = "colour4", "fill" = "fill4"))) + 
  scale_color_manual(aesthetics = "colour4", values = c("#EC8E00", "#EC8E00"), guide = "legend", 
                     name = str_wrap("Final Outcome of Disease-Exposed", 20), breaks = c("DD_minus_DH_sig", "DD_minus_DH_nonsig"), labels = c("Significant", "Non-Significant")) +
  scale_fill_manual(aesthetics = "fill4", values = c("#EC8E00", "transparent"), guide = "legend", 
                    name = str_wrap("Final Outcome of Disease-Exposed", 20), breaks = c("DD_minus_DH_sig", "DD_minus_DH_nonsig"), labels = c("Significant", "Non-Significant")) +
  coord_flip() +
  theme_bw() +
  guides(#outcome
    color1 = guide_legend(order = 2),
    fill1 = guide_legend(order = 2),
    #dose
    color3 = guide_legend(order = 4),
    fill3 = guide_legend(order = 4),
    #DD_DH
    color4 = guide_legend(order = 3),
    fill4 = guide_legend(order = 3))


relayer_logfold_legend <- cowplot::get_legend(relayer_logfold_change_for_legend) 

logfold_change_plot <- ggplot(model_for_plot) +
  geom_hline(yintercept = 0) +
  geom_pointrange(aes(x = asv_id, y = estimate, fill = graph_pval, col = comparison, pch = comparison, ymin = estimate - t_stderr, ymax = estimate + t_stderr), 
                  size = 1, stroke = 1, position = position_dodge(width = 0.65), linewidth = 0.8) + #position = position_dodge(width = 1)
  scale_shape_manual(values = c(24, 22, 23)) +
  scale_color_manual(values = c("outcome" = "#A70000",
                                "dose" = "#732dcf", 
                                "DD_vs_DH_late" = "#EC8E00")) +
  scale_fill_manual(values = c("outcome_nonsig" = "white", 
                               "outcome_sig" = "#A70000",
                               "exposure_nonsig" = "white", 
                               "dose_nonsig" = "white", 
                               "dose_sig" = "#732dcf", #DF369B
                               "DD_vs_DH_late_nonsig" = "white", 
                               "DD_vs_DH_late_sig" = "#EC8E00"), guide = "none") +
  coord_flip() +
  theme_bw() +
  theme(#strip.background = element_blank(),
        #strip.text.y = element_blank(),
        legend.position = "none") +
  scale_x_discrete(breaks = for_plot$asv_id, labels = for_plot$prep) +
  facet_grid(rows = vars(grouped_signatures), space = "free", scales = "free") +
  xlab(NULL) +
  ylab("Logfold Change (Diseased vs. Healthy)")
  
(logfold_change_plot | relayer_logfold_legend) + plot_layout(widths = c(5, 1))
