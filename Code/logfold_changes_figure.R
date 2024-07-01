asvs_for_logfold <- bacterial_signature_asv %>%
  #filter(asv_id == "ASV_65") %>%
  group_by(asv_id) %>%
  reframe(signatures = str_c(signatures, collapse = ', ')) %>%
  mutate(grouped_signatures = case_when(str_detect(signatures, "Aquaria") ~ "Tank Effect",
                                        str_detect(signatures, "Field") ~ "Tank Effect",
                                        str_detect(signatures, "DiseaseOutcome") & str_detect(signatures, "DD_early") & str_detect(signatures, "DD_vs_DH_early") & asv_id %in% sig_homog_asv_list ~ "Putative Early Pathogens",
                                        str_detect(signatures, "DiseaseOutcome") & str_detect(signatures, "DD_late") & str_detect(signatures, "DD_vs_DH_late") & asv_id %in% sig_homog_asv_list ~ "Putative Late Pathogens",
                                        str_detect(signatures, "DiseaseOutcome") & str_detect(signatures, "DD_late") & str_detect(signatures, "DD_vs_DH_late") & !asv_id %in% sig_homog_asv_list~ "Specialized Opportunists",
                                        str_detect(signatures, "DiseaseOutcome") ~ "Unlikely Pathogens",
                                        str_detect(signatures, "HealthyOutcome") ~ "Healthy Associated",
                                        TRUE ~ signatures)) %>%
  #don't bother plotting tank effect, etc...
  filter(grouped_signatures %in% c("Putative Early Pathogens", "Putative Late Pathogens", "Specialized Opportunists", "Healthy Associated")) 


#%>%
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
                       
  )
  ) %>%
  rowwise %>%
  mutate(plot_info = list(emmeans(model, ~treatment) %>%
                            broom::tidy(conf.int = TRUE) %>%
                            separate(treatment, into = c('time', 'exposure', 'final_disease_state')) %>%
                            mutate(graph_cat = ifelse(time == "T0", NA, 
                                                      paste(exposure, final_disease_state, sep = "_"))) %>%
                            {. ->> intermed } %>%
                            mutate(graph_cat = ifelse(time == "T0", "D_D", 
                                                      graph_cat)) %>%
                            dplyr::slice(1) %>%
                            rbind(intermed) %>%
                            mutate(graph_cat = ifelse(is.na(graph_cat), "D_H", 
                                                      graph_cat)) %>%
                            dplyr::slice(rep(1:2, 1)) %>%
                            rbind(intermed) %>%
                            mutate(graph_cat = ifelse(is.na(graph_cat), "H_H", 
                                                      graph_cat)) %>%
                            mutate(c_time = parse_number(time)) %>%
                            mutate(facet_lab = "Experimental") %>%
                            mutate(c_time = ifelse(time == "T0", c_time,
                                                   case_when(graph_cat == "D_D" ~ c_time + 0.30,
                                                             graph_cat == "D_H" ~ c_time,
                                                             graph_cat == "H_H" ~ c_time - 0.30))) %>%
                            mutate(graph_cat = factor(graph_cat, levels = c("D_D", "D_H", "H_H"), labels = c("D_D", "D_H", "H_H"))))) %>%
  left_join(homogenate_models, by = join_by(asv_id))

hmm



#compare
  #disease outcome
  #disease exposure
  #dose


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
  mutate(across(ends_with('t_test'), map, broom::tidy)) %>%
  unnest(dose_t_test) %>%
  rename("D" = "estimate1", "H" = "estimate2", "dose_D_minus_H" = "estimate") %>%
  select(-c(parameter, conf.low, conf.high, method, alternative)) %>%
  mutate(fdr_dose = p.adjust(p.value)) %>%
  select(asv_id, dose_D_minus_H, fdr_dose) %>%
  pivot_longer(dose_D_minus_H, names_to = "comparison", values_to = "estimate") %>%
  rename("fdr_pval" = "fdr_dose")

logfold_dose <- logfold_dose_1 %>% ungroup() %>% 
  dplyr::slice(1:3) %>% mutate(across(everything(), ~NA)) %>% 
  mutate(asv_id = c("ASV_301", "ASV_372", "ASV_458"), fdr_pval = c(1, 1, 1), 
         comparison = c("dose_D_minus_H", "dose_D_minus_H", "dose_D_minus_H"), estimate = c(0, 0, 0)) %>%
  rbind(logfold_dose)
  
#outcome

logfold_outcome <- significant_models %>%
  filter(asv_id %in% asvs_for_logfold$asv_id) %>%
  unnest(data) %>%
  filter(time == "T7") %>%
  group_by(asv_id) %>%
  summarise_at(vars(log2_cpm), list(outcome_t_test = ~ list(t.test(. ~ final_disease_state)))) %>%
  mutate(across(ends_with('t_test'), map, broom::tidy)) %>%
  unnest(outcome_t_test) %>%
  rename("D" = "estimate1", "H" = "estimate2", "outcome_D_minus_H" = "estimate") %>%
  select(-c(parameter, conf.low, conf.high, method, alternative)) %>%
  mutate(fdr_outcome = p.adjust(p.value)) %>%
  select(asv_id, outcome_D_minus_H, fdr_outcome) %>%
  pivot_longer(outcome_D_minus_H, names_to = "comparison", values_to = "estimate") %>%
  rename("fdr_pval" = "fdr_outcome")
  
#exposure

logfold_exposure <- significant_models %>%
  filter(asv_id %in% asvs_for_logfold$asv_id) %>%
  unnest(data) %>%
  filter(time == "T7") %>%
  group_by(asv_id) %>%
  summarise_at(vars(log2_cpm), list(exposure_t_test = ~ list(t.test(. ~ exposure)))) %>%
  mutate(across(ends_with('t_test'), map, broom::tidy)) %>%
  unnest(exposure_t_test) %>%
  rename("D" = "estimate1", "H" = "estimate2", "exposure_D_minus_H" = "estimate") %>%
  select(-c(parameter, conf.low, conf.high, method, alternative)) %>%
  mutate(fdr_exposure = p.adjust(p.value)) %>%
  select(asv_id, exposure_D_minus_H, fdr_exposure) %>%
  pivot_longer(exposure_D_minus_H, names_to = "comparison", values_to = "estimate") %>%
  rename("fdr_pval" = "fdr_exposure")

#exposure on outcome

logfold_DD_DH <- significant_models %>%
  filter(asv_id %in% asvs_for_logfold$asv_id) %>%
  unnest(data) %>%
  filter(treatment %in% c("T7_D_D", "T7_D_H")) %>%
  group_by(asv_id) %>%
  summarise_at(vars(log2_cpm), list(DD_DH_t_test = ~ list(t.test(. ~ treatment)))) %>%
  mutate(across(ends_with('t_test'), map, broom::tidy)) %>%
  unnest(DD_DH_t_test) %>%
  rename("D_D" = "estimate1", "D_H" = "estimate2", "DD_minus_DH" = "estimate") %>%
  select(-c(parameter, conf.low, conf.high, method, alternative)) %>%
  mutate(fdr_DD_DH = p.adjust(p.value)) %>%
  select(asv_id, DD_minus_DH, fdr_DD_DH) %>%
  pivot_longer(DD_minus_DH, names_to = "comparison", values_to = "estimate") %>%
  rename("fdr_pval" = "fdr_DD_DH")
   
#combine
for_plot <- logfold_outcome %>%
  rbind(logfold_exposure, logfold_DD_DH, logfold_dose) %>%
  mutate(graph_pval = str_c(comparison, ifelse(fdr_pval < 0.05, "sig", "nonsig"), sep = "_")) %>%
  left_join(asvs_for_logfold %>% select(asv_id, grouped_signatures), by = join_by(asv_id)) %>%
  left_join(the_plots1 %>% ungroup() %>% select(asv_id, prep), by = join_by(asv_id)) %>%
  mutate(grouped_signatures = factor(grouped_signatures, levels = c("Putative Late Pathogens",
                                            "Specialized Opportunists", "Healthy Associated"))) %>%
  arrange(grouped_signatures) %>%
  mutate(asv_id = factor(asv_id), ordered = TRUE)
  #left_join(significant_models %>% select(asv_id, Domain:Species), by = join_by(asv_id)) %>%


  ggplot(for_plot) +
  geom_hline(yintercept = 0) +
  geom_point(aes(x = asv_id, y = estimate, fill = graph_pval, col = comparison, pch = comparison), 
             size = 3, stroke = 1, position = position_dodge(width = 1)) + #position = position_dodge(width = 1)
  scale_shape_manual(values = c(21, 22, 23, 24)) +
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
          strip.text.y = element_blank()) +
  scale_x_discrete(breaks = for_plot$asv_id, labels = for_plot$prep) +
  facet_grid(rows = vars(grouped_signatures), space = "free", scales = "free") +
  xlab(NULL) +
  ylab("Logfold Difference (D-H)")


  








