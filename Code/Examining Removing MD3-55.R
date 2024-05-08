# alpha div after removing asv 1

everything_except_asv1 <- subset(x = otu_table(altered_microbiome_data), select = colnames(otu_table(altered_microbiome_data))[2:length(colnames(otu_table(altered_microbiome_data)))])
altered_microbiome_data_no_asv1 <- merge_phyloseq(everything_except_asv1, tax_table(altered_microbiome_data), sample_data(altered_microbiome_data), phy_tree(altered_microbiome_data), refseq(altered_microbiome_data))

md355 <- c("ASV_1", "ASV_13", "ASV_122", "ASV_168", "ASV_345", "ASV_1063")


#looks exact same as only removing asv 1
non_md355_asvs <- colnames(otu_table(altered_microbiome_data)) %>%
  as_tibble() %>%
  filter(!value %in% md355) %>%
  pull(value)

everything_except_md355 <- subset(x = otu_table(altered_microbiome_data), select = non_md355_asvs)
altered_microbiome_data_except_md355 <- merge_phyloseq(everything_except_md355, tax_table(altered_microbiome_data), sample_data(altered_microbiome_data), phy_tree(altered_microbiome_data), refseq(altered_microbiome_data))


alpha_table_except_md355 <- microbiome::alpha(altered_microbiome_data_except_md355, index = "all") %>%
  as_tibble(rownames = 'sample_id') %>%
  inner_join(metadata, by = 'sample_id') %>%
  mutate(fragment_id = str_c(str_replace_na(exposure, 'NA'), tank, genotype, sep = '_'))

mod_alpha_tab_except_md355 <- alpha_table_except_md355 %>%
  filter(!tank %in% c("HOMO", "homogenate_fragment")) %>%
  mutate(final_disease_state = ifelse(exposure == "F", "F", final_disease_state)) %>%
  mutate(treatment = str_c(time, exposure, final_disease_state, sep = '_')) %>%
  pivot_longer(cols = !c(colnames(metadata), "fragment_id", "treatment"),
               names_to = 'metric',
               values_to = 'alpha_div_value') %>%
  select(-c(susceptability, resistance, clone_group)) %>%
  mutate(tank_field = if_else(str_detect(treatment, 'F'), 'field', 'tank'), .after = final_disease_state) %>%
  nest_by(metric) %>%
  summarise(alpha_model_except_md355 = list(lmer(alpha_div_value ~ treatment + 
                                      (1 | genotype) + (0 + dummy(tank_field, c("tank")) | tank),
                                    data = data))) %>% 
  rowwise() %>% 
  mutate(p_value = anova(alpha_model_except_md355) %>% 
           rownames_to_column(var = "sig_term") %>% 
           as_tibble() %>% 
           dplyr::rename("p_val" = `Pr(>F)`) %>%
           pull(p_val)) %>%
  ungroup() %>%
  mutate(fdr_p_val = p.adjust(p_value, method = 'fdr')) %>%
  filter(fdr_p_val < 0.05) %>%
  mutate(alpha_type = ifelse(metric %in% c("chao1", "observed"), "richness", str_extract(metric, "[^_]+")))

## figuring out the fig

alpha_graphs_manuscript_except_md355 <- mod_alpha_tab_except_md355 %>%
  mutate(metric = ifelse(str_detect(metric, "chao1"), "richness_chao1", metric)) %>%
  mutate(for_manuscript = ifelse(metric %in% c("diversity_shannon", "dominance_core_abundance", "evenness_camargo", "richness_chao1"), TRUE, FALSE)) %>%
  rowwise() %>%
  mutate(plot_info = list(emmeans(alpha_model_except_md355, ~treatment) %>%
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
  rowwise() %>%
  mutate(plot = list(
    ggplot(data = plot_info, aes(x = c_time, y = estimate, ymin = conf.low, ymax = conf.high)) +
      (geom_line(data = (plot_info %>% filter(graph_cat %in% c("D_H", "D_D"))), aes(colour1 = graph_cat, linetype = graph_cat)) %>%
         rename_geom_aes(new_aes = c("colour" = "colour1"))) + 
      (geom_line(data = (plot_info %>% filter(graph_cat %in% c("H_H"))), aes(colour2 = graph_cat, linetype = graph_cat)) %>%
         rename_geom_aes(new_aes = c("colour" = "colour2"))) +
      (geom_errorbar(data = (plot_info %>% filter(graph_cat %in% c("D_H", "D_D"))), width = 0, aes(colour1 = graph_cat)) %>%
         rename_geom_aes(new_aes = c("colour" = "colour1"))) + 
      (geom_errorbar(data = (plot_info %>% filter(graph_cat %in% c("H_H"))), width = 0, aes(colour2 = graph_cat)) %>%
         rename_geom_aes(new_aes = c("colour" = "colour2"))) +
      (geom_point(data = (plot_info %>% filter(graph_cat %in% c("D_H", "D_D"))), size = 3, aes(colour1 = graph_cat, pch = graph_cat)) %>%
         rename_geom_aes(new_aes = c("colour" = "colour1"))) + 
      (geom_point(data = (plot_info %>% filter(graph_cat %in% c("H_H"))), size = 3, aes(colour2 = graph_cat, pch = graph_cat)) %>%
         rename_geom_aes(new_aes = c("colour" = "colour2"))) +
      
      
      (geom_point(data = (plot_info %>% filter(graph_cat == "dose" & final_disease_state == "na")), size = 3.7, aes(colour3 = exposure), shape = "diamond") %>%
         rename_geom_aes(new_aes = c("colour" = "colour3"))) +
      (geom_errorbar(data = (plot_info %>% filter(graph_cat == "dose" & final_disease_state == "na")), width = 0, aes(colour3 = exposure)) %>%
         rename_geom_aes(new_aes = c("colour" = "colour3"))) + 
      (geom_point(data = (plot_info %>% filter(graph_cat == "dose" & is.na(exposure))), size = 3, shape = 1, col = "black", alpha = 0)) +
      
      scale_color_manual(aesthetics = "colour1", values = c("D_H" = "#22A7B6", "D_D" = "#A70000"), guide = "legend", 
                         name = "Disease Exposed", breaks = c("D_H", "D_D"), labels = c("Healthy", "Diseased")) +
      scale_color_manual(aesthetics = "colour3", breaks = c("H", "D"), values = c("H" = "#17D0B4", "D" = "#e30e0e"), guide = "legend", 
                         name = "Doses", labels = c("H" = "Healthy", "D" = "Diseased")) +
      scale_shape_manual(values = c("D_H" = 17, "D_D" = 17, "H_H" = 16), guide = "none") +
      scale_color_manual(aesthetics = "colour2", values = c("H_H" = "#406F23"), guide = "legend", 
                         name = "Healthy Exposed", labels = c("Healthy")) +
      guides(colour1 = guide_legend(
        override.aes=list(linetype = c(1, 1), shape = c(17, 17))),
        colour2 = guide_legend(
          override.aes=list(linetype = c(6), shape = c(16))),
        colour3 = guide_legend(
          override.aes=list(linetype = c(0, 0)))) +
      scale_x_continuous(breaks=c(0, 3, 7)) +
      scale_linetype_manual(values = c("D_H" = 1, "D_D" = 1, "H_H" = 6), guide = "none") +
      theme_bw() +
      xlab("Time") +
      ylab(metric) +
      labs(title = metric)
  )) %>%
  group_by(for_manuscript) %>%
  summarise(combo_plots = list(wrap_plots(plot) + plot_layout(guides = 'collect')))

alpha_graphs_manuscript_except_md355$combo_plots[2]


normalized_asv_counts %>% select(asv_id, Family, Genus, Species) %>%
  distinct() %>%
  filter(Genus == "MD3-55")




alpha_table_except_md355 %>%
  select(sample_id, diversity_gini_simpson, time, final_disease_state, susceptability) %>%
  rename("ginisimpson_no_md355" = "diversity_gini_simpson") %>%
  left_join(alpha_table %>%
              select(sample_id, diversity_gini_simpson, time, final_disease_state, susceptability) %>%
              rename("ginisimpson_full" = "diversity_gini_simpson")) %>%
  filter(!is.na(susceptability)) %>%
  pivot_longer(cols = c(ginisimpson_no_md355, ginisimpson_full), names_to = "filter_type", values_to = "ginisimpson_value") %>%
  ggplot() +
  geom_boxplot(aes(x = filter_type, y = ginisimpson_value, fill = susceptability)) +
  facet_grid(vars(rows = time), vars(cols = final_disease_state))
  


mod_alpha_tab_except_md355_resist <- alpha_table_except_md355 %>%
  mutate(filter_type = "except_md355") %>%
  rbind(alpha_table %>%
          mutate(filter_type = "full")) %>%
  filter(!tank %in% c("HOMO", "homogenate_fragment")) %>%
  mutate(final_disease_state = ifelse(exposure == "F", "F", final_disease_state)) %>%
  filter(time != "T0") %>%
  mutate(treatment = str_c(time, susceptability, final_disease_state, sep = '_')) %>%
  pivot_longer(cols = !c(colnames(metadata), "fragment_id", "treatment", "filter_type"),
               names_to = 'metric',
               values_to = 'alpha_div_value') %>%
  select(-c(exposure, resistance, clone_group)) %>%
  mutate(tank_field = if_else(str_detect(treatment, 'F'), 'field', 'tank'), .after = final_disease_state) %>%
  nest_by(metric) %>%
  summarise(alpha_model = list(lmer(alpha_div_value ~ time*susceptability*final_disease_state*filter_type +
                                                   (1 | genotype) + (1 | tank),
                                                 data = data))) # %>% 
  rowwise() %>% 
  mutate(p_value = anova(alpha_model) %>% 
           rownames_to_column(var = "sig_term") %>% 
           as_tibble() %>% 
           dplyr::rename("p_val" = `Pr(>F)`) %>%
           filter(sig_term == "treatment") %>% #p value is just for treatment not for the interaction
           pull(p_val)) %>%
  ungroup() %>%
  mutate(fdr_p_val = p.adjust(p_value, method = 'fdr')) %>%
  filter(fdr_p_val < 0.05) %>%
  mutate(alpha_type = ifelse(metric %in% c("chao1", "observed"), "richness", str_extract(metric, "[^_]+")))

ginisimpson_model <- mod_alpha_tab_except_md355_resist%>% filter(metric == "diversity_gini_simpson")

emmeans(ginisimpson_model$alpha_model[[1]], ~time*susceptability*final_disease_state*filter_type, type = 'response') %>%
  multcomp::cld(Letters = LETTERS, by = c("time", "final_disease_state")) %>% 
  as_tibble() %>%
  mutate(.group = str_trim(.group)) %>%
  rename(response = emmean,
         asymp.UCL = upper.CL,
         asymp.LCL = lower.CL) %>%
  ggplot(aes(x = filter_type, y = response, ymin = asymp.LCL, ymax = asymp.UCL, colour = susceptability)) +
  geom_pointrange(position = position_dodge(0.5)) +
  geom_text(aes(y = asymp.UCL, label = .group), position = position_dodge(0.5),
            vjust = -1) +
  facet_grid(vars(rows = time), vars(cols = final_disease_state)) +
  theme_bw() +
  scale_color_manual(values = c("R" = "forestgreen", "S" = "red3"))


#T0
mod_alpha_tab_except_md355_resist_t0 <- alpha_table_except_md355 %>%
  mutate(filter_type = "except_md355") %>%
  rbind(alpha_table %>%
          mutate(filter_type = "full")) %>%
  filter(!tank %in% c("HOMO", "homogenate_fragment")) %>%
  mutate(final_disease_state = ifelse(exposure == "F", "F", final_disease_state)) %>%
  filter(time == "T0") %>%
  mutate(treatment = str_c(time, susceptability, final_disease_state, sep = '_')) %>%
  pivot_longer(cols = !c(colnames(metadata), "fragment_id", "treatment", "filter_type"),
               names_to = 'metric',
               values_to = 'alpha_div_value') %>%
  select(-c(exposure, resistance, clone_group)) %>%
  #mutate(tank_field = if_else(str_detect(treatment, 'F'), 'field', 'tank'), .after = final_disease_state) %>%
  nest_by(metric) %>%
  summarise(alpha_model = list(lmer(alpha_div_value ~ susceptability*filter_type +
                                      (1 | genotype),
                                    data = data)))

t0_ginisimpson_model <- mod_alpha_tab_except_md355_resist%>% filter(metric == "diversity_gini_simpson")

emmeans(t0_ginisimpson_model$alpha_model[[1]], ~susceptability*filter_type, type = 'response') %>%
  multcomp::cld(Letters = LETTERS) %>% 
  as_tibble() %>%
  mutate(.group = str_trim(.group)) %>%
  rename(response = emmean,
         asymp.UCL = upper.CL,
         asymp.LCL = lower.CL) %>%
  ggplot(aes(x = susceptability, y = response, ymin = asymp.LCL, ymax = asymp.UCL, colour = susceptability)) +
  geom_pointrange(position = position_dodge(0.5)) +
  geom_text(aes(y = asymp.UCL, label = .group), position = position_dodge(0.5),
            vjust = -1) +
  facet_wrap(~filter_type) +
  theme_bw() +
  scale_color_manual(values = c("R" = "forestgreen", "S" = "red3"))
