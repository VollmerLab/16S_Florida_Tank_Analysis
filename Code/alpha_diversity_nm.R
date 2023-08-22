
#use alpha_table for normal alpha diversity

#use alpha_table_wo_aquaricketts for alpha diversity without ASV_1
edited_microbiome <- microbiome_data

otu_data <- otu_table(edited_microbiome)

otu_data[,"ASV_1"] <- 0

otu_table(edited_microbiome) <- otu_data

alpha_table_wo_aquaricketts <- microbiome::alpha(edited_microbiome, index = "all") %>%
  as_tibble(rownames = 'sample_id') %>%
  inner_join(metadata, by = 'sample_id') %>%
  mutate(fragment_id = str_c(str_replace_na(exposure, 'NA'), tank, genotype, sep = '_'))
#


t3t7_alpha_models <- alpha_table_wo_aquaricketts %>%
  filter(!tank %in% c("HOMO", "homogenate_fragment")) %>%
  mutate(treatment = str_c(time, exposure, susceptability, sep = '_')) %>%
  pivot_longer(cols = !c(colnames(metadata), "fragment_id", "treatment"), #retain_sample
               names_to = 'metric',
               values_to = 'alpha_div_value') %>%
  select(-c(final_disease_state, clone_group)) %>%
  filter(time %in% c("T3", "T7")) %>%
  nest_by(metric) %>%
  summarise(t3t7_model = list(lmer(alpha_div_value ~ time*exposure*susceptability + 
                                     (1 | genotype) + (1 | tank), data = data)))

t0_alpha_models <- alpha_table_wo_aquaricketts %>%
  filter(!tank %in% c("HOMO", "homogenate_fragment")) %>%
  mutate(treatment = str_c(time, exposure, susceptability, sep = '_')) %>%
  pivot_longer(cols = !c(colnames(metadata), "fragment_id", "treatment"), #retain_sample
               names_to = 'metric',
               values_to = 'alpha_div_value') %>%
  select(-c(final_disease_state, clone_group)) %>%
  filter(time %in% c("T0")) %>%
  nest_by(metric) %>%
  summarise(t0_model = list(lm(alpha_div_value ~ susceptability, data = data))) 

t0_plot_prep <- t0_alpha_models %>%
  rowwise() %>%
  mutate(t0_plot_info = list(emmeans(t0_model, ~susceptability) %>%
                               cld(Letters = letters) %>% # , adjust = 'fdr' ??
                               broom::tidy(conf.int = TRUE) %>%
                               mutate(.group = str_trim(.group)) %>%
                               mutate(c_time = 0) %>%
                               mutate(time = NA, exposure = NA, facet_lab = "Field") %>%
                               mutate(facet_lab = factor(facet_lab, levels = c("Field", "Experimental")))
  )) 

t3t7_sig_terms <- t3t7_alpha_models %>%
  rowwise() %>% 
  mutate(sig_terms = list(anova(t3t7_model) %>% 
                            rownames_to_column(var = "sig_term") %>% 
                            as_tibble() %>% 
                            rename("p_val" = `Pr(>F)`) %>%
                            mutate(fdr_p_val = p.adjust(p_val, method = 'fdr')) %>%
                            filter(fdr_p_val < 0.05) %>%
                            filter(sig_term != "time") %>%
                            pull(sig_term))) %>%
  #filter(length(sig_terms) > 0) %>%
  summarise(sig_terms = str_c(sig_terms, collapse = ', '))


t3t7_sig_termswp <- t3t7_alpha_models %>%
  rowwise() %>% 
  mutate(sig_terms = list(anova(t3t7_model) %>% 
                            rownames_to_column(var = "sig_term") %>% 
                            as_tibble() %>% 
                            rename("p_val" = `Pr(>F)`) %>%
                            mutate(fdr_p_val = p.adjust(p_val, method = 'fdr')) %>%
                            filter(fdr_p_val < 0.05) %>%
                            filter(sig_term != "time") %>%
                            select(sig_term, fdr_p_val))) %>%
  filter(length(sig_terms) > 0)

t3t7_sig_termswp %>% unnest(sig_terms)


t3t7_plot_prep <- t3t7_alpha_models %>%
  left_join(t0_plot_prep, by = join_by("metric")) %>%
  left_join(t3t7_sig_terms, by = join_by("metric")) %>%
  filter(!is.na(sig_terms)) %>%
  rowwise() %>%
  mutate(t3t7_plot_info = list(emmeans(t3t7_model, ~time*exposure*susceptability) %>%
                                 cld(Letters = LETTERS) %>% # , adjust = 'fdr' ??
                                 broom::tidy(conf.int = TRUE) %>%
                                 mutate(.group = str_trim(.group)) %>%
                                 mutate(c_time = parse_number(time)) %>%
                                 mutate(facet_lab = "Experimental") %>%
                                 mutate(facet_lab = factor(facet_lab, levels = c("Field", "Experimental")))
  )) %>%
  mutate(both_plot_infos = list(rbind(t0_plot_info, t3t7_plot_info)))



alpha_div_plots_only_four 
alpha_div_plots_no_asv1 <- t3t7_plot_prep %>%
  #filter(metric %in% c("dominance_dbp", "evenness_camargo", "rarity_low_abundance", "diversity_shannon")) %>%
  #mutate(for_manuscript = "yes") %>%
  rowwise() %>%
  mutate(plot_info = list(both_plot_infos %>%
                            mutate(graph_cat = ifelse(c_time == "0", NA, 
                                                      paste(exposure, susceptability, sep = "_"))) %>%
                            {. ->> intermed } %>%
                            mutate(graph_cat = ifelse(c_time == "0", paste("H", susceptability, sep = "_"), 
                                                      graph_cat)) %>%
                            dplyr::slice(rep(1:2, 1)) %>%
                            rbind(intermed) %>%
                            mutate(graph_cat = ifelse(is.na(graph_cat), paste("D", susceptability, sep = "_"), 
                                                      graph_cat)) %>%
                            mutate(c_time = ifelse(c_time == "0", ifelse(susceptability == "S", c_time - 0.1, c_time + 0.1),
                                                   case_when(graph_cat == "D_S" ~ c_time - 0.35,
                                                             graph_cat == "H_S" ~ c_time - 0.15,
                                                             graph_cat == "D_R" ~ c_time + 0.15,
                                                             graph_cat == "H_R" ~ c_time + 0.35))) %>%
                            mutate(graph_cat = factor(graph_cat, levels = c("D_S", "D_R", "H_S", "H_R"), labels = c("D_S", "D_R", "H_S", "H_R"))))) %>%
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
      
      
      (geom_text(data = (plot_info %>% filter(graph_cat %in% c("D_S", "D_R"))), size = 4, vjust = -1, show.legend = FALSE, aes(colour1 = graph_cat, y = conf.high, label = .group)) %>%
         rename_geom_aes(new_aes = c("colour" = "colour1"))) + 
      (geom_text(data = (plot_info %>% filter(graph_cat %in% c("H_S", "H_R"))), size = 4, vjust = -1, show.legend = FALSE, aes(colour2 = graph_cat, y = conf.high, label = .group)) %>%
         rename_geom_aes(new_aes = c("colour" = "colour2"))) +
      
      (geom_point(data = (plot_info %>% filter(graph_cat == "dose" & susceptability == "na")), size = 3.7, aes(colour3 = exposure), shape = "diamond") %>%
         rename_geom_aes(new_aes = c("colour" = "colour3"))) +
      (geom_errorbar(data = (plot_info %>% filter(graph_cat == "dose" & susceptability == "na")), width = 0, aes(colour3 = exposure)) %>%
         rename_geom_aes(new_aes = c("colour" = "colour3"))) + 
      (geom_point(data = (plot_info %>% filter(graph_cat == "dose" & is.na(exposure))), size = 3, shape = 1, col = "black", alpha = 0)) +
      
      scale_color_manual(aesthetics = "colour1", values = c("#F75D5D", "#A70000"), guide = "legend", 
                         name = "Disease Exposed", labels = c("Susceptible", "Resistant")) +
      scale_shape_manual(values = c(17, 16, 17, 16), guide = "none") +
      scale_color_manual(aesthetics = "colour2", values = c("#3DD8EA", "#048291"), guide = "legend", 
                         name = "Healthy Exposed", labels = c("Susceptible", "Resistant")) +
      guides(colour1 = guide_legend(
        override.aes=list(linetype = c(6, 1), shape = c(16, 17))),
        colour2 = guide_legend(
          override.aes=list(linetype = c(6, 1), shape = c(16, 17)))) +
      scale_x_continuous(breaks=c(0, 3, 7)) +
      scale_linetype_manual(values = c(1, 6, 1, 6), guide = "none") +
      theme_bw() +
      xlab("Time") +
      ylab(expression("Normalized log"[2]*" (cpm)")) +
      labs(title = str_c(metric, " - ", sig_terms, sep = ""))
      #facet_grid(cols = vars(facet_lab), space = "free", scales = "free_x")
  )) %>%
  mutate(alpha_type = ifelse(metric %in% c("chao1", "observed"), "richness", str_extract(metric, "[^_]+"))) %>%
  group_by(alpha_type) %>%
  #group_by(for_manuscript) %>%
  summarise(combo_plots = list(wrap_plots(plot) + plot_layout(guides = 'collect') ))#& plot_annotation(title = alpha_type)))

alpha_div_plots_only_four$combo_plots[[1]]

alpha_div_plots$combo_plots[[4]]

alpha_div_plots_no_asv1$combo_plots[[1]]
