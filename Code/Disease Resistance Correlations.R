
#### FUNCTIONS ####

color_facets <- function(plot, fills) {
  # https://github.com/tidyverse/ggplot2/issues/2096
  g <- ggplot_gtable(ggplot_build(plot))
  strip_t <- which(grepl('strip-t', g$layout$name))
  k <- 1
  for (i in strip_t) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
  }
  grid.draw(g)
}

#### ####


dr_corr_all <- full_data %>%
  group_by(asv_id) %>%
  summarize(corr_val = broom::tidy(cor.test(log2_cpm, resistance))) %>%
  left_join(taxonomy_tibble %>% select(asv_names, Family, Genus), by = c("asv_id" = "asv_names")) %>%
  unnest(corr_val) %>%
  select(-c(method, alternative)) %>%
  filter(p.value < 0.05) %>%
  mutate(pos_neg = ifelse(estimate < 0, "Negative Correlation", "Positive Correlation"))

dr_corr_all_plot <- ggplot(dr_corr_all) +
    geom_hline(yintercept = 0) +
    geom_point(aes(x = fct_reorder(paste(Family, asv_id, sep = " "), estimate), y = estimate, col = pos_neg)) +
    coord_flip() +
    facet_wrap(~pos_neg, scales = "free") +
    theme_bw() +
    xlab("ASV") +
    theme(strip.background = element_rect(color="black", linewidth = 1, linetype="solid"),
            legend.position = "none") +
    scale_color_manual(values = c("#EE3434","#34DFAE")) +
    ggtitle("Significant Correlations - All Time Points")
  
color_facets(dr_corr_all_plot, c("#EE3434","#34DFAE"))  
  
#T0 only
  
dr_corr_t0 <- full_data %>%
    filter(time == "T0") %>%
    group_by(asv_id) %>%
    filter(length(unique(zapsmall(log2_cpm))) > 1) %>%
    summarize(corr_val = broom::tidy(cor.test(log2_cpm, resistance))) %>%
    left_join(taxonomy_tibble %>% select(asv_names, Family, Genus), by = c("asv_id" = "asv_names")) %>%
    unnest(corr_val) %>%
    select(-c(method, alternative)) %>%
    filter(p.value < 0.05) %>%
    mutate(pos_neg = ifelse(estimate < 0, "Negative Correlation", "Positive Correlation"))
  
dr_corr_t0_plot <- ggplot(dr_corr_t0) +
    geom_hline(yintercept = 0) +
    geom_point(aes(x = fct_reorder(paste(Family, asv_id, sep = " "), estimate), y = estimate, col = pos_neg)) +
    coord_flip() +
    facet_wrap(~pos_neg, scales = "free") +
    theme_bw() +
    xlab("ASV") +
    theme(strip.background = element_rect(color="black", linewidth = 1, linetype="solid"),
          legend.position = "none") +
    scale_color_manual(values = c("#EE3434","#34DFAE")) +
    ggtitle("Significant Correlations - T7")

color_facets(dr_corr_t0_plot, c("#EE3434","#34DFAE"))  



dr_corr_t0_plot <- ggplot(dr_corr_t0 %>% mutate(fam_col = case_when(Family == "Colwelliaceae" ~ "colwell",
                                                                    Family == "Rhodobacteraceae" ~ "rhodo",
                                                                    Family == "Flavobacteriaceae" ~ "flavo",
                                                                    Family == "Vibrionaceae" ~ "vibrio")) %>%
                            mutate(fam_col = ifelse(is.na(fam_col), "other", fam_col))) +
                            #mutate(colwell_col = ifelse(Family == "Colwelliaceae", "colwell", pos_neg))) +
  geom_hline(yintercept = 0) +
  geom_point(aes(x = fct_reorder(paste(Family, asv_id, sep = " "), estimate), y = estimate, col = fam_col), size = 3) +
  coord_flip() +
  facet_wrap(~pos_neg, scales = "free") +
  theme_bw() +
  xlab("ASV") +
  theme(strip.background = element_rect(color="black", linewidth = 1, linetype="solid")) +
  scale_color_manual(values = c("firebrick1", "orange", "gray", "deepskyblue", "hotpink")) +
  #scale_color_manual(values = c("darkred", "#EE3434","#34DFAE")) +
  ggtitle("Significant Correlations - T0")

color_facets(dr_corr_t0_plot, c("#EE3434","#34DFAE"))  
  
all_corr <- dr_corr_all %>% select(asv_id, pos_neg)  
t0_corr <- dr_corr_t0 %>% select(asv_id, pos_neg)  
  
  ### making plots

disease_corr_plots_all <- asv_models %>% 
  filter(asv_id %in% t0_corr$asv_id) %>% #choose t0 or all
  inner_join(t0_corr, by = 'asv_id') %>% #choose t0 or all
  arrange(pos_neg) %>%
  group_by(asv_id) %>%
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
                            mutate(c_time = ifelse(time == "T0", ifelse(susceptability == "S", c_time - 0.1, c_time + 0.1),
                                                   case_when(graph_cat == "D_S" ~ c_time - 0.35,
                                                             graph_cat == "H_S" ~ c_time - 0.15,
                                                             graph_cat == "D_R" ~ c_time + 0.15,
                                                             graph_cat == "H_R" ~ c_time + 0.35))) %>%
                            mutate(graph_cat = factor(graph_cat, levels = c("D_S", "D_R", "H_S", "H_R"), labels = c("D_S", "D_R", "H_S", "H_R"))))) %>%
  rowwise() %>%
  mutate(plot = list(
    ggplot(data = plot_info, aes(x = c_time, y = estimate, ymin = conf.low, ymax = conf.high, 
                                 shape = graph_cat)) +
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
      scale_color_manual(aesthetics = "colour1", values = c("#F75D5D", "#A70000"), guide = "legend", 
                         name = "Disease Exposed", labels = c("Susceptible", "Resistant")) +
      scale_shape_manual(values = c(16, 17, 16, 17), guide = "none") +
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
      labs(title = str_c(Family, " ", Genus, " (", asv_id, ")", sep = "")) 
  )) %>%
  group_by(pos_neg) %>%
  summarise(combo_plots = list(wrap_plots(plot) + plot_layout(guides = 'collect') & plot_annotation(title = pos_neg)))

  disease_corr_plots_all$combo_plots[[2]] 
  
  disease_corr_plots$combo_plots[[2]]
  
  
  
  
  
