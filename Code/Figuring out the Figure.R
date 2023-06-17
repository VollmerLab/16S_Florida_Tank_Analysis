
#TODO
#comp upset w color for families
#comp upset w directionality?
#examine definitions of signatures


categorized_asvs <- bacterial_signature_asv %>%
  arrange(signatures) %>%
  group_by(asv_id) %>%
  summarise(signatures = str_c(signatures, collapse = ', ')) %>%
  inner_join(significant_models,
             by = 'asv_id') %>%
  rowwise %>%
  mutate(plot = list(emmeans(model, ~treatment) %>%
                       broom::tidy(conf.int = TRUE) %>%
                       separate(treatment, into = c('time', 'exposure', 'susceptability')) %>%
                       ggplot(aes(x = time, y = estimate, ymin = conf.low, ymax = conf.high, 
                                  colour = exposure,
                                  shape = susceptability)) +
                       geom_pointrange(position = position_dodge(0.5)) +
                       labs(title = str_c(Family, Genus, asv_id, sep = ': '))))


categorized_asvs %>% filter(signatures != "probiotic") %>% select(Family, Genus) %>% distinct()


test <- categorized_asvs %>% ungroup() %>% dplyr::slice(1) %>% rowwise() %>%
  mutate(meanvals = list(emmeans(model, ~treatment) %>% broom::tidy(conf.int = TRUE) %>%
                       separate(treatment, into = c('time', 'exposure', 'susceptability')))) %>% select(meanvals) %>% 
  unnest(meanvals)

full_data %>% filter(asv_id %in% categorized_asvs$asv_id) %>% filter(time != "T0") %>% group_by(asv_id, genotype, time) %>%
  filter(asv_id == "ASV_103", genotype == "K1")
  mutate(change = log2_cpm[exposure == "D"] - log2_cpm[exposure == "H"])

mutate(diff = estimate[time == "T7"] - estimate[time == "T3"])


t33 <- categorized_asvs %>% 
  ungroup() %>%
  select(1:14) %>%
  rowwise() %>%
  mutate(test = list(emmeans(model, ~treatment) %>%
                       broom::tidy(conf.int = TRUE) %>%
                       separate(treatment, into = c('time', 'exposure', 'susceptability')) %>%
                       filter(time != "T0") %>% group_by(time, susceptability) %>% mutate(change = estimate[exposure == "D"] - estimate[exposure == "H"]))) %>%
  select(-c(data, model)) %>%
  unnest(test) %>%
  select(-c(estimate, std.error, df, statistic, exposure)) %>%
  group_by(time, susceptability, asv_id) %>%
  unique() %>% 
  group_by(asv_names) %>% 
  mutate(diff = estimate[time == "T7"] - estimate[time == "T3"]) %>%
  ungroup() %>%
  mutate(terms = factor(terms, levels = c("final_disease_state", "time:final_disease_state", "both"),
                        labels = c("FDS", "FDS:Time", "Both"))) %>%
  left_join(taxonomy_tibble, by = join_by(asv_names))


the_plots %>% formattable() %>% as.htmlwidget()

library(scatterD3)

library(scatterD3)
library(formattable)
library(dplyr)

scatter_cell <- function(x,y,...){
  as.character(
    htmltools::as.tags(
      scatterD3(x,y,...)
    )
  )
}

ft <- full_data %>%
  group_by(asv_id) %>%
  #summarize_each(funs(mean)) %>%
  mutate(
    Sepal.Scatter = scatter_cell(
      .$time,
      .$log2_cpm,
      width = 200,
      height = 200
    )
  ) %>%
  formattable() %>%
  as.htmlwidget()

ft$dependencies <- c(
  ft$dependencies,
  htmlwidgets:::widget_dependencies("scatterD3")
)

ft




#
  
  
testtt <- bacterial_signature_asv %>%
  arrange(signatures) %>%
  group_by(asv_id) %>%
  summarise(signatures = str_c(signatures, collapse = ', ')) %>%
  inner_join(significant_models,
             by = 'asv_id') %>%
  rowwise %>%
  mutate(plot = list(emmeans(model, ~treatment) %>%
                       broom::tidy(conf.int = TRUE) %>%
                       separate(treatment, into = c('time', 'exposure', 'susceptability'))))
         
          
prepped_testtt <- testtt%>% unnest(plot)   
t0dat <- subset(prepped_testtt, prepped_testtt$time == "T0")  
prepped_testtt_1 <- prepped_testtt %>% filter(time != "T0") %>% mutate(facet = exposure) %>%
  rbind(t0dat %>% mutate(facet = "D")) %>% rbind(t0dat %>% mutate(facet = "H"))
                     
ggplot(prepped_testtt_1, aes(x = asv_id, y = estimate, ymin = conf.low, ymax = conf.high, 
          colour = time,
          shape = susceptability)) +
geom_pointrange(position = position_dodge(0.5)) +
  coord_flip() +
facet_grid(rows = vars(signatures), cols = vars(facet), scales = "free", space = "free")

## other stuff above

library(relayer)



####
the_plots1 <- bacterial_signature_asv %>%
  #filter(asv_id %in% c("ASV_85", "ASV_72")) %>%
  arrange(signatures) %>%
  group_by(asv_id) %>%
  summarise(signatures = str_c(signatures, collapse = ', ')) %>%
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
       slice(rep(1:2, 1)) %>%
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
      scale_linetype_manual(values = c(6, 1, 6, 1), guide = "none") +
      theme_bw() +
      xlab("Time") +
      ylab(expression("Normalized log"[2]*" (cpm)")) +
      labs(title = str_c(Family, " ", Genus, " (", asv_id, ")", sep = "")) 
  )) %>%
  group_by(signatures) %>%
  summarise(combo_plots = list(wrap_plots(plot) + plot_layout(guides = 'collect') & plot_annotation(title = signatures)))

the_plots1$combo_plots[[1]]
the_plots1$combo_plots[[2]]
the_plots1$combo_plots[[3]]
the_plots1$combo_plots[[4]]
the_plots1$combo_plots[[5]]
the_plots1$combo_plots[[6]]

       


#####

the_plots1 <- bacterial_signature_asv %>%
  filter(asv_id == "ASV_85") %>%
  arrange(signatures) %>%
  group_by(asv_id) %>%
  summarise(signatures = str_c(signatures, collapse = ', ')) %>%
  inner_join(significant_models,
             by = 'asv_id') %>%
  rowwise %>%
  mutate(plot = list(emmeans(model, ~treatment) %>%
     broom::tidy(conf.int = TRUE) %>%
     separate(treatment, into = c('time', 'exposure', 'susceptability')) %>%
     mutate(graph_cat = ifelse(time == "T0", NA, 
                               paste(exposure, susceptability, sep = "_"))) %>%
     {. ->> intermed } %>%
     mutate(graph_cat = ifelse(time == "T0", paste("H", susceptability, sep = "_"), 
                               graph_cat)) %>%
     slice(rep(1:2, 1)) %>%
     rbind(intermed) %>%
     mutate(graph_cat = ifelse(is.na(graph_cat), paste("D", susceptability, sep = "_"), 
                               graph_cat)) %>%
     mutate(c_time = parse_number(time)) %>%
     mutate(c_time = ifelse(time == "T0", ifelse(susceptability == "S", c_time - 0.1, c_time + 0.1),
                            case_when(graph_cat == "D_S" ~ c_time - 0.35,
                                      graph_cat == "H_S" ~ c_time - 0.15,
                                      graph_cat == "D_R" ~ c_time + 0.15,
                                      graph_cat == "H_R" ~ c_time + 0.35))) %>%
     mutate(graph_cat = factor(graph_cat, levels = c("D_S", "D_R", "H_S", "H_R"))) %>%
     ggplot(aes(x = c_time, y = estimate, ymin = conf.low, ymax = conf.high, 
                colour = graph_cat,
                shape = graph_cat)) +
     geom_line(aes(linetype = graph_cat)) +
     geom_errorbar(linewidth = 0.7) +
     geom_point(size = 3) +
     scale_color_manual(values = c("#F75D5D", "#A70000","#3DD8EA", "#048291")) +
     scale_shape_manual(values = c(16, 17, 16, 17)) +
     scale_x_continuous(breaks=c(0, 3, 7)) +
     scale_linetype_manual(values = c("dashed", "solid", "dashed", "solid")) +
     theme_bw() + 
     xlab("Time") +
     labs(title = str_c(Family, " ", Genus, " (", asv_id, ")", sep = "")))) #%>%
the_plots1$plot[1]

  group_by(signatures) %>%
  summarise(plots = list(wrap_plots(plot) + plot_layout(guides = 'collect') & plot_annotation(title = signatures)))


###
  library(stringr)
  rand_forest_asvs <- c(149, 26, 96, 85, 3, 82, 161, 165, 284, 40, 479, 142)

  rand_forest_asvs <- rand_forest_asvs %>%
    as.character() %>%
    as_tibble() %>%
    mutate(value = paste("ASV_", value, sep = "")) %>%
    pull(value)

  sig_asvs_nm
  
  both_models <- c(rand_forest_asvs, sig_asvs_nm) %>% unique() %>% as.tibble() %>%
    mutate(random_forest = ifelse(value %in% rand_forest_asvs, TRUE, FALSE)) %>%
    mutate(planned_comparisons = ifelse(value %in% sig_asvs_nm, TRUE, FALSE))
  
  ggvenn(both_models, fill_color = c("green", "#FF9300"))







