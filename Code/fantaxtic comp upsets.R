
order_placeholders <- fantaxtic_palette %>% mutate(Genus = NA, group_colour = "white", subgroup_colour = "white", group_subgroup = NA) %>%
  distinct() %>% mutate(Genus = sprintf("**%s**", as.character(Order)))

full_palette <- rbind(fantaxtic_palette, order_placeholders)

fantaxtic_testpalette <- full_palette$subgroup_colour
names(fantaxtic_testpalette) <- full_palette$Genus

m_sig_classified_asvs <-  sig_classified_asvs %>% left_join(fantaxtic_genus_colors, by = join_by(Genus)) %>%
  mutate(subgroup_colour = ifelse(is.na(subgroup_colour), "gray75", subgroup_colour))
  


test_list <- c('**Other**' = 'white', 
               'Other' = 'gray75',
               'Thalassotalea' = '#6A3E15',
               'Alteromonas' = '#8B613B',
               'Pseudoalteromonas' = '#AC8661',
               'Algicola' = '#CDAA87',
               'Other Alteromonadales' = '#EFCEAE',
               'Seonamhaeicola' = '#F986C7',
               'Dokdonia' = '#FA99D0',
               'Allomuricauda' = '#FBACD9',
               'Flavobacterium' = '#FCBFE2',
               'Other Flavobacteriales' = '#FDD3EB',
               'Oleiphilus' = '#0AAFC5',
               'Oceanobacter' = '#34C0D2',
               'Neptuniibacter' = '#5ED1E0',
               'Alcanivorax' = '#88E2ED',
               'Other Oceanospirillales' = '#B2F3FB',
               'Coraliomargarita' = '#7E0277',
               'Verruc-01' = "#9D2997",
               'Lentimonas' = "#BD51B7",
               'Cerasicoccus' = "#DD78D7",
               'Other Puniceicoccales' = '#FDA0F8',
               'Ruegeria' = '#F97900',
               'Thalassovita' = '#FA912D',
               'Shimia' = '#FCA95B',
               'Sulfitobacter' = '#FDC188',
               'Other Rhodobacterales' = '#FFDAB6',
               'MD3-55' = '#146600',
               'Rickettsia' = '#3A8C26',
               'Candidatus Cyrtobacter' = '#60B24C',
               'Candidatus Megaira' = '#86D872',
               'Other Rickettsiales' = '#ADFF99',
               'Aureispira' = '#171FB8',
               'Lewinella' = '#3F45C7',
               'Phaeodactylibacter' = '#676CD7',
               'Saprospira' = '#8F93E6',
               'Other Saprospirales' = '#B7BAF6',
               'Spirochaeta 2' = '#8D00D3',
               'Treponema' = '#B245E9',
               'Other Spirochaetales' = '#D88AFF',
               'Cysteiniphilum' = '#970707',
               'Francisella' = '#B02F2F',
               'Methylophaga' = '#C95757',
               'Candidatus Endoecteinascidia' = '#E27F7F',
               'Other Thiotrichales' = '#FBA8A8',
               'Haloferula' = '#A3CE0B',
               'Rubritalea' = '#B5D935',
               'Fucophilus' = '#C7E45F',
               'Roseibacillus' = '#D9EF89',
               'Other Verrucomicrobiales' = '#ECFBB4',
               'Vibrio' = '#C44600',
               'Photobacterium' = '#D76E34',
               'Enterovibrio' = '#EB9768',
               'Other Vibrionales' = '#FFC09D',
               '**Alteromonadales**' = 'white',
               '**Flavobacteriales**' = 'white',
               '**Oceanospirillales**' = 'white',
               '**Puniceicoccales**' = 'white',
               '**Rhodobacterales**' = 'white',
               '**Rickettsiales**' = 'white',
               '**Saprospirales**' = 'white',
               '**Spirochaetales**' = 'white',
               '**Thiotrichales**' = 'white',
               '**Verrucomicrobiales**' = 'white',
               '**Vibrionales**' = 'white')

               
               



m_sig_classified_asvs <-  sig_classified_asvs %>% 
  left_join(fantaxtic_genus_colors, by = join_by(Genus)) %>%
  mutate(Genus = ifelse(is.na(subgroup_colour), "Other", Genus),
         subgroup_colour = ifelse(is.na(subgroup_colour), "gray75", subgroup_colour))

ordered_fantaxtic_factoring <- fantaxtic_palette %>%
  group_by(Order) %>% 
  group_modify(~add_row(.x, .before = 0)) %>%
  ungroup() %>%
  mutate(subgroup_colour = ifelse(is.na(subgroup_colour), "#FFFFFF", subgroup_colour),
         Genus = ifelse(is.na(Genus), sprintf("**%s**", as.character(Order)), as.character(Genus)),
         group_subgroup = ifelse(is.na(group_subgroup), Genus, group_subgroup))

test123 <- ordered_fantaxtic_factoring$subgroup_colour
names(test123) <- ordered_fantaxtic_factoring$Genus

m_sig_classified_asvs <- sig_classified_asvs %>% 
  left_join(fantaxtic_palette %>% select(-group_colour), by = join_by(Order, Genus)) %>%
  mutate(plot_genus = Genus) %>%
  mutate(plot_genus = case_when(is.na(subgroup_colour) & Order == "Alteromonadales" ~ "Other Alteromonadales",
                                is.na(subgroup_colour) & Order == "Flavobacteriales" ~ "Other Flavobacteriales",
                                is.na(subgroup_colour) & Order == "Oceanospirillales" ~ "Other Oceanospirillales",
                                is.na(subgroup_colour) & Order == "Rhodobacterales" ~ "Other Rhodobacterales",
                                is.na(subgroup_colour) & Order == "Saprospirales" ~ "Other Saprospirales",
                                is.na(subgroup_colour) & Order == "Spirochaetales" ~ "Other Spirochaetales",
                                is.na(subgroup_colour) & Order == "Thiotrichales" ~ "Other Thiotrichales",
                                is.na(subgroup_colour) & Order == "Verrucomicrobiales" ~ "Other Verrucomicrobiales",
                                is.na(subgroup_colour) & Order == "Vibrionales" ~ "Other Vibrionales",
                                is.na(subgroup_colour) & Order == "Rickettsiales" ~ "Other Rickettsiales",
                                TRUE ~ plot_genus)) %>%
  mutate(plot_genus = ifelse(plot_genus %in% fantaxtic_genus_colors$Genus, plot_genus, "Other")) %>%
  select(-subgroup_colour) %>%
  left_join(fantaxtic_genus_colors %>% rename("plot_genus" = Genus), by = join_by(plot_genus)) %>%
  relocate(c(plot_genus, subgroup_colour), .after = asv_id) %>%
  mutate(Order = ifelse(Order %in% fantaxtic_palette$Order, Order, "Other")) %>%
  group_by(Order) %>% 
  group_modify(~add_row(.x, .before = 0)) %>%
  ungroup() %>%
  mutate(subgroup_colour = ifelse(is.na(subgroup_colour), "#FFFFFF", subgroup_colour),
         plot_genus = ifelse(is.na(plot_genus), sprintf("**%s**", as.character(Order)), as.character(plot_genus)),
         group_subgroup = ifelse(is.na(group_subgroup), plot_genus, group_subgroup)) %>%
  mutate(outline = ifelse(str_detect(plot_genus, "\\*"), "no", "yes")) %>%
  mutate(across(DD_early:TankEffect, ~ifelse(is.na(.), FALSE, .))) %>%
  mutate(group_subgroup = ifelse(str_detect(group_subgroup, "Other") & !str_detect(group_subgroup, "\\*"), str_c(Order, " - ", plot_genus), group_subgroup)) %>%
  mutate(group_subgroup = factor(group_subgroup, ordered = T, levels = ordered_fantaxtic_factoring$group_subgroup),
         plot_genus = factor(plot_genus, ordered = T, levels = ordered_fantaxtic_factoring$Genus),
         Order = factor(Order, ordered = T)) %>%
  arrange(group_subgroup) %>%
  rename("Healthy Outcome" = HealthyOutcome, "Diseased Outcome" = DiseaseOutcome, 
         "Tank Effect" = TankEffect,
         "Healthy Exposure" = HealthyExposed, "Diseased Exposure" = DiseaseExposed)





upset(m_sig_classified_asvs %>% select(-c(contains("DD"))),
      colnames(m_sig_classified_asvs %>% select(-c(asv_id, subgroup_colour, plot_genus, outline, group_subgroup, colnames(taxonomy_tibble), contains("DD"))) %>%
                 relocate(`Tank Effect`, contains("Healthy"), contains("Outcome"))),
      base_annotations=list(
        'Intersection size'=intersection_size(counts=T, text = aes(size = 6.5),
                                              bar_number_threshold = 25,
                                              mapping=aes(fill=plot_genus, col = outline) #, col = "gray10"
         ) +
          #scale_fill_manual(values = test_list, limits = m_sig_classified_asvs$plot_genus %>% unique()) +
          scale_fill_manual(values = test123, limits = m_sig_classified_asvs$plot_genus %>% unique()) +
          scale_color_manual(values = c("yes" = "red", "no" = "white")) +
          theme_nested(legend.position=c(1.06,0.2)) +
          guides(fill=guide_legend(title=substitute(bold(bd)~nb, list(bd = "Order", nb = "/ Genus")),
                                   ncol = 1))
      ),
      
      matrix=(
        intersection_matrix(geom=geom_point(shape = "circle filled", size=4, stroke = 0.5, color = "gray10"))
        + scale_color_manual(
          values=c(
            #"Field" = "#70d134",
            #"HealthyExposed" = "#78acff",
            "Healthy Outcome" = "#0d50ba",
            "Tank Effect"= "#c389e0",
            #"DiseaseExposed"= "#ff7878",
            "Diseased Outcome"= "#ba0d0d"
          )
        )
      ),
      queries=list(
        #upset_query(set = "Field", fill = "#70d134"),
        #upset_query(set = "HealthyExposed", fill = "#78acff"),
        upset_query(set = "Healthy Outcome", fill = "#0d50ba"),
        upset_query(set = "Tank Effect", fill = "#c389e0"),
        #upset_query(set = "DiseaseExposed", fill = "#ff7878"),
        upset_query(set = "Diseased Outcome", fill = "#ba0d0d")
      ),
      name='ASVs', width_ratio=0.1, sort_sets = FALSE, min_degree = 1, sort_intersections=FALSE, 
      set_sizes=FALSE,
      #set_sizes=upset_set_size(geom = geom_point(stat = "count", shape = "square", size = 3, color = c("#ba0d0d", "#c389e0", "#0d50ba"))) + theme(axis.ticks.x=element_line()),
      intersections=list(
        'Diseased Outcome',
        'Healthy Outcome',
        c('Diseased Outcome', 'Tank Effect'),
        c('Healthy Outcome', 'Tank Effect'),
        'Tank Effect'
      )) +
  ggtitle("Significant Main Effects") +
  theme(plot.margin = margin(1, 5.5, 1, 1, "cm"))


#no outline
no_outline_fantaxtic_upset <- upset(m_sig_classified_asvs %>% select(-c(contains("DD"))),
      colnames(m_sig_classified_asvs %>% select(-c(asv_id, subgroup_colour, plot_genus, outline, group_subgroup, colnames(taxonomy_tibble %>% select(-asv_names)), contains("DD"))) %>%
                 relocate(`Tank Effect`, `Diseased Exposure`, `Healthy Exposure`, 
                          `Diseased Outcome`, `Healthy Outcome`)),
      base_annotations=list(
        'Intersection size'=intersection_size(counts=T, text = aes(size = 6.5),
                                              bar_number_threshold = 25,
                                              mapping=aes(fill=plot_genus) #, col = "gray10"
        ) +
          scale_fill_manual(values = test_list, limits = m_sig_classified_asvs$plot_genus %>% unique()) + 
          #scale_color_manual(values = c("yes" = "gray10", "no" = "white")) +
          theme_nested() + #legend.position=c(1.06,0.2)
          guides(fill=guide_legend(title=substitute(bold(bd)~nb, list(bd = "Order", nb = "/ Genus")),
                                   ncol = 1)) + 
          theme(legend.position = "none",
                axis.text = element_text(size=12),
                axis.title = element_text(size=13))
      ),
      
      matrix=(
        intersection_matrix(geom=geom_point(shape = "circle filled", size=4, stroke = 0.5, color = "gray10"))
        + scale_color_manual(
          values=c(
            #"Field" = "#70d134",
            "Healthy Exposure" = "#78acff",
            "Healthy Outcome" = "#0d50ba",
            "Tank Effect"= "#c389e0",
            "Disease Exposure"= "#ff7878",
            "Diseased Outcome"= "#ba0d0d"
          )
        )
      ),
      queries=list(
        #upset_query(set = "Field", fill = "#70d134"),
        upset_query(set = "Healthy Exposure", fill = "#78acff"),
        upset_query(set = "Healthy Outcome", fill = "#0d50ba"),
        upset_query(set = "Tank Effect", fill = "#c389e0"),
        upset_query(set = "Disease Exposure", fill = "#ff7878"),
        upset_query(set = "Diseased Outcome", fill = "#ba0d0d")
      ),
      name='ASVs', width_ratio=0.1, sort_sets = FALSE, min_degree = 1, sort_intersections=FALSE, 
      set_sizes=FALSE, keep_empty_groups=TRUE,
      #set_sizes=upset_set_size(geom = geom_point(stat = "count", shape = "square", size = 3, color = c("#ba0d0d", "#c389e0", "#0d50ba"))) + theme(axis.ticks.x=element_line()),
      intersections=list(
        'Diseased Outcome',
        'Healthy Outcome',
        c('Diseased Outcome', 'Tank Effect'),
        c('Healthy Outcome', 'Tank Effect'),
        'Tank Effect'
      )) +
  ggtitle("Significant Main Effects") +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12))

no_outline_fantaxtic_upset
  #theme(plot.margin = margin(1, 5.5, 1, 1, "cm"))

#export 2000x1100

legend_only_cu <- ggplot(m_sig_classified_asvs %>% select(-c(contains("DD")))) +
  geom_bar(mapping=aes(`Diseased Outcome`, fill=plot_genus), stat='count', position='fill') +
  scale_fill_manual(values = test_list, limits = m_sig_classified_asvs$plot_genus %>% unique()) + 
  #scale_color_manual(values = c("yes" = "gray10", "no" = "white")) +
  theme_nested(legend.text = element_text(size=10),
               legend.title = element_text(size=12)) +
  guides(fill=guide_legend(title=substitute(bold(bd)~nb, list(bd = "Order", nb = "/ Genus")),
                           ncol = 2))

savethis <- ggpubr::get_legend(legend_only_cu)


(no_outline_fantaxtic_upset | savethis) + plot_layout(widths = c(6, 2))
#export 1500x900


#### no tank effect ####

add_dose_data <- homogenate_models %>%
  mutate(dose = ifelse(adj_homog_p_val < 0.05, TRUE, FALSE)) %>%
  select(asv_id, dose)

disease_criteria_upset_prelim <- m_sig_classified_asvs %>%
  left_join(add_dose_data, by = join_by("asv_id")) %>%
  filter(!`Tank Effect`) %>%
  select(-`Tank Effect`) %>%
  rename("Early Responder" = "DD_early", "Late Responder" = "DD_late", "Dose" = "dose",
         "Outcome Posthoc" = "DD_vs_DH_late", "Continuous Responder" = "DD_continuous")

tank_only_labels <- disease_criteria_upset %>% group_by(Order) %>% reframe(n = n()) %>% filter(n == 1) %>% pull(Order)

disease_criteria_upset <- disease_criteria_upset_prelim %>%
  filter(!Order %in% tank_only_labels) %>%
  mutate(`Late Responder` = ifelse(`Healthy Outcome`, FALSE, `Late Responder`)) #we only intended to apply this filter to disease-associated ASVs

get_cu_legend <- ggplot(disease_criteria_upset) +
  geom_bar(aes(x = Order, fill=plot_genus))+
  scale_fill_manual(values = test_list, limits = disease_criteria_upset$plot_genus %>% unique()) + 
  #scale_color_manual(values = c("yes" = "gray10", "no" = "white")) +
  theme_nested() +
  guides(fill=guide_legend(title=substitute(bold(bd)~nb, list(bd = "Order", nb = "/ Genus")),
                           ncol = 1))

fantax_cu_legend <- cowplot::get_legend(get_cu_legend)



fantax_cu <- upset(disease_criteria_upset,
      colnames(disease_criteria_upset %>% select(-c(asv_id, subgroup_colour, plot_genus, outline, `Continuous Responder`, group_subgroup, colnames(taxonomy_tibble %>% select(-asv_names)))) %>%
                 relocate(`Dose`,
                           `Outcome Posthoc`, `Diseased Outcome`, `Late Responder`,
                          `Early Responder`,
                          `Healthy Outcome`,
                          `Diseased Exposure`, `Healthy Exposure`)),
      base_annotations=list(
        'Intersection size'=intersection_size(counts=T, text = aes(size = 6.5),
                                              bar_number_threshold = 25,
                                              mapping=aes(fill=plot_genus) #, col = "gray10"
        ) +
          scale_fill_manual(values = test_list, limits = disease_criteria_upset$plot_genus %>% unique()) + 
          #scale_color_manual(values = c("yes" = "gray10", "no" = "white")) +
          theme_nested(legend.position=c(1.06,0.2)) +
          guides(fill = "none")
          #guides(fill=guide_legend(title=substitute(bold(bd)~nb, list(bd = "Order", nb = "/ Genus")),
          #                         ncol = 1))
      ),
      
      matrix=(
        intersection_matrix(geom=geom_point(shape = "circle filled", size=4, stroke = 0.5, color = "gray10"))
        + scale_color_manual(
          values=c(
            #"Healthy Exposure" = "#78acff",
            "Healthy Outcome" = "#0d50ba",
            #"Disease Exposure"= "#FFB37A",
            "Diseased Outcome"= "#FF7D1A",
            "Early Responder" = "#B9E989",
            #"Continuous Responder" = "#7CBF3A",
            "Late Responder" = "#4B870E",
            "Outcome Posthoc" = "#A70000",
            "Dose" = "#732dcf"
          )
        )
      ),
      queries=list(
       # upset_query(set = "Healthy Exposure", fill = "#78acff"),
        upset_query(set = "Healthy Outcome", fill = "#0d50ba"),
        #upset_query(set = "Disease Exposure", fill = "#FFB37A"),
        upset_query(set = "Diseased Outcome", fill = "#FF7D1A"),
        upset_query(set = "Outcome Posthoc", fill = "#BA0D0D"),
        upset_query(set = "Early Responder", fill = "#FFEFA5"),
        #upset_query(set = "Continuous Responder", fill = "#7CBF3A"),
        upset_query(set = "Late Responder", fill = "#FED20D"),
        upset_query(set = "Dose", fill = "#732dcf")
      ),
      name='ASVs', width_ratio=0.1, sort_sets = FALSE, min_degree = 1, sort_intersections=FALSE, 
      set_sizes=FALSE, keep_empty_groups=TRUE,
      intersections = list(c("Late Responder", "Diseased Outcome", "Outcome Posthoc", "Dose"),
                           c("Late Responder", "Diseased Outcome", "Outcome Posthoc"),
                           c("Early Responder", "Diseased Outcome", "Dose"),
                           c("Late Responder", "Diseased Outcome"),
                           c("Early Responder", "Diseased Outcome"),
                           "Healthy Outcome")
      #set_sizes=upset_set_size(geom = geom_point(stat = "count", shape = "square", size = 3, color = c("#ba0d0d", "#c389e0", "#0d50ba"))) + theme(axis.ticks.x=element_line()),
      #intersections=list('Diseased Outcome','Healthy Outcome',c('Diseased Outcome', 'Tank Effect'))
      ) +
  ggtitle("Significant Effects") #+
  theme(plot.margin = margin(1, 5.5, 1, 1, "cm"))

fantax_cu | fantax_cu_legend  
  

design = "
A#
AB
A#
"

list(
  fantax_cu, # A
  fantax_cu_legend # B
) %>% 
  wrap_plots() + 
  plot_layout(heights = c(0.125, 1, 0.125), widths = c(200, 50), design = design)
#export 1300x800
  
  




