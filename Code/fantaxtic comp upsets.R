fantaxtic_genus_colors <- fantaxtic_palette %>% 
  ungroup() %>%
  select(Genus, subgroup_colour)

test <- fantaxtic_palette %>% 
  ungroup() %>%
  select(Genus, subgroup_colour) %>%
  mutate(combo = str_c("'", Genus, "' = '", subgroup_colour, "'")) %>%
  pull(combo)




m_sig_classified_asvs <-  sig_classified_asvs %>% left_join(fantaxtic_genus_colors, by = join_by(Genus)) %>%
  mutate(subgroup_colour = ifelse(is.na(subgroup_colour), "gray75", subgroup_colour))
  


test_list <- c('**Other**' = 'white', 
               'Other' = 'gray75',
               'Thalassotalea' = '#6A3E15',
               'Alteromonas' = '#8B613B',
               'Pseudoalteromonas' = '#AC8661',
               'Psychrosphaera' = '#CDAA87',
               'Other Alteromonadales' = '#EFCEAE',
               'Seonamhaeicola' = '#F986C7',
               'Dokdonia' = '#FA99D0',
               'Allomuricauda' = '#FBACD9',
               'Tenacibaculum' = '#FCBFE2',
               'Other Flavobacteriales' = '#FDD3EB',
               'Oceanobacter' = '#0AAFC5',
               'Neptuniibacter' = '#34C0D2',
               'Alcanivorax' = '#5ED1E0',
               'Marinicella' = '#88E2ED',
               'Other Oceanospirillales' = '#B2F3FB',
               'Coraliomargarita' = '#E60ADA',
               'Ruficoccus' = '#F971F1',
               'Ruegeria' = '#F97900',
               'Thalassovita' = '#FA912D',
               'Shimia' = '#FCA95B',
               'Sulfitobacter' = '#FDC188',
               'Other Rhodobacterales' = '#FFDAB6',
               'MD3-55' = '#177500',
               'Candidatus Cyrtobacter' = '#45A32D',
               'Candidatus Megaira' = '#72D15B',
               'Other Rickettsiales' = '#A1FF8A',
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
               'Thiomicrorhabdus' = '#C95757',
               'Methylophaga' = '#E27F7F',
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
  rename("Healthy Outcome" = HealthyOutcome, "Diseased Outcome" = DiseaseOutcome, "Tank Effect" = TankEffect)


plot_colors <- m_sig_classified_asvs %>% select(plot_genus, subgroup_colour, outline) %>% distinct()


upset(m_sig_classified_asvs %>% select(-c(contains("DD"))),
      colnames(m_sig_classified_asvs %>% select(-c(asv_id, subgroup_colour, plot_genus, outline, group_subgroup, colnames(taxonomy_tibble %>% select(-asv_names)), contains("DD"))) %>%
                 relocate(`Tank Effect`, contains("Healthy"), contains("Outcome"))),
      base_annotations=list(
        'Intersection size'=intersection_size(counts=T, text = aes(size = 6.5),
                                              bar_number_threshold = 25,
                                              mapping=aes(fill=plot_genus, col = outline) #, col = "gray10"
         ) +
          scale_fill_manual(values = test_list, limits = m_sig_classified_asvs$plot_genus %>% unique()) + #plot_colors$subgroup_colour, labels = plot_colors$plot_genus
          scale_color_manual(values = c("yes" = "gray10", "no" = "white")) +
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
upset(m_sig_classified_asvs %>% select(-c(contains("DD"))),
      colnames(m_sig_classified_asvs %>% select(-c(asv_id, subgroup_colour, plot_genus, outline, group_subgroup, colnames(taxonomy_tibble %>% select(-asv_names)), contains("DD"))) %>%
                 relocate(`Tank Effect`, contains("Healthy"), contains("Outcome"))),
      base_annotations=list(
        'Intersection size'=intersection_size(counts=T, text = aes(size = 6.5),
                                              bar_number_threshold = 25,
                                              mapping=aes(fill=plot_genus) #, col = "gray10"
        ) +
          scale_fill_manual(values = test_list, limits = m_sig_classified_asvs$plot_genus %>% unique()) + #plot_colors$subgroup_colour, labels = plot_colors$plot_genus
          #scale_color_manual(values = c("yes" = "gray10", "no" = "white")) +
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

