
clean_bacstrats <- bacterial_signature_asv %>% 
  group_by(asv_id) %>% 
  reframe(bacstrat = paste(list(signatures))) %>%
  mutate(bacstrat = case_when(bacstrat == "c(\"early_pathogen\", \"continuous_pathogen\", \"late_pathogen\")" ~ "continuous_pathogen",
                              bacstrat == 'c("crasher_t3", "crasher_t7")' ~ "continuous_crasher",
                              bacstrat == 'c("early_pathogen", "early_opportunist")' ~ "early_pathogen",
                              bacstrat == 'c("late_pathogen", "late_opportunist")' ~ "late_pathogen",
                              bacstrat == "c(\"probiotic_t7_strict\", \"late_opportunist\")" ~ "late_probiotic",
                              bacstrat == "crasher_t7" ~ "late_crasher",
                              bacstrat == "probiotic_t7_strict" ~ "late_probiotic",
                              bacstrat == 'c("early_pathogen", "late_opportunist")' ~ "early_pathogen",
                              TRUE ~ bacstrat))


#make plot
scores(test_nmds)$sites %>%
  as_tibble(rownames = 'asv_names') %>%
  #left_join(metadata, by = 'asv_names') %>%
  left_join(clean_bacstrats, by = c('asv_names' = 'asv_id')) %>%
  filter(!is.na(bacstrat)) %>%
  #mutate(bacstrat = ifelse(is.na(bacstrat), "disease_unrelated", bacstrat)) %>%
  mutate(shape_determination = ifelse(str_detect(bacstrat, "disease"), 
        "disease_unrelated", ifelse(str_detect(bacstrat, "path"), "pathogen", 
        "opportunist"))) %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_point(data = as_tibble(scores(test_nmds)$species, rownames = aggregation_level) %>% mutate(timepoint = str_extract_numbers(none)),
             size = 1, shape = 1, color = "gray90") +
  geom_text(data = as_tibble(scores(test_nmds)$species, rownames = aggregation_level) %>% mutate(timepoint = str_extract_numbers(none)),
              aes(label = none), size = 3) +
  geom_point(aes(col = bacstrat, shape = shape_determination), size = 3) + #alpha = shape_determination
  #scale_alpha_manual(values = c(0.3, 1,1)) +
  scale_color_manual(values = c(
    "early_pathogen" = "#FF7575",
    "late_pathogen" = "#CF0C0C",
    "early_opportunist" = "lightskyblue",
    "late_opportunist" = "deepskyblue3",
    "late_opportunist&pathogen" = "purple",
    "continuous_opportunist" = "royalblue4",
    "early_path_cont_opp" = "hotpink3",
    "early_path_late_opp" = "orange"
    #"disease_unrelated" = "gray"
  )) +
  theme_bw()

#plot species or site alone
  plot(test_nmds, "species")
  orditorp(test_nmds, "species")
  
  
  
  asv_nmds <- full_data %>% 
    select(asv_id, sample_id, log2_cpm) %>%
    pivot_wider(names_from = sample_id, values_from = log2_cpm) %>%
    left_join(clean_bacstrats, by = join_by("asv_id")) %>%
    filter(!is.na(bacstrat)) %>%
    select(-bacstrat) %>%
    column_to_rownames('asv_id') %>%
    t() %>%
    as.matrix()
  
  test_nmds <- metaMDS(t(asv_nmds), distance = 'bray', k = 2, trymax = 100, autotransform = FALSE, verbose = TRUE)
  sppscores(test_nmds) <- t(asv_nmds)
  
  #shepard plot
  plot(test_nmds$diss, test_nmds$dist)
  
suscep <- full_data %>% 
    select(asv_id, sample_id, resistance) %>%
    pivot_wider(names_from = sample_id, values_from = resistance) %>%
    left_join(clean_bacstrats, by = join_by("asv_id")) %>%
    filter(!is.na(bacstrat)) %>%
    select(-bacstrat) %>%
    column_to_rownames('asv_id') %>%
    t() %>%
    as.matrix()
  
  #adonis2 fails but adonis works
  perm.results <- vegan::adonis(asv_nmds ~ suscep, method="bray",perm=999)
  perm.results$aov.tab #p value less than 0.05
  
  
  
  
  not0_asv_nmds <- full_data %>% 
    select(asv_id, sample_id, log2_cpm) %>%
    #mutate(log2_cpm = str_trim(as.numeric(log2_cpm))) %>%
    pivot_wider(names_from = sample_id, values_from = log2_cpm) %>%
    left_join(clean_bacstrats, by = join_by("asv_id")) %>%
    filter(!is.na(bacstrat)) %>%
    select(-bacstrat) %>%
    select(!contains("T0")) %>%
    column_to_rownames('asv_id') %>%
    t() %>%
    as.matrix()
  
  
  not0_suscep <- full_data %>% 
    select(asv_id, sample_id, resistance) %>%
    #mutate(log2_cpm = str_trim(as.numeric(log2_cpm))) %>%
    pivot_wider(names_from = sample_id, values_from = resistance) %>%
    left_join(clean_bacstrats, by = join_by("asv_id")) %>%
    filter(!is.na(bacstrat)) %>%
    select(-bacstrat) %>%
    select(!contains("T0")) %>%
    column_to_rownames('asv_id') %>%
    t() %>%
    as.matrix()  
  
  not0_perm.results <- vegan::adonis(not0_asv_nmds ~ not0_suscep, method="bray",perm=999)
  not0_perm.results$aov.tab #p value less than 0.05
  
  
  not0_test_nmds <- metaMDS(t(not0_asv_nmds), distance = 'bray', k = 2, trymax = 100, autotransform = FALSE, verbose = TRUE)
  
  
  plot(not0_test_nmds)
  orditorp(not0_test_nmds, "site")
  
   scores(test_nmds)$site %>%
    as_tibble(rownames = 'asv_names') %>%
    left_join(clean_bacstrats, by = c('asv_names' = 'asv_id')) %>%
    filter(!is.na(bacstrat)) %>%
    mutate(bacstrat = case_when(bacstrat == "c(\"early_pathogen\", \"continuous_pathogen\", \"late_pathogen\")" ~ "continuous_pathogen",
                                  bacstrat == 'c("crasher_t3", "crasher_t7")' ~ "continuous_crasher",
                                  bacstrat == 'c("early_pathogen", "early_opportunist")' ~ "early_pathogen",
                                  bacstrat == 'c("late_pathogen", "late_opportunist")' ~ "late_pathogen",
                                  bacstrat == "c(\"probiotic_t7_strict\", \"late_opportunist\")" ~ "late_probiotic",
                                  bacstrat == "crasher_t7" ~ "late_crasher",
                                  bacstrat == "probiotic_t7_strict" ~ "late_probiotic",
                                  bacstrat == 'c("early_pathogen", "late_opportunist")' ~ "early_pathogen",
                                  TRUE ~ bacstrat)) %>%
    mutate(shape_determination = case_when(
     str_detect(bacstrat, "path")  ~ "pathogen",
     str_detect(bacstrat, "opp")  ~ "opportunist",
     str_detect(bacstrat, "crash")  ~ "crasher",
     str_detect(bacstrat, "pro")  ~ "probiotic"
   )) %>%
   ggplot(aes(x = NMDS1, y = NMDS2)) +
    stat_ellipse(data = as_tibble(scores(test_nmds)$species, rownames = aggregation_level) %>% 
                   mutate(exp_conditions = str_split(none, "_")) %>%
                   rowwise() %>%
                   mutate(time = exp_conditions[1], exposure = exp_conditions[2], resist = exp_conditions[3]), 
                 geom = "polygon", alpha = 0.1, level = 0.95, aes(col = time, fill = time)) +
    geom_point(data = as_tibble(scores(test_nmds)$species, rownames = aggregation_level) %>% 
                 mutate(exp_conditions = str_split(none, "_")) %>%
                 rowwise() %>%
                 mutate(time = exp_conditions[1], exposure = exp_conditions[2], resist = exp_conditions[3]),
               size = 1, aes(color = time, shape = exposure)) +
    geom_point(aes(col = bacstrat, shape = shape_determination), size = 4) + #alpha = shape_determination
    geom_text(aes(label = parse_number(asv_names)), size = 3) +
    scale_shape_manual(values = c(
      "opportunist" = 16,
      "pathogen" = 17,
      "crasher" = 15,
      "probiotic" = 8,
      "H" = 6,
      "D" = 4,
      "REP1" = 3
    )) +
    scale_fill_manual(values = c(
      "T0" = "gray80",
      "T3" = "gray40",
      "T7" = "gray10"
    ), guide = "none") +
    scale_color_manual(values = c(
      "continuous_crasher" = "sandybrown",
      "continuous_pathogen" = "firebrick4",
      "late_crasher" = "darkorange2",
      "late_probiotic" = "aquamarine",
      "early_pathogen" = "hotpink",
      "late_pathogen" = "firebrick1",
      "early_opportunist" = "deepskyblue3",
      "late_opportunist" = "royalblue2",
      "T0" = "gray80",
      "T3" = "gray40",
      "T7" = "gray10"
    )) +
    guides(color = guide_legend(override.aes = list(fill = NA))) +
    theme_bw()
 
  
   
   
   scores(test_nmds)$site %>%
     as_tibble(rownames = 'asv_names') %>%
     left_join(clean_bacstrats, by = c('asv_names' = 'asv_id')) %>%
     filter(!is.na(bacstrat)) %>%
     mutate(bacstrat = case_when(bacstrat == "c(\"early_pathogen\", \"continuous_pathogen\", \"late_pathogen\")" ~ "continuous_pathogen",
                                 bacstrat == 'c("crasher_t3", "crasher_t7")' ~ "continuous_crasher",
                                 bacstrat == 'c("early_pathogen", "early_opportunist")' ~ "early_pathogen",
                                 bacstrat == 'c("late_pathogen", "late_opportunist")' ~ "late_pathogen",
                                 bacstrat == "c(\"probiotic_t7_strict\", \"late_opportunist\")" ~ "late_probiotic",
                                 bacstrat == "crasher_t7" ~ "late_crasher",
                                 bacstrat == "probiotic_t7_strict" ~ "late_probiotic",
                                 bacstrat == 'c("early_pathogen", "late_opportunist")' ~ "early_pathogen",
                                 TRUE ~ bacstrat)) %>%
     mutate(shape_determination = case_when(
       str_detect(bacstrat, "path")  ~ "pathogen",
       str_detect(bacstrat, "opp")  ~ "opportunist",
       str_detect(bacstrat, "crash")  ~ "crasher",
       str_detect(bacstrat, "pro")  ~ "probiotic"
     )) %>%
     ggplot(aes(x = NMDS1, y = NMDS2)) +
     stat_ellipse(geom = "polygon", alpha = 0.1, level = 0.95, aes(col = bacstrat, fill = bacstrat)) +
     geom_point(data = as_tibble(scores(test_nmds)$species, rownames = aggregation_level) %>% 
                  mutate(exp_conditions = str_split(none, "_")) %>%
                  rowwise() %>%
                  mutate(time = exp_conditions[1], exposure = exp_conditions[2], resist = exp_conditions[3]),
                size = 1, aes(color = time, shape = exposure)) +
     geom_point(aes(col = bacstrat, shape = shape_determination), size = 4) + #alpha = shape_determination
     geom_text(aes(label = parse_number(asv_names)), size = 3) +
     scale_shape_manual(values = c(
       "opportunist" = 16,
       "pathogen" = 17,
       "crasher" = 15,
       "probiotic" = 8,
       "H" = 6,
       "D" = 4,
       "REP1" = 3
     )) +
     scale_fill_manual(values = c(
       "continuous_crasher" = "sandybrown",
       "continuous_pathogen" = "firebrick4",
       "late_crasher" = "darkorange2",
       "late_probiotic" = "aquamarine",
       "early_pathogen" = "hotpink",
       "late_pathogen" = "firebrick1",
       "early_opportunist" = "deepskyblue3",
       "late_opportunist" = "royalblue2"
     ), guide = "none") +
     scale_color_manual(values = c(
       "continuous_crasher" = "sandybrown",
       "continuous_pathogen" = "firebrick4",
       "late_crasher" = "darkorange2",
       "late_probiotic" = "aquamarine",
       "early_pathogen" = "hotpink",
       "late_pathogen" = "firebrick1",
       "early_opportunist" = "deepskyblue3",
       "late_opportunist" = "royalblue2",
       "T0" = "gray80",
       "T3" = "gray40",
       "T7" = "gray10"
     )) +
     guides(color = guide_legend(override.aes = list(fill = NA))) +
     theme_bw() +
     ylim(-0.3, 0.3) +
     xlim(-0.3, 0.25)
   

   
   
   
   
   
   

   
   

   
   fhf <- vegdist(t(asv_nmds))
   pcoa <- cmdscale (fhf, eig = TRUE)
   ordiplot (pcoa, display = 'sites', type = 'text')
   
   
   
   pcoa$points %>%
     as_tibble(rownames = 'asv_names') %>%
     left_join(clean_bacstrats, by = c('asv_names' = 'asv_id')) %>%
     filter(!is.na(bacstrat)) %>%
     mutate(bacstrat = case_when(bacstrat == "c(\"early_pathogen\", \"continuous_pathogen\", \"late_pathogen\")" ~ "continuous_pathogen",
                                 bacstrat == 'c("crasher_t3", "crasher_t7")' ~ "continuous_crasher",
                                 bacstrat == 'c("early_pathogen", "early_opportunist")' ~ "early_pathogen",
                                 bacstrat == 'c("late_pathogen", "late_opportunist")' ~ "late_pathogen",
                                 bacstrat == "c(\"probiotic_t7_strict\", \"late_opportunist\")" ~ "late_probiotic",
                                 bacstrat == "crasher_t7" ~ "late_crasher",
                                 bacstrat == "probiotic_t7_strict" ~ "late_probiotic",
                                 bacstrat == 'c("early_pathogen", "late_opportunist")' ~ "early_pathogen",
                                 TRUE ~ bacstrat)) %>%
     mutate(shape_determination = case_when(
       str_detect(bacstrat, "path")  ~ "pathogen",
       str_detect(bacstrat, "opp")  ~ "opportunist",
       str_detect(bacstrat, "crash")  ~ "crasher",
       str_detect(bacstrat, "pro")  ~ "probiotic"
     )) %>%
     ggplot(aes(x = V1, y = V2)) +
     stat_ellipse(geom = "polygon", alpha = 0.1, level = 0.95, aes(col = bacstrat, fill = bacstrat)) +
     # geom_point(data = as_tibble(scores(test_nmds)$species, rownames = aggregation_level) %>% 
     #              mutate(exp_conditions = str_split(none, "_")) %>%
     #              rowwise() %>%
     #              mutate(time = exp_conditions[1], exposure = exp_conditions[2], resist = exp_conditions[3]),
     #            size = 1, aes(color = time, shape = exposure)) +
     geom_point(aes(col = bacstrat, shape = shape_determination), size = 4) + #alpha = shape_determination
     geom_text(aes(label = parse_number(asv_names)), size = 3) +
     scale_shape_manual(values = c(
       "opportunist" = 16,
       "pathogen" = 17,
       "crasher" = 15,
       "probiotic" = 8,
       "H" = 6,
       "D" = 4,
       "REP1" = 3
     )) +
     scale_fill_manual(values = c(
       "continuous_crasher" = "sandybrown",
       "continuous_pathogen" = "firebrick4",
       "late_crasher" = "darkorange2",
       "late_probiotic" = "aquamarine",
       "early_pathogen" = "hotpink",
       "late_pathogen" = "firebrick1",
       "early_opportunist" = "deepskyblue3",
       "late_opportunist" = "royalblue2"
     ), guide = "none") +
     scale_color_manual(values = c(
       "continuous_crasher" = "sandybrown",
       "continuous_pathogen" = "firebrick4",
       "late_crasher" = "darkorange2",
       "late_probiotic" = "aquamarine",
       "early_pathogen" = "hotpink",
       "late_pathogen" = "firebrick1",
       "early_opportunist" = "deepskyblue3",
       "late_opportunist" = "royalblue2",
       "T0" = "gray80",
       "T3" = "gray40",
       "T7" = "gray10"
     )) +
     guides(color = guide_legend(override.aes = list(fill = NA))) +
     theme_bw() +
     ylim(-0.3, 0.3) +
     xlim(-0.3, 0.25)
   