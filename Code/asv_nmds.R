
tgfff <- bacterial_signature_asv %>% 
  group_by(asv_id) %>% 
  reframe(bacstrat = paste(list(signatures))) %>%
  mutate(bacstrat = ifelse(bacstrat == 
                             'c("late_pathogen", "late_opportunist")', "late_opportunist&pathogen", 
                           bacstrat)) %>%
  mutate(bacstrat = ifelse(bacstrat == 
                             'c("early_opportunist", "continuous_opportunist", "late_opportunist")', 
                           "continuous_opportunist", bacstrat)) %>%
  mutate(bacstrat = ifelse(bacstrat == 
                             'c("early_pathogen", "late_opportunist")', 
                           "early_path_late_opp", bacstrat)) %>%
  mutate(bacstrat = ifelse(bacstrat == 
                             'c("early_pathogen", "early_opportunist", "continuous_opportunist", "late_opportunist")', 
                           "early_path_cont_opp", bacstrat))


testnmds <- full_data %>% 
  mutate(category = paste(time, exposure, susceptability, sep = "_"), .after = asv_id) %>%
  select(asv_id, category, log2_cpm) %>%
  group_by(asv_id, category) %>%
  summarize(ave_cpm = mean(log2_cpm)) %>%
  mutate(ave_cpm = str_trim(as.numeric(ave_cpm))) %>%
  pivot_wider(names_from = category, values_from = ave_cpm) %>%
  left_join(tgfff, by = join_by("asv_id")) %>%
  filter(!is.na(bacstrat)) %>%
  select(-bacstrat) %>%
  as.matrix()

mat_nmds <- as.matrix(testnmds[, -1])
rownames(mat_nmds) <- testnmds[,1]

test_num_mat <- matrix(as.numeric(mat_nmds),    # Convert to numeric matrix
       ncol = ncol(mat_nmds))

rownames(test_num_mat) <- rownames(mat_nmds)
colnames(test_num_mat) <- colnames(mat_nmds)

dist.f <- vegdist(test_num_mat, method = "bray")

test_nmds <- metaMDS(dist.f, distance = 'bray', k = 2, trymax = 100, autotransform = FALSE, verbose = TRUE)
sppscores(test_nmds) <- test_num_mat



#make plot
scores(test_nmds)$sites %>%
  as_tibble(rownames = 'asv_names') %>%
  #left_join(metadata, by = 'asv_names') %>%
  left_join(tgfff, by = c('asv_names' = 'asv_id')) %>%
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
