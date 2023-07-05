

testnmds <- full_data %>% 
  mutate(category = paste(time, exposure, susceptability, sep = "_"), .after = asv_id) %>%
  select(asv_id, category, log2_cpm) %>%
  group_by(asv_id, category) %>%
  summarize(ave_cpm = mean(log2_cpm)) %>%
  mutate(ave_cpm = str_trim(as.numeric(ave_cpm))) %>%
  pivot_wider(names_from = category, values_from = ave_cpm) %>%
  as.matrix()

mat_nmds <- as.matrix(testnmds[, -1])
rownames(mat_nmds) <- testnmds[,1]

test_num_mat <- matrix(as.numeric(mat_nmds),    # Convert to numeric matrix
       ncol = ncol(mat_nmds))

rownames(test_num_mat) <- rownames(mat_nmds)
colnames(test_num_mat) <- colnames(mat_nmds)

dist.f <- vegdist(test_num_mat, method = "bray")

test_nmds <- metaMDS(dist.f, distance = 'bray', k = 2, trymax = 100, autotransform = FALSE, verbose = TRUE)
sppscores(test_nmds) <- test5

tgfff <- bacterial_signature_asv %>% 
  group_by(asv_id) %>% 
  reframe(bacstrat = paste(list(signatures))) %>%
  mutate(bacstrat = ifelse(bacstrat == 
        'c("late_opportunist", "late_pathogen")', "late_opportunist&pathogen", 
        bacstrat)) %>%
  mutate(bacstrat = ifelse(bacstrat == 
        'c("early_opportunist", "continuous_opportunist", "late_opportunist")', 
        "continuous_opportunist", bacstrat))
 
#make plot
  scores(test_nmds)$sites %>%
  as_tibble(rownames = 'asv_names') %>%
  #left_join(metadata, by = 'asv_names') %>%
  left_join(tgfff, by = c('asv_names' = 'asv_id')) %>%
  mutate(bacstrat = ifelse(is.na(bacstrat), "disease_unrelated", bacstrat)) %>%
  mutate(shape_determination = ifelse(str_detect(bacstrat, "disease"), 
        "disease_unrelated", ifelse(str_detect(bacstrat, "pathogen"), "pathogen", 
        "opportunist"))) %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_point(data = as_tibble(scores(test_nmds)$species, rownames = aggregation_level) %>% mutate(timepoint = str_extract_numbers(none)),
             size = 1, shape = 1) +
  geom_point(aes(alpha = shape_determination, col = bacstrat, shape = shape_determination)) +
  scale_alpha_manual(values = c(0.3, 1,1)) +
  scale_color_manual(values = c(
    "early_pathogen" = "pink2",
    "late_pathogen" = "red",
    "early_opportunist" = "lightskyblue1",
    "late_opportunist" = "mediumblue",
    "late_opportunist&pathogen" = "purple",
    "continuous_opportunist" = "deepskyblue",
    "disease_unrelated" = "gray"
  )) +
  theme_bw()

#plot species or site alone
  plot(test_nmds, "species")
  orditorp(test_nmds, "species")
