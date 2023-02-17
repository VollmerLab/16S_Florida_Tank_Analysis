cpm(gene_normalize_factors, log = TRUE, prior.count = 2) %>%
  rowMeans %>%
  quantile(0.05)


cpm(gene_normalize_factors, log = TRUE, prior.count = 2) %>%
  rowMeans %>%
  tibble(x = .) %>%
  ggplot(aes(x = x)) +
  geom_histogram(fill = "grey", bins = 100) +
  theme_classic() +
  labs(y = "Density", x = "Filtered read counts (logCPM)",
       title = "Distribution of normalized, filtered read counts")

plot_pcoa <- function(cpm_counts){
  filtered_pcoa <- t(cpm_counts) %>%
    vegdist(method = 'euclidean') %>%
    divide_by(1000) %>%
    pcoa()
  
  percent_variance <- filtered_pcoa$values$Eigenvalues / sum(filtered_pcoa$values$Eigenvalue)
  
  filtered_pcoa$vectors %>%
    as_tibble(rownames = 'sequence_id') %>%
    dplyr::select(sequence_id, Axis.1, Axis.2) %>%
    inner_join(sample_metadata,
               by = 'sequence_id') %>%
    ggplot(aes(x = Axis.1, y = Axis.2, colour = treat_outcome, 
               shape = timepoint, group = fragment_id)) +
    geom_point() +
    geom_path() +
    labs(x = str_c('PCoA 1 (', scales::percent(percent_variance[1]), ')'),
         y = str_c('PCoA 2 (', scales::percent(percent_variance[2]), ')')) +
    theme_classic()
  
}

cpm(gene_normalize_factors, log = TRUE, prior.count = 2) %>% plot_pcoa()
