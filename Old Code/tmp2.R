microbiome_data <- read_rds("../intermediate_files/preprocess_microbiome.rds") %>%
  subset_samples(tank %in% c('HOMO'))

tmp <- taxon_abundances %>%
  summarise(broom::tidy(t.test(value ~ exposure, data = data)),
            .groups = 'drop') %>%
  mutate(p.adj = p.adjust(p.value, 'fdr')) %>%
  arrange(p.adj)
tmp %>%
  filter(p.adj < 0.05) %>%
  pull(taxon)
tmp$taxon


microbiome_data <- read_rds("../intermediate_files/preprocess_microbiome.rds") %>%
  subset_samples(tank %in% c('homogenate_fragment'))

tmp2 <- taxon_abundances %>%
  summarise(broom::tidy(t.test(value ~ exposure, data = data)),
            .groups = 'drop') %>%
  mutate(p.adj = p.adjust(p.value, 'fdr')) %>%
  arrange(p.adj)
tmp2 %>%
  filter(p.adj < 0.05)
