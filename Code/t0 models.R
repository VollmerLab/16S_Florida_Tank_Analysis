normalized_asv_counts

t0_only_data <- normalized_asv_counts %>% 
  filter(time == "T0") %>%
  group_by(asv_id)


##By ASV

# no ASVs correlated w continuous resistance
t0_only_data %>%
  group_by(asv_id) %>%
  summarize(corr_val = broom::tidy(cor.test(log2_cpm, resistance))) %>%
  unnest(corr_val) %>%
  filter(!is.na(estimate)) %>%
  mutate(p.value = p.adjust(p.value, method = "fdr")) %>%
  arrange(p.value) %>%
  ggplot() +
  geom_point(aes(x = fct_reorder(asv_id, p.value), y = p.value)) +
  geom_hline(yintercept = 0.05, col = "deepskyblue2") +
  coord_flip() +
  labs(title = "Genus, continuous")

# no ASVs correlated w discrete resistance
t0_only_data %>%
  group_by(asv_id) %>%
  reframe(broom::glance(lm(log2_cpm~susceptability))) 

t0_only_data %>%
  group_by(asv_id) %>%
  reframe(aov_val = broom::tidy(aov(log2_cpm~resistance))) %>%
  unnest(aov_val) %>%
  filter(term == "resistance") %>%
  mutate(p.value = p.adjust(p.value, method = "fdr")) %>%
  arrange(p.value) %>%
  ggplot() +
  geom_point(aes(x = fct_reorder(asv_id, p.value), y = p.value)) +
  geom_hline(yintercept = 0.05, col = "deepskyblue2") +
  coord_flip() +
  labs(title = "Genus, discrete")

##By Genus

# no genera correlated w continuous resistance
# Sulfurovum has p val of 0.0637
t0_only_data %>%
  group_by(Genus) %>%
  filter(!is.na(Genus)) %>%
  reframe(corr_val = broom::tidy(cor.test(log2_cpm, resistance))) %>%
  unnest(corr_val) %>%
  filter(!is.na(estimate)) %>%
  mutate(p.value = p.adjust(p.value, method = "fdr")) %>%
  arrange(p.value) %>%
  ggplot() +
  geom_point(aes(x = fct_reorder(Genus, p.value), y = p.value)) +
  geom_hline(yintercept = 0.05, col = "deepskyblue2") +
  coord_flip() +
  labs(title = "Genus, continuous")

t0_only_data %>%
  filter(Genus == "Sulfurovum") %>%
  ggplot() +
  geom_point(aes(x = log2_cpm, y = resistance)) +
  labs(title = "Sulfurovum (p = 0.0637)")

# no genera correlated w discrete resistance
# Sulfurovum has p val of 0.0725
t0_only_data %>%
  filter(!is.na(Genus)) %>%
  group_by(Genus) %>%
  reframe(aov_val = broom::tidy(aov(log2_cpm~resistance))) %>%
  unnest(aov_val) %>%
  filter(term == "resistance") %>%
  mutate(p.value = p.adjust(p.value, method = "fdr")) %>%
  arrange(p.value) %>%
  ggplot() +
  geom_point(aes(x = fct_reorder(Genus, p.value), y = p.value)) +
  geom_hline(yintercept = 0.05, col = "deepskyblue2") +
  coord_flip() +
  labs(title = "Genus, discrete")

#By Family

# continuous
# Roseobacteraceae 0.0119
# Sulfurovaceae 0.0176
# Alteromonadaceae 0.0644

t0_only_data %>%
  group_by(Family) %>%
  filter(!is.na(Family)) %>%
  reframe(corr_val = broom::tidy(cor.test(log2_cpm, resistance))) %>%
  unnest(corr_val) %>%
  filter(!is.na(estimate)) %>%
  mutate(p.value = p.adjust(p.value, method = "fdr")) %>%
  arrange(p.value) %>%
  ggplot() +
  geom_point(aes(x = fct_reorder(Family, p.value), y = p.value)) +
  geom_hline(yintercept = 0.05, col = "deepskyblue2") +
  coord_flip() +
  labs(title = "Family, continuous")

# discrete
# Roseobacteraceae 0.0132
# Sulfurovaceae 0.0194
# Alteromonadaceae 0.0711

t0_only_data %>%
  filter(!is.na(Family)) %>%
  group_by(Family) %>%
  reframe(aov_val = broom::tidy(aov(log2_cpm~resistance))) %>%
  unnest(aov_val) %>%
  filter(term == "resistance") %>%
  mutate(p.value = p.adjust(p.value, method = "fdr")) %>%
  arrange(p.value) %>%
  ggplot() +
  geom_point(aes(x = fct_reorder(Family, p.value), y = p.value)) +
  geom_hline(yintercept = 0.05, col = "deepskyblue2") +
  coord_flip() +
  labs(title = "Family, discrete")

# sig correlations? really? ....

t0_only_data %>%
  filter(Family == "Roseobacteraceae") %>%
  ggplot() +
  geom_point(aes(x = log2_cpm, y = resistance)) +
  labs(title = "Roseobacteraceae (p = 0.0119)")

t0_only_data %>%
  filter(Family == "Sulfurovaceae") %>%
  ggplot() +
  geom_point(aes(x = log2_cpm, y = resistance)) +
  labs(title = "Sulfurovaceae (p = 0.0176)")

t0_only_data %>%
  filter(Family == "Alteromonadaceae") %>%
  ggplot() +
  geom_point(aes(x = log2_cpm, y = resistance)) +
  labs(title = "Alteromonadaceae (p = 0.0644)")





t0_only_data %>%
  group_by(Family) %>%
  reframe(aov_val = list(lm(log2_cpm~susceptability))) %>%
  rowwise() %>%
  mutate(emmeansvals = as_tibble(emmeans(aov_val, ~susceptability) %>%
                                   contrast('pairwise'))) %>%
  unnest(emmeansvals) %>%
  ungroup() %>%
  mutate(p.value = p.adjust(p.value, method = "fdr")) %>%
  arrange(p.value)

#pairwise emmeans: nothing sig


testing22 <- t0_only_data %>%
  group_by(Family) %>%
  reframe(aov_val = list(lm(log2_cpm~susceptability))) %>%
  rowwise() %>%
  mutate(emmeansvals = list(emmeans(aov_val, ~susceptability, type = 'response') %>%
                            multcomp::cld(Letters = LETTERS) %>% 
                            as_tibble %>%
                            mutate(.group = str_trim(.group))))

testing22 %>%
  tidyr::unnest(cols = c(emmeansvals)) %>%
  group_by(Family) %>%
  filter(!c(emmean[susceptability == "R"] < 5.30 & emmean[susceptability == "S"] < 5.30)) %>%
  ungroup() %>%
  ggplot(aes(x = Family, y = emmean, ymin = lower.CL, ymax = upper.CL, colour = susceptability)) +
  geom_pointrange(position = position_dodge(0.5)) +
  geom_text(aes(y = upper.CL, label = .group), position = position_dodge(0.5),
            vjust = -1) +
  coord_flip()


testing22 %>%
  tidyr::unnest(cols = c(emmeansvals)) %>%
  filter(Family %in% c(Family[.group == "B"])) %>%
  ggplot(aes(x = Family, y = emmean, ymin = lower.CL, ymax = upper.CL, colour = susceptability)) +
  geom_pointrange(position = position_dodge(0.5)) +
  geom_text(aes(y = upper.CL, label = .group), position = position_dodge(0.5),
            vjust = -1) +
  coord_flip()


boop <- testing22 %>%
  mutate(lm_table = list(broom::tidy(aov_val))) %>%
  unnest(cols = c(lm_table))

boop %>%
  filter(term == "susceptabilityS") %>%
  mutate(p.value = p.adjust(p.value, method = "fdr")) %>%
  filter(p.value < 0.05)

#after FDR correcting, no sig Families by discrete susceptibility !

t0_only_data %>%
  group_by(Family) %>%
  reframe(aov_val = list(lm(log2_cpm~susceptability))) %>%
  rowwise() %>%
  mutate(emmeansvals = list(emmeans(aov_val, ~susceptability, type = 'response') %>%
                              multcomp::cld(Letters = LETTERS) %>% 
                              as_tibble %>%
                              mutate(.group = str_trim(.group)))) %>%
  mutate(lm_table = list(broom::tidy(aov_val))) %>%
  unnest(cols = c(lm_table)) %>%
  filter(term == "susceptabilityS") %>%
  mutate(p.value = p.adjust(p.value, method = "fdr")) %>%
  filter(p.value < 0.05)


##
testing %>%
  unnest(emmeansvals) %>%
  mutate(colors = ifelse(Family %in% c("Roseobacteraceae", "Sulfurovaceae", "Alteromonadaceae"), "red", "black")) %>%
  ggplot() +
  geom_point(aes(x = Family, y = emmean, col = colors)) +
  geom_errorbar(aes(x = Family, ymin = lower.CL, ymax = upper.CL))

t0_only_data %>%
  filter(Family %in% c("Roseobacteraceae", "Sulfurovaceae", "Alteromonadaceae")) %>%
  group_by(asv_id, susceptability) %>%
  ggplot() +
  geom_point(aes(x = asv_id, y = log2_cpm, col = susceptability), alpha = 0.7) +
  coord_flip() +
  facet_grid(rows = vars(Family), scales = "free", space = "free") +
  theme_bw()
  

t0_only_data %>%
  filter(Family %in% c("Roseobacteraceae", "Sulfurovaceae", "Alteromonadaceae")) %>%
  group_by(asv_id, susceptability) %>%
  ggplot() +
  geom_hline(yintercept = 5.29) +
  geom_boxplot(aes(x = asv_id, y = log2_cpm, fill = susceptability, col = susceptability)) +
  coord_flip() +
  facet_grid(rows = vars(Family), scales = "free", space = "free") +
  scale_fill_manual(values = c("R" = "#3DD8EA", "S" = "#F75D5D")) +
  scale_color_manual(values = c("R" = "#048291", "S" = "#A70000")) +
  theme_bw()

t0_only_data %>%
  filter(Family %in% c("Roseobacteraceae", "Sulfurovaceae", "Alteromonadaceae")) %>%
  group_by(Family, susceptability) %>%
  ggplot() +
  geom_hline(yintercept = 5.29) +
  geom_boxplot(aes(x = Family, y = log2_cpm, fill = susceptability, col = susceptability)) +
  coord_flip() +
  #facet_grid(rows = vars(Family), scales = "free", space = "free") +
  scale_fill_manual(values = c("R" = "#3DD8EA", "S" = "#F75D5D")) +
  scale_color_manual(values = c("R" = "#048291", "S" = "#A70000")) +
  theme_bw()
