library(ggvenn)

#setting up data frame including observed richness and # of reads
richness_data <- left_join(select(alpha_table, sample_id, observed),
          otu_table(microbiome_data) %>%
            rowSums() %>%
            enframe(name = 'sample_id',
                    'n_reads'),
          by = 'sample_id') %>%
  left_join(sample_data(microbiome_data) %>%
              as_tibble(rownames = 'sample_id'),
            by = 'sample_id')

#graph relationship between # of reads and observed species richness
richness_data %>%
  filter(n_reads > 10000) %>%
  ggplot(aes(x = n_reads, y = observed, colour = exposure, shape = time, linetype = time)) +
  geom_vline(xintercept = 10000) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE)

#graph boxplot of richness by exposure and time
richness_data %>%
  ggplot(aes(x = interaction(exposure, time), y = observed)) +
  geom_boxplot() +
  geom_jitter()

#graph boxplot of # of reads by exposure and time
richness_data %>%
  ggplot(aes(x = interaction(exposure, time), y = n_reads)) +
  geom_boxplot() +
  geom_jitter()

#T3 has lowest richness but highest # of reads - strange


#convert phyloseq data to tibble (metadata + abundance + taxonomy)
the_samples <- psmelt(microbiome_data) %>%
  as_tibble() %>%
  filter(Abundance > 0) %>%
  left_join(otu_table(microbiome_data) %>% #add column for # of reads
              rowSums() %>%
              enframe(name = 'sample_id',
                      'n_reads'),
            by = c('Sample' = 'sample_id'))

#filter for # of reads and remove acerv fragments that got homogenized
the_samples %>%
  filter(n_reads > 10000) %>%
  filter(tank != 'homogenate_fragment') %>%
  group_by(exposure, time, tank) %>%
  #summarize data for each exposure and timepoint to a single tank value and plot
  summarise(richness = n_distinct(OTU),
            n_reads = sum(Abundance),
            n_samples = n_distinct(Sample)) %>%
  ggplot(aes(x = interaction(exposure, time), y = richness)) +
  geom_boxplot() +
  geom_jitter()

#combine all T0 data to one metric and plot richness per tank
the_samples %>%
  # filter(n_reads > 10000) %>%
  filter(tank != 'homogenate_fragment') %>%
  mutate(exposure = if_else(time == 'T0', 'T0', exposure),
         tank = if_else(time == 'T0', 'T0', tank)) %>%
  group_by(exposure, tank, time) %>%
  summarise(richness = n_distinct(OTU),
            n_reads = sum(Abundance),
            n_samples = n_distinct(Sample)) %>%
  ggplot(aes(x = interaction(exposure, time), y = richness)) +
  geom_hline(yintercept = n_distinct(the_samples$OTU)) +
  geom_boxplot() +
  geom_jitter()
  
### Things in healthy/disease/field pool 

#abundances of each species in T0 pools(D or H) and field
the_samples %>%
  # filter(n_reads > 10000) %>%
  filter(tank != 'homogenate_fragment') %>%
  filter(time == 'T0') %>%
  group_by(exposure, OTU, Species, Genus, Family, Order) %>%
  summarise(abundance = sum(Abundance)) %>%
  pivot_wider(names_from = exposure, values_from = abundance)

#4 way venn diagram with ~10% prevalence filter
the_samples %>%
  # filter(n_reads > 10000) %>%
  filter(tank != 'homogenate_fragment') %>% 
  filter(exposure != 'Field') %>%
  group_by(OTU) %>%
  filter(n_distinct(Sample) > 21) %>%
  ungroup %>%
  mutate(time = if_else(time == 'T0', exposure, time)) %>%
  count(time, OTU) %>%
  pivot_wider(names_from = time, values_from = n, values_fill = 0L) %>%
  mutate(across(-OTU, ~. > 0)) %>%
  ggvenn(c('D', 'H', 'T3', 'T7'))

#same set up as venn but added in taxonomy table and saved to variable
target_upset_otus <- the_samples %>%
  # filter(n_reads > 10000) %>%
  filter(tank != 'homogenate_fragment') %>% 
  # filter(exposure != 'Field') %>%
  group_by(OTU) %>%
  filter(n_distinct(Sample) > 21) %>%
  ungroup %>%
  mutate(time = if_else(time == 'T0', exposure, time)) %>%
  count(time, OTU) %>%
  pivot_wider(names_from = time, values_from = n, values_fill = 0L) %>%
  mutate(across(-OTU, ~. > 0)) %>%
  left_join(taxonomy_tibble, by = c('OTU' = 'asv_names')) %>%
  select(-c(Kingdom, Phylum)) 


target_upset_otus %>%
  filter(T3 & T7)

#make complex upset plot of likely suspects for primary pathogen and opportunistic pathogen
upset(filter(target_upset_otus,
             # T3 & T7,
             (D & T7 & T3)), 
      c('D', 'H', 'T3', 'T7'), 
      annotations = list(
        # 2nd method - using ggplot
        'Order'=(
          ggplot(mapping=aes(fill=Order)) 
          + geom_bar(stat = 'count', position = 'fill') 
          + scale_y_continuous(labels = scales::percent_format())
        ) +
          ylab('Order') +
          theme(legend.position = 'top')
      ),
      name='asv_names', width_ratio=0.1, min_size = 1)


#likely suspects
target_otus <- filter(target_upset_otus,
                      # T3 & T7,
                      (D & T7 & T3))


target_otus %>%
  count(Order, Family) %>%
  arrange(-n)

#set up raw microbiome data
raw_target_data <- cpm(otu_tmm, log = TRUE, prior.count = 2) %>%
  t %>%
  as_tibble(rownames = "sample_id") %>%
  full_join(metadata, by = "sample_id") %>%
  pivot_longer(cols = -any_of(colnames(metadata)), 
               names_to = "asv_names", values_to = "value") %>%
  mutate(across(c(exposure, final_disease_state), factor)) %>%
  filter(time %in% c('T3', 'T7')) %>%
  mutate(fragment_id = str_c(exposure, tank, genotype, final_disease_state))
  
#only keep ASVs in list of likely suspects and plot
raw_target_data %>%
  inner_join(target_otus,
            by = c('asv_names' = 'OTU')) %>%
  filter(time == 'T7') %>%
  group_by(sample_id, final_disease_state, Order, Genus) %>%
  filter(!is.na(Genus)) %>%
  summarise(value = sum(value)) %>%
  ungroup %>%
  #dot plot
  ggplot(aes(y = Genus, x = sample_id, size = value, colour = Order)) + 
  geom_point() +
  
  # ggplot(aes(x = sample_id, y = value, fill = Order)) +
  # geom_col() +
  facet_wrap(~final_disease_state, scales = 'free_x') +
  theme(axis.text.x = element_blank())

#make repeated measures model for all asvs
all_models <- raw_target_data %>%
  nest_by(asv_names) %>%
  ungroup %>%
  slice(-971) %>% 
  rowwise(asv_names) %>%
  summarise(model = list(aov_4(value ~ time * (exposure + final_disease_state) + 
                        (time | fragment_id),
                      data = data)))

#make repeated measures model for likely suspects
subset_models <- raw_target_data %>%
  inner_join(target_otus,
             by = c('asv_names' = 'OTU')) %>%
  nest_by(asv_names) %>%
  mutate(model = list(aov_4(value ~ time * (exposure + final_disease_state) + 
                                 (time | fragment_id),
                               data = data)))

#process data for performing complex upset - all ASVs
asv_comp_upset_all <- all_models %>%
  ungroup %>%
  # slice(1:10) %>%
  rowwise(asv_names) %>%
  reframe(as_tibble(model$anova_table, rownames = 'param') %>%
            select(param, `Pr(>F)`)) %>%
  rename(p = `Pr(>F)`) %>%
  group_by(param) %>% #param = term or interaction
  mutate(p = p.adjust(p, 'fdr')) %>%
  ungroup %>%
  mutate(p = p < 0.05) %>% #give true/false
  pivot_wider(names_from = 'param', values_from = p, names_prefix = 'p_') %>% #pivot each term into its own column
  left_join(taxonomy_tibble, by = join_by(asv_names)) %>%
  select(-c(Kingdom, Phylum)) %>%
  filter(!if_all(starts_with('p_'), ~!.)) #remove ASVs where all p vals are not significant

#process data for performing complex upset - likely suspects
asv_comp_upset_subset <- subset_models %>%
  ungroup %>%
  # slice(1:10) %>%
  rowwise(asv_names) %>%
  reframe(as_tibble(model$anova_table, rownames = 'param'),
          d_v_h = emmeans(model, ~final_disease_state) %>%
            as_tibble %>%
            select(emmean) %>%
            pull(1) %>%
            diff) %>%
  rename(p = `Pr(>F)`) %>%
  
  group_by(param) %>%
  mutate(p = p.adjust(p, 'fdr')) %>%
  ungroup %>%
  mutate(p = p < 0.05) %>%
  select(asv_names, param, p, d_v_h) %>%
  pivot_wider(names_from = 'param', values_from = p, names_prefix = 'p_', values_fill = FALSE) %>%
  left_join(taxonomy_tibble, by = join_by(asv_names)) %>%
  select(-c(Kingdom, Phylum)) %>%
  select(-p_time) %>%
  filter(!if_all(starts_with('p_'), ~!.)) %>%
  filter(p_final_disease_state | `p_final_disease_state:time`) #must have significant disease state term

#complex upset for all ASVs
all_upset <- upset(asv_comp_upset_all, 
      colnames(select(asv_comp_upset_all, starts_with('p_'))), 
      
      annotations = list(
        # 2nd method - using ggplot
        'Order'=(
          ggplot(mapping=aes(fill=Order)) 
          + geom_bar(stat = 'count', position = 'fill') 
          + scale_y_continuous(labels = scales::percent_format())
        ) +
          ylab('Order') +
          theme(legend.position = 'top')
      ),
      
      name='asv_names', width_ratio=0.1, min_size = 1)

#complex upset for likely suspect ASVs that are more abundant in Diseased (not necessarily significantly)
subset_upset_moreDisease <- upset(filter(asv_comp_upset_subset, d_v_h < 0) %>%
                        select(-d_v_h),  #negative = more in disease final state, positive = more in healthy disease state
      colnames(select(asv_comp_upset_subset, starts_with('p_'))), 
      
      annotations = list(
        # 2nd method - using ggplot
        'Genus'=(
          ggplot(mapping=aes(fill=Genus)) 
          + geom_bar(stat = 'count', position = 'fill') 
          + scale_y_continuous(labels = scales::percent_format())
        ) +
          ylab('Order') +
          theme(legend.position = 'top')
      ),
      
      name='asv_names', width_ratio=0.1, min_size = 1)

#complex upset for likely suspect ASVs that are more abundant in Healthy (not necessarily significantly)
subset_upset_moreHealthy <- upset(filter(asv_comp_upset_subset, d_v_h >0) %>%
                                    select(-d_v_h),  #negative = more in disease final state, positive = more in healthy disease state
                                  colnames(select(asv_comp_upset_subset, starts_with('p_'))), 
                                  
                                  annotations = list(
                                    # 2nd method - using ggplot
                                    'Genus'=(
                                      ggplot(mapping=aes(fill=Genus)) 
                                      + geom_bar(stat = 'count', position = 'fill') 
                                      + scale_y_continuous(labels = scales::percent_format())
                                    ) +
                                      ylab('Order') +
                                      theme(legend.position = 'top')
                                  ),
                                  
                                  name='asv_names', width_ratio=0.1, min_size = 1)
#add titles
subset_upset_moreDisease &
  labs(title = 'More in Disease')
subset_upset_moreHealthy &
  labs(title = 'More in Healthy')


#what are the distinct genera that are more abundant in diseased in likely suspects list
asv_comp_upset_subset %>%
  filter(d_v_h < 0) %>% #more abundant in diseased
  select(Genus) %>%
  filter(!is.na(Genus)) %>%
  distinct

#dot plot that wasn't very helpful
raw_target_data %>%
  inner_join(target_otus,
             by = c('asv_names' = 'OTU')) %>%
  filter(time == 'T7') %>%
  group_by(sample_id, final_disease_state, Order, Genus) %>%
  filter(!is.na(Genus)) %>%
  summarise(value = sum(value)) %>%
  ungroup %>%
  
  inner_join(asv_comp_upset_subset %>%
               filter(d_v_h < 0) %>%
               select(Genus) %>%
               filter(!is.na(Genus)) %>%
               distinct,
             by = 'Genus') %>%
  
  ggplot(aes(y = Genus, x = sample_id, size = value, colour = Order)) +
  geom_point() +
  
  # ggplot(aes(x = sample_id, y = value, fill = Order)) +
  # geom_col() +
  facet_wrap(~final_disease_state, scales = 'free_x') +
  theme(axis.text.x = element_blank())


#genera that are more abundant in diseased in likely suspects list by final disease state, exposure, and time
raw_target_data %>% 
  left_join(tax_table(microbiome_data) %>%
              as.data.frame() %>%
              as_tibble(rownames = 'asv_names'),
            by = 'asv_names') %>%
  inner_join(asv_comp_upset_subset %>%
               filter(d_v_h < 0) %>%
               select(Genus) %>%
               filter(!is.na(Genus)) %>%
               distinct,
             by = 'Genus') %>%
  
  group_by(time, exposure, final_disease_state, Genus, sample_id) %>%
  summarise(value = sum(value),
            .groups = 'drop_last') %>%
  summarise(mean_cpm = mean(value),
            se_cpm = sd(value) / sqrt(n())) %>%
  
  ggplot(aes(x = time, y = mean_cpm, shape = exposure, colour = final_disease_state)) +
  geom_pointrange(aes(ymin = mean_cpm - se_cpm, ymax = mean_cpm + se_cpm),
                  position = position_dodge(0.5)) +
  facet_wrap(~Genus, scales = 'free_y')
  