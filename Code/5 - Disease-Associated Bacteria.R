#data to analyze differential abundances of disease-associated bacteria

setwd("~/Desktop/Screenshots/Career/Vollmer Lab/GitHub/16S_Florida_Tank_Analysis/Code")

#### Packages ####
library(ggvenn)
library(multcomp)
library(phyloseq)
library(microbiome)
library(vegan)
library(lme4)
library(afex)
library(emmeans)
library(car)
library(edgeR)
library(metagMisc)
library(ape)
library(ggdist)
library(gghalves)
library(patchwork)
library(magrittr)
library(ComplexUpset)
library(strex)
library(forcats)
library(wesanderson)
library(msa)
library(Biostrings)
library(phangorn)
library(tidyverse)

#### Read in Data ####

#convert phyloseq data to tibble (metadata + abundance + taxonomy)
aggregation_level <- 'none' #or none

microbiome_data <- read_rds("../intermediate_files/preprocess_microbiome.rds")
metadata <- sample_data(microbiome_data) %>%
  as_tibble(rownames = 'sample_id') %>%
  select(-retain_sample)

if(aggregation_level != 'none'){
  microbiome_data <- aggregate_taxa(microbiome_data, aggregation_level)
  taxa_names(microbiome_data) <- str_replace_all(taxa_names(microbiome_data), ' |-', '_')
} else {
  #sequences <- taxa_names(microbiome_data)
  taxa_names(microbiome_data) <- str_c('ASV', 1:length(taxa_names(microbiome_data)), sep = '_')
  #names(sequences) <- taxa_names(microbiome_data)
}

the_samples <- psmelt(microbiome_data) %>%
  as_tibble() %>%
  filter(Abundance > 0) %>%
  left_join(otu_table(microbiome_data) %>% #add column for # of reads
              rowSums() %>%
              enframe(name = 'sample_id',
                      'n_reads'),
            by = c('Sample' = 'sample_id'))


alpha_table <- microbiome::alpha(microbiome_data, index = "all") %>%
  as_tibble(rownames = 'sample_id') %>%
  inner_join(metadata, by = 'sample_id') %>%
  mutate(fragment_id = str_c(str_replace_na(exposure, 'NA'), tank, genotype, sep = '_'))

#otu_tmm
otu_tmm <- microbiome_data %>%
  phyloseq_filter_prevalence(prev.trh = 0.1) %>%
  otu_table() %>% 
  t %>% #NOTE: *genus and family do not need the t but ASVs need the t*
  as.data.frame %>%
  as.matrix %>% 
  DGEList(remove.zeros = TRUE) %>%
  edgeR::calcNormFactors(method = 'TMMwsp')

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

#### Analysis of Species Richness ####

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
  
#### Differential Abundance Analysis ####

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

taxonomy_tibble <- tax_table(microbiome_data) %>% 
  as.data.frame %>%
  as_tibble(rownames = "asv_names")

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

#### All Bacteria ####

#make repeated measures model for all asvs
all_models <- raw_target_data %>%
  nest_by(asv_names) %>%
  ungroup %>%
  slice(-971) %>% 
  rowwise(asv_names) %>%
  summarise(model = list(aov_4(value ~ time * (exposure + final_disease_state) + 
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


#complex upset for all ASVs by interaction type
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



#### Likely Suspects ####

#likely suspects
target_otus <- filter(target_upset_otus,
                      (D & T7 & T3))

#make complex upset plot of likely suspects for primary pathogen and opportunistic pathogen
upset(filter(target_upset_otus,
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

#make repeated measures model for likely suspects
subset_models <- raw_target_data %>%
  inner_join(target_otus,
             by = c('asv_names' = 'OTU')) %>%
  nest_by(asv_names) %>%
  mutate(model = list(aov_4(value ~ time * (exposure + final_disease_state) + 
                              (time | fragment_id),
                            data = data)))

#process data for performing complex upset - likely suspects
asv_comp_upset_subset <- subset_models %>%
  ungroup %>%
  # slice(1:10) %>%
  rowwise(asv_names) %>%
  reframe(as_tibble(model$anova_table, rownames = 'param'),
          d_v_h =  emmeans(model, ~final_disease_state) %>%
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

#### emmeans analysis ####

##more in disease - very likely suspects

vl_suspects_list <- filter(asv_comp_upset_subset, d_v_h < 0)$asv_names # which ASVs are more in Disease

v_likely_suspects <- raw_target_data %>%
  filter(asv_names %in% vl_suspects_list) #only select asvs more in disease

vl_suspects_aov <- lmer(value ~ asv_names*final_disease_state*time + (1 | genotype) + (1 | tank), 
                         data = v_likely_suspects)
anova(vl_suspects_aov)


emmeans(vl_suspects_aov, ~final_disease_state*time | asv_names, type = 'response') %>%
  cld(Letters = LETTERS, adjust = 'fdr') %>%
  as_tibble() %>%
  left_join(taxonomy_tibble, by = join_by(asv_names)) %>%
  mutate(graph_color = paste(time, final_disease_state, sep = "_")) %>%
  mutate(.group = str_trim(.group)) %>%
  #rename(emmean = response) %>%
  ggplot(aes(x = paste(Family, " ", Genus, "(", asv_names, ")", sep = ""), y = emmean, ymin = emmean - SE, ymax = emmean + SE,
             colour = graph_color, pch = time)) +
  geom_pointrange(position = position_dodge(0.5)) +
  geom_text(aes(y = (emmean + SE), label = .group),
            position = position_dodge(0.5), vjust = -1) +
  scale_color_manual(values = c("hotpink1", "deepskyblue", "firebrick1", "dodgerblue3")) +
  coord_flip() +
  labs(title = "Very Likely Suspects") +
  xlab("ASV")


#likely suspects

likely_suspects <- raw_target_data %>%
  filter(asv_names %in% asv_comp_upset_subset$asv_names)

likely_suspects_list <- unique(likely_suspects$asv_names)


#for_rds <- list(likely_suspects_list, vl_suspects_list)
#write_rds(for_rds, "likely_suspects_list.rds")

l_suspects_aov <- lmer(value ~ asv_names*final_disease_state*time + (1 | genotype) + (1 | tank), 
                        data = likely_suspects)
anova(l_suspects_aov)

emmeans(l_suspects_aov, ~final_disease_state*time | asv_names, type = 'response') %>%
  cld(Letters = LETTERS, adjust = 'fdr') %>%
  as_tibble() %>%
  left_join(taxonomy_tibble, by = join_by(asv_names)) %>%
  mutate(graph_color = paste(time, final_disease_state, sep = "_")) %>%
  mutate(.group = str_trim(.group)) %>%
  #rename(emmean = response) %>%
  ggplot(aes(x = paste(Family, " ", Genus, "(", asv_names, ")", sep = ""), y = emmean, ymin = emmean - SE, ymax = emmean + SE,
             colour = graph_color, pch = time)) +
  geom_pointrange(position = position_dodge(0.5)) +
  geom_text(aes(y = (emmean + SE), label = .group),
            position = position_dodge(0.5), vjust = -1) +
  scale_color_manual(values = c("hotpink1", "deepskyblue", "firebrick1", "dodgerblue3")) +
  coord_flip() +
  labs(title = "Likely Suspects") +
  xlab("ASV")

#LS not in VLS

only_likely_suspects <- likely_suspects %>%
  filter(!asv_names %in% vl_suspects_list)

only_l_suspects_aov <- lmer(value ~ asv_names*final_disease_state*time + (1 | genotype) + (1 | tank), 
                       data = only_likely_suspects)
anova(only_l_suspects_aov)

emmeans(only_l_suspects_aov, ~final_disease_state*time | asv_names, type = 'response') %>%
  cld(Letters = LETTERS, adjust = 'fdr') %>%
  as_tibble() %>%
  left_join(taxonomy_tibble, by = join_by(asv_names)) %>%
  mutate(graph_color = paste(time, final_disease_state, sep = "_")) %>%
  mutate(.group = str_trim(.group)) %>%
  #rename(emmean = response) %>%
  ggplot(aes(x = paste(Family, " ", Genus, "(", asv_names, ")", sep = ""), y = emmean, ymin = emmean - SE, ymax = emmean + SE,
             colour = graph_color, pch = time)) +
  geom_pointrange(position = position_dodge(0.5)) +
  geom_text(aes(y = (emmean + SE), label = .group),
            position = position_dodge(0.5), vjust = -1) +
  scale_color_manual(values = c("hotpink1", "deepskyblue", "firebrick1", "dodgerblue3")) +
  coord_flip() +
  labs(title = "Only Likely Suspects not in Very Likely Suspects List") +
  xlab("ASV")


# most sig

most_sig_diffs <- subset_models %>%
  ungroup %>%
  # slice(1:10) %>%
  rowwise(asv_names) %>%
  reframe(as_tibble(model$anova_table, rownames = 'param'),
          d_v_h = emmeans(model, ~final_disease_state) %>%
            as_tibble %>%
            select(emmean) %>%
            pull(1) %>%
            diff) %>%
  dplyr::rename(p = `Pr(>F)`) %>%
  group_by(param) %>%
  mutate(p = p.adjust(p, 'fdr')) %>%
  ungroup %>%
  filter(param %in% c("final_disease_state", "final_disease_state:time"))

#likely suspects

interaction_types <- most_sig_diffs %>%
  filter(p < 0.05) %>%
  nest_by(asv_names) %>%
  rowwise() %>%
  mutate(type = ifelse(nrow(data) == 2, "both", NA)) %>%
  unnest(cols = c(data)) %>%
  mutate(type = ifelse(is.na(type), param, type)) %>%
  mutate(d_v_h_sign = ifelse(d_v_h < 0, "More in Disease", "More in Healthy")) %>%
  select(asv_names, type, d_v_h, d_v_h_sign)

both_list <- interaction_types %>%
  filter(type == "both") %>%
  .$asv_names %>%
  unique()


## likely suspects - all together in one model
all_types <- raw_target_data %>%
  filter(asv_names %in% interaction_types$asv_names)

all_types_model <- lmer(value ~ asv_names*final_disease_state*time + (1 | genotype) + (1 | tank), 
     data = all_types)

emmeans(all_types_model, ~final_disease_state*time | asv_names, type = 'response') %>%
  cld(Letters = LETTERS) %>%
  as_tibble() %>%
  left_join(taxonomy_tibble, by = join_by(asv_names)) %>%
  left_join(interaction_types, by = join_by(asv_names), multiple = "all") %>%
  mutate(graph_color = paste(time, final_disease_state, sep = "_")) %>%
  mutate(.group = str_trim(.group)) %>%
  #rename(emmean = response) %>%
  ggplot(aes(x = fct_reorder(paste(Family, " ", Genus, " (", asv_names, ") - ", sep = ""), d_v_h, .desc = TRUE), y = emmean, ymin = emmean - SE, ymax = emmean + SE,
             colour = graph_color, pch = time)) +
  geom_pointrange(position = position_dodge(0.5)) +
  geom_text(aes(y = (emmean + SE), label = .group),
            position = position_dodge(0.5), vjust = -1) +
  scale_color_manual(values = c("hotpink1", "deepskyblue", "firebrick1", "dodgerblue3")) +
  coord_flip() +
  labs(title = "Likely Suspects") +
  xlab("ASV Name") +
  facet_wrap(~type, scales = "free", nrow = 3)


suspect_group <- raw_target_data %>%
  filter(asv_names %in% vl_suspects_list)
  
### FDS

fds_facet <- suspect_group %>%
  filter(asv_names %in% interaction_types$asv_names[interaction_types$type == "final_disease_state"])

fds_facet_aov <- lmer(value ~ asv_names*final_disease_state + (1 | genotype) + (1 | tank), 
                        data = fds_facet)

fds_facet_graph <- emmeans(fds_facet_aov, ~final_disease_state | asv_names, type = 'response') %>%
  cld(Letters = LETTERS) %>%
  as_tibble() %>%
  left_join(taxonomy_tibble, by = join_by(asv_names)) %>%
  left_join(interaction_types, by = join_by(asv_names), multiple = "all") %>%
  mutate(.group = str_trim(.group)) %>%
  #rename(emmean = response) %>%
  ggplot(aes(x = fct_reorder(paste(Family, " ", Genus, " (", asv_names, ")", sep = ""), d_v_h, .desc = TRUE), y = emmean, ymin = emmean - SE, ymax = emmean + SE,
             colour = final_disease_state)) +
  geom_pointrange(position = position_dodge(0.5)) +
  geom_text(aes(y = (emmean + SE), label = .group),
            position = position_dodge(0.5), vjust = -1) +
  scale_color_manual(values = c("firebrick1", "dodgerblue3")) +
  coord_flip() +
  labs(title = "Final Disease State") +
  xlab("ASV Name")

### FDS:time

fds_t_facet <- suspect_group %>%
  filter(asv_names %in% interaction_types$asv_names[interaction_types$type == "final_disease_state:time"])

fds_t_facet_aov <- lmer(value ~ asv_names*final_disease_state*time + (1 | genotype) + (1 | tank), 
                        data = fds_t_facet)

fds_t_facet_graph <- emmeans(fds_t_facet_aov, ~final_disease_state*time | asv_names, type = 'response') %>%
  cld(Letters = LETTERS) %>%
  as_tibble() %>%
  left_join(taxonomy_tibble, by = join_by(asv_names)) %>%
  left_join(interaction_types, by = join_by(asv_names), multiple = "all") %>%
  mutate(graph_color = paste(time, final_disease_state, sep = "_")) %>%
  mutate(.group = str_trim(.group)) %>%
  #rename(emmean = response) %>%
  ggplot(aes(x = fct_reorder(paste(Family, " ", Genus, " (", asv_names, ")", sep = ""), d_v_h, .desc = TRUE), y = emmean, ymin = emmean - SE, ymax = emmean + SE,
             colour = graph_color, pch = time)) +
  geom_pointrange(position = position_dodge(0.5)) +
  geom_text(aes(y = (emmean + SE), label = .group),
            position = position_dodge(0.5), vjust = -1) +
  scale_color_manual(values = c("hotpink1", "deepskyblue", "firebrick1", "dodgerblue3")) +
  coord_flip() +
  labs(title = "Final Disease State:Time") +
  xlab("ASV Name")

### both

both_facet <- suspect_group %>%
  filter(asv_names %in% interaction_types$asv_names[interaction_types$type == "both"])

both_facet_aov <- lmer(value ~ asv_names*final_disease_state*time + (1 | genotype) + (1 | tank), 
                        data = both_facet)

both_facet_graph <- emmeans(both_facet_aov, ~final_disease_state*time | asv_names, type = 'response') %>%
  cld(Letters = LETTERS) %>%
  as_tibble() %>%
  left_join(taxonomy_tibble, by = join_by(asv_names)) %>%
  left_join(interaction_types, by = join_by(asv_names), multiple = "all") %>%
  mutate(graph_color = paste(time, final_disease_state, sep = "_")) %>%
  mutate(.group = str_trim(.group)) %>%
  #rename(emmean = response) %>%
  ggplot(aes(x = fct_reorder(paste(Family, " ", Genus, " (", asv_names, ")", sep = ""), d_v_h, .desc = TRUE), y = emmean, ymin = emmean - SE, ymax = emmean + SE,
             colour = graph_color, pch = time)) +
  geom_pointrange(position = position_dodge(0.5)) +
  geom_text(aes(y = (emmean + SE), label = .group),
            position = position_dodge(0.5), vjust = -1) +
  scale_color_manual(values = c("hotpink1", "deepskyblue", "firebrick1", "dodgerblue3")) +
  coord_flip() +
  labs(title = "Both FDS and FDS:Time") +
  xlab("ASV Name")

#all together

fds_facet_graph / fds_t_facet_graph / both_facet_graph + plot_annotation(title = "Very Likely Suspects")



## very likely suspects - all together in one model
vls_all_types <- raw_target_data %>%
  filter(asv_names %in% interaction_types$asv_names) %>%
  filter(asv_names %in% vl_suspects_list)

vls_all_types_model <- lmer(value ~ asv_names*final_disease_state*time + (1 | genotype) + (1 | tank), 
                        data = vls_all_types)

emmeans(vls_all_types_model, ~final_disease_state*time | asv_names, type = 'response') %>%
  cld(Letters = LETTERS) %>%
  as_tibble() %>%
  left_join(taxonomy_tibble, by = join_by(asv_names)) %>%
  left_join(interaction_types, by = join_by(asv_names), multiple = "all") %>%
  mutate(graph_color = paste(time, final_disease_state, sep = "_")) %>%
  mutate(.group = str_trim(.group)) %>%
  #rename(emmean = response) %>%
  ggplot(aes(x = fct_reorder(paste(Family, " ", Genus, " (", asv_names, ") - ", sep = ""), d_v_h, .desc = TRUE), y = emmean, ymin = emmean - SE, ymax = emmean + SE,
             colour = graph_color, pch = time)) +
  geom_pointrange(position = position_dodge(0.5)) +
  geom_text(aes(y = (emmean + SE), label = .group),
            position = position_dodge(0.5), vjust = -1) +
  scale_color_manual(values = c("hotpink1", "deepskyblue", "firebrick1", "dodgerblue3")) +
  coord_flip() +
  labs(title = "Very Likely Suspects") +
  xlab("ASV Name") +
  facet_wrap(~type, scales = "free", nrow = 3)





#### previous attempt ####
##fds time interaction

more_disease_fds_t <- most_sig_diffs %>% filter(param %in% c("final_disease_state:time")) %>% arrange(d_v_h) %>% 
  filter(p < 0.05) %>% filter(d_v_h < 0)
  
fds_t_md <- raw_target_data %>%
  filter(asv_names %in% more_disease_fds_t$asv_names) %>%
  filter(!asv_names %in% both_list)

fds_t_md_aov <- lmer(value ~ asv_names*final_disease_state*time + (1 | genotype) + (1 | tank), 
                       data = fds_t_md)
anova(fds_t_md_aov)

emmeans(fds_t_md_aov, ~final_disease_state*time | asv_names, type = 'response') %>%
  cld(Letters = LETTERS, adjust = 'fdr') %>%
  as_tibble() %>%
  left_join(taxonomy_tibble, by = join_by(asv_names)) %>%
  mutate(graph_color = paste(time, final_disease_state, sep = "_")) %>%
  mutate(.group = str_trim(.group)) %>%
  left_join(more_disease_fds_t, by = join_by(asv_names)) %>%
  filter(d_v_h < -0.40) %>%
  #rename(emmean = response) %>%
  ggplot(aes(x = fct_reorder(paste(Family, " ", Genus, " (", asv_names, ") - ", sep = ""), d_v_h, .desc = TRUE), y = emmean, ymin = emmean - SE, ymax = emmean + SE,
             colour = graph_color, pch = time)) +
  geom_pointrange(position = position_dodge(0.5)) +
  geom_text(aes(y = (emmean + SE), label = .group),
            position = position_dodge(0.5), vjust = -1) +
  scale_color_manual(values = c("hotpink1", "deepskyblue", "firebrick1", "dodgerblue3")) +
  #scale_color_manual(values = wes_palette("Zissou1", 4, type = "discrete")) +
  coord_flip() +
  labs(title = "Very Likely Suspects - Interaction Only") +
  xlab("ASV Name")

## final disease state only

more_disease_fds <- most_sig_diffs %>% filter(param %in% c("final_disease_state")) %>% arrange(d_v_h) %>% 
  filter(p < 0.05) %>% filter(d_v_h < 0)

fds_md <- raw_target_data %>%
  filter(asv_names %in% more_disease_fds$asv_names) %>%
  filter(!asv_names %in% both_list)

fds_md_aov <- lmer(value ~ asv_names*final_disease_state + (1 | genotype) + (1 | tank), 
                     data = fds_md)
anova(fds_md_aov)

emmeans(fds_md_aov, ~final_disease_state | asv_names, type = 'response') %>%
  cld(Letters = LETTERS, adjust = 'fdr') %>%
  as_tibble() %>%
  left_join(taxonomy_tibble, by = join_by(asv_names)) %>%
  mutate(.group = str_trim(.group)) %>%
  left_join(more_disease_fds, by = join_by(asv_names)) %>%
  #rename(emmean = response) %>%
  ggplot(aes(x = fct_reorder(paste(Family, " ", Genus, " (", asv_names, ")", sep = ""), d_v_h, .desc = TRUE), y = emmean, ymin = emmean - SE, ymax = emmean + SE,
             colour = final_disease_state)) +
  geom_pointrange(position = position_dodge(0.5)) +
  #geom_text(aes(y = (emmean + SE), label = .group),
  #          position = position_dodge(0.5), vjust = -1) +
  scale_color_manual(values = c("firebrick1", "dodgerblue3")) +
  #scale_color_manual(values = wes_palette("Zissou1", 2, type = "discrete")) +
  coord_flip() +
  labs(title = "Very Likely Suspects - Final Disease Only") +
  xlab("ASV Name")


## both final disease state and interaction

more_both_fds <- most_sig_diffs %>% filter(p < 0.05) %>% nest_by(asv_names) %>% 
  rowwise() %>% filter(nrow(data) == 2) %>% unnest(cols = c(data)) %>% arrange(d_v_h) %>% filter(d_v_h < 0)

both_md <- raw_target_data %>%
  filter(asv_names %in% both_list)

both_md_aov <- lmer(value ~ asv_names*final_disease_state*time + (1 | genotype) + (1 | tank), 
                   data = both_md)
anova(both_md_aov)

emmeans(both_md_aov, ~final_disease_state*time | asv_names, type = 'response') %>%
  cld(Letters = LETTERS, adjust = 'fdr') %>%
  as_tibble() %>%
  left_join(taxonomy_tibble, by = join_by(asv_names)) %>%
  mutate(graph_color = paste(time, final_disease_state, sep = "_")) %>%
  mutate(.group = str_trim(.group)) %>%
  left_join(more_disease_fds, by = join_by(asv_names), multiple = "all") %>%
  #rename(emmean = response) %>%
  ggplot(aes(x = fct_reorder(paste(Family, " ", Genus, " (", asv_names, ")", sep = ""), d_v_h, .desc = TRUE), y = emmean, ymin = emmean - SE, ymax = emmean + SE,
             colour = graph_color, pch = time)) +
  geom_pointrange(position = position_dodge(0.5)) +
  geom_text(aes(y = (emmean + SE), label = .group),
            position = position_dodge(0.5), vjust = -1) +
  geom_text(aes(y = 11, label = d_v_h)) +
  scale_color_manual(values = c("hotpink1", "deepskyblue", "firebrick1", "dodgerblue3")) +
  #scale_color_manual(values = wes_palette("Zissou1", 2, type = "discrete")) +
  coord_flip() +
  labs(title = "Very Likely Suspects - Both Final Disease State and FDS:time Interaction") +
  xlab("ASV Name")




#buffer
#### Misc. Plots ####

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
  #ggplot(aes(y = Genus, x = sample_id, size = value, colour = Order)) + 
  #geom_point() +
  
   ggplot(aes(x = sample_id, y = value, fill = Order)) +
   geom_col() +
  facet_wrap(~final_disease_state, scales = 'free_x') +
  theme(axis.text.x = element_blank())

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


#WIP genera that are more abundant in diseased in likely suspects list by final disease state, exposure, and time
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
  



#### Trees ####


tmp <- enframe(sequences, name = 'asv_names', value = 'sequence') %>%
  left_join(taxonomy_tibble, by = 'asv_names') %>%
  nest(data = c(Genus, Species, asv_names, sequence))

tmp %>%
  rowwise %>%
  mutate(n_asv = nrow(data)) %>%
  ungroup %>%
  arrange(n_asv) %>%
  filter(n_asv > 10)

tmp_dat <- tmp %>%
  filter(Family == 'Omnitrophaceae') %>%
  select(data) %>%
  unnest(data) %$%
  set_names(sequence, asv_names) %>%
  DNAStringSet()

?msa
msa_out <- msa(tmp_dat)
msa_phy <- msaConvert(msa_out, 'phangorn::phyDat')



mt <- modelTest(msa_phy)
ml_tree <- pml_bb(mt)
boot_ml <- bootstrap.pml(ml_tree, bs = 100, optNni = TRUE,
                         control = pml.control(trace = 1))
plotBS(midpoint(ml_tree$tree))
plotBS(midpoint(ml_tree$tree), boot_ml)
