multiple_aggregation_levels <- c("none", "Genus", "Family", "Order")

all_corr_results <- c()

#aggregation_level <-  "Genus"

agg_microbiome_data <- read_rds("../intermediate_files/updated_microbiome_data.rds")
for(var in multiple_aggregation_levels) {

aggregation_level <- var

if(aggregation_level != 'none'){
  agg_microbiome_data <- aggregate_taxa(agg_microbiome_data, aggregation_level)
  taxa_names(agg_microbiome_data) <- str_replace_all(taxa_names(agg_microbiome_data), ' |-', '_')
} else {
  sequences <- taxa_names(agg_microbiome_data)
  taxa_names(agg_microbiome_data) <- str_c('ASV', 1:length(taxa_names(agg_microbiome_data)), sep = '_')
  names(sequences) <- taxa_names(agg_microbiome_data)
}

agg_metadata <- sample_data(agg_microbiome_data) %>%
  as_tibble(rownames = 'sample_id') %>%
  dplyr::select(-retain_sample) %>%
  mutate(fragment_id = str_c(str_replace_na(exposure, 'NA'), tank, genotype, sep = '_'),
         .after = sample_id)

agg_melted_ps <- phyloseq_filter_prevalence(agg_microbiome_data, 
                                        prev.trh = 0.2) %>%
  psmelt() %>%
  as_tibble()

agg_longitudinal_genos <- agg_melted_ps %>%
  select(time, genotype) %>%
  group_by(genotype) %>%
  distinct() %>%
  reframe(all_timepoints = str_c(time, collapse = "_")) %>%
  filter(all_timepoints != "T0") %>% # remove 5 fragments used for making doses
  filter(str_detect(all_timepoints, "T3") & str_detect(all_timepoints, "T7")) %>% #removes 2 T7 only and 2 T0/T7
  pull(genotype)

agg_model_samples <- filter(agg_metadata, !str_detect(tank, 'homo|HOMO')) %>%
  filter(genotype %in% agg_longitudinal_genos) %>% #genos present in T3 and T7
  filter(!(exposure == "H" & final_disease_state == "D")) %>%
  pull(sample_id)


#### Make Venn showing ASVs to keep ####
agg_otu_timepoint_presence <- agg_melted_ps %>%
  mutate(across(c(exposure, final_disease_state), factor)) %>%
  filter(time %in% c('T3', 'T7') | (time == "T0" & tank == "HOMO")) %>%
  filter(Abundance > 0) %>%
  mutate(time = if_else(time == 'T0', exposure, time)) %>%
  #filter(OTU %in% bacterial_signature_asv$taxa_id) %>%
  group_by(time, OTU) %>%
  summarise(n = sum(Abundance),
            .groups = 'drop') %>%
  pivot_wider(names_from = time, values_from = n, values_fill = 0L) %>%
  mutate(across(-OTU, ~. > 0))

#ggvenn(otu_timepoint_presence, c('D', 'H', 'T3', 'T7')) + ggtitle("ASV Presence")

agg_venn_all_times_and_doses <- agg_melted_ps %>%
  mutate(across(c(exposure, final_disease_state), factor)) %>%
  filter(tank != "homogenate_fragment") %>%
  filter(Abundance > 0) %>%
  mutate(time = if_else(time == 'T0' & tank == "HOMO", exposure, time)) %>%
  #filter(OTU %in% venn_group) %>%
  group_by(time, OTU) %>%
  summarise(n = sum(Abundance),
            .groups = 'drop') %>%
  pivot_wider(names_from = time, values_from = n, values_fill = 0L) %>%
  mutate(across(-OTU, ~. > 0)) %>% 
  mutate(T0_H = ifelse(H | T0, TRUE, FALSE))

#ggvenn(venn_all_times_and_doses, c('D', 'T0_H', 'T3', 'T7')) + ggtitle("Minimal Filtering")

agg_not_t3t7_only <- filter(agg_venn_all_times_and_doses, 
                        !c(!D & !H & !T0 & T3 & T7 & !T0_H)) %>%
  filter(!c(!D & !H & !T0 & !T3 & T7 & !T0_H)) %>%
  filter(!c(!D & !H & !T0 & T3 & !T7 & !T0_H)) %>%
  pull(OTU)

#### Normalize based on all samples & any other filtering ####
if(aggregation_level == "none"){
  agg_otu_tmm <- agg_microbiome_data %>%
    phyloseq_filter_prevalence(prev.trh = 0.2) %>%
    otu_table() %>% 
    t %>% #NOTE: *genus and family do not need the t but ASVs need the t*
    as.data.frame %>%
    as.matrix %>% 
    DGEList(remove.zeros = TRUE) %>%
    
    #Add any other filtering here
    filter_missingness(agg_model_samples, 0.9) %>%
    filter_missing_groups(agg_metadata, 1) %>%
    
    edgeR::calcNormFactors(method = 'TMMwsp') %>%
    filter_venn(agg_not_t3t7_only) %>% #remove things that are only in T3 and T7
    #filter_venn(otus_to_analyze) %>% #only things in D and T3 and T7, aka potential pathogens
    filter_samples(agg_model_samples) #remove samples not to be analyzed
  
} else {
  agg_otu_tmm <- agg_microbiome_data %>%
    phyloseq_filter_prevalence(prev.trh = 0.2) %>%
    otu_table() %>% 
    #t %>% #NOTE: *genus and family do not need the t but ASVs need the t*
    as.data.frame %>%
    as.matrix %>% 
    DGEList(remove.zeros = TRUE) %>%
    
    #Add any other filtering here
    filter_missingness(agg_model_samples, 0.9) %>%
    filter_missing_groups(agg_metadata, 1) %>%
    
    edgeR::calcNormFactors(method = 'TMMwsp') %>%
    filter_venn(agg_not_t3t7_only) %>% #remove things that are only in T3 and T7
    #filter_venn(otus_to_analyze) %>% #only things in D and T3 and T7, aka potential pathogens
    filter_samples(agg_model_samples) #remove samples not to be analyzed
  
}

agg_venn_group <- agg_otu_tmm %>%
  cpm(log = TRUE, prior.count = 0.5,
      normalized.lib.sizes = TRUE) %>%
  as_tibble(rownames = 'taxa_id') %>% pull(taxa_id)

param <- SnowParam(parallel::detectCores() - 1, "SOCK", progressbar = TRUE)
agg_dream_weights_fullInteraction <- voomWithDreamWeights(counts = agg_otu_tmm,
                                                      formula = ~ model_comp + (1 | genotype) + (1 | tank),
                                                      
                                                      data = filter(agg_metadata, !str_detect(tank, 'homo|HOMO')) %>%
                                                        filter(genotype %in% agg_longitudinal_genos) %>% #genos present in T3 and T7
                                                        #filter(!(exposure == "H" & final_disease_state == "D")) %>%
                                                        arrange(sample_id) %>%
                                                        mutate(model_comp = str_c(time, exposure, susceptability)) %>%
                                                        column_to_rownames('sample_id'),
                                                      BPPARAM = param,
                                                      plot = TRUE)


#### ASV Modelling ####
agg_full_data <- agg_otu_tmm %>%
  cpm(log = TRUE, prior.count = 0.5,
      normalized.lib.sizes = TRUE) %>%
  as_tibble(rownames = 'taxa_id') %>%
  pivot_longer(cols = -taxa_id,
               names_to = 'sample_id',
               values_to = 'log2_cpm') %>%
  left_join(agg_dream_weights_fullInteraction$weights %>%
              set_colnames(colnames(agg_dream_weights_fullInteraction$E)) %>%
              set_rownames(rownames(agg_dream_weights_fullInteraction$E)) %>%
              as_tibble(rownames = 'taxa_id') %>%
              pivot_longer(cols = -taxa_id,
                           names_to = 'sample_id',
                           values_to = 'weight'),
            by = c('taxa_id', 'sample_id')) %>%
  left_join(as_tibble(agg_otu_tmm$counts, rownames = 'taxa_id') %>%
              pivot_longer(cols = -taxa_id,
                           names_to = 'sample_id',
                           values_to = 'read_count'),
            by = c('taxa_id', 'sample_id')) %>%
  
  left_join(agg_metadata, 
            by = 'sample_id') %>%
  left_join(as_tibble(agg_otu_tmm$samples, rownames = 'sample_id') %>%
              select(-group),
            by = 'sample_id') %>%
  left_join(tax_table(agg_microbiome_data) %>%
              as.data.frame() %>%
              as_tibble(rownames = 'taxa_id'),
            by = c('taxa_id'))

agg_normalized_asv_counts <- agg_full_data %>%
  mutate(time = factor(time, ordered = TRUE)) %>%
  mutate(final_disease_state = ifelse(time == "T0", "F", final_disease_state)) %>%
  filter(!c(exposure == "H" & final_disease_state == "D")) %>%
  mutate(treatment = str_c(time, exposure, final_disease_state, sep = '_'),
         time_exposure = str_c(time, exposure, sep = '_'),
         timeC = str_extract(time, '[0-9]+') %>% as.numeric) %>%#,
         #across(Domain:Species, str_replace_na)) %>%
  mutate(asv_number = str_extract(taxa_id, '[0-9]+') %>% as.integer)


t0_corr_test <- agg_normalized_asv_counts %>%
  filter(time == "T0") %>%
  group_by(taxa_id) %>%
  reframe(corr_val = broom::tidy(cor.test(log2_cpm, resistance))) %>%
  unnest(corr_val) %>%
  mutate(fdr_p.value = p.adjust(p.value, method = "fdr")) %>%
  arrange(fdr_p.value) %>%
  mutate(agg_level = aggregation_level)

all_corr_results <- rbind(all_corr_results, t0_corr_test)
}

all_corr_results %>% filter(fdr_p.value < 0.05)

#previous probiotic associations
all_corr_results %>% filter(taxa_id %in% c("MD3_55", "Endozoicomonas", "Myxococcales"))

asv_corr <- all_corr_results %>% 
  filter(agg_level == "none") %>% 
  left_join(normalized_asv_counts %>% select(asv_id, Order:Genus) %>% distinct(), 
            by = join_by("taxa_id" == "asv_id"))

genus_corr <- all_corr_results %>% 
  filter(agg_level == "Genus") %>% 
  left_join(normalized_asv_counts %>% select(Order:Genus) %>% distinct(), 
            by = join_by("taxa_id" == "Genus")) %>%
  mutate(Genus = taxa_id)

family_corr <- all_corr_results %>% 
  filter(agg_level == "Family") %>% 
  left_join(normalized_asv_counts %>% select(Order:Genus) %>% distinct(), 
            by = join_by("taxa_id" == "Family")) %>%
  mutate(Family = taxa_id, Genus = NA) %>%
  distinct() %>%
  filter(is.na(Order) | Order != "Chitinophagales") %>% #duplicate
  mutate(estimate = ifelse(is.na(estimate), 0, estimate)) %>%
  arrange(estimate) %>%
  mutate(Family = factor(Family, ordered = T, levels = .$Family))

corr_plot <- family_corr %>%
  rbind(genus_corr, asv_corr) %>%
  arrange(Family) %>%
  mutate(taxa_id = factor(taxa_id, ordered = T)) %>%
  filter(!is.na(Family) & Family != "NA")

ggplot(corr_plot) +
  geom_hline(yintercept = 0) +
  geom_point(aes(x = Family, y = estimate, col = agg_level), position = position_dodge(0.75)) + #data = (corr_plot %>% filter(agg_level == "none")),
  scale_color_manual(values = c("none" = "#CB3309", "Genus" = "#E6AC0E", "Family" = "#397DBB")) +
  coord_flip() +
  theme_bw()
#export 1200x900


ggplot(corr_plot, aes(x = Family, y = estimate, col = agg_level)) +
  geom_hline(yintercept = 0) +
  geom_boxplot(position = position_dodge(0.75)) +
  geom_point(data = corr_plot %>% filter(agg_level == "Family"), fill = "#397DBB", 
             position = position_dodge(0.75), pch = 21, size = 2) + #data = (corr_plot %>% filter(agg_level == "none")),
  scale_color_manual(values = c("none" = "#CB3309", "Genus" = "#E6AC0E", "Family" = "transparent"),
                     name = "Aggregation Level", breaks = c("none", "Genus", "Family"), 
                     labels = c("ASV", "Genus", "Family")) +
  coord_flip() +
  theme_bw() +
  ylab("Correlation Coefficient (r)")
#export 1200x900
