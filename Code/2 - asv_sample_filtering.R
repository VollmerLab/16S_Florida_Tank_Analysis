setwd("/Users/emilytrytten/Desktop/Screenshots/Career/Vollmer Lab/GitHub/16S_Florida_Tank_Analysis/Code")

#### Libraries ####
library(tidyverse)
library(magrittr)
library(phyloseq)
library(microbiome)
library(metagMisc)
library(edgeR)
library(variancePartition)
library(ggvenn)
library(cowplot)
library(ComplexUpset)
library(microshades)

#### Functions ####
# data <- otu_tmm; prop_missing <- 0.25
filter_missingness <- function(data, model_samples, prop_missing){
  #data = a DGElist object, model_sample = samples to use for the calculation of % missingness
  #prop_missing = maximum percentage of samples which can be 0 and still keep ASV in dataset
  keep <- rowMeans(data$counts[,colnames(data$counts) %in% model_samples] == 0) <= prop_missing
  
  message('\n')
  message('ASVs Removed for not being expressed in enough samples: ', scales::comma(table(keep)[1]))
  message('ASVs Kept by Filter: ', scales::comma(table(keep)[2]))
  message('\n')
  data[keep, keep.lib.sizes = FALSE]
}

# max_missing_group <- 0.5
filter_missing_groups <- function(data, meta, max_missing_group){
  meta <- mutate(metadata, group_var = str_c(time, exposure, susceptability)) %>%
    dplyr::select(sample_id, group_var) %>%
    filter(!is.na(group_var))
  
  keep <- as_tibble(data$counts, rownames = 'gene_id') %>%
    pivot_longer(cols = -gene_id,
                 names_to = 'sample_id') %>%
    inner_join(meta, by = 'sample_id') %>%
    group_by(gene_id, group_var) %>%
    summarise(prop_missing = mean(value == 0),
              .groups = 'drop') %>%
    pivot_wider(names_from = group_var,
                values_from = prop_missing) %>%
    filter(!if_any(where(is.numeric), ~. > max_missing_group)) %>%
    pull(gene_id)
  keep <- rownames(data) %in% keep
  
  message('\n')
  message('ASVs Removed for having more than ', scales::percent(max_missing_group), ' missing data in at least one group: ', scales::comma(sum(!keep)[1]))
  message('ASVs Kept by Filter: ', scales::comma(sum(keep)))
  message('\n')
  data[keep, keep.lib.sizes = FALSE]
}

filter_samples <- function(data, model_samples){
  #data is a DGEList
  keep <- colnames(data) %in% model_samples
  data[,keep, keep.lib.sizes = TRUE]
}

# venn_asv <- otus_to_analyze
filter_venn <- function(data, venn_asv){
  #venn should be a vector of asv's which are kept by the venn diagram
  keep <- rownames(data$counts) %in% venn_asv
  message('Retained ', sum(keep), ' ASVs found in overlap of Venn diagram')
  data[keep, keep.lib.sizes = TRUE]
}

# data <- otu_tmm; min_count <- 100; meta <- metadata
filter_asv_meanCount <- function(data, meta, min_count){
  mean_counts <- cpm(data, normalized.lib.sizes = FALSE, log = FALSE,
                     prior.count = 0) %>%
    as_tibble(rownames = 'asv_id') %>%
    pivot_longer(cols = -asv_id,
                 names_to = 'sample_id',
                 values_to = 'n') %>%
    left_join(metadata, by = 'sample_id') %>%
    group_by(asv_id) %>%
    summarise(mean_count = mean(n)) %>%
    pull(mean_count)
  
  keep <- mean_counts >= min_count
  message('\n')
  message('ASVs Removed for having less than ', scales::comma(min_count), ' mean number of reads across all samples')
  message('ASVs Kept by Filter: ', scales::comma(sum(keep)))
  message('\n')
  data[keep, keep.lib.sizes = TRUE]
}

#### Data ####
aggregation_level <- 'none' #or none

microbiome_data <- read_rds("../intermediate_files/preprocess_microbiome.rds")
metadata <- sample_data(microbiome_data) %>%
  as_tibble(rownames = 'sample_id') %>%
  select(-retain_sample)

model_samples <- filter(metadata, !str_detect(tank, 'homo|HOMO')) %>% pull(sample_id)

#### Aggregate Samples ####
if(aggregation_level != 'none'){
  microbiome_data <- aggregate_taxa(microbiome_data, aggregation_level)
  taxa_names(microbiome_data) <- str_replace_all(taxa_names(microbiome_data), ' |-', '_')
} else {
  sequences <- taxa_names(microbiome_data)
  taxa_names(microbiome_data) <- str_c('ASV', 1:length(taxa_names(microbiome_data)), sep = '_')
  names(sequences) <- taxa_names(microbiome_data)
}


#### Microshades Microbe Abundance ####

mdf_prep_test1 <- microbiome_data %>%
  tax_glom("Genus") %>%
  psmelt() 

mdf_processed_data <- mdf_prep_test1 %>%
  filter(Abundance > 0) %>%
  mutate(category = paste(time, exposure, final_disease_state, 
                          sep = "_")) %>%
  mutate(category = ifelse(str_detect(category, "F"), "T0", category)) %>%
  filter(!category %in% c("T3_H_D", "T7_H_D")) %>%
  mutate(time = ifelse(category %in% c("T0_D_D", "T0_H_H"), "Doses", time)) %>%
  mutate(category = ifelse(category == "T0_D_D", "Diseased", ifelse(category == "T0_H_H", "Healthy", category))) %>%
  mutate(category = ifelse(time == "T3", paste(time, exposure, sep = "_"), category)) %>%
  group_by(category) %>%
  mutate(total = sum(Abundance)) %>%
  ungroup() %>%
  select(-c(retain_sample)) %>%
  group_by(category, Genus) %>%
  reframe(Kingdom, Phylum, Class, Order, Family, time, total, rel_abun = sum(Abundance)/total) %>%
  distinct() %>%
  rename(Sample = category, Abundance = rel_abun) %>%
  mutate(Sample = factor(Sample, levels = c("Healthy", "Diseased", "T0", "T3_H", "T3_D", "T7_H_H", "T7_D_H", "T7_D_D"))) %>%
  as.data.frame()

## ORDER GENUS

color_objs_ordergenus <- create_color_dfs(mdf_processed_data, group_level = "Order", 
                                          selected_groups = c("Rickettsiales", "Enterobacterales", "Flavobacteriales", 
                                                              "Pseudomonadales",  "Rhodobacterales"), cvd = TRUE)
mdf_ordergenus <- color_objs_ordergenus$mdf
cdf_ordergenus <- color_objs_ordergenus$cdf

legend_ordergenus <-custom_legend(mdf_ordergenus, cdf_ordergenus, group_level = "Order")

plot_ordergenus_prelim <- plot_microshades(mdf_ordergenus, cdf_ordergenus) + 
  scale_y_continuous(labels = scales::percent, expand = expansion(0)) +
  facet_grid(cols = vars(time), scales = "free", space = "free") +
  theme_bw() +
  theme(legend.position = "none", plot.margin = margin(6,20,6,6)) +
  #labs(title = "Order Genus") +
  scale_x_discrete(name = element_blank(), labels = c("T0" = "Field", "T3_H" = "Healthy", "T3_D" =
                                "Disease", "T7_H_H" = "Healthy (Healthy)",
                              "T7_D_H" = "Disease (Healthy)", "T7_D_D" = "Disease (Diseased)"))

plot_grid(plot_ordergenus_prelim, legend_ordergenus,  rel_widths = c(1, .25))


#ORDER GENUS - SUSCEPTIBILITY

mdf_processed_suscep_data <- mdf_prep_test1 %>%
  filter(Abundance > 0) %>%
  mutate(category = paste(time, exposure, susceptability, 
                          sep = "_")) %>%
  mutate(category = ifelse(str_detect(category, "F"), paste("T0", susceptability, sep = "_"), category)) %>%
  mutate(time = ifelse(category %in% c("T0_D_NA", "T0_H_NA"), "Doses", time)) %>%
  mutate(category = ifelse(category == "T0_D_NA", "Diseased", ifelse(category == "T0_H_NA", "Healthy", category))) %>%
  group_by(category) %>%
  mutate(total = sum(Abundance)) %>%
  ungroup() %>%
  select(-c(retain_sample)) %>%
  group_by(category, Genus) %>%
  reframe(Domain, Phylum, Class, Order, Family, time, total, rel_abun = sum(Abundance)/total) %>%
  distinct() %>%
  rename(Sample = category, Abundance = rel_abun) %>%
  mutate(Sample = factor(Sample, levels = c("Healthy", "Diseased", "T0_S", "T0_R", "T3_H_S", "T3_H_R", "T3_D_S", "T3_D_R", "T7_H_S", "T7_H_R", "T7_D_S", "T7_D_R"))) %>%
  as.data.frame()

color_objs_suscep <- create_color_dfs(mdf_processed_suscep_data, group_level = "Order", 
                                          selected_groups = c("Rickettsiales", "Alteromonadales", "Oceanospirillales", 
                                                              "Francisellales",  "Flavobacteriales"), cvd = TRUE)
mdf_suscep <- color_objs_suscep$mdf
cdf_suscep <- color_objs_suscep$cdf

legend_suscep <-custom_legend(mdf_suscep, cdf_suscep, group_level = "Order")

plot_suscep_prelim <- plot_microshades(mdf_suscep, cdf_suscep) + 
  scale_y_continuous(labels = scales::percent, expand = expansion(0)) +
  facet_grid(cols = vars(time), scales = "free", space = "free") +
  theme_bw() +
  theme(legend.position = "none", plot.margin = margin(6,20,6,6)) +
  labs(title = "Order Genus - Susceptibility")

plot_grid(plot_suscep_prelim, legend_suscep,  rel_widths = c(1, .25))


#### Make Venn showing ASVs to keep ####
otu_timepoint_presence <- phyloseq_filter_prevalence(microbiome_data, 
                                                     prev.trh = 0.2) %>%
  psmelt() %>%
  as_tibble() %>%
  mutate(across(c(exposure, final_disease_state), factor)) %>%
  filter(time %in% c('T3', 'T7') | (time == "T0" & tank == "HOMO")) %>%
  filter(Abundance > 0) %>%
  mutate(time = if_else(time == 'T0', exposure, time)) %>%
  #filter(OTU %in% bacterial_signature_asv$asv_id) %>%
  group_by(time, OTU) %>%
  summarise(n = sum(Abundance),
            .groups = 'drop') %>%
  pivot_wider(names_from = time, values_from = n, values_fill = 0L) %>%
  mutate(across(-OTU, ~. > 0))

ggvenn(otu_timepoint_presence, c('D', 'H', 'T3', 'T7')) + ggtitle("ASV Presence")

venn_all_times_and_doses <- phyloseq_filter_prevalence(microbiome_data, 
                                         prev.trh = 0.2) %>%
  psmelt() %>%
  as_tibble  %>%
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

ggvenn(venn_all_times_and_doses, c('D', 'T0_H', 'T3', 'T7')) + ggtitle("Minimal Filtering")

upset(venn_all_times_and_doses %>% left_join(taxonomy_tibble, by = c("OTU" = "asv_names")), 
      c("T7", "T3", "T0", "D", "H"), 
      
      base_annotations=list(
        'Intersection size'=intersection_size(counts=T, text = aes(size = 6), fill = "slategray4")
      ),
      
      queries=list(upset_query(set='D', color="#DF0000", fill = "#DF0000"),
                   upset_query(set='T0', color="#D98EFF", fill = "#D98EFF"),
                   upset_query(set='T3', color="#B21BFF", fill = "#B21BFF"),
                   upset_query(set='T7', color="#650197", fill = "#650197"),
                   upset_query(set='H', color="#0FAB02", fill = "#0FAB02")),
      
      name='asv_names', width_ratio=0.1, min_size = 1, sort_sets = FALSE) + 
  ggtitle("Minimal Filtering")

otus_to_analyze <- filter(otu_timepoint_presence, 
                          (D & T3 & T7)) %>%
  pull(OTU)

not_t3t7_only <- filter(venn_all_times_and_doses, 
                    !c(!D & !H & !T0 & T3 & T7 & !T0_H)) %>%
  filter(!c(!D & !H & !T0 & !T3 & T7 & !T0_H)) %>%
  filter(!c(!D & !H & !T0 & T3 & !T7 & !T0_H)) %>%
  pull(OTU)
  
#### Normalize based on all samples & any other filtering ####
otu_tmm <- microbiome_data %>%
  phyloseq_filter_prevalence(prev.trh = 0.2) %>%
  otu_table() %>% 
  t %>% #NOTE: *genus and family do not need the t but ASVs need the t*
  as.data.frame %>%
  as.matrix %>% 
  DGEList(remove.zeros = TRUE) %>%
  
  #Add any other filtering here
  filter_missingness(model_samples, 0.9) %>%
  filter_missing_groups(metadata, 1) %>%
 
  edgeR::calcNormFactors(method = 'TMMwsp') %>%
  filter_venn(not_t3t7_only) %>% #remove things that are only in T3 and T7
  #filter_venn(otus_to_analyze) %>% #only things in D and T3 and T7, aka potential pathogens
  filter_samples(model_samples) %>% #remove samples not to be analyzed
  filter_asv_meanCount(metadata, 100) #Remove ASVs with an average of less than N CPM per sample


venn_group <- otu_tmm %>%
  cpm(log = TRUE, prior.count = 0.5,
      normalized.lib.sizes = TRUE) %>%
  as_tibble(rownames = 'asv_id') %>% pull(asv_id)

#### Plot Group by number of ASVs = 0 ####
otu_tmm$counts %>%
  as_tibble(rownames = 'asv_id') %>%
  pivot_longer(cols = -asv_id,
               names_to = 'sample_id',
               values_to = 'n') %>%
  left_join(metadata, by = 'sample_id') %>%
  group_by(time, exposure, susceptability, asv_id) %>%
  summarise(pct_missing = sum(n == 0) / n(),
            .groups = 'drop_last') %>%
  summarise(mean_missing = mean(pct_missing),
            sd_missing = sd(pct_missing),
            .groups = 'drop') %>%
  ggplot(aes(x = time, y = mean_missing, ymin = mean_missing - sd_missing, 
             ymax = mean_missing + sd_missing, 
             colour = interaction(exposure, susceptability))) +
  geom_pointrange(position = position_dodge(0.5)) +
  labs(x = NULL,
       y = 'Average Percent of Samples ASVs are 0 counts from')


#### Variance weighting ####
param <- SnowParam(parallel::detectCores() - 1, "SOCK", progressbar = TRUE)
dream_weights_fullInteraction <- voomWithDreamWeights(counts = otu_tmm, 
                                  formula = ~ model_comp + (1 | genotype) + (1 | tank),
                                  
                                  data = filter(metadata, !str_detect(tank, 'homo|HOMO')) %>%
                                    arrange(sample_id) %>%
                                    mutate(model_comp = str_c(time, exposure, susceptability)) %>%
                                    column_to_rownames('sample_id'),
                                  BPPARAM = param, 
                                  plot = TRUE)
#87917 vs 36251
#### ASV Modelling ####
full_data <- otu_tmm %>%
  cpm(log = TRUE, prior.count = 0.5,
      normalized.lib.sizes = TRUE) %>%
  as_tibble(rownames = 'asv_id') %>%
  pivot_longer(cols = -asv_id,
               names_to = 'sample_id',
               values_to = 'log2_cpm') %>%
  left_join(dream_weights_fullInteraction$weights %>% 
              set_colnames(colnames(dream_weights_fullInteraction$E)) %>%
              set_rownames(rownames(dream_weights_fullInteraction$E)) %>%
              as_tibble(rownames = 'asv_id') %>%
              pivot_longer(cols = -asv_id,
                           names_to = 'sample_id',
                           values_to = 'weight'),
            by = c('asv_id', 'sample_id')) %>%
  left_join(as_tibble(otu_tmm$counts, rownames = 'asv_id') %>%
              pivot_longer(cols = -asv_id,
                           names_to = 'sample_id',
                           values_to = 'read_count'),
            by = c('asv_id', 'sample_id')) %>%
  
  left_join(metadata, 
            by = 'sample_id') %>%
  left_join(as_tibble(otu_tmm$samples, rownames = 'sample_id') %>%
              select(-group),
            by = 'sample_id') %>%
  left_join(tax_table(microbiome_data) %>%
              as.data.frame() %>%
              as_tibble(rownames = 'asv_id'),
            by = c('asv_id'))


full_data %>% group_by(Order) %>% summarize(counts = sum(log2_cpm)) %>% arrange(desc(counts))


#filtered for the 
write_csv(full_data, '../intermediate_files/fully_preprocessed_samples.csv.gz')


### Output ASV CPMs for doses ####
homogenate_data <- microbiome_data %>%
  subset_samples(str_detect(tank, 'HOMO')) %>%
  otu_table %>%
  t() %>%
  as.data.frame %>%
  as.matrix %>% 
  DGEList(remove.zeros = FALSE) %>%
  cpm(log = TRUE, prior.count = 0.5,
      normalized.lib.sizes = TRUE) %>%
  as_tibble(rownames = 'asv_id') %>%
  pivot_longer(cols = -asv_id,
               names_to = 'sample_id',
               values_to = 'log2_cpm') %>%
  filter(asv_id %in% unique(full_data$asv_id)) %>%
  mutate(exposure = str_extract(sample_id, '[DH]'))
write_csv(homogenate_data, '../intermediate_files/homogenate_cpm.csv')

#### Plot ASV number vs mean cpm ####
otu_tmm %>%
  cpm(log = TRUE, prior.count = 0.5,
      normalized.lib.sizes = TRUE) %>%
  as_tibble(rownames = 'asv_id') %>%
  pivot_longer(cols = -asv_id,
               names_to = 'sample_id',
               values_to = 'log2_cpm') %>%
  mutate(asv_number = str_extract(asv_id, '[0-9]+') %>% as.integer()) %>%
  group_by(asv_number) %>%
  summarise(log2_cpm_mean = mean(log2_cpm),
            log2_cpm_se = sd(log2_cpm) / sqrt(n())) %>%
  ggplot(aes(x = asv_number, y = log2_cpm_mean, 
             ymin = log2_cpm_mean - log2_cpm_se, ymax = log2_cpm_mean + log2_cpm_se)) +
  geom_pointrange() +
  geom_vline(xintercept = 200, colour = 'red')


otu_tmm %>%
  cpm(log = TRUE, prior.count = 0.5,
      normalized.lib.sizes = TRUE) %>%
  as_tibble(rownames = 'asv_id') %>%
  pivot_longer(cols = -asv_id,
               names_to = 'sample_id',
               values_to = 'log2_cpm') %>%
  mutate(asv_number = str_extract(asv_id, '[0-9]+') %>% as.integer()) %>%
  group_by(asv_number) %>%
  summarise(log2_cpm_mean = mean(log2_cpm),
            log2_cpm_se = sd(log2_cpm) / sqrt(n()),
            sum_log2_cpm = sum(log2_cpm)) %>%
  ggplot(aes(x = log2_cpm_mean, y = sum_log2_cpm)) +
  geom_point() 

otu_tmm %>%
  cpm(log = TRUE, prior.count = 0.5,
      normalized.lib.sizes = TRUE) %>%
  as_tibble(rownames = 'asv_id') %>%
  pivot_longer(cols = -asv_id,
               names_to = 'sample_id',
               values_to = 'log2_cpm') %>%
  mutate(asv_number = str_extract(asv_id, '[0-9]+') %>% as.integer()) %>%
  group_by(asv_number) %>%
  summarise(log2_cpm_mean = mean(log2_cpm),
            log2_cpm_se = sd(log2_cpm) / sqrt(n())) %>%
  filter(asv_number <= 200)

otu_tmm %>%
  cpm(log = FALSE, prior.count = 0,
      normalized.lib.sizes = FALSE) %>%
  as_tibble(rownames = 'asv_id') %>%
  pivot_longer(cols = -asv_id,
               names_to = 'sample_id',
               values_to = 'log2_cpm') %>%
  mutate(asv_number = str_extract(asv_id, '[0-9]+') %>% as.integer()) %>%
  group_by(asv_number) %>%
  summarise(mean_count = mean(log2_cpm),
            sqrt_sd = sqrt(sd(log2_cpm)),
            sum_count = sum(log2_cpm)) %>%
  ggplot(aes(x = mean_count, y = sqrt_sd, colour = asv_number < 200)) +
  geom_point() +
  geom_vline(xintercept = 100, colour = 'darkgreen') +
  scale_x_continuous(limits = c(0, 1000))



otu_tmm %>%
  cpm(log = FALSE, prior.count = 0.5,
      normalized.lib.sizes = FALSE) %>%
  as_tibble(rownames = 'asv_id') %>%
  pivot_longer(cols = -asv_id,
               names_to = 'sample_id',
               values_to = 'log2_cpm') %>%
  mutate(asv_number = str_extract(asv_id, '[0-9]+') %>% as.integer()) %>%
  group_by(asv_number) %>%
  summarise(log2_cpm_mean = mean(log2_cpm),
            sqrt_sd = sqrt(sd(log2_cpm))) %>%
  filter(log2_cpm_mean >= 100)
