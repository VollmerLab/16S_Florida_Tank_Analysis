
#### Libraries ####
library(tidyverse)
library(magrittr)
library(phyloseq)
library(microbiome)
library(metagMisc)
library(edgeR)
library(variancePartition)
library(ggvenn)

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
    select(sample_id, group_var) %>%
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

#### Agregate Samples ####
if(aggregation_level != 'none'){
  microbiome_data <- aggregate_taxa(microbiome_data, aggregation_level)
  taxa_names(microbiome_data) <- str_replace_all(taxa_names(microbiome_data), ' |-', '_')
} else {
  #sequences <- taxa_names(microbiome_data)
  taxa_names(microbiome_data) <- str_c('ASV', 1:length(taxa_names(microbiome_data)), sep = '_')
  #names(sequences) <- taxa_names(microbiome_data)
}


#### Microshades Microbe Abundance ####


#### Make Venn showing ASVs to keep ####
otu_timepoint_presence <- phyloseq_filter_prevalence(microbiome_data, 
                                                     prev.trh = 0.1) %>%
  psmelt() %>%
  as_tibble  %>%
  mutate(across(c(exposure, final_disease_state), factor)) %>%
  filter(time %in% c('T3', 'T7') | (time == "T0" & tank == "HOMO")) %>%
  filter(Abundance > 0) %>%
  mutate(time = if_else(time == 'T0', exposure, time)) %>%
  group_by(time, OTU) %>%
  summarise(n = sum(Abundance),
            .groups = 'drop') %>%
  pivot_wider(names_from = time, values_from = n, values_fill = 0L) %>%
  mutate(across(-OTU, ~. > 0))

ggvenn(otu_timepoint_presence, c('D', 'H', 'T3', 'T7'))
  
otus_to_analyze <- filter(otu_timepoint_presence, 
                          (D & T3 & T7)) %>%
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
  filter_venn(otus_to_analyze) %>%
  filter_samples(model_samples) %>% #remove samples not to be analyzed
  filter_asv_meanCount(metadata, 100) #Remove ASVs with an average of less than N CPM per sample


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
  
  left_join(metadata, 
            by = 'sample_id') %>%
  left_join(tax_table(microbiome_data) %>%
              as.data.frame() %>%
              as_tibble(rownames = 'asv_id'),
            by = c('asv_id'))

write_csv(full_data, '../intermediate_files/fully_preprocessed_samples.csv.gz')
n_distinct(full_data$asv_id)

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
