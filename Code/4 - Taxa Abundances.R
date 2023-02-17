#code to examine changes to taxa abundances

setwd("~/Documents/GitHub/16S_Florida_Tank_Analysis/Code")

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
library(tidyverse)

select <- dplyr::select

aggregation_level <- 'Family' #or none

#### Read in Data ####
microbiome_data <- read_rds("../intermediate_files/preprocess_microbiome.rds") %>%
  subset_samples(time %in% c('T3', 'T7'))
metadata <- sample_data(microbiome_data) %>%
  as_tibble(rownames = 'sample_id') %>%
  select(-retain_sample) %>%
  mutate(fragment_id = str_c(str_replace_na(exposure, 'NA'), tank, genotype, sep = '_'),
         .after = sample_id)

if(aggregation_level != 'none'){
  microbiome_data <- aggregate_taxa(microbiome_data, aggregation_level)
  taxa_names(microbiome_data) <- str_replace_all(taxa_names(microbiome_data), ' |-', '_')
} else {
  taxa_names(microbiome_data) <- str_c('ASV', 1:length(taxa_names(microbiome_data)), sep = '_')
}

#### Filtering and normalizing data ####
otu_tmm <- microbiome_data %>%
  phyloseq_filter_prevalence(prev.trh = 0.1) %>%
  otu_table() %>% 
  as.data.frame %>%
  as.matrix %>% 
  DGEList(remove.zeros = TRUE) %>%
  edgeR::calcNormFactors(method = 'TMMwsp') #change method to TMM


cpm(otu_tmm, log = TRUE, prior.count = 2) %>%
  rowMeans %>%
  quantile(0.05)


cpm(otu_tmm, log = TRUE, prior.count = 2) %>%
  rowMeans %>%
  tibble(x = .) %>%
  ggplot(aes(x = x)) +
  geom_histogram(fill = "grey", bins = 100) +
  theme_classic() +
  labs(y = "Density", x = "Filtered number of taxa (logCPM)",
       title = "Distribution of normalized, filtered taxa")

plot_pcoa <- function(cpm_counts){
  filtered_pcoa <- t(cpm_counts) %>%
    vegdist(method = 'euclidean') %>%
    magrittr::divide_by(1000) %>%
    pcoa()
  
  percent_variance <- filtered_pcoa$values$Eigenvalues / sum(filtered_pcoa$values$Eigenvalue)
  
  filtered_pcoa$vectors %>%
    as_tibble(rownames = 'sample_id') %>%
    dplyr::select(sample_id, Axis.1, Axis.2) %>%
    inner_join(metadata,
               by = 'sample_id') %>%
    ggplot(aes(x = Axis.1, y = Axis.2, colour = final_disease_state, 
               shape = time, group = fragment_id)) +
    geom_point() +
    geom_path() +
    labs(x = str_c('PCoA 1 (', scales::percent(percent_variance[1]), ')'),
         y = str_c('PCoA 2 (', scales::percent(percent_variance[2]), ')')) +
    theme_classic()
  
}

cpm(otu_tmm, log = TRUE, prior.count = 2) %>% plot_pcoa() #x axis explains greatest amount of variance
  #y axis is next largest amount of variance
  #cpm, when given a DGEList, defaults to applying normalization factors for us

#nesting data by taxon abundances
taxon_abundances <- cpm(otu_tmm, log = TRUE, prior.count = 2) %>%
  t %>%
  as_tibble(rownames = "sample_id") %>%
  full_join(metadata, by = "sample_id") %>%
  pivot_longer(cols = -any_of(colnames(metadata)), 
               names_to = "taxon", values_to = "value") %>%
  mutate(across(c(exposure, final_disease_state, time), factor)) %>%
  nest_by(taxon)

taxon_abundances$taxon
data <- taxon_abundances$data[[31]]


norm_mod <- lmer(value ~ time * (exposure + final_disease_state) + (1 | fragment_id), data = data) 

rm_mod <- aov_4(value ~ time * (exposure + final_disease_state) + (time | fragment_id), data = data)

gamma_mod <- glmer(value ~ time * (exposure + final_disease_state) + (1 | fragment_id), data = data, family = Gamma(link = log))

bind_rows(
  gamma = as_tibble(emmeans(gamma_mod, ~time:final_disease_state, type = 'response')) %>%
    rename(emmean = response,
           lower.CL = asymp.LCL,
           upper.CL = asymp.UCL) ,
  normal = as_tibble(emmeans(norm_mod, ~time:final_disease_state, type = 'response')),
  rm = as_tibble(emmeans(rm_mod, ~time:final_disease_state, type = 'response')),
  .id = 'model'
) %>%
  ggplot(aes(x = interaction(time, final_disease_state), y = emmean, ymin = lower.CL, ymax = upper.CL, colour = model)) +
  geom_pointrange(position = position_dodge(0.5))
  
  
#### Model Each Taxon Independently ####
all_models <- taxon_abundances %>%
  mutate(model = list(aov_4(value ~ time * (exposure + final_disease_state) + (time | fragment_id), data = data)))


all_models %>%
  reframe(model = list(model), 
          sig_terms = find_unique_significant_terms(model, 0.05))

#### Functions #####

#needs to be tidied

find_unique_significant_terms <- function(model, alpha){
  significant_terms <- model$anova_table %>%
    as_tibble(rownames = 'param') %>%
    janitor::clean_names() %>%
    filter(pr_f < alpha) %>%
    pull(param)
  
  unique_values <- outer(significant_terms, significant_terms, str_count) %>%
    colSums() %>%
    equals(1)
  
  str_replace(significant_terms[unique_values], ':', '*')
}

make_emmean_model <- function(model, form, alpha){
  emmeans(model, specs = form, type = 'response') %>%
    cld(Letters = LETTERS, reversed = TRUE, alpha = alpha) 
}

make_model_plot <- function(predOut, rawData, var_inclusion){
  
  all_vars <- str_extract_all(var_inclusion, 'time|final_disease_state|exposure') %>%
    unlist
  
  if(any(str_detect(all_vars, 'time'))){
    x_var <- str_subset(all_vars, 'time.*')
  } else {
    x_var <- all_vars[1]
  }
  
  if(length(all_vars) > 1){
    colour_var <- str_subset(all_vars, x_var, negate = TRUE)
  } else {
    colour_var <- NULL
  }
  
  y_var <- str_subset(colnames(predOut), 'emmean|response')
  y_data <- str_subset(colnames(rawData), 'value')
  
  if(is.null(colour_var)){
    as_tibble(predOut) %>%
      mutate(.group = str_trim(.group)) %>%
      ggplot(aes(x = !!sym(x_var))) +
      
      stat_halfeye(data = rawData, aes(y = !!sym(y_data)),
                   adjust = 0.5, width = 0.6, .width = 0, 
                   alpha = 0.5, show.legend = FALSE,
                   fatten_point = 0, justification = -0.25, 
                   position = position_dodge(0.5), size = 0) +
      
      geom_half_point(data = rawData, aes(y = !!sym(y_data)),
                      side = 'l', range_scale = 0.1, alpha = 1,
                      position = position_dodge(width = 0.5),
                      show.legend = FALSE,
                      transformation = position_jitter(height = 0, width = 0.05)) +
      
      geom_pointrange(aes(y = !!sym(y_var),
                          ymin = lower.CL, ymax = upper.CL),
                      position = position_dodge(0.5)) +
      geom_text(aes(y = (upper.CL), label = .group),
                position = position_dodge(0.5), vjust = -0.1, 
                show.legend = FALSE) +
      labs(x = case_when(x_var %in% c('time', 'timepoint') ~ 'Time (d)',
                         x_var == 'exposure' ~ 'Exposure',
                         x_var == 'final_disease_state' ~ 'Disease State'))
    
  } else {
    as_tibble(predOut) %>%
      mutate(.group = str_trim(.group)) %>%
      ggplot(aes(x = !!sym(x_var), colour = !!sym(colour_var))) +
      
      stat_halfeye(data = rawData, aes(y = !!sym(y_data), fill = !!sym(colour_var)),
                   adjust = 0.5, width = 0.6, .width = 0, 
                   alpha = 0.5, show.legend = FALSE,
                   fatten_point = 0, justification = -0.25, 
                   position = position_dodge(0.5), size = 0) +
      
      geom_half_point(data = rawData, aes(y = !!sym(y_data), colour = !!sym(colour_var)),
                      side = 'l', range_scale = 0.1, alpha = 1,
                      position = position_dodge(width = 0.5),
                      show.legend = FALSE,
                      transformation = position_jitter(height = 0, width = 0.05)) +
      
      geom_pointrange(aes(y = !!sym(y_var),
                          ymin = lower.CL, ymax = upper.CL),
                      position = position_dodge(0.5)) +
      geom_text(aes(y = (upper.CL), label = .group),
                position = position_dodge(0.5), vjust = -0.1,
                show.legend = FALSE) +
      labs(x = case_when(x_var %in% c('time', 'timepoint') ~ 'Time (d)',
                         x_var == 'exposure' ~ 'Exposure',
                         x_var == 'final_disease_state' ~ 'Disease State'),
           colour = case_when(colour_var %in% c('time', 'timepoint') ~ 'Time (d)',
                              colour_var == 'exposure' ~ 'Exposure',
                              colour_var == 'final_disease_state' ~ 'Disease State'))
  }
  
}

make_aov_summary <- function(model){
  model$anova_table %>%
    as_tibble(rownames = 'effect') %>%
    janitor::clean_names() %>%
    mutate(across(c(num_df, den_df), round, digits = 3),
           df = str_c(num_df, den_df, sep = ', '),
           p = pr_f) %>%
    dplyr::select(effect, df, f, p) %>%
    pivot_wider(names_from = 'effect',
                values_from = c('df', 'f', 'p'), 
                names_vary = 'slowest')
}

make_model_plot(tmp$em_out[[1]], tmp$data[[1]], tmp$terms[[1]])

tmp <- all_models %>%
  ungroup %>%
  rowwise %>%
  mutate(terms = list(find_unique_significant_terms(model, 0.05))) %>%
  unnest(terms, keep_empty = TRUE) %>%
  rowwise %>%
  mutate(em_out = list(possibly(make_emmean_model, otherwise = NULL)(model, 
                                                                     as.formula(str_c('~', terms)), 
                                                                     0.05))) %>%
  mutate(plot = list(possibly(make_model_plot, otherwise = NULL, quiet = TRUE)(em_out, data, terms))) %>%
  group_by(taxon, data, model) %>%
  reframe(plot = ifelse(any(is.na(terms)),
                        list(NULL),
                        list(wrap_plots(plot) & 
                               labs(y = 'log2(CPM)') &
                               plot_annotation(title = taxon) & 
                               theme_classic() &
                               theme(panel.background = element_rect(colour = 'black', fill = NA),
                                     axis.text = element_text(colour = 'black', size = 12),
                                     axis.title = element_text(colour = 'black', size = 16))))) %>%
  rowwise %>%
  mutate(possibly(make_aov_summary, otherwise = NULL)(model)) %>%
  ungroup 

tmp$model[[32]]
tmp$plot[[31]]


tmp %>% 
  filter(`p_final_disease_state:time` < 0.05,
         p_exposure < 0.05,
         `p_exposure:time` > 0.05) %>%
  slice(2) %>%
  pull(model) %>%
  pluck(1)
