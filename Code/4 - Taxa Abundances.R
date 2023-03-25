#code to examine changes to taxa abundances

setwd("~/Desktop/Screenshots/Career/Vollmer Lab/GitHub/16S_Florida_Tank_Analysis/Code")

#### Packages ####

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
library(tidyverse)

select <- dplyr::select


#### Functions ####
find_unique_significant_terms_rmANOVA <- function(model, alpha){
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

find_unique_significant_terms <- function(model, alpha){
  significant_terms <- car::Anova(model) %>%
    as_tibble(rownames = 'param') %>%
    janitor::clean_names() %>%
    filter(pr_chisq < alpha) %>%
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

make_lmer_aov_summary <- function(model){
  car::Anova(model) %>%
    as_tibble(rownames = 'effect') %>%
    janitor::clean_names() %>%
    mutate(df = round(df, digits = 3), p = pr_chisq) %>%
    dplyr::select(effect, df, chisq, p) %>%
    pivot_wider(names_from = 'effect',
                values_from = c('df', 'chisq', 'p'), 
                names_vary = 'slowest')
}

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

plot_agg_pcoa <- function(cpm_counts){
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
         y = str_c('PCoA 2 (', scales::percent(percent_variance[2]), ')'),
         title = agg_title) +
    theme_classic()
  
}

get_tidy_p_values <- function(model){
  #for tidy(anova(aov(model))), extract the p-values only
  p_vals <- model %>% 
    dplyr::select(p.value) %>%
    filter(!is.na(p.value)) %>%
    pull()
  
  return(p_vals)
}
#### Patchwork PCoA Plot for ASVs, Genus, Family ####
microbiome_data <- read_rds("../intermediate_files/preprocess_microbiome.rds") %>%
  subset_samples(time %in% c('T3', 'T7'))
metadata <- sample_data(microbiome_data) %>%
  as_tibble(rownames = 'sample_id') %>%
  dplyr::select(-retain_sample) %>%
  mutate(fragment_id = str_c(str_replace_na(exposure, 'NA'), tank, genotype, sep = '_'),
         .after = sample_id)

agg_levels <- c("none", "Genus", "Family")

for(i in 1:length(agg_levels)){
  aggregation_level = agg_levels[i]
  
  if(aggregation_level != 'none'){
    microbiome_data <- aggregate_taxa(microbiome_data, aggregation_level)
    taxa_names(microbiome_data) <- str_replace_all(taxa_names(microbiome_data), ' |-', '_')
    agg_title <-  aggregation_level
  } else {
    taxa_names(microbiome_data) <- str_c('ASV', 1:length(taxa_names(microbiome_data)), sep = '_')
    agg_title <-  "ASV"
  }
  
  if(i == 1){
    otu_tmm <- microbiome_data %>%
      phyloseq_filter_prevalence(prev.trh = 0.1) %>%
      otu_table() %>% 
      t %>% #NOTE: *genus and family do not need the t but ASVs need the t*
      as.data.frame %>%
      as.matrix %>% 
      DGEList(remove.zeros = TRUE) %>%
      edgeR::calcNormFactors(method = 'TMMwsp')
    plota <- plot_agg_pcoa(cpm(otu_tmm, log = TRUE, prior.count = 2))
    
  }else{ 
    otu_tmm <- microbiome_data %>%
      phyloseq_filter_prevalence(prev.trh = 0.1) %>%
      otu_table() %>% 
      as.data.frame %>%
      as.matrix %>% 
      DGEList(remove.zeros = TRUE) %>%
      edgeR::calcNormFactors(method = 'TMMwsp')
    if(i == 2){
      plotb <- plot_agg_pcoa(cpm(otu_tmm, log = TRUE, prior.count = 2)) +
          scale_x_reverse() +
          scale_y_reverse()
    }else if(i == 3){
      plotc <- plot_agg_pcoa(cpm(otu_tmm, log = TRUE, prior.count = 2))
    }
  }
}
plota + plotb + plotc

#### Read in Data ####

aggregation_level <- 'Family' #or none

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
#FIXME NOTE: set # below for whether working with ASV level or not
otu_tmm <- microbiome_data %>%
  phyloseq_filter_prevalence(prev.trh = 0.1) %>%
  otu_table() %>% 
  #t %>% #NOTE: *genus and family do not need the t but ASVs need the t*
  as.data.frame %>%
  as.matrix %>% 
  DGEList(remove.zeros = TRUE) %>%
  edgeR::calcNormFactors(method = 'TMMwsp') #TMMwsp is for high prevalence of 0s
  
#cpm is counts per million, can be used as a descriptive measure for the expression level of a gene
cpm(otu_tmm, log = TRUE, prior.count = 2) %>% #prior.count is # to add to each value so not log(0)
  rowMeans %>%
  quantile(0.05) 


cpm(otu_tmm, log = TRUE, prior.count = 2) %>%
  rowMeans %>% #get ave value for each group in aggregation level (ex: family)
  tibble(x = .) %>% #remove labels, only keep ave values
  ggplot(aes(x = x)) +
  geom_histogram(fill = "grey", bins = 100) +
  theme_classic() +
  labs(y = "Density", x = "Filtered number of taxa (logCPM)",
       title = "Distribution of filtered and normalized taxa")

cpm(otu_tmm, log = TRUE, prior.count = 2) %>% plot_pcoa() 

  #x axis explains greatest amount of variance, y axis is next largest amount of variance
  #cpm, when given a DGEList, defaults to applying normalization factors for us
  #time explains 29% of variation and final disease state explains 9%

  #PCA focuses on shared variance: it tries to summarize multiple variables in the minimum number 
  #of components so that each component explains the most variance. 
  #PCoA on the other hand focuses on distances, and it tries to extract the dimensions that account 
  #for the maximum distances

taxonomy_tibble <- tax_table(microbiome_data) %>% 
  as.data.frame %>%
  as_tibble(rownames = "asv_names")

#nesting by ASV
taxon_abundances <- cpm(otu_tmm, log = TRUE, prior.count = 2) %>%
  t %>%
  as_tibble(rownames = "sample_id") %>%
  full_join(metadata, by = "sample_id") %>%
  pivot_longer(cols = -any_of(colnames(metadata)), 
               names_to = "asv_names", values_to = "value") %>%
  mutate(across(c(exposure, final_disease_state, time), factor)) %>%
  nest_by(asv_names)
  

#nesting data by taxon abundances (aggregated data)
taxon_abundances <- cpm(otu_tmm, log = TRUE, prior.count = 2) %>%
  t %>%
  as_tibble(rownames = "sample_id") %>%
  full_join(metadata, by = "sample_id") %>%
  pivot_longer(cols = -any_of(colnames(metadata)), 
               names_to = "taxon", values_to = "value") %>%
  mutate(across(c(exposure, final_disease_state, time), factor)) %>% #making these into factors
  nest_by(taxon)


#### What Changed ####
#to get this, read in the data w/o filtering for T3 and T7

all_timepoints <- cpm(otu_tmm, log = TRUE, prior.count = 2) %>%
  t %>%
  as_tibble(rownames = "sample_id") %>%
  full_join(metadata, by = "sample_id") %>%
  pivot_longer(cols = -any_of(colnames(metadata)), 
               names_to = "asv_names", values_to = "value") %>%
  mutate(across(c(exposure, final_disease_state), factor)) %>%
  mutate(time = readr::parse_number(time)) %>%
  left_join(taxonomy_tibble, by = join_by(asv_names))

homog_frags_asv <- all_timepoints %>%
  filter(tank == "homogenate_fragment")

homog_asv <- all_timepoints %>%
  filter(tank == "HOMO") %>%
  select(-disease_state)

field_frags_asv <- all_timepoints %>%
  filter(exposure == "Field")

clean_asv_data <- all_timepoints %>%
  filter(!tank %in% c("homogenate_fragment", "HOMO")) %>%
  filter(exposure != "Field")

#selecting significant ASVs

#get p values comparing D to H in homogenate to determine interesting ASVs
sig_homog_asv1 <- homog_asv %>% 
  nest_by(asv_names) %>%
  mutate(model = list(generics::tidy(anova(aov(value~final_disease_state, data = data))))) %>%
  rowwise %>%
  mutate(p_disease_state = get_tidy_p_values(model))

sig_homog_asv1_1 <- sig_homog_asv1 %>% #109 asvs
  filter(p_disease_state < 0.05)

list_of_asvs <- sig_homog_asv1_1$asv_names

#examine those ASVs in the T3-T7 data and keep ASVs w a significant interaction time*final_disease
homog_asv_changes <- clean_asv_data %>%
  filter(asv_names %in% list_of_asvs) %>%
  nest_by(asv_names) %>%
  mutate(model = list(generics::tidy(anova(lmer(value ~ time*final_disease_state - time 
                                                - final_disease_state+ (1 | fragment_id), 
                           data = data))))) %>%
  mutate(p_interaction = get_tidy_p_values(model)) %>%
  filter(p_interaction < 0.05)

sig_interaction_asvs <- homog_asv_changes$asv_names

#time series ASVs w sig interaction:
interaction_asv_changes <- clean_asv_data %>%
  filter(asv_names %in% sig_interaction_asvs)

#aov model for time series ASVs w sig interaction
homog_changes_aov <- lmer(value ~ time*final_disease_state*asv_names + (1 | fragment_id), 
                          data = interaction_asv_changes)
anova(homog_changes_aov)

#random effect of higher taxonomy types

emmeans(homog_changes_aov, ~time*final_disease_state | asv_names, type = 'response')%>%
  cld(Letters = LETTERS, adjust = 'fdr') %>%
  as_tibble() %>%
  mutate(.group = str_trim(.group)) %>%
  #rename(emmean = response) %>%
  ggplot(aes(x = asv_names, y = emmean, ymin = emmean - SE, ymax = emmean + SE,
             colour = time, pch = final_disease_state)) +
  geom_pointrange(position = position_dodge(0.5)) +
  geom_text(aes(y = (emmean + SE), label = .group),
            position = position_dodge(0.5), vjust = -1) +
  coord_flip()

# T7 values D > H

more_diseased <- interaction_asv_changes %>%
  filter(time == 7) #%>%
  #pivot_wider(names_from = final_disease_state, values_from = value) 

#TODO keep working on pivot
  


ave_homog_changes <- homog_asv_changes %>%
  group_by(asv_names, time, final_disease_state) %>%
  summarize(ave_abun = mean(value)) %>%
  ungroup()

ggplot(data = ave_homog_changes) +
  geom_line(aes(x = time, y = ave_abun, col = asv_names), alpha = 0.5) +
  geom_point(aes(x = time, y = ave_abun, col = asv_names), alpha = 0.5) +
  facet_wrap(~final_disease_state) +
  theme(legend.position = "none")

homog_changes_aov <- lmer(value ~ time*final_disease_state*asv_names + (1 | fragment_id), 
                          data = homog_asv_changes)

anova(homog_changes_aov)

emmeans(homog_changes_aov, ~time*final_disease_state*asv_names, type = 'response')%>%
  cld(Letters = LETTERS) %>%
  as_tibble() %>%
  mutate(.group = str_trim(.group)) %>%
  #rename(emmean = response) %>%
  ggplot(aes(x = asv_names, y = emmean, ymin = emmean - SE, ymax = emmean + SE,
             colour = time, pch = final_disease_state)) +
  geom_pointrange(position = position_dodge(0.5)) +
  geom_text(aes(y = (emmean + SE), label = .group),
            position = position_dodge(0.5), vjust = -1) +
  coord_flip()


#old method
top_homog_asv <- homog_asv %>%
  arrange(desc(ave_abundance)) %>%
  group_by(final_disease_state) %>%
  #filter(ave_abundance > 7.05) %>%
  slice(1:20) %>% #keep only n most abundant ASVs for healthy and diseased
  ungroup()

top_homog_list <- unique(dplyr::pull(top_homog_asv, asv_names))

ggplot(top_homog_asv) +
  geom_col(aes(x = fct_reorder(asv_names, ave_abundance), y = ave_abundance, 
               fill = final_disease_state), position = "dodge") +
  coord_flip()


homog_asv_changes <- clean_asv_data %>%
  filter(asv_names %in% top_homog_list) 

ave_homog_changes <- homog_asv_changes %>%
  group_by(asv_names, time, final_disease_state) %>%
  summarize(ave_abun = mean(value)) %>%
  ungroup()

ggplot(data = ave_homog_changes) +
  geom_line(aes(x = time, y = ave_abun, col = asv_names), alpha = 0.5) +
  geom_point(aes(x = time, y = ave_abun, col = asv_names), alpha = 0.5) +
  facet_wrap(~final_disease_state) +
  theme(legend.position = "none")


homog_changes_aov <- lmer(value ~ time*final_disease_state*asv_names + (1 | fragment_id), 
                          data = homog_asv_changes)
anova(homog_changes_aov)

emmeans(homog_changes_aov, ~time*final_disease_state | asv_names, type = 'link') %>%
  contrast('pairwise', adjust = "fdr")

emmeans(homog_changes_aov, ~final_disease_state * time | asv_names, type = 'link') %>%
  contrast('pairwise', adjust = "fdr")

tidy(anova(aov(y ~ x)))
tidy(t.test(y ~ x)) #apply the t test to the homogenates and THEN look at those ASVs

emmeans(homog_changes_aov, ~time*final_disease_state*asv_names, type = 'response') 

emmeans(homog_changes_aov, ~time*final_disease_state | asv_names, type = 'response')%>%
  cld(Letters = LETTERS, adjust = 'fdr') %>%
  as_tibble() %>%
  mutate(.group = str_trim(.group)) %>%
  #rename(emmean = response) %>%
  ggplot(aes(x = asv_names, y = emmean, ymin = emmean - SE, ymax = emmean + SE,
             colour = final_disease_state)) +
  geom_pointrange(position = position_dodge(0.5)) +
  geom_text(aes(y = (emmean + SE), label = .group),
            position = position_dodge(0.5), vjust = -1) +
  coord_flip()

emmeans(homog_changes_aov, ~time*final_disease_state*asv_names, type = 'response') %>%
  cld(Letters = LETTERS) %>%
  as_tibble() %>%
  mutate(.group = str_trim(.group)) %>%
  #rename(emmean = response) %>%
  ggplot(aes(x = asv_names, y = emmean, ymin = emmean - SE, ymax = emmean + SE,
             colour = time, pch = final_disease_state)) +
  geom_pointrange(position = position_dodge(0.5)) +
  geom_text(aes(y = (emmean + SE), label = .group),
            position = position_dodge(0.5), vjust = -1) +
  coord_flip()


#### Comparing Models ####
taxon_abundances$taxon
data <- taxon_abundances$data[[31]] #pick taxa of interest


norm_mod <- lmer(value ~ time * (exposure + final_disease_state) + (1 | fragment_id), data = data)
#lmer is for mixed-effect models, assumes that the residual error has a Gaussian distribution

rm_mod <- aov_4(value ~ time * (exposure + final_disease_state) + (time | fragment_id), data = data)
#repeated measures, only one random effect at a time

gamma_mod <- glmer(value ~ time * (exposure + final_disease_state) + (1 | fragment_id), data = data, 
                   family = Gamma(link = log))
#glmer allows you to predict response variables with non-Gaussian distributions


#graph to compare results from linear mixed-effects, repeated measures, and generalized lmer models
bind_rows(
  gamma = as_tibble(emmeans(gamma_mod, ~time:final_disease_state, type = 'response')) %>%
    rename(emmean = response,
           lower.CL = asymp.LCL,
           upper.CL = asymp.UCL) ,
  normal = as_tibble(emmeans(norm_mod, ~time:final_disease_state, type = 'response')),
  rm = as_tibble(emmeans(rm_mod, ~time:final_disease_state, type = 'response')),
  .id = 'model'
) %>%
  ggplot(aes(x = interaction(time, final_disease_state), y = emmean, ymin = lower.CL, 
             ymax = upper.CL, colour = model)) +
  geom_pointrange(position = position_dodge(0.5))
  
  
#### Model Each Taxon Independently ####
library(multidplyr)
cluster <- new_cluster(parallel::detectCores() - 1)
cluster_library(cluster, c('lme4', 'dplyr', 'tidyr', 'magrittr', 'stringr'))
cluster_copy(cluster, c('find_unique_significant_terms'))

all_models <- taxon_abundances %>%
  #partition(cluster) %>%
  mutate(model = list(lmer(value ~ time * (exposure + final_disease_state) + (1 | genotype) + 
                             (1 | tank), data = data))) #%>% #adds a column w anova results
  collect()

#for aggregated samples
aov_and_graphs <- all_models %>%
  ungroup %>%
  rowwise %>% #computes on a data frame one row at a time
  #partition(cluster) %>%
  mutate(terms = list(find_unique_significant_terms(model, 0.05))) %>%
  #collect %>%
  unnest(terms, keep_empty = TRUE) %>%
  rowwise %>%
  mutate(em_out = list(possibly(make_emmean_model, otherwise = NULL)(model, 
                                                                     as.formula(str_c('~', terms)), 
                                                                     0.05))) %>%
  mutate(plot = list(possibly(make_model_plot, otherwise = NULL, quiet = TRUE)
                     (em_out, data, terms))) %>%
  group_by(asv_names, data, model) %>%
  reframe(plot = ifelse(any(is.na(terms)), #if NAs in data, then no plot.  else, plot
                        list(NULL),
                        list(wrap_plots(plot) & 
                               labs(y = 'log2(CPM)') &
                               plot_annotation(title = asv_names) & 
                               theme_classic() &
                               theme(panel.background = element_rect(colour = 'black', fill = NA),
                                     axis.text = element_text(colour = 'black', size = 12),
                                     axis.title = element_text(colour = 'black', size = 16))))) %>%
  rowwise %>% 
  mutate(possibly(make_lmer_aov_summary, otherwise = NULL)(model)) %>%
  ungroup

#### Rickettsias - less good attempt ####

asv_level_aov_graphs <- all_models %>%
  ungroup %>%
  rowwise %>% #computes on a data frame one row at a time
  mutate(terms = list(find_unique_significant_terms(model, 0.05))) %>%
  unnest(terms, keep_empty = TRUE) %>%
  rowwise %>%
  mutate(em_out = list(possibly(make_emmean_model, otherwise = NULL)(model, 
                                                                     as.formula(str_c('~', terms)), 
                                                                     0.05))) %>%
  mutate(plot = list(possibly(make_model_plot, otherwise = NULL, quiet = TRUE)
                     (em_out, data, terms))) %>%
  mutate(plot = ifelse(any(is.na(terms)), #if NAs in data, then no plot.  else, plot
                        list(NULL),
                        list(wrap_plots(plot) & 
                               labs(y = 'log2(CPM)') &
                               plot_annotation(title = asv_names) & 
                               theme_classic() &
                               theme(panel.background = element_rect(colour = 'black', fill = NA),
                                     axis.text = element_text(colour = 'black', size = 12),
                                     axis.title = element_text(colour = 'black', size = 16))))) %>%
  rowwise %>% 
  mutate(possibly(make_lmer_aov_summary, otherwise = NULL)(model)) %>%
  ungroup %>%
  select(-c(em_out, terms))


#write_rds(asv_level_aov_graphs, "asv_aov_and_plots.rds")


asv_comp_upset <- asv_level_aov_graphs %>%
  select(-c(contains("df"), contains("chisq"), plot, data, model)) %>%
  full_join(taxonomy_tibble, by = join_by(asv_names)) %>%
  select(-c(Kingdom, Phylum)) %>%
  mutate(across(starts_with("p_"), ~p.adjust(.x, method = "fdr")))

p_vals_cu <- asv_comp_upset %>%
  select(c(asv_names, starts_with("p_"))) %>%
  mutate(across(starts_with("p_"), ~. < 0.05))




upset(mutate(asv_comp_upset, across(starts_with('p_'), ~. < 0.05)) %>%
        filter(!if_all(starts_with('p_'), ~!.)), 
      colnames(select(asv_comp_upset, starts_with('p_'))), 
      
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


#how to check for scenarios that fulfill specific combos of significant conditions

aov_and_graphs %>% 
  filter(`p_final_disease_state:time` < 0.05,
         p_exposure < 0.05,
         `p_exposure:time` > 0.05) %>%
  slice(2) %>%
  pull(model) %>%
  pluck(1)


aov_and_graphs %>%
  filter(str_detect(taxon, '[Gg]eitlerinema')) %>%
  pull(plot) %>%
  pluck(1)

ungroup(all_models) %>% slice(24) %>% pull(model) %>% pluck(1) %>% anova



left_join(asv_level_aov_graphs, taxonomy_tibble, by = join_by(asv_names)) %>%
  filter(Order == 'Rickettsiales') %>%
  sample_n(5) %>%
  pull(plot) %>%
  wrap_plots()


left_join(asv_level_aov_graphs, taxonomy_tibble, by = join_by(asv_names)) %>%
  filter(Order == 'Rickettsiales') %>%
  mutate(across(starts_with('p_'), ~. < 0.05)) %>%
  summarise(across(starts_with('p_'), sum),
            n = n())



tmp <- left_join(asv_level_aov_graphs, taxonomy_tibble, by = join_by(asv_names)) %>%
  filter(Order == 'Rickettsiales') %>%
  dplyr::select(asv_names, model) %>%
  rowwise(asv_names) %>%
  reframe(emmeans(model, ~final_disease_state) %>%
         as_tibble())


tmp %>%
  ggplot(aes(x = final_disease_state, y = log(emmean), group = asv_names)) +
  geom_path() +
  geom_point()


left_join(asv_level_aov_graphs, taxonomy_tibble, by = join_by(asv_names)) %>%
  filter(Order == 'Rickettsiales') %>%
  dplyr::select(asv_names, data) %>%
  unnest(data) %>%
  group_by(asv_names, final_disease_state) %>%
  summarise(value = mean(value),
            .groups = 'drop') %>%
  ggplot(aes(x = final_disease_state, y = value)) +
  geom_jitter()




rick_data <- left_join(asv_level_aov_graphs, taxonomy_tibble, by = join_by(asv_names)) %>%
  filter(Order == 'Rickettsiales') %>%
  dplyr::select(asv_names, data) %>%
  unnest(data) %>%
  left_join(taxonomy_tibble, by = join_by(asv_names))


rick_model <- lmer(value ~ time * (exposure + final_disease_state) + 
       (1 | genotype) + (1 | tank) + 
       (1 | Genus) + (1 | Family),
     data = rick_data)

summary(rick_model)
car::Anova(rick_model)

emmeans(rick_model, ~exposure) %>%
  as_tibble %>%
  ggplot(aes(x = exposure, y = emmean, ymin = emmean - SE, ymax = emmean + SE, colour = exposure)) +
  geom_pointrange(position = position_dodge(0.5))

emmeans(rick_model, ~final_disease_state) %>%
  contrast('pairwise')


rick_data %>%
  group_by(genotype, tank, time, exposure, final_disease_state) %>%
  summarise(value = sum(value)) %>%
  
  lmer(value ~ time * (exposure + final_disease_state) + 
         (1 | genotype) + (1 | tank),
       data = .) %>%
  summary()

#### Abundant Families ####

taxonomy_tibble1 <- tax_table(microbiome_data) %>% 
  as.data.frame %>%
  as_tibble(rownames = "family_names")

abundance_data <- cpm(otu_tmm, log = TRUE, prior.count = 2) %>%
  t %>%
  as_tibble(rownames = "sample_id") %>%
  full_join(metadata, by = "sample_id") %>%
  pivot_longer(cols = -any_of(colnames(metadata)), 
               names_to = "family_names", values_to = "value") %>%
  mutate(across(c(exposure, final_disease_state, time), factor)) %>%
  select(-c(sample_id, tank, genotype, disease_state)) %>%
  inner_join(taxonomy_tibble1, by = join_by(family_names)) %>%
  mutate(graph_category = paste(time, final_disease_state, sep = "_"))

abundance_data$graph_category <- factor(abundance_data$graph_category, levels = c("T3_H", "T7_H", "T3_D", "T7_D"))
  

abundance_data_summed <- abundance_data %>%
  group_by(graph_category, family_names) %>%
  summarize(abundance = sum(value)) %>%
  ungroup()

ranks <- abundance_data_summed %>%
  arrange(desc(abundance)) %>%
  group_by(graph_category) %>%
  slice(1:26)

ggplot(ranks, aes(graph_category, abundance, fill = family_names, label = family_names)) +
  geom_col(position = "fill") +
  geom_text(size = 3, position = position_fill(vjust = 0.5)) +
  theme_bw() +
  theme(legend.position = "none")
  

top20 <- unique(ranks$family_names)

abundance_data_summed$family_names <- ifelse(abundance_data_summed$family_names %in% top20, abundance_data_summed$family_names, "Other")

ggplot(abundance_data_summed) +
  geom_col(aes(graph_category, abundance, fill = family_names))
  







