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
library(wesanderson)
library(tidytext)
library(strex)
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

#EXAMPLE of emmeans within ASVs

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

#### Rickettsias - Jason's graphs ####

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
  filter(time %in% c("T3", "T7") | (time == "T0" & tank == "HOMO")) %>%
  select(-c(sample_id, tank, genotype, disease_state)) %>%
  inner_join(taxonomy_tibble1, by = join_by(family_names)) %>%
  mutate(graph_category = paste(time, final_disease_state, sep = "_"))

abundance_data$graph_category <- factor(abundance_data$graph_category, 
                                        levels = c("T0_H", "T3_H", "T7_H", "T0_D", "T3_D", "T7_D"))
  

abundance_data_summed <- abundance_data %>%
  group_by(graph_category, family_names, final_disease_state) %>%
  summarize(abundance = sum(value)) %>%
  ungroup()

ranks <- abundance_data_summed %>%
  arrange(desc(abundance)) %>%
  group_by(graph_category) %>%
  slice(1:26)

total_ranks <- ranks %>%
  summarize(totals = sum(abundance))

ranks <- ranks %>%
  full_join(total_ranks, by = join_by(graph_category)) %>%
  mutate(rel_abund = 100*abundance/totals) %>%
  mutate(family_names = factor(family_names))
  
ranks %>%
  mutate(family_sort = fct_reorder(reorder_within(graph_category, family_names, rel_abund), rel_abund, .desc = TRUE)) 

num_of_fams <- length(unique(ranks$family_names))
getPalette <- colorRampPalette(brewer.pal(8, "Viridis"))

library(RColorBrewer)
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,num_of_fams), col=sample(col_vector, num_of_fams))

set.seed(1245)
ggplot(ranks, aes(graph_category, abundance, group = reorder_within(family_names, rel_abund, graph_category) %>% fct_rev(), 
                  label = paste(family_names, paste(round(rel_abund, digits = 2), "%", sep = ""), sep = " - "),
                  fill = family_names)) +
  geom_col(position = "fill", col = "black") +
  geom_text(size = 3, position = position_fill(vjust = 0.5)) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(title = "Relative Abundances of 26 Most Abundant Families Per Sample") +
  ylab("Relative Abundance") +
  xlab("Sample") +
  scale_x_discrete(labels = c("Healthy Homogenate", "Healthy T3", "Healthy T7", 
                              "Diseased Homogenate", "Diseased T3", "Diseased T7")) +
  scale_fill_manual(values = sample(col_vector, num_of_fams))




#fill = fct_reorder(reorder_within(graph_category, family_names, rel_abund), rel_abund, .desc = TRUE)


#bait abundances

aggregation_level <- 'none' #or none

microbiome_data <- read_rds("../intermediate_files/preprocess_microbiome.rds") %>%
  subset_samples(time %in% 'T0')
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

#species

otu_tmm <- microbiome_data %>%
  phyloseq_filter_prevalence(prev.trh = 0.1) %>%
  otu_table() %>% 
  t %>% #NOTE: *genus and family do not need the t but ASVs need the t*
  as.data.frame %>%
  as.matrix %>% 
  DGEList(remove.zeros = TRUE) %>%
  edgeR::calcNormFactors(method = 'TMMwsp') #TMMwsp is for high prevalence of 0s

bait_data <- cpm(otu_tmm, log = TRUE, prior.count = 2) %>%
  t %>%
  as_tibble(rownames = "sample_id") %>%
  full_join(metadata, by = "sample_id") %>%
  pivot_longer(cols = -any_of(colnames(metadata)), 
               names_to = "asv_names", values_to = "value") %>%
  mutate(across(c(exposure, final_disease_state, time), factor)) %>%
  filter(tank == "HOMO") %>%
  select(-c(disease_state, time, sample_id, tank)) %>%
  group_by(final_disease_state, asv_names) %>%
  summarize(abundance = sum(value)) %>%
  ungroup() %>%
  left_join(taxonomy_tibble, by = join_by(asv_names)) %>%
  arrange(desc(abundance)) %>%
  group_by(final_disease_state) %>%
  slice(1:30)

bait_totals <- bait_data %>%
  summarize(totals = sum(abundance))

family_totals <- bait_data %>%
  group_by(final_disease_state, Family) %>%
  summarize(fam_totals = sum(abundance)) %>%
  arrange(desc(fam_totals))

all_bait_data <- bait_data %>%
  full_join(bait_totals, by = join_by(final_disease_state)) %>%
  full_join(family_totals, by = join_by(final_disease_state, Family)) %>%
  ungroup() %>%
  group_by(final_disease_state) %>%
  mutate(fam_order = 10*rank(desc(fam_totals))) %>%
  group_by(final_disease_state, Family) %>%
  mutate(spec_order = rank(desc(abundance))) %>%
  mutate(test = as.character(log(spec_order+1)/3)) %>%
  mutate(test1 = str_after_last(test, "\\.")) %>%
  mutate(overall_order = as.numeric(paste(fam_order, test1, sep = "."))) %>%
  mutate(rel_abund = 100*abundance/totals)

num_of_asvs <- length(unique(all_bait_data$asv_names))

set.seed(032897)
ggplot(all_bait_data, aes(final_disease_state, abundance, fill = asv_names, group = reorder_within(asv_names, rel_abund, final_disease_state) %>% fct_rev(), 
        label = paste(Family, Genus, Species, paste("(", asv_names, ")", " - ", round(rel_abund, digits = 2),"%", 
                                                    sep = ""),  sep = " "))) +
  geom_col(position = "fill", col = "black") +
  geom_text(size = 3, position = position_fill(vjust = 0.5)) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(title = "30 Most Abundant Species") +
  scale_fill_manual(values = sample(col_vector, num_of_asvs))


#families

otu_tmm_f <- microbiome_data %>%
  phyloseq_filter_prevalence(prev.trh = 0.1) %>%
  otu_table() %>% 
  #t %>% #NOTE: *genus and family do not need the t but ASVs need the t*
  as.data.frame %>%
  as.matrix %>% 
  DGEList(remove.zeros = TRUE) %>%
  edgeR::calcNormFactors(method = 'TMMwsp') #TMMwsp is for high prevalence of 0s


bait_data_f <- cpm(otu_tmm_f, log = TRUE, prior.count = 2) %>%
  t %>%
  as_tibble(rownames = "sample_id") %>%
  full_join(metadata, by = "sample_id") %>%
  pivot_longer(cols = -any_of(colnames(metadata)), 
               names_to = "family_names", values_to = "value") %>%
  mutate(across(c(exposure, final_disease_state, time), factor)) %>%
  filter(tank == "HOMO") %>%
  select(-c(disease_state, time, sample_id, tank)) %>%
  group_by(final_disease_state, family_names) %>%
  summarize(abundance = sum(value)) %>%
  ungroup() %>%
  left_join(taxonomy_tibble1, by = join_by(family_names)) %>%
  arrange(desc(abundance)) %>%
  group_by(final_disease_state) %>%
  slice(1:30)

bait_data_f$final_disease_state <- factor(bait_data_f$final_disease_state, levels = c("H", "D"))


ggplot(bait_data_f, aes(final_disease_state, abundance, fill = fct_reorder(family_names, abundance, .desc = TRUE), 
                        label = paste(family_names, round(abundance, digits = 1), sep = " - "))) +
  geom_col(position = "fill") +
  geom_text(size = 3, position = position_fill(vjust = 0.5)) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(title = "30 Most Abundant Families in Baits")

  
### Complex Upset for Baits

cu_bait_data <- cpm(otu_tmm, log = TRUE, prior.count = 2) %>%
  t %>%
  as_tibble(rownames = "sample_id") %>%
  full_join(metadata, by = "sample_id") %>%
  pivot_longer(cols = -any_of(colnames(metadata)), 
               names_to = "asv_names", values_to = "value") %>%
  mutate(across(c(exposure, final_disease_state, time), factor)) %>%
  filter(tank == "HOMO") %>%
  select(-c(disease_state, time, sample_id, tank)) %>%
  group_by(asv_names, final_disease_state) %>%
  summarize(ave_val = mean(value)) %>%
  rowwise() %>%
  mutate(present = ifelse(ave_val > 7.04, TRUE, FALSE)) %>%
  select(-ave_val) %>%
  pivot_wider(names_from = final_disease_state, values_from = present)

disease_states <- c("H", "D")

upset(cu_bait_data, disease_states, name="Bait Type", width_ratio=0.1) +
  labs(title = "ASVs in Baits")

  
#### Logfold Changes ####
"red" #get otu_tmm w/o filtering for T3 and T7 for this chunk:

log2_data <- cpm(otu_tmm, log = TRUE, prior.count = 2) %>%
  t %>%
  as_tibble(rownames = "sample_id") %>%
  full_join(metadata, by = "sample_id") %>%
  pivot_longer(cols = -any_of(colnames(metadata)), 
               names_to = "asv_names", values_to = "value") %>%
  mutate(across(c(exposure, final_disease_state, time), factor)) %>%
  #mutate(value = log2(value)) %>%
  filter(time %in% c("T3", "T7") | (time == "T0" & tank == "HOMO")) %>%
  group_by(time, final_disease_state, asv_names) %>%
  summarize(ave_val = mean(value)) %>%
  ungroup() %>%
  nest_by(time, asv_names) %>%
  rowwise() %>%
  filter(nrow(data) > 1) %>%
  unnest(cols = c(data)) %>%
  group_by(time, asv_names) %>%  
  summarize(logfold = ave_val[final_disease_state == "D"] - ave_val[final_disease_state == "H"])

from_rds <- read_rds("likely_suspects_list.rds")
likely_suspects_list <- from_rds[[1]]
v_likely_suspects_list <- from_rds[[2]]

ls_log2_data <- log2_data %>%
  filter(asv_names %in% likely_suspects_list)

ggplot(ls_log2_data) +
  geom_hline(yintercept = 0, col = "black") +
  geom_point(aes(x = asv_names, y = logfold, col = time)) +
  coord_flip() +
  theme_bw() +
  labs(title = "Likely Suspects Logfold Change")

vls_log2_data <- log2_data %>%
  filter(asv_names %in% v_likely_suspects_list)

ggplot(vls_log2_data) +
  geom_hline(yintercept = 0, col = "black") +
  geom_point(aes(x = asv_names, y = logfold, col = time)) +
  coord_flip() +
  theme_bw() +
  labs(title = "Very Likely Suspects Logfold Change")


## boxplot version

log2_data_eb <- cpm(otu_tmm, log = TRUE, prior.count = 2) %>%
  t %>%
  as_tibble(rownames = "sample_id") %>%
  full_join(metadata, by = "sample_id") %>%
  pivot_longer(cols = -any_of(colnames(metadata)), 
               names_to = "asv_names", values_to = "value") %>%
  mutate(across(c(exposure, final_disease_state, time), factor)) %>%
  mutate(value = log2(value)) %>%
  filter(time %in% c("T3", "T7") | (time == "T0" & tank == "HOMO")) %>%
  group_by(time, asv_names, genotype) %>%
  summarize(logfold = value[final_disease_state == "D"] - value[final_disease_state == "H"])

#likely suspects
ls_log2_data_eb <- log2_data_eb %>%
  filter(asv_names %in% likely_suspects_list)

ggplot(ls_log2_data_eb) +
  geom_hline(yintercept = 0, col = "black") +
  geom_boxplot(aes(x = asv_names, y = logfold, fill = time, col = time)) +
  coord_flip() +
  theme_bw() +
  scale_fill_manual(values = c("indianred1", "seagreen2", "deepskyblue1")) +
  scale_color_manual(values = c("indianred4", "seagreen4", "royalblue2")) +
  labs(title = "Likely Suspects Logfold Change")

#very likely suspects
vls_log2_data_eb <- log2_data_eb %>%
  filter(asv_names %in% v_likely_suspects_list)

ggplot(vls_log2_data_eb) +
  geom_hline(yintercept = 0, col = "black") +
  geom_boxplot(aes(x = asv_names, y = logfold, fill = time, col = time)) +
  coord_flip() +
  theme_bw() +
  scale_fill_manual(values = c("indianred1", "seagreen2", "deepskyblue1")) +
  scale_color_manual(values = c("indianred4", "seagreen4", "royalblue2")) +
  labs(title = "Very Likely Suspects Logfold Change")

### logfold emmeans

log_emmeans_data <- cpm(otu_tmm, log = TRUE, prior.count = 2) %>%
  t %>%
  as_tibble(rownames = "sample_id") %>%
  full_join(metadata, by = "sample_id") %>%
  pivot_longer(cols = -any_of(colnames(metadata)), 
               names_to = "asv_names", values_to = "value") %>%
  filter(time %in% c("T3", "T7") | (time == "T0" & tank == "HOMO")) %>%
  mutate(across(c(exposure, final_disease_state, time), factor)) %>%
  select(value, time, final_disease_state, asv_names, tank, genotype) %>%
  filter(asv_names %in% likely_suspects_list) %>%
  mutate(log_value = log(value))

log_emmeans_aov <- lmer(log_value ~time*final_disease_state*asv_names + (1 | tank) + (1 | genotype), 
                        data = log_emmeans_data)

logfold_emmeans_contrasts <- emmeans(log_emmeans_aov, ~final_disease_state | asv_names*time) %>%
  contrast('pairwise', adjust = 'fdr')

logfold_emmeans_contrasts %>%
  as_tibble() %>%
  filter(time == "T0") %>%
  filter(estimate > 0) %>%
  {.$asv_names ->> more_in_disease_logfold}


logfold_emmeans_contrasts %>%
  as_tibble() %>%
  filter(asv_names %in% more_in_disease_logfold) %>%
  group_by(asv_names) %>%
  mutate(whats_more = ifelse(estimate[time == "T7"] > estimate[time == "T3"], "More in T7", "More in T3")) %>%
  filter(estimate[time == "T7"] > 0 | estimate[time == "T3"] > 0) %>%
  mutate(growth = estimate[time == "T7"] - estimate[time == "T3"]) %>%
  arrange(asv_names) %>%
  ungroup() %>%
  mutate(alpha_val = ifelse(time == "T0", "less","more")) %>%
  left_join(taxonomy_tibble, by = join_by(asv_names)) %>%
  ggplot(aes(x = fct_reorder(paste(Family, " ", Genus, " (", asv_names, ")", sep = ""), growth), y = estimate, 
             ymin = estimate - SE, ymax = estimate + SE, colour = time, pch = time)) +
  geom_hline(yintercept = 0) +
  geom_pointrange(position = position_dodge(0.5)) +
  labs(title = "Likely Suspects Logfold Change") +
  coord_flip() +
  facet_grid(rows = vars(whats_more), scales = "free", space = "free") +
  scale_color_manual(values = c("goldenrod1", "darkorange1", "firebrick1")) +
  xlab("ASV Name")
  #scale_alpha_discrete(range = c(0.7,1))
  
#VLS

vls_log_emmeans_data <- log_emmeans_data %>%
  filter(asv_names %in% v_likely_suspects_list)

vls_log_emmeans_aov <- lmer(log_value ~time*final_disease_state*asv_names + (1 | tank) + (1 | genotype), 
                        data = vls_log_emmeans_data)

logfold_vls_emmeans_contrasts <- emmeans(vls_log_emmeans_aov, ~final_disease_state | asv_names*time) %>%
  contrast('pairwise', adjust = 'fdr')

logfold_vls_emmeans_contrasts %>%
  as_tibble() %>%
  filter(time == "T0") %>%
  filter(estimate > 0) %>%
  {.$asv_names ->> vls_more_in_disease_logfold}


logfold_vls_emmeans_contrasts %>%
  as_tibble() %>%
  filter(asv_names %in% vls_more_in_disease_logfold) %>%
  group_by(asv_names) %>%
  mutate(whats_more = ifelse(estimate[time == "T7"] > estimate[time == "T3"], "More in T7", "More in T3")) %>%
  filter(estimate[time == "T7"] > 0 | estimate[time == "T3"] > 0) %>%
  mutate(growth = estimate[time == "T7"] - estimate[time == "T3"]) %>%
  arrange(asv_names) %>%
  ungroup() %>%
  mutate(alpha_val = ifelse(time == "T0", "less","more")) %>%
  left_join(taxonomy_tibble, by = join_by(asv_names)) %>%
  ggplot(aes(x = fct_reorder(paste(Family, " ", Genus, " (", asv_names, ")", sep = ""), growth), y = estimate, 
             ymin = estimate - SE, ymax = estimate + SE, colour = time, pch = time)) +
  geom_hline(yintercept = 0) +
  geom_pointrange(position = position_dodge(0.5)) +
  labs(title = "Very Likely Suspects Logfold Change") +
  coord_flip() +
  facet_grid(rows = vars(whats_more), scales = "free", space = "free") +
  scale_color_manual(values = c("firebrick1", "goldenrod1", "darkorange1")) +
  xlab("ASV Name")

#### several models ####

log_emmeans_data_37 <- cpm(otu_tmm, log = TRUE, prior.count = 2) %>%
  t %>%
  as_tibble(rownames = "sample_id") %>%
  full_join(metadata, by = "sample_id") %>%
  pivot_longer(cols = -any_of(colnames(metadata)), 
               names_to = "asv_names", values_to = "value") %>%
  filter(time %in% c("T3", "T7")) %>%
  mutate(across(c(exposure, final_disease_state, time), factor)) %>%
  filter(asv_names %in% likely_suspects_list) %>% #v_likely_suspects_list or likely_suspects_list
  mutate(log_value = log(value)) %>%
  mutate(frag_ASV_id = paste(fragment_id, asv_names, sep = "_")) %>%
  select(-c(sample_id, disease_state))

log_emmeans_data_0 <- cpm(otu_tmm, log = TRUE, prior.count = 2) %>%
  t %>%
  as_tibble(rownames = "sample_id") %>%
  full_join(metadata, by = "sample_id") %>%
  pivot_longer(cols = -any_of(colnames(metadata)), 
               names_to = "asv_names", values_to = "value") %>%
  filter(time == "T0" & tank == "HOMO") %>%
  mutate(across(c(exposure, final_disease_state, time), factor)) %>%
  select(value, time, final_disease_state, asv_names, tank, genotype) %>%
  filter(asv_names %in% likely_suspects_list) %>% #v_likely_suspects_list or likely_suspects_list
  mutate(log_value = log(value)) %>%
  mutate(genotype = ifelse(final_disease_state == "H", as.numeric(genotype) +5, as.numeric(genotype)))

log_emmeans_aov_0 <- lmer(log_value ~final_disease_state*asv_names + (1 | genotype), 
                        data = log_emmeans_data_0)

logfold_emmeans_contrasts_0 <- emmeans(log_emmeans_aov_0, ~final_disease_state | asv_names) %>%
  contrast('pairwise', adjust = 'fdr')

log_emmeans_aov_37 <- aov_4(log_value ~time*final_disease_state*asv_names + (1 + time | frag_ASV_id), 
                          data = log_emmeans_data_37)

logfold_emmeans_contrasts_37 <- emmeans(log_emmeans_aov_37, ~final_disease_state | asv_names*time) %>%
  contrast('pairwise', adjust = 'fdr')

logfold_emmeans_contrasts_0 %>%
  as_tibble() %>%
  filter(estimate > 0) %>%
  {.$asv_names ->> more_in_disease_logfold_0}


logfold_emmeans_contrasts_37 %>%
  as_tibble() %>%
  rbind(as_tibble(logfold_emmeans_contrasts_0) %>% mutate(time = "T0")) %>%
  filter(asv_names %in% more_in_disease_logfold_0) %>%
  group_by(asv_names) %>%
  mutate(whats_more = ifelse(estimate[time == "T7"] > estimate[time == "T3"], "More in T7", "More in T3")) %>%
  filter(estimate[time == "T7"] > 0 | estimate[time == "T3"] > 0) %>%
  mutate(growth = estimate[time == "T7"] - estimate[time == "T3"]) %>%
  arrange(asv_names) %>%
  ungroup() %>%
  mutate(time = factor(time, levels = c("T0", "T3", "T7"))) %>%
  mutate(alpha_val = ifelse(time == "T0", "less","more")) %>%
  left_join(interaction_types, by = join_by(asv_names), multiple = "all") %>%
  left_join(taxonomy_tibble, by = join_by(asv_names)) %>%
  ggplot(aes(x = fct_reorder(paste(Family, " ", Genus, " (", asv_names, ")", sep = ""), growth), y = estimate, 
             ymin = estimate - SE, ymax = estimate + SE, colour = time, pch = time, alpha = whats_more)) +
  geom_hline(yintercept = 0) +
  geom_pointrange(position = position_dodge(0.5)) +
  labs(title = "Likely Suspects Logfold Change - Combined Models") +
  coord_flip() +
  facet_grid(rows = vars(type), scales = "free", space = "free") +
  scale_color_manual(values = c("goldenrod1", "darkorange1", "firebrick1")) +
  xlab("ASV Name") +
  scale_alpha_discrete(range = c(0.4, 1), guide = "none")


## TEST

log2_data$value <- log(log2_data$value)

test <- log2_data %>% left_join(taxonomy_tibble) %>% filter(Genus %in% "Pseudoalteromonas")

test_pseudo_aov <- lmer(value ~time*final_disease_state*asv_names + (1 | tank) + (1 | genotype), 
                          data = test)

test_pseudo <- emmeans(test_pseudo_aov, ~final_disease_state | time*asv_names) %>%
  contrast('pairwise', adjust = 'fdr')

test_pseudo %>%
  as_tibble() %>%
  left_join(taxonomy_tibble, by = join_by(asv_names)) %>%
  group_by(asv_names) %>%
  mutate(whats_more = ifelse(estimate[time == "T7"] > estimate[time == "T3"], "More in T7", "More in T3")) %>%
  mutate(estimate_val = estimate[time == "T0"]) %>%
  arrange(asv_names) %>%
  ungroup() %>%
  ggplot(aes(x = fct_reorder(paste(Genus, " ", Species, "(", asv_names, ")", sep = ""), estimate_val), y = estimate, ymin = estimate - SE, ymax = estimate + SE,
             colour = time, pch = time)) +
  geom_hline(yintercept = 0) +
  geom_pointrange(position = position_dodge(0.5)) +
  coord_flip() +
  facet_grid(rows = vars(whats_more), scales = "free", space = "free") +
  scale_color_manual(values = c("firebrick1", "goldenrod1", "darkorange1")) +
  labs(title = "Pseudoalteromonas") +
  xlab("ASV")

#### Tank Effect ####

tank_data <- cpm(otu_tmm, log = TRUE, prior.count = 2) %>%
  t %>%
  as_tibble(rownames = "sample_id") %>%
  full_join(metadata, by = "sample_id") %>%
  pivot_longer(cols = -any_of(colnames(metadata)), 
               names_to = "asv_names", values_to = "value") %>%
  filter(time %in% c("T3", "T7")) %>%
  mutate(across(c(exposure, final_disease_state, time), factor))

model1 <- lmer(value ~ final_disease_state + (1 | time) + (1 | tank) + (1 | asv_names), data = tank_data)
    #w/o ASVs -  tank = 0.01625, time = 0.05286 
    #w/ ASVs - asv_names = 0.6951, tank = 0.0177, time = 0.0529 

tank_model <- lmer(value ~ time + (1 | tank), 
              data = filter(tank_data, final_disease_state == "D"))
#H: ~ time, tank is 0.0287, ~ tank, time is 0.08549
#D: ~ time, tank is 0.006784, ~ tank, time is 0.002925


tank_exposure_model <- lmer(value ~ final_disease_state + (1 | tank) + (1 | time) + (1 | asv_names), 
                   data = filter(tank_data, exposure == "H"))
#H: asv_names - 0.71949, tank - 0.02200, time - 0.07243 
#D: asv_names - 0.744188, tank - 0.003628, time - 0.034126

emmeans(tank_exposure_model, ~final_disease_state, type = 'response') %>%
  cld(Letters = LETTERS, adjust = 'fdr') %>%
  as_tibble() %>%
  mutate(.group = str_trim(.group)) %>%
  #rename(emmean = response) %>%
  ggplot(aes(x = final_disease_state, y = emmean, ymin = emmean - SE, ymax = emmean + SE,
             colour = final_disease_state)) +
  geom_pointrange(position = position_dodge(0.5)) +
  geom_text(aes(y = (emmean + SE), label = .group),
            position = position_dodge(0.5), vjust = -1) +
  labs(title = "Tank Effect Model")

