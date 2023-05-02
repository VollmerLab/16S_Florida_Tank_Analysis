#Code that streamlines the process from Script #5: Disease-Associated Bacteria

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

select <- dplyr::select

#### Functions ####
find_unique_significant_terms_mixed <- function(model, alpha){
  significant_terms <- model$anova_table %>%
    as_tibble(rownames = 'param') %>%
    #janitor::clean_names() %>%
    filter(`Pr(>F)` < alpha) %>%
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

#otu_tmm
otu_tmm <- microbiome_data %>%
  phyloseq_filter_prevalence(prev.trh = 0.1) %>%
  otu_table() %>% 
  t %>% #NOTE: *genus and family do not need the t but ASVs need the t*
  as.data.frame %>%
  as.matrix %>% 
  DGEList(remove.zeros = TRUE) %>%
  edgeR::calcNormFactors(method = 'TMMwsp')

#raw model data - contains all 0 values
model_data <- cpm(otu_tmm, log = TRUE, prior.count = 2) %>%
  t %>%
  as_tibble(rownames = "sample_id") %>%
  full_join(metadata, by = "sample_id") %>%
  pivot_longer(cols = -any_of(colnames(metadata)), 
               names_to = "asv_names", values_to = "value") %>%
  mutate(across(c(exposure, final_disease_state), factor)) %>%
  filter(time %in% c('T3', 'T7') | (time == "T0" & tank == "HOMO")) %>%
  mutate(fragment_id = str_c(exposure, tank, genotype, final_disease_state)) %>%
  mutate(fragment_id = if_else(time == 'T0', 'homogenate', fragment_id))

#target microbiome data - zeroes removed and each ASV must be in 10+% of individuals
target_data <- cpm(otu_tmm, log = TRUE, prior.count = 2) %>%
  t %>%
  as_tibble(rownames = "sample_id") %>%
  full_join(metadata, by = "sample_id") %>%
  pivot_longer(cols = -any_of(colnames(metadata)), 
               names_to = "asv_names", values_to = "value") %>%
  mutate(across(c(exposure, final_disease_state), factor)) %>%
  filter(time %in% c('T3', 'T7') | (time == "T0" & tank == "HOMO")) %>%
  mutate(fragment_id = str_c(exposure, tank, genotype, final_disease_state)) %>%
  group_by(asv_names) %>%
  filter(n_distinct(sample_id[value > 7.04]) > 21) %>%
  ungroup()

#likely suspects upset prep
target_upset_data <- target_data %>%
  mutate(time = if_else(time == 'T0', exposure, time)) %>%
  group_by(time, asv_names) %>%
  summarise(n = sum(value > 7.04)) %>%
  pivot_wider(names_from = time, values_from = n, values_fill = 0L) %>%
  mutate(across(-asv_names, ~. > 0)) %>%
  left_join(taxonomy_tibble, by = "asv_names") %>%
  select(-c(Kingdom, Phylum)) %>%
  filter(D & T7 & T3)

taxonomy_tibble <- tax_table(microbiome_data) %>% 
  as.data.frame %>%
  as_tibble(rownames = "asv_names")


#### Rowwise Model ####

ls_model <- model_data %>% 
  filter(time != 'T0') %>%
  inner_join(target_upset_data) %>%
  nest_by(asv_names) %>%
  # mutate(model = list(mixed(value ~time * (final_disease_state + exposure) + 
  #                        (0 + dummy(time, c('T3', 'T7')) | fragment_id), data = data, method = 'KR',
  #                        control = variancePartition:::vpcontrol)))
  mutate(model = list(mixed(value ~time * (final_disease_state + exposure) + 
                              (1 | fragment_id), data = data, method = 'KR',
                            control = variancePartition:::vpcontrol)))
  

ls_model_w_terms <- ls_model %>%
  rowwise %>% 
  mutate(terms = list(find_unique_significant_terms_mixed(model, 0.05))) %>%
  rowwise(everything()) %>%
  summarise(enframe(terms) %>%
           mutate(name = TRUE) %>%
           pivot_wider(names_from = 'value',
                       values_from = 'name'),
           .groups = 'drop') %>%
  mutate(across(where(is.logical), ~if_else(is.na(.), FALSE, .)))

#### Prep for Comp Upsets ####

#process data for performing complex upset - LIKELY SUSPECTS
subset_asv_comp_upset <- ls_model %>% #reduces from 305 to 249 bc 56 are significant for nothing, 133 of these are ~time
  ungroup() %>%
  rowwise(asv_names) %>%
  reframe(as_tibble(model$anova_table, rownames = 'param'),
          d_v_h =  emmeans(model, ~final_disease_state) %>%
            as_tibble %>%
            select(emmean) %>%
            pull(1) %>%
            diff) %>%
  dplyr::rename(p = `Pr(>F)`) %>%
  group_by(param) %>%
  mutate(p = p.adjust(p, 'fdr')) %>%
  ungroup %>%
  mutate(p = p < 0.05) %>%
  select(asv_names, param, p, d_v_h) %>%
  pivot_wider(names_from = 'param', values_from = p, names_prefix = 'p_', values_fill = FALSE) %>%
  left_join(taxonomy_tibble, by = join_by(asv_names)) %>%
  select(-c(Kingdom, Phylum)) %>%
  #select(-p_time) %>% #TODO decide whether to include time in upsets here
  filter(!if_all(starts_with('p_'), ~!.))
  #filter(p_final_disease_state | `p_final_disease_state:time`)


#### Complex Upsets ####
  
#Visualizing Important Intersections in Venn Diagram
  upset(target_upset_data, 
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

  
#Likely Suspects, 249 total; 116 w/out time
upset(subset_asv_comp_upset,
        colnames(select(subset_asv_comp_upset, starts_with('p_'))), 
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
        name='asv_names', width_ratio=0.1, min_size = 0) +
  ggtitle("Likely Suspects")

upset(subset_asv_comp_upset,
                    colnames(select(subset_asv_comp_upset, starts_with('p_'))),
                    name='asv_names', width_ratio=0.1, min_size = 0) +
  ggtitle("Likely Suspects")
  
# All likely suspects - more abundant in Diseased
upset(filter(subset_asv_comp_upset, d_v_h < 0) %>%
        select(-d_v_h),  #negative = more in disease final state, positive = more in healthy disease state
      colnames(select(subset_asv_comp_upset, starts_with('p_'))), 
      
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
      
      name='asv_names', width_ratio=0.1, min_size = 1) +
  ggtitle("More in Diseased")


# All likely suspects - more abundant in Healthy
upset(filter(subset_asv_comp_upset, d_v_h > 0) %>%
        select(-d_v_h),  #negative = more in disease final state, positive = more in healthy disease state
      colnames(select(subset_asv_comp_upset, starts_with('p_'))), 
      
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
      
      name='asv_names', width_ratio=0.1, min_size = 1) +
  ggtitle("More in Healthy")

#### Graphing ####

#venn diagram
target_data %>%
  mutate(time = if_else(time == 'T0', exposure, time)) %>%
  # count(time, asv_names) %>%
  group_by(time, asv_names) %>%
  summarise(n = sum(value > 7.04)) %>%
  pivot_wider(names_from = time, values_from = n, values_fill = 0L) %>%
  mutate(across(-asv_names, ~. > 0)) %>%
  ggvenn(c('D', 'H', 'T3', 'T7'))



#FIXME combine these pipes, fix the graph
test <- ls_model_w_terms %>%
  filter(final_disease_state) %>%
  ungroup() %>%
  dplyr::slice(1:2) %>%
  rowwise() %>%
  mutate(test = list(make_emmean_model(model, as.formula("~final_disease_state"), 0.05)))

ttt <- test %>%
  group_by(asv_names) %>%
  reframe(test =  as_tibble(test, .name_repair = "unique")) %>%
  unnest(cols = c(test)) %>%
  unnest(cols = `...1`) %>%
  mutate(.group = str_trim(.group)) %>%
  ggplot(aes(x = asv_names, y = emmean, ymin = emmean - SE, ymax = emmean + SE,
             col = final_disease_state)) +
  geom_pointrange(position = position_dodge(0.5)) +
  geom_text(aes(y = (emmean + SE), label = .group),
            position = position_dodge(0.5), vjust = -1) +
  #scale_color_manual(values = c("hotpink1", "deepskyblue", "firebrick1", "dodgerblue3")) +
  coord_flip()





#### emmeans - currently bad ####

emm_dat <- model_data %>%
  inner_join(target_upset_data, by = c("asv_names")) %>%
  mutate(frag_ASV_id = paste(fragment_id, asv_names, sep = "_"))

#FIXME which data to include in this emmeans model

emm_aov <- aov_4(value ~time*final_disease_state*asv_names + (1 + time | frag_ASV_id), 
      data = emm_dat)

fds_only_list <- emmeans(emm_aov, ~final_disease_state | asv_names) %>%
  contrast('pairwise', adjust = 'fdr') %>%
  as_tibble() %>%
  filter(p.value < 0.05) %>%
  .$asv_names


emmeans(emm_aov, ~final_disease_state | asv_names, type = 'response') %>%
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




