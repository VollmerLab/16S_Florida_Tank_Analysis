#Code for making phylogenetic trees
setwd("~/Desktop/Screenshots/Career/Vollmer Lab/GitHub/16S_Florida_Tank_Analysis/Code")

#TODO phylolm

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
library(treedataverse)
library(multidplyr)
library(tidyverse)

select <- dplyr::select

#### Functions ####
make_tree <- function(asv_data, nboot = 100, make_plot = FALSE){
  prelim_tree <- asv_data %$%
    set_names(sequence, asv_names) %>%
    DNAStringSet() %>%
    msa() %>%
    msaConvert(., 'phangorn::phyDat') %>%
    modelTest() %>% 
    pml_bb()
  
  if(!(is.na(nboot) | is.null(nboot))){
    bootstrapping <- bootstrap.pml(prelim_tree, bs = nboot, optNni = TRUE,
                                   control = pml.control(trace = 1))
  } else {
    bootstrapping <- NULL
  }
  
  if(make_plot){
    if(!(is.na(nboot) | is.null(nboot))){
      plot(midpoint(prelim_tree$tree))
    } else {
      plotBS(midpoint(prelim_tree$tree), bootstrapping)
    }
  }
  return(list(prelim_tree = prelim_tree, bootstraps = bootstrapping))
}


make_tree_plot <- function(tree_output, taxonomy_tibble, fam_name){
  
  # plotBS(midpoint(tree_output$prelim_tree$tree), tree_output$bootstraps)
  transferBootstrap(tree = midpoint(tree_output$prelim_tree$tree), 
                    BStrees = tree_output$bootstraps) %>%
    as.treedata() %>%
    as_tibble %>%
    mutate(tip.label = str_extract(label, 'ASV_[0-9]+'),
           boot_support = str_extract(label, '^[0-9]+$') %>% as.integer) %>%
    #left_join(metadata,
    #          by = c('tip.label' = 'asv_names')) %>%
    # filter(!str_detect(label, 'ASV'))
    as.treedata() %>%
    #mutate(likely_color = ifelse(tip.label %in% asvs_opp, "opportunist", 
    #                             ifelse(tip.label %in% asvs_path, "pathogen", "other"))) %>%
    left_join(clean_bacstrats %>% filter(!bacstrat %in% c("late_crasher","continuous_crasher")), by = c('tip.label' = 'asv_id')) %>%
    mutate(bacstrat = ifelse(is.na(bacstrat), "other", bacstrat)) %>%
    left_join(taxonomy_tibble %>% select(asv_names, Family, Genus, Species), by = c('tip.label' = 'asv_names')) %>%
    ggtree(ladderize = TRUE, layout = 'rectangular') +
    geom_text(aes(label = ifelse(str_detect(tip.label, 'ASV'), paste(Genus, " (", parse_number(tip.label), ")", 
                                                                     sep = ""),tip.label), colour = bacstrat), hjust = 0) +
    geom_text(aes(label = boot_support), hjust = 0) + 
    #scale_color_manual(values = c("very likely" = "red", "likely" = "orange", "unlikely" = "black")) +
    scale_color_manual(values = c(
        "early_pathogen" = "hotpink",
        "late_opportunist" = "royalblue2",
        "late_probiotic" = "aquamarine3",
        "early_opportunist" = "deepskyblue3",
        "late_pathogen" = "firebrick1",
        "continuous_pathogen" = "firebrick4",
        "other" = "black"
        )) +
    theme(legend.position='top', 
          legend.justification='left',
          legend.direction='horizontal') +
    ggtitle(fam_name)
}

#### Read in Data ####

#only ASVs we care about - need to fix implementation
relevant_asvs <- read_csv('../intermediate_files/nonvenn_filtered_fully_preprocessed_samples.csv.gz') %>%
  pull(asv_id) %>% unique()

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
  sequences <- taxa_names(microbiome_data)
  taxa_names(microbiome_data) <- str_c('ASV', 1:length(taxa_names(microbiome_data)), sep = '_')
  names(sequences) <- taxa_names(microbiome_data)
}

otu_tmm <- microbiome_data %>%
  phyloseq_filter_prevalence(prev.trh = 0.1) %>%
  otu_table() %>% 
  t %>% #NOTE: *genus and family do not need the t but ASVs need the t*
  as.data.frame %>%
  as.matrix %>% 
  DGEList(remove.zeros = TRUE) %>%
  edgeR::calcNormFactors(method = 'TMMwsp')

#raw model data - contains all 0 values
full_metadata <- cpm(otu_tmm, log = TRUE, prior.count = 2) %>%
  t %>%
  as_tibble(rownames = "sample_id") %>%
  full_join(metadata, by = "sample_id") %>%
  pivot_longer(cols = -any_of(colnames(metadata)), 
               names_to = "asv_names", values_to = "value") %>%
  mutate(across(c(exposure, final_disease_state), factor)) %>%
  filter(time %in% c('T3', 'T7') | (time == "T0" & tank == "HOMO")) %>%
  select(-c(value, clone_group, resistance, susceptability, reads))

taxonomy_tibble <- tax_table(microbiome_data) %>% 
  as.data.frame %>%
  as_tibble(rownames = "asv_names")



#### Trees ####

kept_asvs <- full_data %>% pull(asv_id) %>% unique()

important_asvs <- clean_bacstrats %>% pull(asv_id) %>% unique()

important_families <- clean_bacstrats %>% 
  left_join(taxonomy_tibble, by = join_by("asv_id" == "asv_names")) %>% 
  pull(Family) %>% unique()


#format data for becoming a tree
current_fam <- enframe(sequences, name = 'asv_names', value = 'sequence') %>%
  left_join(taxonomy_tibble, by = 'asv_names') %>%
  filter(asv_names %in% kept_asvs) %>% #only make trees w ASVs that passed filtering
  rowwise() %>%
  nest(data = c(Genus, Species, asv_names, sequence)) %>%
  filter(Family %in% important_families) %>%
  rowwise %>%
  mutate(n_asv = nrow(data)) %>%
  ungroup %>%
  filter(n_asv > 3)

current_fam <- current_fam %>% filter(Family %in% c("Francisellaceae"))

  # cluster <- new_cluster(parallel::detectCores() - 1)
  # cluster_library(cluster, c('dplyr', 'msa', 'Biostrings', 'ape', 'phangorn', 'magrittr', 'treedataverse'))
  # cluster_copy(cluster, c('make_tree', 'make_tree_plot'))

plot_list <- current_fam %>%
  unnest(data) %>%
  left_join(full_metadata, by = join_by(asv_names), multiple = "all") %>%
  mutate(asv_names = paste(asv_names, " (", Genus, " ", Species, ")", sep = "")) %>%
  rowwise() %>%
  nest(data = c(Genus, Species, asv_names, sequence)) %>%
  filter(!is.na(final_disease_state)) %>% ### the metadata for the ones w all NAs were removed by prev filter
  rowwise() %>%
  nest(metadata = c(sample_id, time, exposure, tank, genotype, final_disease_state)) %>%
  #partition(cluster) %>%
  rowwise() %>%
  mutate(forest = list(make_tree(data, nboot = 100, make_plot = FALSE)),
         tree_plot = list(possibly(make_tree_plot, otherwise = NULL)(forest, taxonomy_tibble, Family))) %>%
  #collect %>%
  identity()


#write_rds(plot_list, "../intermediate_files/family_trees.rds")

#### Jason's tutorial ####

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


#### ####



