#Code for making phylogenetic trees
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
  sequences <- taxa_names(microbiome_data)
  taxa_names(microbiome_data) <- str_c('ASV', 1:length(taxa_names(microbiome_data)), sep = '_')
  names(sequences) <- taxa_names(microbiome_data)
}



#### Trees ####

"red"
#TODO work in progress, adjust based on new model


# interesting families
likely_families <- likely_suspects_list %>% 
  as_tibble() %>% 
  reframe(asv_names = value) %>% 
  left_join(taxonomy_tibble, by = join_by(asv_names)) %>%
  .$Family %>%
  unique() %>%
  sort()

very_likely_families <- vl_suspects_list %>% 
  as_tibble() %>% 
  reframe(asv_names = value) %>% 
  left_join(taxonomy_tibble, by = join_by(asv_names)) %>%
  .$Family %>%
  unique() %>%
  sort()

likely_families <- filter(subset_asv_comp_upset, d_v_h < 0) %>% .$Family %>% unique()


#tree stuff
tmp <- enframe(sequences, name = 'asv_names', value = 'sequence') %>%
  left_join(taxonomy_tibble, by = 'asv_names') %>%
  rowwise() %>%
  #mutate(asv_names = ifelse(asv_names %in% vl_suspects_list, paste("**", asv_names, sep = ""), asv_names)) %>%
  #mutate(asv_names = ifelse(asv_names %in% likely_suspects_list, paste("*", asv_names, sep = ""), asv_names)) %>%
  nest(data = c(Genus, Species, asv_names, sequence)) %>%
  filter(Family %in% likely_families) %>%
  rowwise %>%
  mutate(n_asv = nrow(data)) %>%
  ungroup
  #mutate(likelihood = ifelse(Family %in% very_likely_families, "very likely", "likely")) 

tmp %>%
  arrange(n_asv)
current_fam <- tmp %>% filter(Order %in% c("Campylobacterales", "Myxococcales")) %>%
  filter(n_asv < 30)

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
# asv_data <- current_fam$data[[1]]
# make_tree(current_fam$data[[1]], nboot = NA, make_plot = TRUE)

make_tree_plot <- function(tree_output, metadata){
  plotBS(midpoint(tree_output$prelim_tree$tree), tree_output$bootstraps) %>%
    as.treedata() %>%
    as_tibble %>%
    mutate(tip.label = str_extract(label, 'ASV_[0-9]+'),
           boot_support = str_extract(label, '^[0-9]+$') %>% as.integer) %>%
    left_join(metadata,
              by = c('tip.label' = 'ASV_name')) %>%
    # filter(!str_detect(label, 'ASV'))
    as.treedata() %>%
    ggtree(ladderize = TRUE, layout = 'rectangular') +
    geom_text(aes(label = tip.label, colour = INTERESTING), hjust = 0) +
    geom_text(aes(label = boot_support), hjust = 0)
}

library(treedataverse)
library(multidplyr)
cluster <- new_cluster(parallel::detectCores() - 1)
cluster_library(cluster, c('dplyr', 'msa', 'Biostrings', 'ape', 'phangorn', 'magrittr', 'treedataverse'))
cluster_copy(cluster, c('make_tree', 'make_tree_plot'))

plot_list <- current_fam %>%
  left_join(metadata = ) %>%
  rowwise() %>%
  # partition(cluster) %>%
  mutate(forest = list(make_tree(data, nboot = 10, make_plot = FALSE)),
         tree_plot = list(make_tree_plot(forest, metadata))) %>%
  # collect %>%
  identity()




tree_output<- plot_list$forest[[1]]


plot_list$forest[[1]]$prelim_tree$tree %>%
  as.treedata() %>% 
  # as_tibble
  ggtree() +
  geom_text(aes(x = branch, label = label))

#Jason's tutorial

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

