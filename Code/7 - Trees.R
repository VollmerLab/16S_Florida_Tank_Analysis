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

#tree stuff
tmp <- enframe(sequences, name = 'asv_names', value = 'sequence') %>%
  left_join(taxonomy_tibble, by = 'asv_names') %>%
  rowwise() %>%
  mutate(asv_names = ifelse(asv_names %in% vl_suspects_list, paste("**", asv_names, sep = ""), asv_names)) %>%
  mutate(asv_names = ifelse(asv_names %in% likely_suspects_list, paste("*", asv_names, sep = ""), asv_names)) %>%
  nest(data = c(Genus, Species, asv_names, sequence)) %>%
  filter(Family %in% likely_families) %>%
  mutate(likelihood = ifelse(Family %in% very_likely_families, "very likely", "likely"))

current_fam <- tmp %>% filter(Family %in% c("Thiomicrospiraceae", "Arenicellaceae"))

make_tree <- function(current_fam){
  prelim_tree <- current_fam %>%
    select(data) %>%
    unnest(data) %$%
    set_names(sequence, asv_names) %>%
    DNAStringSet() %>%
    msa() %>%
    msaConvert(., 'phangorn::phyDat') %>%
    modelTest() %>% 
    pml_bb()
  
  bootstrapping <- bootstrap.pml(prelim_tree, bs = 100, optNni = TRUE,
                                 control = pml.control(trace = 1))
  tree_plot <- plotBS(midpoint(prelim_tree$tree), bootstrapping)
  
  return(tree_plot)
}

make_tree(current_fam)

plot_list <- current_fam %>%
  rowwise() %>%
  mutate(plots = list(make_tree(.)))


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

