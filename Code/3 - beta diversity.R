#code to calculate the beta diversity
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
library(permute)
library(tidyverse)

#### Read in Data ####
aggregation_level <- 'none' #or none

microbiome_data <- read_rds("../intermediate_files/preprocess_microbiome.rds")
metadata <- sample_data(microbiome_data) %>%
  as_tibble(rownames = 'sample_id') %>%
  dplyr::select(-retain_sample)

if(aggregation_level != 'none'){
  microbiome_data <- aggregate_taxa(microbiome_data, aggregation_level)
  taxa_names(microbiome_data) <- str_replace_all(taxa_names(microbiome_data), ' |-', '_')
} else {
  taxa_names(microbiome_data) <- str_c('ASV', 1:length(taxa_names(microbiome_data)), sep = '_')
}

mb_data <- microbiome_data %>%
  #subset_samples(time %in% c('T3', 'T7')) %>%
  phyloseq_transform_css %>% #normalizing by column and then log transform
  phyloseq_filter_prevalence(prev.trh = 0.1) %>% #filter for only 10%+ prevalence
  otu_table %>%
  t #transposes the data (flips the axes)

#### Analysis ####

#sequential permanova 

#envfit(ord = nmds, env = condensed otu table)

#normalize and do log2 cpm 

#nmds for each time point separately
otu_nmds <- metaMDS(mb_data, distance = 'mountford', k = 2, trymax = 100, autotransform = FALSE, verbose = TRUE)
#plot(otu_nmds)

nmds_plot <- scores(otu_nmds)$sites %>%
  as_tibble(rownames = 'sample_id') %>%
  left_join(metadata, by = 'sample_id') %>%
  
  ggplot(aes(x = NMDS1, y = NMDS2, colour = final_disease_state, shape = time)) +
  geom_point(data = as_tibble(scores(otu_nmds)$species, rownames = aggregation_level),
             colour = 'gray60', size = 0.1, shape = 'circle') +
  
  geom_point() +
  theme_classic()

env_dat <- t(mb_data) %>%
  as.data.frame %>%
  as_tibble(rownames = 'asv') %>%
  left_join(as_tibble(as.data.frame(tax_table(microbiome_data)), rownames = 'asv'), 
            by = 'asv') %>%
  group_by(Order) %>%
  summarise(across(where(is.numeric), sum)) %>%
  pivot_longer(cols = -Order) %>%
  pivot_wider(names_from = 'Order',
              values_from = 'value') %>%
  column_to_rownames('name')
 
env_model <- envfit(otu_nmds, env_dat)


env_arrows <- env_model$vectors$arrows %>%
  as_tibble(rownames = 'Order') %>%
  mutate(p = env_model$vectors$pvals) %>%
  filter(p < 0.05) 

nmds_plot + 
  geom_segment(data = env_arrows, aes(xend = 0, yend = 0, x = NMDS1 / 5, y = NMDS2 / 5),
               inherit.aes = FALSE) +
  geom_text(data = filter(env_arrows, str_detect(Order, 'Rick')), aes(x = NMDS1 / 5, y = NMDS2 / 5, label = Order),
            inherit.aes = FALSE)

scores(otu_nmds)$sites %>%
  as_tibble(rownames = 'sample_id') %>%
  filter(NMDS1 < -2)

adonis2(mb_data ~ time + exposure + final_disease_state, data = filter(metadata, time %in% c('T3', 'T7')), 
        permutations = 999, method = 'bray', by = NULL)

adonis2(mb_data ~ time + exposure + final_disease_state, data = filter(metadata, time %in% c('T3', 'T7')), 
        permutations = 999, method = 'bray', by = 'margin')

adonis2(mb_data ~ time * exposure * final_disease_state - time:exposure:final_disease_state - exposure:final_disease_state, data = filter(metadata, time %in% c('T3', 'T7')), 
        permutations = 999, method = 'bray', by = 'margin')

adonis2(mb_data ~ time * exposure * final_disease_state - time:exposure:final_disease_state - exposure:final_disease_state, data = filter(metadata, time %in% c('T3', 'T7')), 
        permutations = 999, method = 'bray', by = 'terms')


tmp_data <- filter(metadata, time %in% c('T3', 'T7')) %>%
  mutate(fragment_id = str_c(exposure, tank, genotype, sep = '_')) 
rda1 <- rda(mb_data ~ time * (exposure + final_disease_state) + Condition(fragment_id), data = tmp_data) 

h <- how(within = Within(type = 'none'), plots = Plots(strata = tmp_data$fragment_id, type = 'free'))
anova(rda1, permutations = h, model = 'reduced', by = NULL)
anova(rda1, permutations = h, model = 'reduced', by = 'term')
anova(rda1, permutations = h, model = 'reduced', by = 'margin')

plot(rda1)

adonis2(mb_data ~ time * (exposure + final_disease_state) , data = tmp_data,
        permutations = h, by = 'margin')


h2 <- how(within = Within(type = 'series', constant = TRUE), plots = Plots(strata = tmp_data$fragment_id, type = 'free'))
h2

anova(rda1, permutations = h2, model = 'reduced', by = 'term')

#### Notes on metaMDS distance metrics ####

#view info and equations for calculating under ?vegdist
#can also find some good explanations at 
      #https://r.qcbs.ca/workshop09/book-en/types-of-distance-coefficients.html  

#good in detecting underlying ecological gradients:
  #"gower" - divides all distances by the number of observations (rows) and scales each column to 
      #unit range, differs from "altGower" by the scaling, how different two records are from 0-1
      #can include categorical variables
  #"bray" - semimetric, the default
      #semimetric means it violates the triangle inequality property
  #"jaccard" - metric, calculated using the Bray–Curtis dissimilarity
      #they think this should be preferred over the default bray
  #"kulczynski" - can be substitute for jaccard when some species have small ranges that are 
      #subsets of larger ranges

#able to handle different sample sizes:
  #"morisita" - can be used with genuine count data (integers) only
      #a statistical measure of dispersion of individuals in a population. 
      #It is used to compare overlap among samples. This formula is based on the assumption that 
      #increasing the size of the samples will increase the diversity 
      #because it will include different habitats
  #"horn" (Horn–Morisita) - variant of morisita, able to handle any abundance data
  #"binomial" - derived from Binomial deviance under null hypothesis that the two compared communities 
      #are equal
  #"chao" - tries to take into account the number of unseen species pairs
  #"cao" - a minimally biased index for high beta diversity and variable sampling intensity
      #intended for count (integer) data, and it is undefined for zero abundance

#presence–absence data (should be able to handle unknown (and variable) sample sizes):
  #"mountford" - inverse of Fisher's alpha, takes into consideration the # of species in each separate
      #community and the # of species that are present in both
      #far less dependent on sample size than jaccards, sorenson, and kulczynski
  #"raup" (Raup–Crick) - based on the probability of observing at least j species in shared 
      #in compared communities

#for compositional data
  #"aitchison" - equivalent to Euclidean distance between CLR-transformed samples (centered log ratio) 
      #and deals with positive compositional data
      #aitchison thinks it's better than bray bc it has a better stability to subsetting and aggregation
      #and is a proper distance
  #"robust.aitchison" - uses robust CLR ("rlcr"), making it applicable to non-negative data 
      #including zeroes (unlike the standard Aitchison)

#not good in gradient separation without proper standardization:
  #"manhattan"
  #"euclidean"

#miscellaneous:
  
  #"canberra" - weighted version of the Manhattan distance, a numerical measure of the distance between 
      #pairs of points in a vector space
      #has also been used to analyze the gut microbiome in different disease states
  #"altGower" - how different two records are from 0-1
      #omits double-zeros and divides by the number of pairs with at least one above-zero value
      #does not scale columns
  #"mahalanobis" - Euclidean distances of a matrix where columns are centred, have unit variance, and 
      #are uncorrelated. The index is not commonly used for community data, but it is sometimes used for 
      #environmental variables
  #"chisq" - Euclidean distances of Chi-square transformed data
  #"chord" - Euclidean distance of a matrix where rows are standardized to unit norm 
      #(their sums of squares are 1) 
  #"hellinger" - used to quantify the similarity between two probability distributions
      #it's a type of f-divergence
  #"clark" - ?

#### ####
