#code to calculate the beta diversity
setwd("~/Documents/GitHub/16S_Florida_Tank_Analysis/Code")

library(tidyverse)
library(phyloseq)
library(microbiome)
library(vegan)
library(metagMisc)

aggregation_level <- 'Genus' #or none

#### Read in Data ####
microbiome_data <- read_rds("../intermediate_files/preprocess_microbiome.rds")
metadata <- sample_data(microbiome_data) %>%
  as_tibble(rownames = 'sample_id') %>%
  select(-retain_sample)

if(aggregation_level != 'none'){
  microbiome_data <- aggregate_taxa(microbiome_data, aggregation_level)
  taxa_names(microbiome_data) <- str_replace_all(taxa_names(microbiome_data), ' |-', '_')
} else {
  taxa_names(microbiome_data) <- str_c('ASV', 1:length(taxa_names(microbiome_data)), sep = '_')
}

mb_data <- microbiome_data %>%
  subset_samples(time %in% c('T3', 'T7')) %>%
  phyloseq_transform_css %>% #normalizing by column and then log transform
  phyloseq_filter_prevalence(prev.trh = 0.1) %>% #filter for only 10%+ prevalence
  otu_table %>%
  t #transposes the data (flips the axes)


#look into distance metrics

otu_nmds <- metaMDS(mb_data, distance = 'bray', k = 3, trymax = 100, autotransform = FALSE, verbose = TRUE)
plot(otu_nmds)

scores(otu_nmds)$sites %>%
  as_tibble(rownames = 'sample_id') %>%
  left_join(metadata, by = 'sample_id') %>%
  
  ggplot(aes(x = NMDS1, y = NMDS2, colour = final_disease_state)) +
  
  geom_point(data = as_tibble(scores(otu_nmds)$species, rownames = aggregation_level),
             colour = 'gray50', size = 0.1) +
  
  geom_point() +
  theme_classic()


adonis2(mb_data ~ time + exposure + final_disease_state, data = filter(metadata, time %in% c('T3', 'T7')), 
        permutations = 999, method = 'bray', by = NULL)

adonis2(mb_data ~ time + exposure + final_disease_state, data = filter(metadata, time %in% c('T3', 'T7')), 
        permutations = 999, method = 'bray', by = 'margin')

adonis2(mb_data ~ time * exposure * final_disease_state - time:exposure:final_disease_state - exposure:final_disease_state, data = filter(metadata, time %in% c('T3', 'T7')), 
        permutations = 999, method = 'bray', by = 'margin')

adonis2(mb_data ~ time * exposure * final_disease_state - time:exposure:final_disease_state - exposure:final_disease_state, data = filter(metadata, time %in% c('T3', 'T7')), 
        permutations = 999, method = 'bray', by = 'terms')

library(permute)
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


