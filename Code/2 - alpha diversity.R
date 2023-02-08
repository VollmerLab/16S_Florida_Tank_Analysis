#code to calculate the alpha diversity
setwd("~/Documents/GitHub/16S_Florida_Tank_Analysis/Code")
library(tidyverse)
library(phyloseq)
library(microbiome)
library(vegan)

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

#data filtering step?

#### Alpha Diversity ####
alpha_table <- microbiome::alpha(microbiome_data, index = "all") %>%
  as_tibble(rownames = 'sample_id') %>%
  inner_join(metadata, by = 'sample_id') %>%
  mutate(fragment_id = str_c(str_replace_na(exposure, 'NA'), tank, genotype, sep = '_'))

alpha_table

#only looking at time points 3 and 7
timepoint_data <- alpha_table %>%
  filter(time %in% c('T3', 'T7'))

lm(observed ~ (exposure + final_disease_state) * time, data = timepoint_data) %>% car::Anova(type = 2)
glm(observed ~ (exposure + final_disease_state) * time, data = timepoint_data, family = 'poisson')

library(lme4)
linear_model <- lmer(observed ~ (exposure + final_disease_state) * time + (1 | fragment_id) + 
                       (1 | tank), data = timepoint_data)
count_model_1 <- glmer(observed ~ (exposure + final_disease_state) * time + (1 | fragment_id) + (1 | tank), 
      data = timepoint_data, family = 'poisson') 
car::Anova(count_model_1, type = 2) #poisson is for count data
summary(count_model_1)

count_model_nb <- glmer.nb(observed ~ (exposure + final_disease_state) * time + (1 | fragment_id) + (1 | tank), 
         data = timepoint_data)
summary(count_model_nb)
car::Anova(count_model_nb, type = 2)

AIC(count_model_nb, count_model_1, linear_model) #better model has smaller AIC value

library(afex)

mixed(observed ~ (exposure + final_disease_state) * time + (1 | fragment_id) + (1 | tank), data = timepoint_data)
      #does essentially the same thing as repeated measures but can do multiple random effects
aov_4(observed ~ exposure + final_disease_state + time + (1 + time | fragment_id), data = timepoint_data)
      #only one random effect at a time

library(emmeans)
select(timepoint_data, exposure, final_di)

emmeans(count_model_nb, ~final_disease_state * time, type = 'link') %>% contrast('pairwise')
emmeans(count_model_nb, ~final_disease_state * time, type = 'response') %>%
  multcomp::cld(Letters = LETTERS) %>% 
  as_tibble %>%
  mutate(.group = str_trim(.group)) %>%
  # rename(response = emmean,
  #        asymp.UCL = upper.CL, 
  #        asymp.LCL = lower.CL) %>%
  ggplot(aes(x = time, y = response, ymin = asymp.LCL, ymax = asymp.UCL, colour = final_disease_state)) +
  geom_pointrange(position = position_dodge(0.5)) +
  geom_text(aes(y = asymp.UCL, label = .group), position = position_dodge(0.5),
            vjust = -1) +
  scale_colour_manual(values = c('D' = 'red', 'H' = 'blue'))


select(alpha_table, sample_id, observed, any_of(colnames(metadata)))

colnames(timepoint_data)

tmp <- timepoint_data %>%
  pivot_longer(cols = observed:rarity_rare_abundance,
               names_to = 'metric',
               values_to = 'value') %>%
  nest_by(metric) %>%
  summarise(model = list(lmer(value ~ (exposure + final_disease_state) * time + 
                                (1 | fragment_id) + (1 | tank), data = data)))

tmp$data[[1]]

tmp

#### Alpha Diversity Index Notes ####

#richness 

#observed is observed species richness
#chao1 is nonparametric method for estimating the number of species in a community, 
    #based on the concept that rare species infer the most information about the number of missing species
    #particularly useful for data sets skewed toward the low-abundance species

#dominance

#dbp is the relative abundance of the most abundant species of the sample. 
    #Index gives values in interval 0 to 1, where bigger value represent greater dominance
#dmn is the sum of relative abundances of the two most abundant species of the sample
#absolute index equals to the absolute abundance of the most dominant n species of the sample 
    #(specify the number with the argument ntaxa)
#relative index equals to the relative abundance of the most dominant n species of the sample 
    #(specify the number with the argument ntaxa). This index gives values in interval 0 to 1
#simpson's lambda is  the probability that two randomly chosen individuals belongs to the same species. 
    #The higher the probability, the greater the dominance
#core_abundance is Core abundance index is sum of relative abundances of core species in the sample. 
    #Index gives values in interval 0 to 1, where bigger value represent greater dominance
    #Core species are species that are most abundant in all samples
#gini measures how unevenly abundances are distributed
    #If there is small group of species that represent large portion of total abundance of microbes, 
    #the inequality is large and Gini index closer to 1. If all species has equally large abundances, 
    #the equality is perfect and Gini index equals 0

#Rarity and low abundance

#log_modulo_skewness is a rarity index that characterizes the concentration of species at low abundance. 
    #It uses the skewness of the frequency distribution of arithmetic abundance classes
#low abundance gives the concentration of species at low abundance, or the relative proportion of rare species in [0,1].
    #The species that are below the indicated detection threshold are considered rare. #use "detection = " (ex: 0.2/100)
    #Note that population prevalence is not considered. If the detection argument is a vector, 
    #then a data.frame is returned, one column for each detection threshold.
#rare abundance gives the relative proportion of rare species (ie. those that are not part of the core microbiota)
    #in the interval [0,1]. This is the complement (1-x) of the core abundance. The rarity function provides the 
    #abundance of the least abundant taxa within each sample, regardless of the population prevalence.

#Coverage

#coverage gives the number of groups needed to have a given proportion of the ecosystem occupied (by default is 0.5 ie 50%)

#Diversity

#inverse_simpson is an indication of the richness in a community with uniform evenness that would have the 
    #same level of diversity, calculated as 1/lambda where lambda is the simpson index #difficult to find info about
#gini-simpson measures the probability that two randomly selected individuals belong to different species
    #1 - lambda where lambda is the simpson index






#### exploratory plots ####

timepoint_data %>%
  ggplot(aes(x = observed)) +
  geom_histogram(bins = 100)

alpha_table %>%
  arrange(exposure) %>%
ggplot() +
  geom_point(aes(x = sample_id, y = observed, col = final_disease_state))

ggplot(data = alpha_table) +
  geom_boxplot(aes(x = disease_state, y = observed, col = final_disease_state)) +
  scale_color_manual(values = c("tomato", "skyblue")) +
  facet_wrap(~exposure)

ggplot(data = timepoint_data) +
  geom_histogram(aes(x = disease_state, fill = final_disease_state), stat="count") +
  scale_fill_manual(values = c("tomato", "skyblue")) +
  facet_wrap(~exposure)

ggplot(data = alpha_table) +
  geom_boxplot(aes(x = disease_state, y = observed, col = final_disease_state)) +
  scale_color_manual(values = c("tomato", "skyblue")) +
  facet_wrap(~time)

ggplot(data = timepoint_data) +
  geom_histogram(aes(x = disease_state, fill = final_disease_state), stat="count") +
  scale_fill_manual(values = c("tomato", "skyblue")) +
  facet_wrap(~time)

ggplot(data = timepoint_data) +
  geom_jitter(aes(x = final_disease_state, y = diversity_shannon, col = disease_state))

timepoint_data %>%
  filter(time == "T7") %>%
ggplot() +
  geom_jitter(aes(x = final_disease_state, y = diversity_shannon, col = disease_state))

timepoint_data %>%
  filter(time == "T7") %>%
  ggplot() +
  geom_jitter(aes(x = final_disease_state, y = observed, col = disease_state))

ggplot(data = timepoint_data) +
  geom_point(aes(x = observed, y = diversity_shannon, col = final_disease_state)) #+
  geom_smooth(method = "lm", se = TRUE)


#### ####
















