#code to calculate the alpha diversity
setwd("~/Desktop/Screenshots/Career/Vollmer Lab/GitHub/16S_Florida_Tank_Analysis/Code")

#### Packages ####
library(phyloseq)
library(microbiome)
library(vegan)
library(lme4)
library(afex)
library(emmeans)
library(car)
library(emmeans)
library(multcomp)
library(tidyverse)

#### Read in Data ####
aggregation_level <- 'Genus' #or none

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
    #must use untrimmed data set for alpha diversity measures of richness to get meaningful results
    #functions depend heavily on singletons
#TODO add read abundance sort of 1000 reads or more
  
    #there are 2 samples with 5 or less species observed, might be skewing the data

#### Alpha Diversity ####
alpha_table <- microbiome::alpha(microbiome_data, index = "all") %>%
  as_tibble(rownames = 'sample_id') %>%
  inner_join(metadata, by = 'sample_id') %>%
  mutate(fragment_id = str_c(str_replace_na(exposure, 'NA'), tank, genotype, sep = '_'))

alpha_table

#only looking at time points 3 and 7
timepoint_data <- alpha_table %>%
  filter(time %in% c('T3', 'T7'))


#### Various Models ####
lm(observed ~ (exposure + final_disease_state) * time, data = timepoint_data) %>% car::Anova(type = 2)
glm(observed ~ (exposure + final_disease_state) * time, data = timepoint_data, family = 'poisson')


linear_model <- lmer(observed ~ (exposure + final_disease_state) * time + (1 | fragment_id) + 
                       (1 | tank), data = timepoint_data)
count_model_1 <- glmer(observed ~ (exposure + final_disease_state) * time + 
                         (1 | fragment_id) + (1 | tank), data = timepoint_data, family = 'poisson') 
car::Anova(count_model_1, type = 2) #poisson is for count data
summary(count_model_1)

count_model_nb <- glmer.nb(observed ~ (exposure + final_disease_state) * time + 
                             (1 | fragment_id) + (1 | tank), 
         data = timepoint_data)
summary(count_model_nb)
car::Anova(count_model_nb, type = 2)

AIC(count_model_nb, count_model_1, linear_model) #better model has smaller AIC value

mixed(observed ~ (exposure + final_disease_state) * time + 
        (1 | fragment_id) + (1 | tank), data = timepoint_data)
      #does essentially the same thing as repeated measures but can do multiple random effects
aov_4(observed ~ exposure + final_disease_state + time + (1 + time | fragment_id), data = timepoint_data)
      #only one random effect at a time

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

#### Analyzing Alpha Table ####

#all diversity metrics in table
select(alpha_table, sample_id, observed, any_of(colnames(metadata)))

colnames(timepoint_data)

tmp <- timepoint_data %>%
  pivot_longer(cols = observed:rarity_rare_abundance,
               names_to = 'metric',
               values_to = 'value') %>%
  nest_by(metric) %>%
  summarise(model = list(lmer(value ~ (exposure + final_disease_state) * time + 
                                (1 | fragment_id) + (1 | tank), data = data)))

#Selected measures:

#richness - observed; dominance - dbp, gini; rarity - rare abundance
#diversity - shannon; evenness - camargo, bulla

kept_columns <- c("sample_id", "fragment_id", "observed", "dominance_dbp", 
                  "dominance_gini", "rarity_rare_abundance",
                  "diversity_shannon", "evenness_camargo", "evenness_bulla")
alpha_kept_columns <- paste("alpha_table$", kept_columns, sep = "")

reduced_alpha_table <- timepoint_data %>% #only T3 and T7
  select(all_of(kept_columns), any_of(colnames(metadata)))

r_alpha_models <- reduced_alpha_table %>%
  pivot_longer(cols = observed:evenness_bulla,
               names_to = 'metric',
               values_to = 'value') %>%
  nest_by(metric) %>%
  summarise(model = list(lmer(value ~ (exposure + final_disease_state) * time + 
                                (1 | fragment_id) + (1 | tank), data = data)))


#core abundance and relative abundance

#broom , broom mixed, fixef  -> adapt make_aov_summary function

#just fragment id as random effect

#what happens to things without genus and in general #look at documentation

#family genus and ASVs
    #extra credit if i do all three at once

#collapse ASVs and plot NAs as gray

#of the things that we added, how did they change

#mds plot

#make filter and add to data preprocessing 1 file and recode NAs

#add T0 to alpha?

#### Shannon Diversity ####

#how diverse the species in a given community are

#greater diversity at T7 than T3

anova(r_alpha_models$model[[1]])
#time ***

#exposure:time *
emmeans(r_alpha_models$model[[1]], ~(exposure) * time, 
                         type = 'response') %>%
  cld(Letters = LETTERS) %>%
  as_tibble() %>%
  mutate(.group = str_trim(.group)) %>%
  #rename(emmean = response) %>%
  ggplot(aes(x = time, y = emmean, ymin = emmean - SE, ymax = emmean + SE, 
             col = exposure)) +
  geom_pointrange(position = position_dodge(0.5)) +
  geom_text(aes(y = (emmean + SE), label = .group),
            position = position_dodge(0.5), vjust = -1) +
  ylab(r_alpha_models$metric[1])

#final_disease_state:time *
emmeans(r_alpha_models$model[[1]], ~(final_disease_state) * time, 
        type = 'response') %>%
  cld(Letters = LETTERS) %>%
  as_tibble() %>%
  mutate(.group = str_trim(.group)) %>%
  #rename(emmean = response) %>%
  ggplot(aes(x = time, y = emmean, ymin = emmean - SE, ymax = emmean + SE, 
             col = final_disease_state)) +
  geom_pointrange(position = position_dodge(0.5)) +
  geom_text(aes(y = (emmean + SE), label = .group),
            position = position_dodge(0.5), vjust = -1) +
  ylab(r_alpha_models$metric[1])

#all
emmeans(r_alpha_models$model[[1]], ~(exposure + final_disease_state) * time, 
        type = 'response') %>%
  cld(Letters = LETTERS) %>%
  as_tibble() %>%
  mutate(.group = str_trim(.group)) %>%
  #rename(emmean = response) %>%
  ggplot(aes(x = time, y = emmean, ymin = emmean - SE, ymax = emmean + SE, 
             col = exposure, shape = final_disease_state)) +
  facet_wrap(~final_disease_state) +
  geom_pointrange(position = position_dodge(0.5)) +
  geom_text(aes(y = (emmean + SE), label = .group),
            position = position_dodge(0.5), vjust = -1) +
  ylab(r_alpha_models$metric[1])

#### DBP Dominance ####

#relative abundance of most abundant species, 0-1 & bigger #s means more dominant

#greater dominance at T3 than T7, matches the increase in shannon diversity at T7

anova(r_alpha_models$model[[2]])
#time ***
#final_disease_state:time .

#exposure:time *
emmeans(r_alpha_models$model[[2]], ~(exposure) * time, 
        type = 'response') %>%
  cld(Letters = LETTERS) %>%
  as_tibble() %>%
  mutate(.group = str_trim(.group)) %>%
  #rename(emmean = response) %>%
  ggplot(aes(x = time, y = emmean, ymin = emmean - SE, ymax = emmean + SE, 
             col = exposure)) +
  geom_pointrange(position = position_dodge(0.5)) +
  geom_text(aes(y = (emmean + SE), label = .group),
            position = position_dodge(0.5), vjust = -1) +
  ylab(r_alpha_models$metric[2])

#all
emmeans(r_alpha_models$model[[2]], ~(exposure + final_disease_state) * time, 
        type = 'response') %>%
  cld(Letters = LETTERS) %>%
  as_tibble() %>%
  mutate(.group = str_trim(.group)) %>%
  #rename(emmean = response) %>%
  ggplot(aes(x = time, y = emmean, ymin = emmean - SE, ymax = emmean + SE, 
             col = exposure, shape = final_disease_state)) +
  facet_wrap(~final_disease_state) +
  geom_pointrange(position = position_dodge(0.5)) +
  geom_text(aes(y = (emmean + SE), label = .group),
            position = position_dodge(0.5), vjust = -1) +
  ylab(r_alpha_models$metric[2])

#### Gini Dominance ####

#how unevenly abundances are distributed, 0-1 & perfect equality is 0

#slight increase in evenness at T7, matches the decreased dominance
#scale of y-axis is very slight

anova(r_alpha_models$model[[3]])
#time ***
#exposure:time .

#final_disease_state:time *
emmeans(r_alpha_models$model[[3]], ~(final_disease_state) * time, 
        type = 'response') %>%
  cld(Letters = LETTERS) %>%
  as_tibble() %>%
  mutate(.group = str_trim(.group)) %>%
  #rename(emmean = response) %>%
  ggplot(aes(x = time, y = emmean, ymin = emmean - SE, ymax = emmean + SE, 
             col = final_disease_state)) +
  geom_pointrange(position = position_dodge(0.5)) +
  geom_text(aes(y = (emmean + SE), label = .group),
            position = position_dodge(0.5), vjust = -1) +
  ylab(r_alpha_models$metric[3])

#all
emmeans(r_alpha_models$model[[3]], ~(exposure + final_disease_state) * time, 
        type = 'response') %>%
  cld(Letters = LETTERS) %>%
  as_tibble() %>%
  mutate(.group = str_trim(.group)) %>%
  #rename(emmean = response) %>%
  ggplot(aes(x = time, y = emmean, ymin = emmean - SE, ymax = emmean + SE, 
             col = exposure, shape = final_disease_state)) +
  facet_wrap(~final_disease_state) +
  geom_pointrange(position = position_dodge(0.5)) +
  geom_text(aes(y = (emmean + SE), label = .group),
            position = position_dodge(0.5), vjust = -1) +
  ylab(r_alpha_models$metric[3])

#### Bulla Evenness ####

#evenness w/ equal weight to all species, sensitive to rare species

#nothing significant

anova(r_alpha_models$model[[4]])
#exposure:time .

      #not worth plotting (plot showed nothing)

#### Camargo's Evenness ####

#proportions of indivs between sites, 0-1 & 1 is even, 0 is patchy

#T3 is really even and T7 is patchier

anova(r_alpha_models$model[[5]])
#final_disease_state **
#time ***

#final_disease_state:time **
emmeans(r_alpha_models$model[[5]], ~(final_disease_state) * time, 
        type = 'response') %>%
  cld(Letters = LETTERS) %>%
  as_tibble() %>%
  mutate(.group = str_trim(.group)) %>%
  #rename(emmean = response) %>%
  ggplot(aes(x = time, y = emmean, ymin = emmean - SE, ymax = emmean + SE, 
             col = final_disease_state)) +
  geom_pointrange(position = position_dodge(0.5)) +
  geom_text(aes(y = (emmean + SE), label = .group),
            position = position_dodge(0.5), vjust = -1) +
  ylab(r_alpha_models$metric[5])

#### Observed Species Richness ####

#number of species observed

#more species richness at T7 than T3, matches shannon diversity

anova(r_alpha_models$model[[6]])
#time ***

#final_disease_state:time **
emmeans(r_alpha_models$model[[6]], ~(final_disease_state) * time, 
        type = 'response') %>%
  cld(Letters = LETTERS) %>%
  as_tibble() %>%
  mutate(.group = str_trim(.group)) %>%
  #rename(emmean = response) %>%
  ggplot(aes(x = time, y = emmean, ymin = emmean - SE, ymax = emmean + SE, 
             col = final_disease_state)) +
  geom_pointrange(position = position_dodge(0.5)) +
  geom_text(aes(y = (emmean + SE), label = .group),
            position = position_dodge(0.5), vjust = -1) +
  ylab(r_alpha_models$metric[6])

#### Rare Abundance ####

#relative proportion of rare species(i.e. not most abundant in all sites), 0-1

#for healthy samples, rare species increased in abundance at T7

anova(r_alpha_models$model[[7]])
#final_disease_state ***
#time ***

#exposure:time *
emmeans(r_alpha_models$model[[7]], ~(exposure) * time, 
        type = 'response') %>%
  cld(Letters = LETTERS) %>%
  as_tibble() %>%
  mutate(.group = str_trim(.group)) %>%
  #rename(emmean = response) %>%
  ggplot(aes(x = time, y = emmean, ymin = emmean - SE, ymax = emmean + SE, 
             col = exposure)) +
  geom_pointrange(position = position_dodge(0.5)) +
  geom_text(aes(y = (emmean + SE), label = .group),
            position = position_dodge(0.5), vjust = -1) +
  ylab(r_alpha_models$metric[7])

#final_disease_state:time ***
emmeans(r_alpha_models$model[[7]], ~(final_disease_state) * time, 
        type = 'response') %>%
  cld(Letters = LETTERS) %>%
  as_tibble() %>%
  mutate(.group = str_trim(.group)) %>%
  #rename(emmean = response) %>%
  ggplot(aes(x = time, y = emmean, ymin = emmean - SE, ymax = emmean + SE, 
             col = final_disease_state)) +
  geom_pointrange(position = position_dodge(0.5)) +
  geom_text(aes(y = (emmean + SE), label = .group),
            position = position_dodge(0.5), vjust = -1) +
  ylab(r_alpha_models$metric[7])

#all
emmeans(r_alpha_models$model[[7]], ~(exposure + final_disease_state) * time, 
        type = 'response') %>%
  cld(Letters = LETTERS) %>%
  as_tibble() %>%
  mutate(.group = str_trim(.group)) %>%
  #rename(emmean = response) %>%
  ggplot(aes(x = time, y = emmean, ymin = emmean - SE, ymax = emmean + SE, 
             col = exposure, shape = final_disease_state)) +
  facet_wrap(~final_disease_state) +
  geom_pointrange(position = position_dodge(0.5)) +
  geom_text(aes(y = (emmean + SE), label = .group),
            position = position_dodge(0.5), vjust = -1) +
  ylab(r_alpha_models$metric[7])

#buffer line 

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
#simpson's lambda is the probability that two randomly chosen individuals belongs to the same species. 
    #The higher the probability, the greater the dominance
#core_abundance is the sum of relative abundances of core species in the sample. 
    #Index gives values in interval 0 to 1, where bigger value represent greater dominance
    #Core species are species that are most abundant in all samples
#gini measures how unevenly abundances are distributed
    #If there is small group of species that represent large portion of total abundance of microbes, 
    #the inequality is large and Gini index closer to 1. If all species has equally large abundances, 
    #the equality is perfect and Gini index equals 0

#Rarity and low abundance

#log_modulo_skewness is a rarity index that characterizes the concentration of species at low abundance. 
    #It uses the skewness of the frequency distribution of arithmetic abundance classes
#low abundance gives the concentration of species at low abundance, or the relative proportion of rare 
    #species in [0,1].The species that are below the indicated detection threshold are considered rare.
    # #use "detection = " (ex: detection = 0.2/100)
    #Note that population prevalence is not considered. If the detection argument is a vector, 
    #then a data.frame is returned, one column for each detection threshold.
#rare abundance gives the relative proportion of rare species in the interval [0,1].
    #(rare = those that are not part of the core microbiota)
    #This is the complement (1-x) of the core abundance. The rarity function provides the 
    #abundance of the least abundant taxa within each sample, regardless of the population prevalence.

#Diversity

#inverse_simpson is an indication of the richness in a community with uniform evenness that would have
    #the same level of diversity, calculated as 1/lambda where lambda is the simpson index 
#gini-simpson measures the probability that two randomly selected individuals belong to different species
    #1 - lambda where lambda is the simpson index
#shannon shows how diverse the species in a given community are
    #It rises with the number of species and the evenness of their abundance
#fisher's alpha describes mathematically the relationship between the number of species 
    #and the number of individuals in those species
#coverage gives the number of groups needed to have a given proportion of the ecosystem occupied 
    #(by default is 0.5 ie 50%)

#Evenness

#camargo's compares proportions of individuals between sites, with 1 being even and 0 being patchy
    #relatively unaffected by sites with very few organisms, and is unaffected by site richness
#simpson's evenness is a variant of the reciprocal Simpson index
    #Index values range from near 0 (1/s) (patchy or skewed) to 1 (even), and the index 
    #is relatively unaffected by sites with very few individuals.
#pielou is shannon diversity index value divided by the maximum possible shannon diversity index given
    #complete evenness (proportion 0 to 1)
#evar is based on the variance in abundance over the species taken over log abundance
    #so proportional differences are compared, then converted to a 0-1 scale by arctan
#bulla gives equal weight to all species regardless of abundance so it's 
    #sensitive to the presence of rare species


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
















