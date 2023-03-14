#code to calculate the alpha diversity
setwd("~/Desktop/Screenshots/Career/Vollmer Lab/GitHub/16S_Florida_Tank_Analysis/Code")

#TODO:
#what happens to things without genus and in general #look at documentation
#of the things that we added, how did they change
#mds plot
#T0 alpha
#complex upset

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
library(strex)
library(tidyverse)

#### Functions ####
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
                names_vary = 'fastest')
}

#model <- all_metrics_tp$model[[1]]

get_fixed_effects <- function(model){
  temp_fe <- fixef(model) %>%
  as_tibble(rownames = 'component')
  
  temp_fe$component <- temp_fe$component %>%
  str_replace_all(":",".") %>%
  str_remove_all("[()]") %>%
  str_remove_all("final_") %>%
  str_remove_all("_state")
  
  temp_fe <- temp_fe %>%
  pivot_wider(names_from = 'component',
              values_from = c('value'))
  temp_fe
  
  return(temp_fe)
}

get_aov_p_values <- function(model){
  tibble_results <- car::Anova(model) %>%
  as_tibble(rownames = 'effect') %>%
  janitor::clean_names() %>%
  mutate(p = pr_chisq) %>%
  dplyr::select(effect, p) %>%
  pivot_wider(names_from = 'effect',
              values_from = c('p'))
  
  tibble_results %>%
  {colnames(.) ->> namevar} %>% #{x ->> y} saves x to y w/o messing up the pipe
  {paste("p", namevar, sep = "_") ->> new_names}
  
  colnames(tibble_results) = new_names
  
  return(tibble_results)
}

#### Read in Data ####
aggregation_level <- 'none' #or none

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

all_metrics_tp <- timepoint_data %>%
  pivot_longer(cols = observed:rarity_rare_abundance,
               names_to = 'metric',
               values_to = 'value') %>%
  nest_by(metric) %>%
  mutate(model = list(lmer(value ~ (exposure + final_disease_state) * time + 
                                (1 | fragment_id), data = data)))

#full model with plots
all_metrics_tp <- all_metrics_tp %>%
  ungroup %>%
  rowwise %>% #computes on a data frame one row at a time
  mutate(terms = list(find_unique_significant_terms(model, 0.05))) %>%
  unnest(terms, keep_empty = TRUE) %>%
  rowwise %>%
  mutate(em_out = list(possibly(make_emmean_model, otherwise = NULL)(model, 
                                                                     as.formula(str_c('~', terms)), 
                                                                     0.05))) %>%
  mutate(plot = list(possibly(make_model_plot, otherwise = NULL, quiet = TRUE)(em_out, data, terms))) %>%
  group_by(metric, data) %>%
  reframe(plot = ifelse(any(is.na(terms)), #if NAs in data, then no plot.  else, plot
                        list(NULL),
                        list(wrap_plots(plot) & 
                               labs(y = 'log2(CPM)') &
                               plot_annotation(title = metric) & 
                               theme_classic() &
                               theme(panel.background = element_rect(colour = 'black', fill = NA),
                                     axis.text = element_text(colour = 'black', size = 12),
                                     axis.title = element_text(colour = 'black', size = 16))))) %>%
  rowwise %>% 
  mutate(possibly(make_aov_summary, otherwise = NULL)(model)) %>%
  ungroup

#all diversity metrics - significance table
all_metrics_tp <- timepoint_data %>%
  pivot_longer(cols = observed:rarity_rare_abundance,
               names_to = 'metric',
               values_to = 'value') %>%
  nest_by(metric) %>%
  summarize(model = list(lmer(value ~ (exposure + final_disease_state) * time + 
                                (1 | fragment_id), data = data))) %>%
  ungroup %>%
  rowwise %>% #computes on a data frame one row at a time
  group_by(metric) %>%
  rowwise %>%
  mutate(possibly(get_aov_p_values, otherwise = NULL)(model)) %>%
  rowwise %>%
  mutate(possibly(get_fixed_effects, otherwise = NULL)(model)) %>%
  select(-model) %>%
  ungroup

#### Adding Stars to Significant Results ####
test_mat <- as.matrix(all_metrics_tp)

for(c in 2:6) {
  for(r in 1:nrow(test_mat)) {
    if(test_mat[r , c] < 0.001){
      test_mat[r , c] <- paste("***", as.character(test_mat[r,c]), sep = "")
    } else if(test_mat[r , c] < 0.01){
      test_mat[r , c] <- paste("**", as.character(test_mat[r,c]), sep = "")
    } else if(test_mat[r , c] < 0.05){
      test_mat[r , c] <- paste("*", as.character(test_mat[r,c]), sep = "")
    }
  }
}

#p values with stars added for all metrics
starred_metric_results <- as_tibble(test_mat)

one_star_metrics <- c()
number1 <- 1
two_star_metrics <- c()
number2 <- 1
three_star_metrics <- c()
number3 <- 1
for(c in 2:6) {
  for(r in 1:nrow(test_mat)) {
    if(grepl("*", test_mat[r,c], fixed=TRUE)){
      one_star_metrics[number1] <- test_mat[r,1]
      number1 <- number1 + 1
    }
    if(grepl("**", test_mat[r,c], fixed=TRUE)){
      two_star_metrics[number2] <- test_mat[r,1]
      number2 <- number2 + 1
    }
    if(grepl("***", test_mat[r,c], fixed=TRUE)){
      three_star_metrics[number3] <- test_mat[r,1]
      number3 <- number3 + 1
    }
  }
}
one_star_metrics <- unique(one_star_metrics)
two_star_metrics <- unique(two_star_metrics)
three_star_metrics <- unique(three_star_metrics)

lengths <- max(c(length(one_star_metrics), length(two_star_metrics), length(three_star_metrics)))
length(one_star_metrics) <- lengths
length(two_star_metrics) <- lengths
length(three_star_metrics) <- lengths

#list of which metrics give significant results in at least one column
sig_metric_list <- as_tibble(cbind(one_star_metrics,two_star_metrics,three_star_metrics))

#outputs with stars:
#list of metrics
sig_metric_list
#p-values with stars
starred_metric_results

#listing significant metrics for different aggregations of the data
significant_metrics_all_aggs <- tibble("aggregation" = "", "metric" = "", "interaction" = "",
                                       "value" = "", "stars" = "")
indexval <-  1
i <- 1
for(c in 5:6) {
  for(r in 1:nrow(test_mat)) {
    if(grepl("*", test_mat[r,c], fixed=TRUE)){
      significant_metrics_all_aggs$aggregation[indexval] = agg_levels[i]
      significant_metrics_all_aggs$metric[indexval] = test_mat[r, 1]
      significant_metrics_all_aggs$interaction[indexval] = dimnames(test_mat)[[2]][c]
      significant_metrics_all_aggs$value[indexval] = str_extract_numbers(test_mat[r,c], 
                                                                          decimals = TRUE)[[1]]
      significant_metrics_all_aggs$stars[indexval] = str_extract_non_numerics(test_mat[r,c], 
                                                                          decimals = TRUE)[[1]]
      indexval <- indexval + 1
      significant_metrics_all_aggs <-  add_row(significant_metrics_all_aggs, "aggregation" = "", 
                                               "metric" = "", "interaction" = "", "value" = "", 
                                               "stars" = "")
    }
    
  }
}

#### PCoA Style Plot ####

pcoa_like_plot <- timepoint_data %>%
  pivot_longer(cols = observed:rarity_rare_abundance,
               names_to = 'metric',
               values_to = 'value') %>%
  mutate(time = readr::parse_number(time)) %>%
  nest_by(metric) %>%
  mutate(model = list(lmer(value ~ (exposure + final_disease_state) * time + 
                             (1 | fragment_id), data = data)))  %>%
  rowwise %>%
  mutate(time_terms = list(find_unique_significant_terms(model, 0.05))) %>%
  mutate(pcoa_style_exposure = "") %>%
  mutate(pcoa_style_disease = "")

for(i in 1:nrow(pcoa_like_plot)){
  if(grepl("exposure\\*time", pcoa_like_plot$time_terms[i])){
    pcoa_like_plot$pcoa_style_exposure[i] = list(ggplot(data = pcoa_like_plot$data[[i]]) +
      geom_line(aes(x = time, y = value, col = exposure, pch = fragment_id)) +
      labs(title = paste(pcoa_like_plot$metric[i], "by exposure", sep = " ")))
  }
  if(grepl("final_disease_state\\*time", pcoa_like_plot$time_terms[i])){
    pcoa_like_plot$pcoa_style_disease[i] = list(ggplot(data = pcoa_like_plot$data[[i]]) +
      geom_line(aes(x = time, y = value, col = final_disease_state, pch = fragment_id)) +
      labs(title = paste(pcoa_like_plot$metric[i], "by final disease state", sep = " ")))
  }
}

#### Multiple Aggregation Levels ####
#read in data
microbiome_data <- read_rds("../intermediate_files/preprocess_microbiome.rds")
metadata <- sample_data(microbiome_data) %>%
  as_tibble(rownames = 'sample_id') %>%
  select(-retain_sample)

agg_levels <- c("none", "Genus", "Family")
significant_metrics_all_aggs <- tibble("aggregation" = "", "metric" = "", "interaction" = "",
                                       "value" = "", "stars" = "")
indexval <-  1

for(i in 1:length(agg_levels)){
  aggregation_level = agg_levels[i]

if(aggregation_level != 'none'){
  microbiome_data <- aggregate_taxa(microbiome_data, aggregation_level)
  taxa_names(microbiome_data) <- str_replace_all(taxa_names(microbiome_data), ' |-', '_')
} else {
  taxa_names(microbiome_data) <- str_c('ASV', 1:length(taxa_names(microbiome_data)), sep = '_')
}

alpha_table <- microbiome::alpha(microbiome_data, index = "all") %>%
  as_tibble(rownames = 'sample_id') %>%
  inner_join(metadata, by = 'sample_id') %>%
  mutate(fragment_id = str_c(str_replace_na(exposure, 'NA'), tank, genotype, sep = '_'))

timepoint_data <- alpha_table %>%
  filter(time %in% c('T3', 'T7'))

all_metrics_tp <- timepoint_data %>%
  pivot_longer(cols = observed:rarity_rare_abundance,
               names_to = 'metric',
               values_to = 'value') %>%
  nest_by(metric) %>%
  mutate(model = list(lmer(value ~ (exposure + final_disease_state) * time + 
                             (1 | fragment_id), data = data)))
##assign(paste0("metric_model_", i), all_metrics_tp, globalenv())

#full model with plots
all_metrics_tp_full <- all_metrics_tp %>%
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
  group_by(metric, data) %>%
  reframe(plot = ifelse(any(is.na(terms)), #if NAs in data, then no plot.  else, plot
                        list(NULL),
                        list(wrap_plots(plot) & 
                              labs(y = 'log2(CPM)') &
                              plot_annotation(title = paste(metric, aggregation_level, sep = ", ")) & 
                              theme_classic() &
                              theme(panel.background = element_rect(colour = 'black', fill = NA),
                                     axis.text = element_text(colour = 'black', size = 12),
                                     axis.title = element_text(colour = 'black', size = 16))))) %>%
  rowwise %>% 
  mutate(possibly(make_aov_summary, otherwise = NULL)(model)) %>%
  ungroup 

assign(paste0("full_metric_model_", i), all_metrics_tp_full, globalenv())

all_metrics_plots <- all_metrics_tp_full %>%
  select(c(metric, plot)) %>%
  mutate("aggregation" = aggregation_level)

assign(paste0("all_metrics_plots_", i), all_metrics_plots, globalenv())

#all diversity metrics - significance table
all_metrics_sig_table <- timepoint_data %>%
  pivot_longer(cols = observed:rarity_rare_abundance,
               names_to = 'metric',
               values_to = 'value') %>%
  nest_by(metric) %>%
  summarize(model = list(lmer(value ~ (exposure + final_disease_state) * time + 
                                (1 | fragment_id), data = data))) %>%
  ungroup %>%
  rowwise %>% #computes on a data frame one row at a time
  group_by(metric) %>%
  rowwise %>%
  mutate(possibly(get_aov_p_values, otherwise = NULL)(model)) %>%
  rowwise %>%
  mutate(possibly(get_fixed_effects, otherwise = NULL)(model)) %>%
  select(-model) %>%
  ungroup

assign(paste0("all_sig_metrics_", i), all_metrics_sig_table, globalenv())

temp_mat <- as.matrix(all_metrics_sig_table)

for(c in 2:6) {
  for(r in 1:nrow(temp_mat)) {
    if(as.numeric(temp_mat[r , c]) < 0.001){
      temp_mat[r , c] <- paste("***", as.character(temp_mat[r,c]), sep = "")
    } else if(as.numeric(temp_mat[r , c]) < 0.01){
      temp_mat[r , c] <- paste("**", as.character(temp_mat[r,c]), sep = "")
    } else if(as.numeric(temp_mat[r , c]) < 0.05){
      temp_mat[r , c] <- paste("*", as.character(temp_mat[r,c]), sep = "")
    }
  }
}

starred_metric_results <- as_tibble(temp_mat)

assign(paste0("starred_metric_results_", i), starred_metric_results, globalenv())

#extracting significant metrics for all aggregation levels

for(r in 1:nrow(temp_mat)) {
  for(c in 5:6) {
    if(grepl("*", temp_mat[r,c], fixed=TRUE)){
      significant_metrics_all_aggs$aggregation[indexval] = agg_levels[i]
      significant_metrics_all_aggs$metric[indexval] = temp_mat[r, 1]
      significant_metrics_all_aggs$interaction[indexval] = dimnames(temp_mat)[[2]][c]
      significant_metrics_all_aggs$value[indexval] = str_replace_all(temp_mat[r,c], "\\*", "")
      significant_metrics_all_aggs$stars[indexval] = str_extract_non_numerics(temp_mat[r,c], 
                                                                              decimals = TRUE)[[1]]
      indexval <- indexval + 1
      significant_metrics_all_aggs <-  add_row(significant_metrics_all_aggs, "aggregation" = "", 
                                               "metric" = "", "interaction" = "", "value" = "", 
                                               "stars" = "")
      }
    }
  }
}
#remove extra empty row
significant_metrics_all_aggs <- significant_metrics_all_aggs[!apply(significant_metrics_all_aggs == 
                                                                      "", 1, all),] 
all_aggs_plots <- rbind(all_metrics_plots_1, all_metrics_plots_2, all_metrics_plots_3)

#final tibble showing the significant metrics and plots for all aggregations:
sig_metrics_final <- left_join(significant_metrics_all_aggs, all_aggs_plots)

str_extract_numbers(temp_mat[r,c], 
                    decimals = TRUE)[[1]]

#### Multiple Aggregations T0 - WIP####

#work in progress

microbiome_data <- read_rds("../intermediate_files/preprocess_microbiome.rds")
metadata <- sample_data(microbiome_data) %>%
  as_tibble(rownames = 'sample_id') %>%
  select(-retain_sample)

agg_levels <- c("none", "Genus", "Family")
significant_metrics_all_aggs_T0 <- tibble("aggregation" = "", "metric" = "", "interaction" = "",
                                       "value" = "", "stars" = "")
indexval <-  1

for(i in 1:length(agg_levels)){
  aggregation_level = agg_levels[i]
  
  if(aggregation_level != 'none'){
    microbiome_data <- aggregate_taxa(microbiome_data, aggregation_level)
    taxa_names(microbiome_data) <- str_replace_all(taxa_names(microbiome_data), ' |-', '_')
  } else {
    taxa_names(microbiome_data) <- str_c('ASV', 1:length(taxa_names(microbiome_data)), sep = '_')
  }
  
  alpha_table <- microbiome::alpha(microbiome_data, index = "all") %>%
    as_tibble(rownames = 'sample_id') %>%
    inner_join(metadata, by = 'sample_id') %>%
    mutate(fragment_id = str_c(str_replace_na(exposure, 'NA'), tank, genotype, sep = '_'))
  
  timepoint_data <- alpha_table %>%
    filter(time %in% c('T0'))
  
  all_metrics_tp <- timepoint_data %>%
    pivot_longer(cols = observed:rarity_rare_abundance,
                 names_to = 'metric',
                 values_to = 'value') %>%
    nest_by(metric) %>%
    mutate(model = list(lmer(value ~ (exposure + final_disease_state) + (1 | tank), data = data)))
  ##assign(paste0("metric_model_", i), all_metrics_tp, globalenv())
  
  #full model with plots
  all_metrics_tp_full <- all_metrics_tp %>%
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
    group_by(metric, data) %>%
    reframe(plot = ifelse(any(is.na(terms)), #if NAs in data, then no plot.  else, plot
                          list(NULL),
                          list(wrap_plots(plot) & 
                                 labs(y = 'log2(CPM)') &
                                 plot_annotation(title = paste(metric, aggregation_level, sep = ", ")) & 
                                 theme_classic() &
                                 theme(panel.background = element_rect(colour = 'black', fill = NA),
                                       axis.text = element_text(colour = 'black', size = 12),
                                       axis.title = element_text(colour = 'black', size = 16))))) %>%
    rowwise %>% 
    mutate(possibly(make_aov_summary, otherwise = NULL)(model)) %>%
    ungroup 
  
  ##assign(paste0("full_metric_model_", i), all_metrics_tp_full, globalenv())
  
  all_metrics_plots <- all_metrics_tp_full %>%
    select(c(metric, plot)) %>%
    mutate("aggregation" = aggregation_level)
  
  assign(paste0("all_metrics_plots_T0_", i), all_metrics_plots, globalenv())
  
  #all diversity metrics - significance table
  all_metrics_sig_table <- timepoint_data %>%
    pivot_longer(cols = observed:rarity_rare_abundance,
                 names_to = 'metric',
                 values_to = 'value') %>%
    nest_by(metric) %>%
    summarize(model = list(lmer(value ~ (exposure + final_disease_state), data = data))) %>%
    ungroup %>%
    rowwise %>% #computes on a data frame one row at a time
    group_by(metric) %>%
    rowwise %>%
    mutate(possibly(get_aov_p_values, otherwise = NULL)(model)) %>%
    rowwise %>%
    mutate(possibly(get_fixed_effects, otherwise = NULL)(model)) %>%
    select(-model) %>%
    ungroup
  
  ##assign(paste0("all_sig_metrics_", i), all_metrics_sig_table, globalenv())
  
  temp_mat <- as.matrix(all_metrics_sig_table)
  
  for(c in 2:6) {
    for(r in 1:nrow(temp_mat)) {
      if(temp_mat[r , c] < 0.001){
        temp_mat[r , c] <- paste("***", as.character(temp_mat[r,c]), sep = "")
      } else if(temp_mat[r , c] < 0.01){
        temp_mat[r , c] <- paste("**", as.character(temp_mat[r,c]), sep = "")
      } else if(temp_mat[r , c] < 0.05){
        temp_mat[r , c] <- paste("*", as.character(temp_mat[r,c]), sep = "")
      }
    }
  }
  
  starred_metric_results <- as_tibble(temp_mat)
  
  ##assign(paste0("starred_metric_results_", i), starred_metric_results, globalenv())
  
  #extracting significant metrics for all aggregation levels
  for(c in 5:6) {
    for(r in 1:nrow(temp_mat)) {
      if(grepl("*", temp_mat[r,c], fixed=TRUE)){
        significant_metrics_all_aggs_T0$aggregation[indexval] = agg_levels[i]
        significant_metrics_all_aggs_T0$metric[indexval] = temp_mat[r, 1]
        significant_metrics_all_aggs_T0$interaction[indexval] = dimnames(temp_mat)[[2]][c]
        significant_metrics_all_aggs_T0$value[indexval] = str_extract_numbers(temp_mat[r,c], 
                                                                           decimals = TRUE)[[1]]
        significant_metrics_all_aggs_T0$stars[indexval] = str_extract_non_numerics(temp_mat[r,c], 
                                                                                decimals = TRUE)[[1]]
        indexval <- indexval + 1
        significant_metrics_all_aggs_T0 <-  add_row(significant_metrics_all_aggs_T0, "aggregation" = "", 
                                                 "metric" = "", "interaction" = "", "value" = "", 
                                                 "stars" = "")
      }
    }
  }
}
#remove extra empty row
significant_metrics_all_aggs_T0 <- significant_metrics_all_aggs_T0[!apply(significant_metrics_all_aggs_T0 == 
                                                                      "", 1, all),] 
all_aggs_plots_T0 <- rbind(all_metrics_plots_T0_1, all_metrics_plots_T0_2, all_metrics_plots_T0_3)

#final tibble showing the significant metrics and plots for all aggregations:
sig_metrics_final_T0 <- left_join(significant_metrics_all_aggs_T0, all_aggs_plots_T0)



#### Selected Metrics ####

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
                                (1 | fragment_id), data = data)))



### Shannon Diversity ###

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

### DBP Dominance ###

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

### Gini Dominance ###

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

### Bulla Evenness ###

#evenness w/ equal weight to all species, sensitive to rare species

#nothing significant

anova(r_alpha_models$model[[4]])
#exposure:time .

      #not worth plotting (plot showed nothing)

### Camargo's Evenness ###

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

### Observed Species Richness ###

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

### Rare Abundance ###

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



#### Melted Phyloseq ####

simple_microbiome_og <- psmelt(microbiome_data) 

  ### Rickettsias
  
simple_microbiome <- simple_microbiome_og %>%
  #filter(time %in% c('T3', 'T7')) %>%
  filter(final_disease_state %in% c('D', 'H')) %>%
  mutate(time = readr::parse_number(time)) %>%
  select(-retain_sample)
  
simple_microbiome <- simple_microbiome %>%
  mutate(fragment_id = str_c(str_replace_na(exposure, 'NA'), tank, genotype, sep = '_'))
  
rick <- simple_microbiome %>%
  filter(Order == "Rickettsiales")

ggplot(data = rick) +
  geom_line(aes(x = time, y = Abundance, col = OTU)) +
  facet_wrap(~final_disease_state) +
  theme(legend.position = "none")

#each line is an asv within an individual
ave_rick <- rick %>%
  group_by(fragment_id, OTU, time, final_disease_state) %>%
  summarize(ave_abun = mean(Abundance)) %>%
  ungroup %>%
  filter(ave_abun > 0)

ggplot(data = ave_rick) +
  geom_line(aes(x = time, y = ave_abun, col = OTU, pch = fragment_id), alpha = 0.5) +
  geom_point(aes(x = time, y = ave_abun, col = OTU), alpha = 0.5) +
  facet_wrap(~final_disease_state) +
  theme(legend.position = "none") #+
  ylim(0, 1000)

#each line is an ASV
  asv_rick <- rick %>%
    group_by(OTU, time, final_disease_state) %>%
    summarize(ave_abun = mean(Abundance)) %>%
    ungroup %>%
    filter(ave_abun > 0)  
  
ggplot(data = asv_rick) +
  geom_line(aes(x = time, y = ave_abun, col = OTU), alpha = 0.5) +
  geom_point(aes(x = time, y = ave_abun, col = OTU), alpha = 0.5) +
  facet_wrap(~final_disease_state) +
  theme(legend.position = "none") #+
  ylim(0, 400)
  
  test <- asv_rick %>%
    filter(ave_abun > 10000) #ASV_1 is the very abundant one
  #ASV_1 is ... Rickettsiales Fokiniaceae MD3-55 <NA>
  
  ### all abundances

  simple_microbiome_og

  
  
  
  
  
  
  
  
  
  
  
  
  