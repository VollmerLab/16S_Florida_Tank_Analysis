setwd("/Users/Emily/Desktop/GitHub/16S_Florida_Tank_Analysis/Code")

#### Libraries ####
library(tidyverse)
library(phyloseq)
library(metagMisc) #devtools::install_github("vmikk/metagMisc")
library(microbiome) #BiocManager::install("microbiome")
library(variancePartition) #BiocManager::install("variancePartition")
library(parallel)
library(multidplyr)
library(lmerTest)
library(emmeans)
library(stringr)



library(magrittr)
library(vegan)
library(edgeR) #BiocManager::install("edgeR")
library(ggvenn)
library(cowplot)
library(ComplexUpset)
library(microshades) #remotes::install_github("KarstensLab/microshades", dependencies = TRUE)
library(relayer)
library(rempsyc)
library(fantaxtic) #devtools::install_github("gmteunisse/fantaxtic")
library(ggnested) #devtools::install_github("gmteunisse/ggnested")

set.seed(68748)

#### Functions ####

#from https://rdrr.io/github/vmikk/metagMisc/src/R/phyloseq_filter.R
#package metagMisc, needed to change the package type used for the function "prevalence" from "data.table" to "base"
#using the base functions does the same thing but is just slower and less efficient 
#but the data.table method was throwing an error after updating R/RStudio
#Gave the error: "Error in prevalence(physeq) : object 'variable' not found"
phyloseq_filter_prevalence <- function(physeq, prev.trh = 0.05, abund.trh = NULL, threshold_condition = "OR", abund.type = "total"){
  
  ## Threshold validation
  if(prev.trh > 1 | prev.trh < 0){ stop("Prevalence threshold should be non-negative value in the range of [0, 1].\n") }
  if(!is.null(abund.trh)){ 
    if(abund.trh <= 0){ stop("Abundance threshold should be non-negative value larger 0.\n") }
  }
  
  # ## Check for the low-prevalence species (compute the total and average prevalences of the features in each phylum)
  # prevdf_smr <- function(prevdf){
  #   plyr::ddply(prevdf, "Phylum", function(df1){ 
  #     data.frame(
  #       Average = mean(df1$Prevalence),
  #       Total = sum(df1$Prevalence))
  #     })
  # }
  # prevdf_smr( prevalence(physeq) )
  
  ## Check the prevalence threshold
  # phyloseq_prevalence_plot(prevdf, physeq)
  
  ## Define prevalence threshold as % of total samples
  ## This function is located in 'phyloseq_prevalence_plot.R' file
  prevalenceThreshold <- prev.trh * phyloseq::nsamples(physeq)
  
  ## Calculate prevalence (number of samples with OTU) and OTU total abundance
  prevdf <- metagMisc::prevalence(physeq, package = "base")
  
  ## Get the abundance type
  if(abund.type == "total") { prevdf$AbundFilt <- prevdf$TotalAbundance }
  if(abund.type == "mean")  { prevdf$AbundFilt <- prevdf$MeanAbundance }
  if(abund.type == "median"){ prevdf$AbundFilt <- prevdf$MedianAbundance }
  
  ## Which taxa to preserve
  if(is.null(abund.trh)) { tt <- prevdf$Prevalence >= prevalenceThreshold }
  if(!is.null(abund.trh)){
    ## Keep OTU if it either occurs in many samples OR it has high abundance
    if(threshold_condition == "OR"){
      tt <- (prevdf$Prevalence >= prevalenceThreshold | prevdf$AbundFilt >= abund.trh)
    }
    
    ## Keep OTU if it occurs in many samples AND it has high abundance
    if(threshold_condition == "AND"){
      tt <- (prevdf$Prevalence >= prevalenceThreshold & prevdf$AbundFilt >= abund.trh)
    }
  }
  
  ## Extract names for the taxa we whant to keep
  keepTaxa <- prevdf$Taxa[ tt ]
  
  ## Execute prevalence filter
  res <- phyloseq::prune_taxa(taxa = keepTaxa, x = physeq)
  return(res)
}

fit_model <- function(formula, data, use_weights = TRUE){
  if(!use_weights){
    data$weight <- 1
  }
  
  full_model <- lmer(formula, 
                     weights = data$weight,
                     data = data, 
                     REML = TRUE,
                     control = variancePartition:::vpcontrol)
  
  
  main_formula <- as.character(formula)
  re_formula <- str_c(main_formula[2], main_formula[1], str_extract_all(main_formula[3], '\\(.*\\)')) %>%
    as.formula()
  
  re_model <- lmer(re_formula, 
                   weights = data$weight,
                   data = data, 
                   REML = TRUE,
                   control = variancePartition:::vpcontrol)
  tibble(model = list(full_model), re_model = list(re_model))
}

process_model <- function(model, re_model, random_anova){
  #create type 3 anova table with KR based p-values, marginal & conditional r2 and eta2 effect size
  #also output as single row all anova results sorted nicely
  #needs fit model and null model with only random effects 
  aov_tab <- anova(model, type = '3', ddf = 'Kenward-Roger')
  
  global_row <- anova(model, re_model) %>% 
    broom::tidy() %>%
    filter(term == 'model') %>%
    select(statistic, df, p.value) %>%
    rename(chisq = statistic,
           pvalue = p.value) %>%
    rename_with(~str_c(., '_global'))
  
  r2_row <- performance::r2(model) %>%
    as_tibble
  
  aov_row <- as_tibble(aov_tab, rownames = 'term') %>%
    mutate(term = str_replace(term, 'time_exposure', 'timeXexposure')) %>%
    rename(ss = 'Sum Sq',
           ms = 'Mean Sq',
           n.DF = NumDF,
           d.DF = DenDF,
           fvalue = 'F value',
           pvalue = 'Pr(>F)') %>%
    rowwise %>%
    mutate(eta2Partial = effectsize::F_to_eta2(fvalue, n.DF, d.DF, ci = NULL)$Eta2_partial) %>%
    tidyr::pivot_wider(names_from = term,
                       values_from = where(is.numeric),
                       names_vary = 'slowest') %>%
    rename_with(~str_replace_all(., ':', 'X'))
  
  varDecomp_row <- VarCorr(model) %>%
    as_tibble() %>%
    mutate(varComp = sdcor^2 / sum(sdcor^2)) %>%
    select(grp, varComp) %>%
    filter(grp != 'Residual') %>%
    left_join(random_anova %>%
                as_tibble(rownames = 'term') %>%
                mutate(term = str_extract(term, '\\| [0-9a-zA-Z]+') %>%
                         str_remove('\\| +')) %>%
                filter(!is.na(term)) %>%
                select(term, Df, LRT, `Pr(>Chisq)`) %>%
                rename(df = Df,
                       chisq = LRT,
                       pvalue = `Pr(>Chisq)`),
              by = c('grp' = 'term')) %>%
    tidyr::pivot_wider(names_from = 'grp',
                       values_from = c('varComp', 'df', 'chisq', 'pvalue'),
                       names_vary = 'slowest') #%>%
  # rename_with(~str_replace_all(., '_', '.'))
  
  tibble::tibble(anova_table = list(aov_tab)) %>%
    bind_cols(global_row, r2_row, varDecomp_row, ., aov_row)
}

run_posthoc <- function(model, contrast_list){
  em_out <- emmeans(model, ~treatment)
  
  contrast_list %>%
    rowwise(direction) %>%
    reframe(emmeans::contrast(em_out,
                              method = contrast$contrasts, 
                              adjust = 'none',
                              side = direction) %>%
              as_tibble)
}

process_postHoc <- function(posthoc){
  post_row <- as_tibble(posthoc) %>%
    dplyr::rename(tvalue = t.ratio,
                  pvalue = p.value) %>%
    mutate(contrast = str_c(contrast, direction, sep = '_'), .keep = 'unused') %>%
    pivot_wider(names_from = c('contrast'),
                values_from = c('estimate', 'SE', 'df', 'tvalue', 'pvalue'),
                names_vary = 'slowest')
  post_row
}

p_adjust <- function(df, exclude_cols = NA_character_){
  exclude_cols <- if_else(is.na(exclude_cols), '@@@', exclude_cols)
  mutate(df, across(c(contains('pvalue'), -contains(exclude_cols)), ~p.adjust(., method = 'fdr'),
                    .names = 'fdr_{.col}')) %>% 
    rename_with(~str_replace_all(., 'fdr_pvalue', 'fdr')) %>% 
    
    mutate(across(c(contains('pvalue'), -contains(exclude_cols)), safe_qvalue,
                  .names = 'qvalue_{.col}')) %>%
    rename_with(~str_replace_all(., 'qvalue_pvalue', 'qvalue')) %>%
    reorder_columns 
}

reorder_columns <- function(df){
  p_cols <- str_subset(colnames(df), 'pvalue')
  fdr_cols <- str_replace(p_cols, 'pvalue', 'fdr')
  q_cols <- str_replace(p_cols, 'pvalue', 'qvalue')
  
  for(col_num in 1:length(p_cols)){
    df <- relocate(df, fdr_cols[col_num], q_cols[col_num], .after = p_cols[col_num])
  }
  df
}

safe_qvalue <- possibly(.f = ~qvalue(.)$qvalues, otherwise = NA_real_)






# data <- otu_tmm; prop_missing <- 0.25
filter_missingness <- function(data, model_samples, prop_missing){
  #data = a DGElist object, model_sample = samples to use for the calculation of % missingness
  #prop_missing = maximum percentage of samples which can be 0 and still keep ASV in dataset
  keep <- rowMeans(data$counts[,colnames(data$counts) %in% model_samples] == 0) <= prop_missing
  
  message('\n')
  message('ASVs Removed for not being expressed in enough samples: ', scales::comma(table(keep)[1]))
  message('ASVs Kept by Filter: ', scales::comma(table(keep)[2]))
  message('\n')
  data[keep, keep.lib.sizes = FALSE]
}

# max_missing_group <- 0.5
filter_missing_groups <- function(data, meta, max_missing_group){
  meta <- mutate(metadata, group_var = str_c(time, exposure, susceptability)) %>%
    dplyr::select(sample_id, group_var) %>%
    filter(!is.na(group_var))
  
  keep <- as_tibble(data$counts, rownames = 'gene_id') %>%
    pivot_longer(cols = -gene_id,
                 names_to = 'sample_id') %>%
    inner_join(meta, by = 'sample_id') %>%
    group_by(gene_id, group_var) %>%
    summarise(prop_missing = mean(value == 0),
              .groups = 'drop') %>%
    pivot_wider(names_from = group_var,
                values_from = prop_missing) %>%
    filter(!if_any(where(is.numeric), ~. > max_missing_group)) %>%
    pull(gene_id)
  keep <- rownames(data) %in% keep
  
  message('\n')
  message('ASVs Removed for having more than ', scales::percent(max_missing_group), ' missing data in at least one group: ', scales::comma(sum(!keep)[1]))
  message('ASVs Kept by Filter: ', scales::comma(sum(keep)))
  message('\n')
  data[keep, keep.lib.sizes = FALSE]
}

filter_samples <- function(data, model_samples){
  #data is a DGEList
  keep <- colnames(data) %in% model_samples
  data[,keep, keep.lib.sizes = TRUE]
}

# venn_asv <- otus_to_analyze
filter_venn <- function(data, venn_asv){
  #venn should be a vector of asv's which are kept by the venn diagram
  keep <- rownames(data$counts) %in% venn_asv
  message('Retained ', sum(keep), ' ASVs found in overlap of Venn diagram')
  data[keep, keep.lib.sizes = TRUE]
}

# data <- otu_tmm; min_count <- 100; meta <- metadata
filter_asv_meanCount <- function(data, meta, min_count){
  mean_counts <- cpm(data, normalized.lib.sizes = FALSE, log = FALSE,
                     prior.count = 0) %>%
    as_tibble(rownames = 'asv_id') %>%
    pivot_longer(cols = -asv_id,
                 names_to = 'sample_id',
                 values_to = 'n') %>%
    left_join(metadata, by = 'sample_id') %>%
    group_by(asv_id) %>%
    summarise(mean_count = mean(n)) %>%
    pull(mean_count)
  
  keep <- mean_counts >= min_count
  message('\n')
  message('ASVs Removed for having less than ', scales::comma(min_count), ' mean number of reads across all samples')
  message('ASVs Kept by Filter: ', scales::comma(sum(keep)))
  message('\n')
  data[keep, keep.lib.sizes = TRUE]
}


#modified version of ggnested that creates the color palette for the ggnested plot based on genus abundance
#so the most abundant overall genus has the darkest color and the least abundant has the lightest
mod_ggnested <- function(data, 
                         mapping = aes(), 
                         ...,
                         legend_labeling = c("sub", "join", "main"),
                         join_str = " - ",
                         legend_title = NULL,
                         main_keys = TRUE,
                         nested_aes = c("fill", "color"),
                         gradient_type = c("both", "shades", "tints"),
                         min_l = 0.05,
                         max_l = 0.95,
                         main_palette = NULL, 
                         base_clr = "#008CF0"
){
  
  # Check if mapping has required args
  aes_args <- names(mapping)
  if (!"main_group" %in% aes_args){
    stop("Error: provide the main_group in the aesthetic mapping argument. For non-nested data, use the regular ggplot2 function.")
  }
  if (!"sub_group" %in% aes_args){
    stop("Error: provide a subgroup in the aesthetic mapping argument. For non-nested data, use the regular ggplot2 function.")
  }
  
  # Show warnings when fill or colour are specified in aes
  if ("fill" %in% aes_args & "fill" %in% nested_aes){
    warning("Warning: fill aesthetics will be ignored in the main ggnested function. Please specify non-nested fill in the geom_* layer. Alternatively,
            remove 'fill' from mapping_aes.")
    mapping$fill <- NULL
  }
  if (("colour" %in% aes_args | "color" %in% aes_args) & ("colour" %in% nested_aes | "color" %in% nested_aes)){
    warning("Warning: colour aesthetics will be ignored in the main ggnested function. Please specify non-nested colour in the geom_* layer. Alternatively,
            remove 'colour' from mapping_aes.")
    mapping$colour <- NULL
    mapping$color <- NULL
  }
  
  # Define group and subgroup
  group <- quo_name(mapping$main_group)
  subgroup <- quo_name(mapping$sub_group)
  
  # Generate the nested palette
  pal0 <- nested_palette(data, group, subgroup, gradient_type, min_l, max_l, 
                         main_palette, base_clr, join_str)
  
  color_order_sub <- pal0$subgroup_colour
  
  pal_step1 <- data %>% select(all_of(subgroup), Abundance) %>%
    group_by(across(all_of(subgroup))) %>%
    reframe(tot_abun = sum(Abundance)) %>%
    arrange(desc(tot_abun))
  
  pal_others <- pal_step1 %>% filter(str_detect(Genus, "Other"))
  
  reordered_pal <- pal_step1 %>%
    filter(!str_detect(Genus, "Other")) %>%
    rbind(., pal_others) %>%
    pull(Genus)
  
  pal <- pal0 %>%
    select(-subgroup_colour) %>%
    group_by(Order) %>%
    slice(match(reordered_pal, Genus)) %>%
    cbind(subgroup_colour = color_order_sub) %>%
    relocate(subgroup_colour, .after = group_colour)
  
  # Extract colours
  colours <- pal %>%
    rename(sublabel = !!subgroup,
           label = !!group) %>%
    as.data.frame()
  
  # Add main_group labels to the legend as extra keys that appear as titles
  if (main_keys){
    colours <- colours %>%
      group_by(label) %>% 
      group_modify(~add_row(.x, .before = 0)) %>%
      ungroup() %>%
      mutate(subgroup_colour = ifelse(is.na(subgroup_colour), "#FFFFFF", subgroup_colour),
             sublabel = ifelse(is.na(sublabel), sprintf("**%s**", as.character(label)), as.character(sublabel)),
             group_subgroup = ifelse(is.na(group_subgroup), sprintf("**%s**", as.character(label)), group_subgroup)) %>%
      as.data.frame() 
  }
  
  # Get the final colours
  vals <- colours$subgroup_colour
  names(vals) <- colours$group_subgroup
  
  # Reorder the data
  df <- left_join(data, pal, by = c(group, subgroup)) %>%
    arrange(group, subgroup) %>%
    mutate(group_subgroup = factor(group_subgroup, ordered = T, levels = colours$group_subgroup),
           !!subgroup := factor(!!sym(subgroup), ordered = T, levels = reordered_pal),
           !!group := factor(!!sym(group), ordered = T)) %>%
    ungroup() %>%
    arrange(group_subgroup)
  
  # Add legend labels and title
  if (legend_labeling[1] == "join"){
    labels <- colours$group_subgroup
    leg_title <- sprintf("%s%s%s", group, join_str, subgroup)
  } else if (legend_labeling[1] == "main"){
    labels <- colours$label
    leg_title <- group
  } else if (legend_labeling[1] == "sub"){
    labels <- colours$sublabel
    leg_title <- subgroup
  } else {
    stop("Invalid option for legend_labeling. Pick one of c('join', 'main', 'sub')")
  }
  
  if (!is.null(legend_title)){
    leg_title <- legend_title
  }
  
  # Generate a scale
  nested_scale <- scale_discrete_manual(..., 
                                        aesthetics = nested_aes, 
                                        name = leg_title, 
                                        values = vals, 
                                        labels = labels, 
                                        drop = F)
  
  # Update mapping
  if ("fill" %in% nested_aes){
    mapping$fill <- quo(group_subgroup)
  }
  if ("colour" %in% nested_aes | "color" %in% nested_aes){
    mapping$colour <- quo(group_subgroup)
  }
  
  # Generate the plot
  p <- ggplot(df, mapping, ...) +
    nested_scale
  if (main_keys){
    p <- p +
      theme_nested(theme)
  }
  return(p)
}

get_mod_ggnested_palette <- function(data, 
                                     mapping = aes(), 
                                     ...,
                                     legend_labeling = c("sub", "join", "main"),
                                     join_str = " - ",
                                     legend_title = NULL,
                                     main_keys = TRUE,
                                     nested_aes = c("fill", "color"),
                                     gradient_type = c("both", "shades", "tints"),
                                     min_l = 0.05,
                                     max_l = 0.95,
                                     main_palette = NULL, 
                                     base_clr = "#008CF0"
){
  # Define group and subgroup
  group <- quo_name(mapping$main_group)
  subgroup <- quo_name(mapping$sub_group)
  
  # Generate the nested palette
  pal0 <- nested_palette(data, group, subgroup, gradient_type, min_l, max_l, 
                         main_palette, base_clr, join_str)
  
  color_order_sub <- pal0$subgroup_colour
  
  pal_step1 <- data %>% select(all_of(subgroup), Abundance) %>%
    group_by(across(all_of(subgroup))) %>%
    reframe(tot_abun = sum(Abundance)) %>%
    arrange(desc(tot_abun))
  
  pal_others <- pal_step1 %>% filter(str_detect(Genus, "Other"))
  
  reordered_pal <- pal_step1 %>%
    filter(!str_detect(Genus, "Other")) %>%
    rbind(., pal_others) %>%
    pull(Genus)
  
  pal <- pal0 %>%
    select(-subgroup_colour) %>%
    group_by(Order) %>%
    slice(match(reordered_pal, Genus)) %>%
    cbind(subgroup_colour = color_order_sub) %>%
    relocate(subgroup_colour, .after = group_colour)
  
  return(pal)
}







#### Read In Data ####
updated_microbiome_data <- read_rds("../intermediate_files/updated_microbiome_data.rds")

metadata <- sample_data(updated_microbiome_data) %>%
  as_tibble(rownames = 'sample_id') %>%
  dplyr::select(-retain_sample) %>%
  mutate(fragment_id = str_c(str_replace_na(exposure, 'NA'), tank, genotype, sep = '_'),
         .after = sample_id) #fragment ID tracks a single fragment regardless of timepoint sampled

#melt phyloseq to examine sample data
melted_ps <- phyloseq_filter_prevalence(updated_microbiome_data, 
                                        prev.trh = 0.2) %>% #prevalence threshold of 20%
  psmelt() %>%
  as_tibble()

#make list of genotypes that are longitudinally present in the tanks (i.e. present at day 3 and at day 7)
longitudinal_genos <- melted_ps %>%
  select(time, genotype) %>%
  group_by(genotype) %>%
  distinct() %>%
  reframe(all_timepoints = str_c(time, collapse = "_")) %>%
  filter(all_timepoints != "T0") %>% # remove 5 fragments used for making doses
  filter(str_detect(all_timepoints, "T3") & str_detect(all_timepoints, "T7")) %>% #remove genotypes that aren't present at both T3 and T7
  pull(genotype)

#make list of samples to be used in downstream analyses
model_samples <- filter(metadata, !str_detect(tank, 'HOMO')) %>% #remove homogenate doses
  filter(genotype %in% longitudinal_genos) %>% #genos present in T3 and T7
  pull(sample_id)

#### Alpha Diversity Post Hoc Contrasts ####

# posthoc_order <- c('T0.F.F', 'T3.D.D', 'T3.D.H', 'T3.H.H', 'T7.D.D', 'T7.D.H', 'T7.H.H')
# order given by emmeans(model, ~treatment)

#compares the average of T3/T7 healthy-exposed healthy corals to T0 field-collected corals
alpha_posthoc_time <- list('aquarium' = c(-1, 0, 0, 1/2, 0, 0, 1/2))

#compares the average of T3/T7 corals exposed to diseased homogenate to the average of T3/T7 corals exposed to healthy homogenate
alpha_posthoc_exp <- list('exposure' = c(0, 1/4, 1/4, -1/2, 1/4, 1/4, -1/2))

#compares the average of T3/T7 corals that contract disease to the average of T3/T7 corals that remain asymptomatic (putatively healthy)
alpha_posthoc_outc <- list('outcome' = c(0, 1/2, -1/4, -1/4, 1/2, -1/4, -1/4))

#compares T3 to T0
alpha_t0_t3 <- list('t0_t3' = c(-1, 1/3, 1/3, 1/3, 0, 0, 0))

#compares T7 to T0
alpha_t0_t7 <- list('t0_t7' = c(-1, 0, 0, 0, 1/3, 1/3, 1/3))

#compares T7 to T3
alpha_t3_t7 <- list('t3_t7' = c(0, -1/3, -1/3, -1/3, 1/3, 1/3, 1/3))


#these comparisons are all two-sided
alpha_two_sided_tests <- tibble(microbial_signature = c('alpha_posthoc_time', 'alpha_posthoc_exp', 'alpha_posthoc_outc',
                                                        'alpha_t0_t3', 'alpha_t0_t7', 'alpha_t3_t7'),
                                contrasts = list(alpha_posthoc_time, alpha_posthoc_exp, alpha_posthoc_outc,
                                                 alpha_t0_t3, alpha_t0_t7, alpha_t3_t7),
                                direction = '=') #=

alpha_posthoc_categories <- alpha_two_sided_tests %>%
  unnest(contrasts) %>%
  mutate(contrast_name = names(contrasts)) %>%
  group_by(contrast_name, contrasts, direction) %>%
  summarise(signatures = list(c(microbial_signature)),
            .groups = 'drop') %>%
  nest(contrast = -direction)


#### Calculate Alpha Diversity ####

#rarefy and calculate alpha diversity metrics
alpha_table <- subset_samples(updated_microbiome_data, sample_names(updated_microbiome_data) %in% model_samples) %>%
  rarefy_even_depth(rngseed = 68748) %>% #rarefying reduces the data down to 4311 taxa
  microbiome::alpha(index = "all") %>% #calculate alpha diversity using the microbiome package
  as_tibble(rownames = 'sample_id') %>%
  inner_join(metadata, by = 'sample_id') %>% #add in metadata
  mutate(fragment_id = str_c(str_replace_na(exposure, 'NA'), tank, genotype, sep = '_'))

#run lmer and posthocs on alpha diversity metrics
mod_alpha_tab <- alpha_table %>%
  mutate(final_disease_state = ifelse(exposure == "F", "F", final_disease_state)) %>%
  mutate(treatment = str_c(time, exposure, final_disease_state, sep = '_')) %>%
  pivot_longer(cols = !c(colnames(metadata), "fragment_id", "treatment"),
               names_to = 'metric',
               values_to = 'alpha_div_value') %>%
  select(-c(susceptability, resistance, clone_group)) %>%
  mutate(tank_field = if_else(str_detect(treatment, 'F'), 'field', 'tank'), .after = final_disease_state,
         weight = 1) %>%
  nest_by(metric) %>%
  mutate(fit_model(alpha_div_value ~ treatment + (1 | genotype) +
                     (0 + dummy(tank_field, c("tank")) | tank),
                   data, 
                   use_weights = FALSE),
         random_anova = list(rand(model)),
         process_model(model, re_model, random_anova),
         posthoc = list(run_posthoc(model, alpha_posthoc_categories))) %>%
  select(-re_model, -ends_with('global')) %>%
  ungroup() %>% 
  p_adjust() %>%
  relocate(anova_table, .after = model) %>% 
  relocate(posthoc, .after = random_anova)

# All 22 metrics are significant for treatment
mod_alpha_tab %>%
  select(metric, starts_with('fdr')) %>% 
  mutate(across(starts_with('fdr'), ~. < 0.05)) %>%
  
  pivot_longer(cols = -metric,
               names_to = c('term'),
               values_to = 'significance',
               names_prefix = 'fdr_') %>%
  filter(significance) %>%
  group_by(metric) %>%
  mutate(sigs = str_c(term, collapse = ", ")) %>%
  select(metric, sigs) %>%
  distinct() %>%
  ungroup() %>%
  group_by(sigs) %>%
  reframe(n = n())

#process the posthocs
alpha_significant_models <- mod_alpha_tab %>%
  filter(fdr_treatment < 0.05) %>%
  rowwise() %>%
  mutate(process_postHoc(posthoc)) %>%
  ungroup() %>%
  p_adjust(exclude_cols = c('treatment', 'tank', 'genotype'))

#p vals for manuscript metrics
alpha_significant_models %>%
  select(metric, `fdr_aquarium_=`) %>% 
  filter(metric %in% c("diversity_shannon", "dominance_core_abundance", "evenness_camargo", "chao1"))

alpha_significant_models %>%
  select(metric, posthoc) %>% 
  filter(metric %in% c("diversity_shannon", "dominance_core_abundance", "evenness_camargo", "chao1")) %>%
  unnest(posthoc) %>%
  filter(contrast == "aquarium")

alpha_significant_models$anova_table[[1]]

formatted_alpha_table <- alpha_significant_models %>%
  select(metric, contains("aquarium"), contains("exposure"), contains("outcome")) %>%
  pivot_longer(-metric, names_to = "value_type", values_to = "value") %>%
  mutate(split_vals = str_split(value_type, "_")) %>%
  rowwise() %>%
  mutate(stat_term = split_vals[[1]][1],
         sig_term = split_vals[[2]][1]) %>%
  select(-c(value_type, split_vals)) %>%
  pivot_wider(names_from = stat_term, values_from = value) %>%
  filter(metric %in% c("chao1", "dominance_core_abundance")) %>%
  select(-c(qvalue, pvalue)) %>%
  rename("pvalue" = "fdr", "t" = "tvalue")

alpha_table <- nice_table(formatted_alpha_table)
print(alpha_table, preview = "docx")  

#which of the posthocs are significant
alpha_metric_signatures <- alpha_significant_models %>%
  select(metric, starts_with('fdr')) %>% 
  select(-contains(c('treatment', 'tank', 'genotype'))) %>%
  mutate(across(starts_with('fdr'), ~. < 0.05)) %>%
  
  pivot_longer(cols = -metric,
               names_to = c('term'),
               values_to = 'significance') %>%
  mutate(term = str_remove(term, 'fdr_')) %>%
  mutate(direction = str_extract(term, '[><=]'),
         term = str_remove(term, '_[><=]')) %>%
  filter(significance) %>%
  ungroup()

# how many sig for each type: 21 for tank/time
alpha_metric_signatures %>%
  group_by(term) %>%
  reframe(n = n())

#rarity_log_modulo_skewness is not sig for time/tank
alpha_significant_models %>%
  select(metric, starts_with('fdr')) %>% 
  select(-contains(c('treatment', 'tank', 'genotype'))) %>%
  mutate(across(starts_with('fdr'), ~. < 0.05)) %>%
  filter(`fdr_aquarium_=` == FALSE)






