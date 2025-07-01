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
library(rempsyc)
library(relayer)
library(cowplot)
library(patchwork)
library(strex)
library(fantaxtic) #devtools::install_github("gmteunisse/fantaxtic")
library(ggnested) #devtools::install_github("gmteunisse/ggnested")
library(ggvenn)
library(edgeR) #BiocManager::install("edgeR")
library(magrittr)
library(seqinr)

#### Functions ####

#from https://rdrr.io/github/vmikk/metagMisc/src/R/phyloseq_filter.R
#package metagMisc, needed to change the package type used for the function "prevalence" from "data.table" to "base"
#using the base functions does the same thing but is just slower and less efficient 
#but the data.table method was throwing an error after updating R/RStudio
#Gave the error: "Error in prevalence(physeq) : object 'variable' not found"
mod_phyloseq_filter_prevalence <- function(physeq, prev.trh = 0.05, abund.trh = NULL, threshold_condition = "OR", abund.type = "total"){
  
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

#get the color palette from the modified ggnested function
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

filter_samples <- function(data, model_samples){
  #data is a DGEList
  keep <- colnames(data) %in% model_samples
  message('Retained ', sum(keep), ' samples')
  data[,keep, keep.lib.sizes = TRUE]
}

# venn_asv <- otus_to_analyze
filter_venn <- function(data, venn_asv){
  #venn should be a vector of ASVs which are kept by the venn diagram
  keep <- rownames(data$counts) %in% venn_asv
  message('Retained ', sum(keep), ' ASVs found in overlap of Venn diagram')
  data[keep, keep.lib.sizes = TRUE]
}

#### Read In Data ####
updated_microbiome_data <- read_rds("../intermediate_files/updated_microbiome_data.rds")

asv_sequences <- tax_table(updated_microbiome_data) %>%
  as.data.frame() %>%
  as_tibble(rownames = "asv_id") %>%
  full_join(refseq(updated_microbiome_data) %>% as.data.frame() %>% as_tibble(rownames = "asv_id"), by = join_by("asv_id")) %>%
  rename("sequence" = "x") %>%
  pivot_longer(Domain:Species, names_to = "taxa_level", values_to = "classification") %>%
  mutate(classification = ifelse(is.na(classification), "NA", classification)) %>%
  mutate(combo_taxa = str_c(taxa_level, classification, sep = " - ")) %>%
  group_by(asv_id, sequence) %>%
  reframe(full_taxonomy = str_c(combo_taxa, collapse = "; ")) %>% 
  mutate(full_taxonomy = str_c(asv_id, full_taxonomy, sep = ": "))
  
#Supplementary Data #1 - 16S Sequences
#write.fasta(as.list(asv_sequences$sequence), names = asv_sequences$full_taxonomy, open = "w", file.out = "../Supplementary Data/ECT2025_supplementary_data_1.fasta")

metadata <- sample_data(updated_microbiome_data) %>%
  as_tibble(rownames = 'sample_id') %>%
  mutate(fragment_id = str_c(str_replace_na(exposure, 'NA'), tank, genotype, sep = '_'),
         .after = sample_id) #fragment ID tracks a single fragment regardless of timepoint sampled

#Supplementary Data #2 - Metadata
#write_csv(metadata, "../Supplementary Data/ECT2025_supplementary_data_2.csv")

#melt phyloseq to examine sample data
melted_ps <- mod_phyloseq_filter_prevalence(updated_microbiome_data, 
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
  filter(genotype %in% longitudinal_genos) %>% #genotypes present in T3 and T7
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
                                direction = '=')

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
  select(-c(susceptability, resistance)) %>%
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

# All 22 alpha diversity metrics are significant for the main effect of treatment
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

#process the post hoc comparisons for significant metrics
alpha_significant_models <- mod_alpha_tab %>%
  filter(fdr_treatment < 0.05) %>%
  rowwise() %>%
  mutate(process_postHoc(posthoc)) %>%
  ungroup() %>%
  p_adjust(exclude_cols = c('treatment', 'tank', 'genotype'))

#get statistics for manuscript metrics
manuscript_alpha_metrics <- c("dominance_core_abundance", "chao1")

alpha_significant_models %>%
  select(metric, `fdr_aquarium_=`, posthoc) %>% #fdr_aquarium_= is the FDR corrected p-value for time/moving corals from the field to aquarium tanks
  filter(metric %in% manuscript_alpha_metrics) %>%
  unnest(posthoc) %>%
  filter(contrast == "aquarium")

#format alpha diversity metrics into an easy to read table
formatted_alpha_table <- alpha_significant_models %>%
  select(metric, contains("aquarium"), contains("exposure"), contains("outcome")) %>%
  pivot_longer(-metric, names_to = "value_type", values_to = "value") %>%
  mutate(split_vals = str_split(value_type, "_")) %>%
  rowwise() %>%
  mutate(stat_term = split_vals[[1]][1],
         sig_term = split_vals[[2]][1]) %>%
  select(-c(value_type, split_vals)) %>%
  pivot_wider(names_from = stat_term, values_from = value) %>%
  filter(metric %in% manuscript_alpha_metrics) %>%
  select(-c(qvalue, pvalue)) %>%
  rename("pvalue" = "fdr", "t" = "tvalue")

#output stats table to Word
alpha_table <- nice_table(formatted_alpha_table)
print(alpha_table, preview = "docx")  

#which of the post hocs are significant
alpha_metric_signatures <- alpha_significant_models %>%
  select(metric, starts_with('fdr')) %>% #only looking at FDR-corrected p-values
  select(-contains(c('treatment', 'tank', 'genotype'))) %>% #not looking at main effects, only at posthocs
  mutate(across(starts_with('fdr'), ~. < 0.05)) %>% #convert significance to TRUE/FALSE
  pivot_longer(cols = -metric,
               names_to = c('term'),
               values_to = 'significance') %>%
  mutate(term = str_remove(term, 'fdr_')) %>%
  mutate(direction = str_extract(term, '[><=]'),
         term = str_remove(term, '_[><=]')) %>% #make directionality its own column
  filter(significance) %>%
  ungroup()

#view how many metrics out of 22 are significant for each post hoc comparison
alpha_metric_signatures %>%
  group_by(term) %>%
  reframe(n = n())

#get statistics for manuscript for whether metrics changed significantly between time points
alpha_significant_models %>%
  select(metric, contains("t0"), contains("t3"), contains("t7")) %>% 
  filter(metric %in% manuscript_alpha_metrics) %>%
  pivot_longer(cols = !metric, names_to = "val", values_to = "num") %>%
  mutate(time = str_after_first(val, "_"), val = str_before_first(val, "_")) %>%
  pivot_wider(names_from = val, values_from = num)

#get intercept (T0) for chao1 
summary((alpha_significant_models %>% filter(metric %in% manuscript_alpha_metrics))$model[[1]])

#get intercept (T0) for core abundance
summary((alpha_significant_models %>% filter(metric %in% manuscript_alpha_metrics))$model[[2]])

#### Alpha Diversity Plots ####

alpha_graphs_manuscript <- alpha_significant_models %>%
  filter(metric %in% manuscript_alpha_metrics) %>%
  rowwise() %>%
  mutate(plot_info = list(emmeans(model, ~treatment) %>%
                            broom::tidy(conf.int = TRUE) %>%
                            separate(treatment, into = c('time', 'exposure', 'final_disease_state')) %>%
                            #the next set of lines is adding in identical rows of the T0 points for each of the 3 combinations
                            #of exposure and outcome.  This is to allow identical T0 points to be plotted so that different
                            #colored lines for the different graph categories can be drawn from T0 to T3
                            mutate(graph_cat = ifelse(time == "T0", NA, 
                                                      paste(exposure, final_disease_state, sep = "_"))) %>% # make graph category variable with exposure and outcome
                            {. ->> intermed } %>% #save this current output to intermed
                            mutate(graph_cat = ifelse(time == "T0", "D_D", 
                                                      graph_cat)) %>% #set this T0 point to D_D
                            dplyr::slice(1) %>% #keep only the first row
                            rbind(intermed) %>% #add the rest of the rows from intermed back in
                            mutate(graph_cat = ifelse(is.na(graph_cat), "D_H", 
                                                      graph_cat)) %>% #set the newly added T0 point to D_H
                            dplyr::slice(rep(1:2, 1)) %>% #keep the top 2 rows
                            rbind(intermed) %>% #add the rest of the rows from intermed back in
                            mutate(graph_cat = ifelse(is.na(graph_cat), "H_H", 
                                                      graph_cat)) %>% #set the newly added T0 point to H_H
                            mutate(c_time = parse_number(time)) %>% #make continuous time variable for plotting
                            mutate(facet_lab = "Experimental") %>%
                            mutate(c_time = ifelse(time == "T0", c_time, #adjust the time based on exposure/outcome so data points aren't overlapping
                                                   case_when(graph_cat == "D_D" ~ c_time + 0.30,
                                                             graph_cat == "D_H" ~ c_time,
                                                             graph_cat == "H_H" ~ c_time - 0.30))) %>%
                            mutate(graph_cat = factor(graph_cat, levels = c("D_D", "D_H", "H_H"), labels = c("D_D", "D_H", "H_H"))))) %>%
  rowwise() %>%
  # I used the relayer package to set up multiple legends
  # colour1 is disease exposed and colour2 is healthy exposed
  mutate(plot = list(
    ggplot(data = plot_info, aes(x = c_time, y = estimate, ymin = estimate - std.error, ymax = estimate + std.error)) +
      (geom_line(data = (plot_info %>% filter(graph_cat %in% c("D_H", "D_D"))), aes(colour1 = graph_cat, linetype = graph_cat)) %>%
         rename_geom_aes(new_aes = c("colour" = "colour1"))) + 
      (geom_line(data = (plot_info %>% filter(graph_cat %in% c("H_H"))), aes(colour2 = graph_cat, linetype = graph_cat)) %>%
         rename_geom_aes(new_aes = c("colour" = "colour2"))) +
      (geom_errorbar(data = (plot_info %>% filter(graph_cat %in% c("D_H", "D_D"))), width = 0, aes(colour1 = graph_cat)) %>%
         rename_geom_aes(new_aes = c("colour" = "colour1"))) + 
      (geom_errorbar(data = (plot_info %>% filter(graph_cat %in% c("H_H"))), width = 0, aes(colour2 = graph_cat)) %>%
         rename_geom_aes(new_aes = c("colour" = "colour2"))) +
      (geom_point(data = (plot_info %>% filter(graph_cat %in% c("D_H", "D_D"))), size = 3, aes(colour1 = graph_cat, pch = graph_cat)) %>%
         rename_geom_aes(new_aes = c("colour" = "colour1"))) + 
      (geom_point(data = (plot_info %>% filter(graph_cat %in% c("H_H"))), size = 3, aes(colour2 = graph_cat, pch = graph_cat)) %>%
         rename_geom_aes(new_aes = c("colour" = "colour2"))) +
      scale_color_manual(aesthetics = "colour1", values = c("D_H" = "#22A7B6", "D_D" = "#A70000"), guide = "legend", 
                         name = "Disease Exposed", breaks = c("D_H", "D_D"), labels = c("Healthy", "Diseased")) +
      scale_shape_manual(values = c("D_H" = 17, "D_D" = 17, "H_H" = 16), guide = "none") +
      scale_color_manual(aesthetics = "colour2", values = c("H_H" = "#406F23"), guide = "legend", 
                         name = "Healthy Exposed", labels = c("Healthy")) +
      guides(colour1 = guide_legend(
        override.aes=list(linetype = c(1, 1), shape = c(17, 17))),
        colour2 = guide_legend(
          override.aes=list(linetype = c(6), shape = c(16)))) +
      scale_x_continuous(breaks=c(0, 3, 7)) +
      scale_linetype_manual(values = c("D_H" = 1, "D_D" = 1, "H_H" = 6), guide = "none") +
      theme_bw() +
      xlab("Time") +
      ylab(metric) +
      labs(title = metric)
  )) #will give warnings about ignoring unknown aesthetics, but still functions as desired

#format data for adding lines showing significant differences between time points
sig_lines <- alpha_significant_models %>%
  select(metric, contains("t0"), contains("t3"), contains("t7")) %>% #get posthoc comparisons between timepoints
  select(metric, contains("fdr")) %>%
  filter(metric %in% manuscript_alpha_metrics) %>%
  pivot_longer(-metric, names_to = c("time1", "time2", "sign"), values_to = "significance", 
               names_sep = "_", names_prefix = "fdr_") %>%
  mutate(pval_label = case_when(significance < 0.001 ~ "***",
                                significance < 0.01 ~ "**",
                                significance < 0.05 ~ "*",
                                TRUE ~ "")) %>% #add stars for significance
  mutate(significance = ifelse(significance < 0.05, "yes", "no")) %>%
  mutate(time1 = parse_number(time1), time2 = parse_number(time2)) %>%
  select(-sign) %>%
  rowwise() %>%
  mutate(conf.low = 0, conf.high = 0, estimate = 0, text_x = mean(c(time1, time2)), #put stars halfway between the timepoints
         line_order = case_when(sum(time1, time2) == 7 ~ 3.5, #adjust height for comparing day 0 to day 7
                                sum(time1, time2) == 10 ~ 2, #adjust height for comparing day 3 to day 7
                                sum(time1, time2) == 3 ~ 1)) #adjust height for comparing day 0 to day 3

#make plot for Chao1
chao1_plot <- alpha_graphs_manuscript$plot[[1]] + theme(legend.position="none") + 
  ylab(str_wrap("Chao1 Richness Index", 15)) + labs(title = NULL) +
  xlab(NULL) +
  #add significance lines
  geom_segment(data = sig_lines %>% filter(metric == "chao1"), 
               aes(x = time1, xend = time2, y = 400 + line_order*26.667, 
                   yend = 400 + line_order*26.667,
                   alpha = significance), col = "gray20") +
  #add significance stars
  geom_text(data = sig_lines %>% filter(metric == "chao1"),
            aes(x = text_x, y = 400 + 10 + line_order*26.667,
                label = pval_label), col = "gray20", size = 6) +
  scale_alpha_manual(values = c("yes" = 1, "no" = 0)) +
  theme(axis.text = element_text(size=11),
        axis.title = element_text(size=12))

#make plot for core abundance
core_abun_plot <- alpha_graphs_manuscript$plot[[2]] + theme(legend.position="none") +
  ylab(str_wrap("Core abundance dominance index", 15)) + labs(title = NULL) +
  #add significance lines
  geom_segment(data = sig_lines %>% filter(metric == "dominance_core_abundance"), 
               aes(x = time1, xend = time2, y = 0.75 + line_order*0.05, 
                   yend = 0.75 + line_order*0.05,
                   alpha = significance), col = "gray20") +
  #add significance stars
  geom_text(data = sig_lines %>% filter(metric == "dominance_core_abundance"),
            aes(x = text_x, y = 0.75 + 0.01 + line_order*0.05,
                label = pval_label), col = "gray20", size = 6) +
  scale_alpha_manual(values = c("yes" = 1, "no" = 0)) +
  ylim(0, 1) +
  theme(axis.text = element_text(size=11),
        axis.title = element_text(size=12))

#extract legend from original plot
alpha_legend <- cowplot::get_legend(alpha_graphs_manuscript$plot[[2]])

#combine the alpha diversity plots
both_alphas <- chao1_plot / core_abun_plot + plot_annotation(tag_levels = "A")

#add legend to plots

tiff("Fig3.tiff", units="in", width=10, height=6, res=300)
(both_alphas | alpha_legend) + 
  plot_layout(widths = c(4, 1)) + plot_annotation(tag_levels = list(c("A", "B")))
dev.off()
#export 1000x600

#### Fantaxtic Relative Abundance ####

#code taken and modified from https://github.com/gmteunisse/fantaxtic

top_level <- "Order"
nested_level <- "Genus"
sample_order <- NULL

#get the 11 top Orders and top 4 genera within each Order
top_asv <- nested_top_taxa(updated_microbiome_data,
                           top_tax_level = "Order",
                           nested_tax_level = "Genus",
                           n_top_taxa = 11, 
                           n_nested_taxa = 4)

# Create names for NA taxa
ps_tmp <- top_asv$ps_obj %>%
  name_na_taxa()

# Add labels to taxa with the same names
ps_tmp <- ps_tmp %>%
  label_duplicate_taxa(tax_level = nested_level)

# Convert physeq to df
psdf <- psmelt(ps_tmp)

#make aesthetic/readability adjustments to labels
psdf <- psdf %>%
  rename("old_sample" = "Sample") %>%
  mutate(final_disease_state = ifelse(tank == "Field", "F", final_disease_state)) %>%
  mutate(Sample = str_c(time, exposure, final_disease_state, sep = "_")) %>%
  filter(!Sample %in% c("T3_H_D", "T7_H_D")) %>%
  mutate(facet_level = case_when(Sample == "T0_F_F" ~ "Day 0",
                                 Sample %in% c("T0_D_D", "T0_H_H") ~ "Doses",
                                 time == "T3" ~ "Day 3",
                                 time == "T7" ~ "Day 7")) %>%
  mutate(Sample = case_when(Sample == "T0_F_F" ~ "Field",
                            Sample == "T0_D_D" ~ "D Dose",
                            Sample == "T0_H_H" ~ "H Dose",
                            Sample %in% c("T3_D_D", "T3_D_H") ~ "T3 Diseased",
                            TRUE ~ Sample)) %>%
  mutate(facet_level = factor(facet_level, levels = c("Doses", "Day 0", "Day 3", "Day 7")),
         Sample = factor(Sample, levels = c("H Dose", "D Dose", "Field", "T3_H_H",
                                            "T3 Diseased", "T7_H_H", "T7_D_H", "T7_D_D")))

# Move the merged labels to the appropriate positions in the plot:
# Top merged labels need to be at the top of the plot,
# nested merged labels at the bottom of each group
psdf <- move_label(psdf = psdf,
                   col_name = top_level,
                   label = "Other",
                   pos = 0)
psdf <- move_nested_labels(psdf,
                           top_level = top_level,
                           nested_level = nested_level,
                           top_merged_label = "Other",
                           nested_label = "Other",
                           pos = Inf)

# Reorder samples
if(!is.null(sample_order)){
  if(all(sample_order %in% unique(psdf$Sample))){
    psdf <- psdf %>%
      mutate(Sample = factor(Sample, levels = sample_order))
  } else {
    stop("Error: not all(sample_order %in% sample_names(ps_obj)).")
  }
  
}

custom_palette <- taxon_colours(top_asv$ps_obj, 
                                tax_level = "Order", 
                                palette = c(Rickettsiales = "#051900", 
                                            Spirochaetales = "#500078",
                                            Alteromonadales = "#2e1b09",
                                            Verrucomicrobiales = "#7B9B08",
                                            Flavobacteriales = "#F873BE",
                                            Rhodobacterales = "#CA6200",
                                            Puniceicoccales = "#380135",
                                            Saprospirales = "#111787",
                                            Oceanospirillales = "#078090",
                                            Vibrionales = "#7C2C00",
                                            Thiotrichales = "#590404",
                                            Other = "gray75"))

# Generate the plot
tiff("Fig2.tiff", units="in", width=15, height=9, res=300)
mod_ggnested(psdf,
             aes_string(main_group = top_level,
                        sub_group = nested_level,
                        x = "Sample",
                        y = "Abundance"),
             main_palette = custom_palette,
             gradient_type = "tints") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_nested(theme_bw) + 
  geom_col(position = position_fill()) + 
  facet_grid(col = vars(facet_level), space = "free", scales = "free") + 
  scale_x_discrete(name = element_blank(), labels = c("Field" = "Field", "D Dose" = "Diseased", "H Dose" =
                                                        "Healthy", "T3 Diseased" = str_wrap("Disease- Exposed", 10), "T3_H_H" = str_wrap("Healthy- Exposed", 10),
                                                      "T7_H_H" = str_wrap("Healthy", 10),
                                                      "T7_D_H" = str_wrap("Disease- Exposed Healthy", 10), "T7_D_D" = "Diseased")) +
  guides(fill=guide_legend(title=substitute(bold(bd)~nb, list(bd = "Order", nb = "/ Genus")), nrow = 35),
         col=guide_legend(title=substitute(bold(bd)~nb, list(bd = "Order", nb = "/ Genus"))), nrow = 35) +
  ylab("Relative Abundance") +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=13),
        legend.title = element_text(size=13),
        strip.text = element_text(size=12))
dev.off()
#export 1500x900

#color palette used in this plot
fantaxtic_palette <- get_mod_ggnested_palette(psdf,
                                              aes_string(main_group = top_level,
                                                         sub_group = nested_level,
                                                         x = "Sample",
                                                         y = "Abundance"),
                                              main_palette = custom_palette,
                                              gradient_type = "tints")

write_csv(fantaxtic_palette, "../intermediate_files/fantaxtic_color_palette.csv")

#### ASV Sample Filtering ####

#where are OTUs present across all times and in the doses
venn_all_times_and_doses <- melted_ps %>%
  mutate(across(c(exposure, final_disease_state), factor)) %>%
  filter(Abundance > 0) %>%
  mutate(time = if_else(time == 'T0' & tank == "HOMO", exposure, time)) %>%
  group_by(time, OTU) %>%
  summarise(n = sum(Abundance),
            .groups = 'drop') %>%
  pivot_wider(names_from = time, values_from = n, values_fill = 0L) %>%
  mutate(across(-OTU, ~. > 0)) %>% 
  mutate(T0_H = ifelse(H | T0, TRUE, FALSE)) #is it in T0 and/or the healthy dose?

ggvenn(venn_all_times_and_doses, c('D', 'T0_H', 'T3', 'T7')) + ggtitle("All Times ASV Presence")

#remove things that are only present in T3 and/or T7 - ASVs reasonably should originate in the native microbiomes (T0) or
#on one of the homogenate doses
not_t3t7_only <- filter(venn_all_times_and_doses, 
                        !c(!D & !H & !T0 & T3 & T7 & !T0_H)) %>% #remove things that are only in T3 and T7
  filter(!c(!D & !H & !T0 & !T3 & T7 & !T0_H)) %>% #remove things that are only in T7
  filter(!c(!D & !H & !T0 & T3 & !T7 & !T0_H)) %>% #remove things that are only in T3
  pull(OTU)

venn_all_times_and_doses %>% filter(OTU %in% not_t3t7_only) %>% ggvenn(c('D', 'T0_H', 'T3', 'T7')) + ggtitle("Removed T3/T7 only")

#normalization based on all samples & implement filtering 
otu_tmm <- updated_microbiome_data %>%
  mod_phyloseq_filter_prevalence(prev.trh = 0.2) %>% #20% prevalence filter
  otu_table() %>% 
  t %>% #NOTE: *genus and family do not need the t but ASVs need the t*
  as.data.frame %>%
  as.matrix %>% 
  DGEList(remove.zeros = TRUE) %>%
  edgeR::calcNormFactors(method = 'TMMwsp') %>% #normalize using TMMwsp method
  filter_venn(not_t3t7_only) %>% #remove things that are only in T3 and T7
  filter_samples(model_samples) #remove samples not to be analyzed

#### Variance Weighting ####
param <- SnowParam(parallel::detectCores() - 1, "SOCK", progressbar = TRUE)
dream_weights_fullInteraction <- voomWithDreamWeights(counts = otu_tmm,
                                                      formula = ~ model_comp + (1 | genotype) + (1 | tank),
                                                      
                                                      data = filter(metadata, !str_detect(tank, 'homo|HOMO')) %>%
                                                        filter(genotype %in% longitudinal_genos) %>% #genotypes present in T3 and T7
                                                        filter(!(exposure == "H" & final_disease_state == "D")) %>%
                                                        arrange(sample_id) %>%
                                                        mutate(model_comp = str_c(time, exposure, susceptability)) %>%
                                                        column_to_rownames('sample_id'),
                                                      BPPARAM = param,
                                                      plot = TRUE)

#### ASV Modelling ####
full_data <- otu_tmm %>%
  cpm(log = TRUE, prior.count = 0.5, normalized.lib.sizes = TRUE) %>% #log count per million
  as_tibble(rownames = 'asv_id') %>%
  pivot_longer(cols = -asv_id, names_to = 'sample_id', values_to = 'log2_cpm') %>%
  left_join(dream_weights_fullInteraction$weights %>%
              set_colnames(colnames(dream_weights_fullInteraction$E)) %>%
              set_rownames(rownames(dream_weights_fullInteraction$E)) %>%
              as_tibble(rownames = 'asv_id') %>%
              pivot_longer(cols = -asv_id,
                           names_to = 'sample_id',
                           values_to = 'weight'),
            by = c('asv_id', 'sample_id')) %>% #add in calculated weights
  left_join(as_tibble(otu_tmm$counts, rownames = 'asv_id') %>%
              pivot_longer(cols = -asv_id,
                           names_to = 'sample_id',
                           values_to = 'read_count'),
            by = c('asv_id', 'sample_id')) %>% #add back in read counts
  left_join(metadata, by = 'sample_id') %>% #add in metadata
  left_join(as_tibble(otu_tmm$samples, rownames = 'sample_id') %>% select(-group), by = 'sample_id') %>% #add library size and normalization factors
  left_join(tax_table(updated_microbiome_data) %>% as.data.frame() %>% as_tibble(rownames = 'asv_id'), by = c('asv_id')) #add taxonomy

#fully filtered and processed data
write_csv(full_data, '../intermediate_files/fully_preprocessed_samples.csv.gz')


#### Process Homogenate Dose Data ####
homogenate_data <- updated_microbiome_data %>%
  subset_samples(str_detect(tank, 'HOMO')) %>%
  otu_table() %>%
  t() %>%
  as.data.frame() %>%
  as.matrix() %>% 
  DGEList(remove.zeros = FALSE) %>%
  cpm(log = TRUE, prior.count = 0.5, normalized.lib.sizes = TRUE) %>%
  as_tibble(rownames = 'asv_id') %>%
  pivot_longer(cols = -asv_id,
               names_to = 'sample_id',
               values_to = 'log2_cpm') %>%
  filter(asv_id %in% unique(full_data$asv_id)) %>%
  mutate(exposure = str_extract(sample_id, '[DH]'))

#fully filtered and processed dose data
write_csv(homogenate_data, '../intermediate_files/homogenate_cpm.csv')
