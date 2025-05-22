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
  message('Retained ', sum(keep), ' taxa found in overlap of Venn diagram')
  data[keep, keep.lib.sizes = TRUE]
}
#### Read In Data ####

agg_microbiome_data <- read_rds("../intermediate_files/updated_microbiome_data.rds")

#what taxonomic levels to aggregate by
multiple_aggregation_levels <- c("none", "Genus", "Family", "Order") #none means the ASV level

#create an empty vector to store the results
all_corr_results <- c()

#### Get Aggregated Taxonomy ####

#ASV Taxonomy
asv_taxonomy_prelim <- agg_microbiome_data
taxa_names(asv_taxonomy_prelim) <- str_c('ASV', 1:length(taxa_names(asv_taxonomy_prelim)), sep = '_')
asv_taxonomy <- tax_table(asv_taxonomy_prelim) %>% as.data.frame() %>% as_tibble(rownames = 'asv_id')

#Genus Taxonomy
genus_taxonomy_prelim <- aggregate_taxa(agg_microbiome_data, "Genus")
taxa_names(genus_taxonomy_prelim) <- str_replace_all(taxa_names(genus_taxonomy_prelim), ' |-', '_')
genus_taxonomy <- tax_table(genus_taxonomy_prelim) %>% as.data.frame() %>% as_tibble() %>%
  select(-unique)

#Family Taxonomy
family_taxonomy_prelim <- aggregate_taxa(agg_microbiome_data, "Family")
taxa_names(family_taxonomy_prelim) <- str_replace_all(taxa_names(family_taxonomy_prelim), ' |-', '_')
family_taxonomy <- tax_table(family_taxonomy_prelim) %>% as.data.frame() %>% as_tibble()%>%
  select(-unique)

#### Correlation Test For Loop ####
for(var in multiple_aggregation_levels) {
  
  #set aggregation level
  aggregation_level <- var
  message('Calculating correlations between each bacterial ', ifelse(aggregation_level == "none", "ASV", aggregation_level), ' and disease resistance')
  
  #set taxa names and aggregate based on aggregation level
  if(aggregation_level != 'none'){
    agg_microbiome_data <- aggregate_taxa(agg_microbiome_data, aggregation_level)
    taxa_names(agg_microbiome_data) <- str_replace_all(taxa_names(agg_microbiome_data), ' |-', '_')
  } else {
    sequences <- taxa_names(agg_microbiome_data)
    taxa_names(agg_microbiome_data) <- str_c('ASV', 1:length(taxa_names(agg_microbiome_data)), sep = '_')
    names(sequences) <- taxa_names(agg_microbiome_data)
  }
  
  agg_metadata <- sample_data(agg_microbiome_data) %>%
    as_tibble(rownames = 'sample_id') %>%
    dplyr::select(-retain_sample) %>%
    mutate(fragment_id = str_c(str_replace_na(exposure, 'NA'), tank, genotype, sep = '_'),
           .after = sample_id) #fragment ID tracks a single fragment regardless of timepoint sampled
  
  #melt phyloseq to examine sample data
  agg_melted_ps <- mod_phyloseq_filter_prevalence(agg_microbiome_data, 
                                              prev.trh = 0.2) %>% #prevalence threshold of 20%
    psmelt() %>%
    as_tibble()
  
  #make list of genotypes that are longitudinally present in the tanks (i.e. present at day 3 and at day 7)
  longitudinal_genos <- agg_melted_ps %>%
    select(time, genotype) %>%
    group_by(genotype) %>%
    distinct() %>%
    reframe(all_timepoints = str_c(time, collapse = "_")) %>%
    filter(all_timepoints != "T0") %>% # remove 5 fragments used for making doses
    filter(str_detect(all_timepoints, "T3") & str_detect(all_timepoints, "T7")) %>% #remove genotypes that aren't present at both T3 and T7
    pull(genotype)
  
  #make list of samples to be used in downstream analyses
  agg_model_samples <- filter(agg_metadata, !str_detect(tank, 'HOMO')) %>% #remove homogenate doses
    filter(genotype %in% longitudinal_genos) %>% #genotypes present in T3 and T7
    pull(sample_id)
  
  #where are OTUs present across all times and in the doses
  agg_venn_all_times_and_doses <- agg_melted_ps %>%
    mutate(across(c(exposure, final_disease_state), factor)) %>%
    filter(Abundance > 0) %>%
    mutate(time = if_else(time == 'T0' & tank == "HOMO", exposure, time)) %>%
    group_by(time, OTU) %>%
    summarise(n = sum(Abundance),
              .groups = 'drop') %>%
    pivot_wider(names_from = time, values_from = n, values_fill = 0L) %>%
    mutate(across(-OTU, ~. > 0)) %>% 
    mutate(T0_H = ifelse(H | T0, TRUE, FALSE)) #is it in T0 and/or the healthy dose?
  
  ggvenn(agg_venn_all_times_and_doses, c('D', 'T0_H', 'T3', 'T7')) + ggtitle("All Times ASV Presence")
  
  #remove things that are only present in T3 and/or T7 - ASVs reasonably should originate in the native microbiomes (T0) or
  #on one of the homogenate doses
  agg_not_t3t7_only <- filter(agg_venn_all_times_and_doses, 
                          !c(!D & !H & !T0 & T3 & T7 & !T0_H)) %>% #remove things that are only in T3 and T7
    filter(!c(!D & !H & !T0 & !T3 & T7 & !T0_H)) %>% #remove things that are only in T7
    filter(!c(!D & !H & !T0 & T3 & !T7 & !T0_H)) %>% #remove things that are only in T3
    pull(OTU)
  
  #ASVs need to be transposed but higher taxa do not, so set up if else statement that will 
  #transpose the data only if aggregated to the ASV level
  if(aggregation_level == "none"){
    agg_otu_tmm <- agg_microbiome_data %>%
      mod_phyloseq_filter_prevalence(prev.trh = 0.2) %>% #20% prevalence filter
      otu_table() %>% 
      t %>% #NOTE: *genus and family do not need the t but ASVs need the t*
      as.data.frame %>%
      as.matrix %>% 
      DGEList(remove.zeros = TRUE) %>%
      edgeR::calcNormFactors(method = 'TMMwsp') %>% #normalize using TMMwsp method
      filter_venn(agg_not_t3t7_only) %>% #remove things that are only in T3 and T7
      filter_samples(agg_model_samples) #remove samples not to be analyzed
  } else {
    agg_otu_tmm <- agg_microbiome_data %>%
      mod_phyloseq_filter_prevalence(prev.trh = 0.2) %>% #20% prevalence filter
      otu_table() %>% 
      #t %>% #NOTE: *genus and family do not need the t but ASVs need the t*
      as.data.frame %>%
      as.matrix %>% 
      DGEList(remove.zeros = TRUE) %>%
      edgeR::calcNormFactors(method = 'TMMwsp') %>% #normalize using TMMwsp method
      filter_venn(agg_not_t3t7_only) %>% #remove things that are only in T3 and T7
      filter_samples(agg_model_samples) #remove samples not to be analyzed
  }
  
  ## Variance Weighting
  param <- SnowParam(parallel::detectCores() - 1, "SOCK", progressbar = TRUE)
  dream_weights_fullInteraction <- voomWithDreamWeights(counts = agg_otu_tmm,
                                                        formula = ~ model_comp + (1 | genotype) + (1 | tank),
                                                        
                                                        data = filter(agg_metadata, !str_detect(tank, 'homo|HOMO')) %>%
                                                          filter(genotype %in% longitudinal_genos) %>% #genotypes present in T3 and T7
                                                          filter(!(exposure == "H" & final_disease_state == "D")) %>%
                                                          arrange(sample_id) %>%
                                                          mutate(model_comp = str_c(time, exposure, susceptability)) %>%
                                                          column_to_rownames('sample_id'),
                                                        BPPARAM = param,
                                                        plot = TRUE)
  
  ## ASV Modelling
  agg_full_data <- agg_otu_tmm %>%
    cpm(log = TRUE, prior.count = 0.5, normalized.lib.sizes = TRUE) %>% #log count per million
    as_tibble(rownames = 'taxa_id') %>%
    pivot_longer(cols = -taxa_id, names_to = 'sample_id', values_to = 'log2_cpm') %>%
    left_join(dream_weights_fullInteraction$weights %>%
                set_colnames(colnames(dream_weights_fullInteraction$E)) %>%
                set_rownames(rownames(dream_weights_fullInteraction$E)) %>%
                as_tibble(rownames = 'taxa_id') %>%
                pivot_longer(cols = -taxa_id,
                             names_to = 'sample_id',
                             values_to = 'weight'),
              by = c('taxa_id', 'sample_id')) %>% #add in calculated weights
    left_join(as_tibble(agg_otu_tmm$counts, rownames = 'taxa_id') %>%
                pivot_longer(cols = -taxa_id,
                             names_to = 'sample_id',
                             values_to = 'read_count'),
              by = c('taxa_id', 'sample_id')) %>% #add back in read counts
    left_join(agg_metadata, by = 'sample_id') %>% #add in metadata
    left_join(as_tibble(agg_otu_tmm$samples, rownames = 'sample_id') %>% select(-group), by = 'sample_id') %>% #add library size and normalization factors
    left_join(tax_table(agg_microbiome_data) %>% as.data.frame() %>% as_tibble(rownames = 'taxa_id'), by = c('taxa_id')) #add taxonomy
  
  #read in and format data
  agg_normalized_asv_counts <- agg_full_data %>%
    mutate(time = factor(time, ordered = TRUE)) %>%
    mutate(final_disease_state = ifelse(time == "T0", "F", final_disease_state)) %>%
    mutate(treatment = str_c(time, exposure, final_disease_state, sep = '_'),
           time_exposure = str_c(time, exposure, sep = '_'),
           timeC = str_extract(time, '[0-9]+') %>% as.numeric) %>%
    mutate(asv_number = str_extract(taxa_id, '[0-9]+') %>% as.integer)
  
  #run correlation test
  t0_corr_test <- agg_normalized_asv_counts %>%
    filter(time == "T0") %>% #only looking at field-collected T0
    group_by(taxa_id) %>%
    reframe(corr_val = broom::tidy(cor.test(log2_cpm, resistance))) %>% #test for correlation between log2 CPM abundance and fragment disease resistance
    unnest(corr_val) %>%
    mutate(fdr_p.value = p.adjust(p.value, method = "fdr")) %>% #FDR correct p-values
    arrange(fdr_p.value) %>%
    mutate(agg_level = aggregation_level) #add column showing aggregation level
  
  all_corr_results <- rbind(all_corr_results, t0_corr_test) #add these results to the table of results
}

#### Examine Results ####

write_rds(all_corr_results, "../intermediate_files/correlation_test_results.rds")

#check if any correlations are significant
all_corr_results %>% filter(fdr_p.value < 0.05)

#check results for previous probiotic associations
all_corr_results %>% filter(taxa_id %in% c("MD3_55", "Endozoicomonas", "Myxococcales")) #Endozoicomonas is not present in this dataset

#ASV-level correlation test results
asv_corr <- all_corr_results %>% 
  filter(agg_level == "none") %>% 
  left_join(asv_taxonomy, by = join_by("taxa_id" == "asv_id")) #add taxonomic info

#Genus-level correlation test results
genus_corr <- all_corr_results %>% 
  filter(agg_level == "Genus") %>% 
  left_join(genus_taxonomy, by = join_by("taxa_id" == "Genus")) %>% #add taxonomic info
  mutate(Genus = taxa_id, Species = NA)

#Family-level correlation test results
family_corr <- all_corr_results %>% 
  filter(agg_level == "Family") %>% 
  left_join(family_taxonomy, by = join_by("taxa_id" == "Family")) %>% #add taxonomic info
  mutate(Family = taxa_id, Genus = NA, Species = NA) %>%
  mutate(Family = factor(Family, ordered = T, levels = .$Family))

#combine ASV, Genus, and Family level correlation results (with taxonomy added back in) into one tibble to make the plot
corr_plot <- family_corr %>%
  rbind(genus_corr, asv_corr) %>%
  arrange(Family) %>%
  mutate(taxa_id = factor(taxa_id, ordered = T)) %>%
  filter(!is.na(Family))

#plot the correlation coefficients by ASV, Genus, and Family
ggplot(corr_plot, aes(x = Family, y = estimate, col = agg_level)) +
  geom_hline(yintercept = 0) +
  geom_boxplot(position = position_dodge(0.75)) +
  geom_point(data = corr_plot %>% filter(agg_level == "Family"), fill = "#397DBB", 
             position = position_dodge(0.75), pch = 21, size = 2) + #data = (corr_plot %>% filter(agg_level == "none")),
  scale_color_manual(values = c("none" = "#CB3309", "Genus" = "#E6AC0E", "Family" = "transparent"),
                     name = "Aggregation Level", breaks = c("none", "Genus", "Family"), 
                     labels = c("ASV", "Genus", "Family")) +
  coord_flip() +
  theme_bw() +
  ylab("Correlation Coefficient (r)")
#export 1200x900
