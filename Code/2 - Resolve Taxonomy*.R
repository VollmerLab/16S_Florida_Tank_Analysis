setwd("/Users/Emily/Desktop/GitHub/16S_Florida_Tank_Analysis/Code")

#### Libraries ####
library(magrittr)
library(phyloseq)
library(microbiome) #BiocManager::install("microbiome")
library(vegan)
library(metagMisc) #devtools::install_github("vmikk/metagMisc")
library(edgeR) #BiocManager::install("edgeR")
library(variancePartition) #BiocManager::install("variancePartition")
library(ggvenn)
library(cowplot)
library(ComplexUpset)
library(microshades) #remotes::install_github("KarstensLab/microshades", dependencies = TRUE)
library(lmerTest)
library(tidyverse)
library(emmeans)
library(relayer)
library(rempsyc)
library(fantaxtic) #devtools::install_github("gmteunisse/fantaxtic")
library(ggnested) #devtools::install_github("gmteunisse/ggnested")

set.seed(68748)

#### Functions ####
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

safe_qvalue <- possibly(.f = ~qvalue(.)$qvalues, otherwise = NA_real_)

# df <- asv_models
reorder_columns <- function(df){
  p_cols <- str_subset(colnames(df), 'pvalue')
  fdr_cols <- str_replace(p_cols, 'pvalue', 'fdr')
  q_cols <- str_replace(p_cols, 'pvalue', 'qvalue')
  
  for(col_num in 1:length(p_cols)){
    df <- relocate(df, fdr_cols[col_num], q_cols[col_num], .after = p_cols[col_num])
  }
  df
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

# model <- asv_models$model[[1]]; re_model <- asv_models$re_model[[1]]; random_anova <- asv_models$random_anova[[1]]
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

# posthoc <- asv_models$posthoc[[1]]
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


#### Data ####
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


#### Update Taxonomy where Possible ####
#old taxonomy info
taxonomy_tibble <- tax_table(microbiome_data) %>% 
  as.data.frame %>%
  as_tibble(rownames = "asv_names")

#read in data
updated_taxonomy <- read_csv('../intermediate_files/updated_taxonomy.csv')

#replace everything less than 80% confidence with NA
updated_taxonomy_above80 <- updated_taxonomy %>%
  mutate(Domain = ifelse(Domain_confidence > 80, Domain, NA),
         Phylum = ifelse(Phylum_confidence > 80, Phylum, NA),
         Class = ifelse(Class_confidence > 80, Class, NA),
         Order = ifelse(Order_confidence > 80, Order, NA),
         Family = ifelse(Family_confidence > 80, Family, NA),
         Genus = ifelse(Genus_confidence > 80, Genus, NA),
         Species = ifelse(Species_confidence > 80, Species, NA)
  )

unclassified_below_class <-  updated_taxonomy_above80 %>%
  filter(!if_all(Domain:Species, ~is.na(.))) %>%
  filter((is.na(Order) & is.na(Family) & is.na(Genus) & is.na(Species))) %>%
  pull(asv_id)

#get old taxonomy info for ASVs that are all NA in updated taxonomy
all_na_in_new_tax <- updated_taxonomy_above80 %>% 
  filter(if_all(Domain:Species, ~is.na(.)) | asv_id %in% unclassified_below_class) %>% 
  mutate(across(contains("confidence"), ~NA)) %>% 
  select(-c(Domain:Species)) %>%
  left_join(taxonomy_tibble, by = join_by("asv_id" == "asv_names"))

#combine the old taxonomy with the updated version
combined_taxonomy <- updated_taxonomy_above80 %>% 
  filter(!(if_all(Domain:Species, ~is.na(.)) | asv_id %in% unclassified_below_class)) %>%
  full_join(all_na_in_new_tax)

#get list of genera that have multiple described taxonomies
multiple_classifications_list <- combined_taxonomy %>%
  select(Domain:Genus) %>%
  group_by(Genus) %>% 
  distinct() %>%
  summarise(n = n()) %>%
  filter(n > 1, !is.na(Genus)) %>%
  plyr::arrange(Genus) %>%
  pull(Genus)

#get the ncbi classifications for the genera with multiple described taxonomies

# ncbi_classifications_by_genus <- tax_name(multiple_classifications_list, db = "ncbi", 
#                                           get = c("Domain", "Phylum", "Class", "Order", "Family")) %>%
#   select(-c(db, Domain)) %>%
#   dplyr::rename("Genus" = "query")
# write_csv(ncbi_classifications_by_genus, "../intermediate_files/ncbi_classifications_for_overlaps.csv")
ncbi_classifications_by_genus <- read_csv("../intermediate_files/new_ncbi_classifications_for_overlaps.csv")





#most updated taxonomy, just phylum:genus
#each genus has only one described classification, is combined with old and new taxonomies
#doesnt contain NA for genus rows
nonoverlapping_taxonomy <- combined_taxonomy %>% 
  select(Phylum:Genus) %>% 
  distinct() %>% 
  filter(!is.na(Genus)) %>%
  filter(!Genus %in% multiple_classifications_list) %>%
  rbind(ncbi_classifications_by_genus) %>%
  distinct()

tax_name("Candidatus Berkiella", db = "ncbi", get = c("Phylum", "Class", "Order", "Family")) %>%
  select(-db) %>%
  dplyr::rename("Genus" = "query")

#our ASVs with the most up to date taxonomy - almost
full_taxonomy <- combined_taxonomy %>%
  select(-c(Phylum:Family)) %>%
  left_join(nonoverlapping_taxonomy, by = join_by(Genus)) %>%
  relocate(Phylum:Family, .after = Domain) %>%
  filter(!asv_id %in% c(combined_taxonomy %>% filter(is.na(Genus)) %>% pull(asv_id))) %>%
  rbind(combined_taxonomy %>% filter(is.na(Genus))) %>%
  arrange(parse_number(asv_id)) %>%
  mutate(Order = ifelse(Family == "Puniceicoccaceae", "Puniceicoccales", Order)) %>%
  mutate(Family = ifelse(Family == "Unknown Family_4", "Coxiellaceae", Family))
#mutate(Class = ifelse(Order == "Verrucomicrobiales", "Verrucomicrobiia", Class), #fixing misclassifications of Class
#      Class = ifelse(Order == "Puniceicoccales", "Opitutia", Class))


#make version of microbiome data ps to update the taxonomy in
updated_microbiome_data <- microbiome_data

tax_table(updated_microbiome_data) <- full_taxonomy %>% 
  select(-contains("confidence")) %>%
  arrange(parse_number(asv_id)) %>%
  column_to_rownames("asv_id") %>%
  as.matrix()

metadata <- sample_data(updated_microbiome_data) %>%
  as_tibble(rownames = 'sample_id') %>%
  dplyr::select(-retain_sample) %>%
  mutate(fragment_id = str_c(str_replace_na(exposure, 'NA'), tank, genotype, sep = '_'),
         .after = sample_id)


#resolve multiple classifications
# i.e. one family is listed as different classes or orders
aggregated_microbiome_data <- aggregate_taxa(updated_microbiome_data, "Family")

nonunique_taxonomy_families <- tax_table(aggregated_microbiome_data) %>%
  as.data.frame() %>%
  as_tibble() %>%
  filter(str_detect(unique, "_")) %>%
  pull(Family) %>%
  unique()

# ncbi_classifications_by_family <- tax_name(nonunique_taxonomy_families, db = "ncbi",
#                                            get = c("Domain", "Phylum", "Class", "Order")) %>%
#   select(-c(db, Domain)) %>%
#   dplyr::rename("Family" = "query") %>%
#   mutate(Class = ifelse(Order == "Nannocystales" & is.na(Class), "Polyangia", Class))
# 
# write_csv(ncbi_classifications_by_family, "../intermediate_files/ncbi_classifications_for_family_overlaps.csv")



ncbi_classifications_by_family <- read_csv("../intermediate_files/ncbi_classifications_for_family_overlaps.csv")

new_tax_just_non_unique <- tax_table(updated_microbiome_data) %>% 
  as.data.frame() %>%
  as_tibble(rownames = "asv_names") %>%
  filter(Family %in% nonunique_taxonomy_families) %>%
  select(-c(Phylum, Class, Order)) %>%
  left_join(ncbi_classifications_by_family, by = join_by("Family")) %>%
  relocate(c(Phylum, Class, Order), .after = Domain)

unique_family_level_taxonomy <-  tax_table(updated_microbiome_data) %>% 
  as.data.frame() %>%
  as_tibble(rownames = "asv_names") %>%
  filter(!Family %in% nonunique_taxonomy_families) %>%
  rbind(new_tax_just_non_unique) %>%
  arrange(parse_number(asv_names)) %>%
  mutate(Class = ifelse(is.na(Class) & Order == "Polyangiales", "Polyangia", Class),
         Class = ifelse(is.na(Class) & Order == "Haliangiales", "Polyangia", Class))

# fix by order now

aggregated_microbiome_data_o <- aggregate_taxa(updated_microbiome_data, "Order")

nonunique_taxonomy_orders <- tax_table(aggregated_microbiome_data_o) %>%
  as.data.frame() %>%
  as_tibble() %>%
  filter(str_detect(unique, "_")) %>%
  pull(Order) %>%
  unique()

# ncbi_classifications_by_order <- tax_name(nonunique_taxonomy_orders, db = "ncbi",
#                                            get = c("Phylum", "Class")) %>%
#   select(-c(db)) %>%
#   dplyr::rename("Order" = "query")
# 
# write_csv(ncbi_classifications_by_order, "../intermediate_files/ncbi_classifications_for_order_overlaps.csv")

ncbi_classifications_by_order <- read_csv("../intermediate_files/ncbi_classifications_for_order_overlaps.csv")

new_tax_just_non_unique_orders <- unique_family_level_taxonomy %>%
  filter(Order %in% nonunique_taxonomy_orders) %>%
  select(-c(Phylum, Class)) %>%
  left_join(ncbi_classifications_by_order, by = join_by("Order")) %>%
  relocate(c(Phylum, Class), .after = Domain)

unique_family_order_level_taxonomy <-  unique_family_level_taxonomy %>%
  filter(!Order %in% nonunique_taxonomy_orders) %>%
  rbind(new_tax_just_non_unique_orders) %>%
  arrange(parse_number(asv_names))

# fix by class now

aggregated_microbiome_data_c <- aggregate_taxa(updated_microbiome_data, "Class")

nonunique_taxonomy_classes <- tax_table(aggregated_microbiome_data_c) %>%
  as.data.frame() %>%
  as_tibble() %>%
  filter(str_detect(unique, "_")) %>%
  pull(Class) %>%
  unique()

# ncbi_classifications_by_class <- tax_name(nonunique_taxonomy_classes, db = "ncbi",
#                                            get = c("Phylum")) %>%
#   select(-c(db)) %>%
#   dplyr::rename("Class" = "query")
# 
# write_csv(ncbi_classifications_by_class, "../intermediate_files/ncbi_classifications_for_class_overlaps.csv")

ncbi_classifications_by_class <- read_csv("../intermediate_files/ncbi_classifications_for_class_overlaps.csv")

new_tax_just_non_unique_classes <- unique_family_order_level_taxonomy %>%
  filter(Class %in% nonunique_taxonomy_classes) %>%
  select(-c(Phylum)) %>%
  left_join(ncbi_classifications_by_class, by = join_by("Class")) %>%
  relocate(c(Phylum), .after = Domain)

unique_family_order_class_level_taxonomy <-  unique_family_order_level_taxonomy %>%
  filter(!Class %in% nonunique_taxonomy_classes) %>%
  rbind(new_tax_just_non_unique_classes) %>%
  arrange(parse_number(asv_names))

#re-add back to phyloseq object
tax_table(updated_microbiome_data) <- unique_family_order_class_level_taxonomy %>%
  column_to_rownames("asv_names") %>%
  as.matrix()

metadata <- sample_data(updated_microbiome_data) %>%
  as_tibble(rownames = 'sample_id') %>%
  dplyr::select(-retain_sample) %>%
  mutate(fragment_id = str_c(str_replace_na(exposure, 'NA'), tank, genotype, sep = '_'),
         .after = sample_id)

#total number of each taxa that were initially found:
tax_table(updated_microbiome_data) %>% 
  as.data.frame() %>% 
  as_tibble() %>%
  select(Family) %>% #change taxa level here
  filter(!is.na(.)) %>%
  distinct() %>%
  nrow()

#write_rds(updated_microbiome_data, "../intermediate_files/updated_microbiome_data.rds")