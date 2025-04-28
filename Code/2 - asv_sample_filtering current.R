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

#### read in updated taxonomy ####
updated_microbiome_data <- read_rds("../intermediate_files/updated_microbiome_data.rds")

metadata <- sample_data(updated_microbiome_data) %>%
  as_tibble(rownames = 'sample_id') %>%
  dplyr::select(-retain_sample) %>%
  mutate(fragment_id = str_c(str_replace_na(exposure, 'NA'), tank, genotype, sep = '_'),
         .after = sample_id)

#melted ps
melted_ps <- phyloseq_filter_prevalence(updated_microbiome_data, 
                                        prev.trh = 0.2) %>%
  psmelt() %>%
  as_tibble()

longitudinal_genos <- melted_ps %>%
  select(time, genotype) %>%
  group_by(genotype) %>%
  distinct() %>%
  reframe(all_timepoints = str_c(time, collapse = "_")) %>%
  filter(all_timepoints != "T0") %>% # remove 5 fragments used for making doses
  filter(str_detect(all_timepoints, "T3") & str_detect(all_timepoints, "T7")) %>% #removes 2 T7 only and 2 T0/T7
  pull(genotype)

model_samples <- filter(metadata, !str_detect(tank, 'homo|HOMO')) %>%
  filter(genotype %in% longitudinal_genos) %>% #genos present in T3 and T7
  filter(!(exposure == "H" & final_disease_state == "D")) %>%
  pull(sample_id)

#roughly check evenness across groups
melted_ps %>% 
  filter(Sample %in% model_samples) %>% 
  select(-c(OTU, Abundance)) %>% 
  select(!colnames(updated_taxonomy %>% select(-c(asv_id, contains("confidence"))))) %>%
  distinct() %>%
  group_by(time, exposure, final_disease_state) %>%
  reframe(n = n())

#### Look at Sequences ####
colwells_of_interest <- c("ASV_144", "ASV_148", "ASV_274", "ASV_88", "ASV_939")

colwell_seqs <- sequences[colwells_of_interest]

outside_colwells <- c("TACGGAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCGTAGGCGGTTTGATAAGCCAGATGTGAAATCCCGGGGCTTAACCTCGGAACTGCATTTGGAACTGTTTGACTAGAGTACTGTAGAGGGTGGTGGAATTTCCAGTGTAGCGGTGAAATGCGTAGAGATTGGAAGGAACATCAGTGGCGAAGGCGGCCACCTGGACAGATACTGACGCTGAGGCACGAAAGCGTGGGGAGCGAACAGG",
                      "TACGGAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCGTAGGCGGATAGTTAAGCGAGATGTGAAATCCCGGGGCTCAACCTCGGAACTGCATTTCGAACTGGCTGTCTAGAGTCTTGTAGAGGGTGGTGGAATTTCCAGTGTAGCGGTGAAATGCGTAGAGATTGGAAGGAACATCAGTGGCGAAGGCGGCCACCTGGACAAAGACTGACGCTGAGGCACGAAAGCGTGGGGAGCGAACAGG")
names(outside_colwells) <- c("schul66", "schul65")
combined_colwells <- c(as.list(colwell_seqs), outside_colwells)

#write.fasta(as.list(combined_colwells), names = names(combined_colwells), open = "w", file.out = "colwells.fa")

cysteiniphilums <- full_taxonomy %>% filter(Genus == "Cysteiniphilum") %>% pull(asv_id)
cysteiniphilum_seqs <- sequences[cysteiniphilums]

sequences["ASV_65"]

#write.fasta(as.list(cysteiniphilum_seqs), names = names(cysteiniphilum_seqs), open = "w", file.out = "cysteiniphilums.fa")
library(seqinr)
write.fasta(as.list(sequences), names = names(sequences), file.out = "florida_ch1_microbes.fa")

temp_tax <- up_melt %>% as_tibble() %>% select(OTU, Domain:Species)

seq_lengths <- sequences %>% stack() %>% as_tibble() %>% 
  rename("asv_id" = "ind", "sequence" = "values") %>% 
  relocate(asv_id, .before = sequence) %>%
  mutate(seq_length = nchar(sequence)) %>%
  left_join(temp_tax, by = join_by("asv_id" == "OTU"))

#how many pass filtering
normalized_asv_counts %>%
  mutate(Genus = ifelse(Genus %in% c("Cysteiniphilum", "Vibrio", "Thalassotalea"), Genus, "Other"),
         Genus = factor(Genus, levels = c("Cysteiniphilum", "Thalassotalea", "Vibrio", "Other"))) %>%
  select(asv_id, Genus, Species) %>%
  distinct() %>% group_by(Genus) %>% reframe(n = n())

#how many initially
seq_lengths %>%
  mutate(Genus = ifelse(Genus %in% c("Cysteiniphilum", "Vibrio", "Thalassotalea"), Genus, "Other"),
         Genus = factor(Genus, levels = c("Cysteiniphilum", "Thalassotalea", "Vibrio", "Other"))) %>%
  distinct() %>% group_by(Genus) %>% reframe(n = n())


fl_seqs <- seq_lengths %>%
  mutate(Genus = ifelse(Genus %in% c("Cysteiniphilum", "Vibrio", "Thalassotalea"), Genus, "Other"),
         Genus = factor(Genus, levels = c("Cysteiniphilum", "Thalassotalea", "Vibrio", "Other"))) %>%
  filter(!is.na(Genus)) %>%
  distinct() %>%
  ggplot() +
  geom_histogram(aes(x = seq_length, fill = Genus), bins = 17) +
  theme_bw() +
  scale_fill_manual(values = c("maroon", "forestgreen", "purple", "gray70")) +
  ggtitle("Florida") + 
  xlim(380, 420)

max(seq_lengths$seq_length)

#### Alpha Posthocs ####
# posthoc_order <- c('T0.F.F', 'T3.D.D', 'T3.D.H', 'T3.H.H', 'T7.D.D', 'T7.D.H', 'T7.H.H')
# emmeans(model, ~treatment)
alpha_posthoc_time <- list('aquarium' = c(-1, 0, 0, 1/2, 0, 0, 1/2))

alpha_posthoc_exp <- list('exposure' = c(0, 1/4, 1/4, -1/2, 1/4, 1/4, -1/2))

alpha_posthoc_outc <- list('outcome' = c(0, 1/2, -1/4, -1/4, 1/2, -1/4, -1/4))

alpha_t0_t3 <- list('t0_t3' = c(-1, 1/3, 1/3, 1/3, 0, 0, 0))

alpha_t0_t7 <- list('t0_t7' = c(-1, 0, 0, 0, 1/3, 1/3, 1/3))

alpha_t3_t7 <- list('t3_t7' = c(0, -1/3, -1/3, -1/3, 1/3, 1/3, 1/3))



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

#### Alpha Diversity ####
# pre-filtering for abundance

alpha_table <- subset_samples(updated_microbiome_data, sample_names(updated_microbiome_data) %in% model_samples) %>%
  rarefy_even_depth(rngseed = 68748) %>% #rarefying reduces the data down to 4311 taxa
  microbiome::alpha(index = "all") %>%
  as_tibble(rownames = 'sample_id') %>%
  inner_join(metadata, by = 'sample_id') %>%
  mutate(fragment_id = str_c(str_replace_na(exposure, 'NA'), tank, genotype, sep = '_'))

mod_alpha_tab <- alpha_table %>%
  mutate(final_disease_state = ifelse(exposure == "F", "F", final_disease_state)) %>%
  mutate(treatment = str_c(time, exposure, final_disease_state, sep = '_')) %>%
  pivot_longer(cols = !c(colnames(metadata), "fragment_id", "treatment"),
               names_to = 'metric',
               values_to = 'alpha_div_value') %>%
  select(-c(susceptability, resistance, clone_group)) %>%
  mutate(tank_field = if_else(str_detect(treatment, 'F'), 'field', 'tank'), .after = final_disease_state) %>%
  nest_by(metric) %>%
  #partition(cluster) %>%
  mutate(fit_model(alpha_div_value ~ treatment + (1 | genotype) + #(1 | tank),
                     (0 + dummy(tank_field, c("tank")) | tank),
                   data, 
                   use_weights = FALSE),
         random_anova = list(rand(model)),
         process_model(model, re_model, random_anova),
         posthoc = list(run_posthoc(model, alpha_posthoc_categories))) %>%
  collect() %>%
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


#old way of doing alpha model, no contrasts integrated
# mod_alpha_tab <- alpha_table %>%
#   mutate(final_disease_state = ifelse(exposure == "F", "F", final_disease_state)) %>%
#   mutate(treatment = str_c(time, exposure, final_disease_state, sep = '_')) %>%
#   pivot_longer(cols = !c(colnames(metadata), "fragment_id", "treatment"),
#                names_to = 'metric',
#                values_to = 'alpha_div_value') %>%
#   select(-c(susceptability, resistance, clone_group)) %>%
#   mutate(tank_field = if_else(str_detect(treatment, 'F'), 'field', 'tank'), .after = final_disease_state) %>%
#   nest_by(metric) %>%
#   rowwise() %>%
#   mutate(alpha_model = list(lmer(alpha_div_value ~ treatment +
#                                       (1 | genotype) + (0 + dummy(tank_field, c("tank")) | tank),
#                                     data = data))) %>%
#   rowwise() %>%
#   mutate(p_value = anova(alpha_model) %>%
#            rownames_to_column(var = "sig_term") %>%
#            as_tibble() %>%
#            dplyr::rename("p_val" = `Pr(>F)`) %>%
#            pull(p_val)) %>%
#   ungroup() %>%
#   mutate(fdr_p_val = p.adjust(p_value, method = 'fdr')) %>%
#   filter(fdr_p_val < 0.05) %>%
#   mutate(alpha_type = ifelse(metric %in% c("chao1", "observed"), "richness", str_extract(metric, "[^_]+")))


#### Alpha Div Plots ####

#alpha div values by time
alpha_significant_models %>%
  mutate(metric = ifelse(str_detect(metric, "chao1"), "richness_chao1", metric)) %>%
  filter(metric %in% c("dominance_core_abundance", "richness_chao1")) %>%
  rowwise() %>%
  mutate(average_values = list(emmeans(model, ~treatment) %>%
                            broom::tidy(conf.int = TRUE) %>%
                            separate(treatment, into = c('time', 'exposure', 'final_disease_state')) %>%
                            group_by(time) %>%
                            reframe(ave_value = mean(estimate))
  )) %>%
  select(metric, average_values) %>%
  unnest(average_values)

alpha_graphs_manuscript <- alpha_significant_models %>%
  filter(metric %in% c("dominance_core_abundance", "chao1")) %>%
  #mutate(metric = ifelse(str_detect(metric, "chao1"), "richness_chao1", metric)) %>%
  #filter(metric %in% c("diversity_shannon", "dominance_core_abundance", "evenness_camargo", "richness_chao1")) %>%
  #mutate(for_manuscript = ifelse(metric %in% c("diversity_shannon", "dominance_core_abundance", "evenness_camargo", "richness_chao1"), TRUE, FALSE)) %>%
  rowwise() %>%
  mutate(plot_info = list(emmeans(model, ~treatment) %>%
                            broom::tidy(conf.int = TRUE) %>%
                            separate(treatment, into = c('time', 'exposure', 'final_disease_state')) %>%
                            mutate(graph_cat = ifelse(time == "T0", NA, 
                                                      paste(exposure, final_disease_state, sep = "_"))) %>%
                            {. ->> intermed } %>%
                            mutate(graph_cat = ifelse(time == "T0", "D_D", 
                                                      graph_cat)) %>%
                            dplyr::slice(1) %>%
                            rbind(intermed) %>%
                            mutate(graph_cat = ifelse(is.na(graph_cat), "D_H", 
                                                      graph_cat)) %>%
                            dplyr::slice(rep(1:2, 1)) %>%
                            rbind(intermed) %>%
                            mutate(graph_cat = ifelse(is.na(graph_cat), "H_H", 
                                                      graph_cat)) %>%
                            mutate(c_time = parse_number(time)) %>%
                            mutate(facet_lab = "Experimental") %>%
                            mutate(c_time = ifelse(time == "T0", c_time,
                                                   case_when(graph_cat == "D_D" ~ c_time + 0.30,
                                                             graph_cat == "D_H" ~ c_time,
                                                             graph_cat == "H_H" ~ c_time - 0.30))) %>%
                            mutate(graph_cat = factor(graph_cat, levels = c("D_D", "D_H", "H_H"), labels = c("D_D", "D_H", "H_H"))))) %>%
  rowwise() %>%
  mutate(plot = list(
    ggplot(data = plot_info, aes(x = c_time, y = estimate, ymin = conf.low, ymax = conf.high)) +
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
      
      
      (geom_point(data = (plot_info %>% filter(graph_cat == "dose" & final_disease_state == "na")), size = 3.7, aes(colour3 = exposure), shape = "diamond") %>%
         rename_geom_aes(new_aes = c("colour" = "colour3"))) +
      (geom_errorbar(data = (plot_info %>% filter(graph_cat == "dose" & final_disease_state == "na")), width = 0, aes(colour3 = exposure)) %>%
         rename_geom_aes(new_aes = c("colour" = "colour3"))) + 
      (geom_point(data = (plot_info %>% filter(graph_cat == "dose" & is.na(exposure))), size = 3, shape = 1, col = "black", alpha = 0)) +
      
      scale_color_manual(aesthetics = "colour1", values = c("D_H" = "#22A7B6", "D_D" = "#A70000"), guide = "legend", 
                         name = "Disease Exposed", breaks = c("D_H", "D_D"), labels = c("Healthy", "Diseased")) +
      scale_color_manual(aesthetics = "colour3", breaks = c("H", "D"), values = c("H" = "#17D0B4", "D" = "#e30e0e"), guide = "legend", 
                         name = "Doses", labels = c("H" = "Healthy", "D" = "Diseased")) +
      scale_shape_manual(values = c("D_H" = 17, "D_D" = 17, "H_H" = 16), guide = "none") +
      scale_color_manual(aesthetics = "colour2", values = c("H_H" = "#406F23"), guide = "legend", 
                         name = "Healthy Exposed", labels = c("Healthy")) +
      guides(colour1 = guide_legend(
        override.aes=list(linetype = c(1, 1), shape = c(17, 17))),
        colour2 = guide_legend(
          override.aes=list(linetype = c(6), shape = c(16))),
        colour3 = guide_legend(
          override.aes=list(linetype = c(0, 0)))) +
      scale_x_continuous(breaks=c(0, 3, 7)) +
      scale_linetype_manual(values = c("D_H" = 1, "D_D" = 1, "H_H" = 6), guide = "none") +
      theme_bw() +
      xlab("Time") +
      ylab(metric) +
      labs(title = metric)
  )) #%>%
  #group_by(for_manuscript) %>%
  #summarise(combo_plots = list(wrap_plots(plot) + plot_layout(guides = 'auto')))

sig_lines <- alpha_significant_models %>%
  select(metric, contains("t0"), contains("t3"), contains("t7")) %>% 
  select(metric, contains("fdr")) %>%
  filter(metric %in% c("dominance_core_abundance", "chao1")) %>%
  pivot_longer(-metric, names_to = c("time1", "time2", "sign"), values_to = "significance", 
               names_sep = "_", names_prefix = "fdr_") %>%
  mutate(pval_label = case_when(significance < 0.001 ~ "***",
                                significance < 0.01 ~ "**",
                                significance < 0.05 ~ "*",
                                TRUE ~ "")) %>%
  mutate(significance = ifelse(significance < 0.05, "yes", "no")) %>%
  mutate(time1 = parse_number(time1), time2 = parse_number(time2)) %>%
  select(-sign) %>%
  rowwise() %>%
  mutate(conf.low = 0, conf.high = 0, text_x = mean(c(time1, time2)),
         line_order = case_when(sum(time1, time2) == 7 ~ 3.5,
                                sum(time1, time2) == 10 ~ 2,
                                sum(time1, time2) == 3 ~ 1))

#CHAO1
chao1_plot <- alpha_graphs_manuscript$plot[[1]] + theme(legend.position="none") +
  ylab(str_wrap("Chao1 Richness Index", 15)) + labs(title = NULL) +
  xlab(NULL) +
  geom_segment(data = sig_lines %>% filter(metric == "chao1"), 
               aes(x = time1, xend = time2, y = 400 + line_order*26.667, 
                   yend = 400 + line_order*26.667,
                   alpha = significance), col = "gray20") +
  geom_text(data = sig_lines %>% filter(metric == "chao1"),
             aes(x = text_x, y = 400 + 10 + line_order*26.667,
                 label = pval_label), col = "gray20", size = 6) +
  scale_alpha_manual(values = c("yes" = 1, "no" = 0)) +
  theme(axis.text = element_text(size=11),
        axis.title = element_text(size=12))

#CORE ABUN
core_abun_plot <- alpha_graphs_manuscript$plot[[2]] + theme(legend.position="none") +
  ylab(str_wrap("Core abundance dominance index", 15)) + labs(title = NULL) +
  geom_segment(data = sig_lines %>% filter(metric == "dominance_core_abundance"), 
               aes(x = time1, xend = time2, y = 0.75 + line_order*0.05, 
                   yend = 0.75 + line_order*0.05,
                   alpha = significance), col = "gray20") +
  geom_text(data = sig_lines %>% filter(metric == "dominance_core_abundance"),
            aes(x = text_x, y = 0.75 + 0.01 + line_order*0.05,
                label = pval_label), col = "gray20", size = 6) +
  scale_alpha_manual(values = c("yes" = 1, "no" = 0)) +
  ylim(0, 1) +
  theme(axis.text = element_text(size=11),
        axis.title = element_text(size=12))


fancy_legend <- cowplot::get_legend(alpha_graphs_manuscript$plot[[2]])

both_alphas <- chao1_plot / core_abun_plot + plot_annotation(tag_levels = "A")

(both_alphas | fancy_legend) + 
  plot_layout(widths = c(4, 1)) + plot_annotation(tag_levels = list(c("A", "B")))
# 1000x600

#outdated
combined_alpha_plots <- ((alpha_graphs_manuscript$plot[[2]] + theme(legend.position="none") + 
       ylab(str_wrap("Shannon-Weiner diversity index", 15)) + labs(title = NULL) +
       xlab(NULL)) / #"Diversity"
    (alpha_graphs_manuscript$plot[[1]] + theme(legend.position="none") + 
       ylab(str_wrap("Chao1 Richness Index", 15)) + labs(title = NULL) +
       xlab(NULL)) / #"Richness"
    (alpha_graphs_manuscript$plot[[4]] + theme(legend.position="none") + 
       ylab(str_wrap("Camargo's evenness index", 15)) + labs(title = NULL) +
       xlab(NULL)) / #"Evenness"
    (alpha_graphs_manuscript$plot[[3]] + theme(legend.position="none") + 
       ylab(str_wrap("Core abundance dominance index", 15)) + labs(title = NULL)) 
    ) + #"Dominance"
  plot_annotation(tag_levels = "A")

#NEW just two metrics
combined_alpha_plots <- (
                           (alpha_graphs_manuscript$plot[[1]] + theme(legend.position="none") + 
                              ylab(str_wrap("Chao1 Richness Index", 15)) + labs(title = NULL) +
                              xlab(NULL)) / #"Richness"
                           (alpha_graphs_manuscript$plot[[3]] + theme(legend.position="none") + 
                              ylab(str_wrap("Core abundance dominance index", 15)) + labs(title = NULL)) 
) + #"Dominance"
  plot_annotation(tag_levels = "A")
    
  

(combined_alpha_plots | fancy_legend) + 
  plot_layout(widths = c(4, 1)) + plot_annotation(tag_levels = list(c("A", "B")))

alpha_graphs_manuscript$combo_plots[[2]]

group_by(alpha_type) %>%
  summarise(combo_plots = list(wrap_plots(plot) + plot_layout(guides = 'collect') & plot_annotation(title = alpha_type)))

alpha_graphs$combo_plots[[2]]


normalized_asv_counts %>%
  select(fragment_id) %>%
  distinct() %>%
  filter(!str_detect(fragment_id, "Field")) %>%
  nrow()
#108 fragments
  
normalized_asv_counts %>%
  select(exposure, fragment_id) %>%
  distinct() %>%
  filter(!str_detect(fragment_id, "Field")) %>%
  group_by(exposure) %>%
  reframe(frag_num = n())
#64 D exposed, 44 H exposed


# get alpha inc/dec stats for manuscript
alpha_significant_models %>%
  select(metric, contains("t0"), contains("t3"), contains("t7")) %>% 
  filter(metric %in% c("dominance_core_abundance", "chao1")) %>%
  pivot_longer(cols = !metric, names_to = "val", values_to = "num") %>%
  mutate(time = str_after_first(val, "_"), val = str_before_first(val, "_")) %>%
  pivot_wider(names_from = val, values_from = num)


#### Microshades Microbe Abundance ####

mdf_prep_test1 <- updated_microbiome_data %>%
  tax_glom("Genus") %>%
  psmelt() 

mdf_prep_test1 %>% select(Abundance, Order) %>% group_by(Order) %>% 
  reframe(tot_abun = sum(Abundance)) %>% arrange(desc(tot_abun))

mdf_processed_data <- mdf_prep_test1 %>%
  filter(Abundance > 0) %>%
  mutate(category = paste(time, exposure, final_disease_state, 
                          sep = "_")) %>%
  mutate(category = ifelse(str_detect(category, "F"), "T0", category)) %>%
  filter(!category %in% c("T3_H_D", "T7_H_D")) %>%
  mutate(time = ifelse(category %in% c("T0_D_D", "T0_H_H"), "Doses", time)) %>%
  mutate(category = ifelse(category == "T0_D_D", "Diseased", ifelse(category == "T0_H_H", "Healthy", category))) %>%
  mutate(category = ifelse(time == "T3", paste(time, exposure, sep = "_"), category)) %>%
  group_by(category) %>%
  mutate(total = sum(Abundance)) %>%
  ungroup() %>%
  select(-c(retain_sample)) %>%
  group_by(category, Genus) %>%
  reframe(Domain, Phylum, Class, Order, Family, time, total, rel_abun = sum(Abundance)/total) %>%
  distinct() %>%
  rename(Sample = category, Abundance = rel_abun) %>%
  mutate(Sample = factor(Sample, levels = c("Healthy", "Diseased", "T0", "T3_H", "T3_D", "T7_H_H", "T7_D_H", "T7_D_D"))) %>%
  as.data.frame()

## ORDER GENUS

color_objs_ordergenus <- create_color_dfs(mdf_processed_data, group_level = "Order", 
                                          selected_groups = c("Rickettsiales", "Alteromonadales", "Flavobacteriales", 
                                                              "Vibrionales",  "Verrucomicrobiales"), cvd = TRUE)
mdf_ordergenus <- color_objs_ordergenus$mdf
cdf_ordergenus <- color_objs_ordergenus$cdf

legend_ordergenus <-custom_legend(mdf_ordergenus, cdf_ordergenus, group_level = "Order")

plot_ordergenus_prelim <- plot_microshades(mdf_ordergenus, cdf_ordergenus) + 
  scale_y_continuous(labels = scales::percent, expand = expansion(0)) +
  facet_grid(cols = vars(time), scales = "free", space = "free") +
  theme_bw() +
  theme(legend.position = "none", plot.margin = margin(6,20,6,6)) +
  #labs(title = "Order Genus") +
  scale_x_discrete(name = element_blank(), labels = c("T0" = "Field", "T3_H" = "Healthy", "T3_D" =
                                "Disease", "T7_H_H" = "Healthy (Healthy)",
                              "T7_D_H" = "Disease (Healthy)", "T7_D_D" = "Disease (Diseased)"))

plot_grid(plot_ordergenus_prelim, legend_ordergenus,  rel_widths = c(1, .25))


#ORDER GENUS - SUSCEPTIBILITY

mdf_processed_suscep_data <- mdf_prep_test1 %>%
  filter(Abundance > 0) %>%
  mutate(category = paste(time, exposure, susceptability, 
                          sep = "_")) %>%
  mutate(category = ifelse(str_detect(category, "F"), paste("T0", susceptability, sep = "_"), category)) %>%
  mutate(time = ifelse(category %in% c("T0_D_NA", "T0_H_NA"), "Doses", time)) %>%
  mutate(category = ifelse(category == "T0_D_NA", "Diseased", ifelse(category == "T0_H_NA", "Healthy", category))) %>%
  group_by(category) %>%
  mutate(total = sum(Abundance)) %>%
  ungroup() %>%
  select(-c(retain_sample)) %>%
  group_by(category, Genus) %>%
  reframe(Domain, Phylum, Class, Order, Family, time, total, rel_abun = sum(Abundance)/total) %>%
  distinct() %>%
  rename(Sample = category, Abundance = rel_abun) %>%
  mutate(Sample = factor(Sample, levels = c("Healthy", "Diseased", "T0_S", "T0_R", "T3_H_S", "T3_H_R", "T3_D_S", "T3_D_R", "T7_H_S", "T7_H_R", "T7_D_S", "T7_D_R"))) %>%
  as.data.frame()

color_objs_suscep <- create_color_dfs(mdf_processed_suscep_data, group_level = "Order", 
                                          selected_groups = c("Rickettsiales", "Alteromonadales", "Oceanospirillales", 
                                                              "Francisellales",  "Flavobacteriales"), cvd = TRUE)
mdf_suscep <- color_objs_suscep$mdf
cdf_suscep <- color_objs_suscep$cdf

legend_suscep <-custom_legend(mdf_suscep, cdf_suscep, group_level = "Order")

plot_suscep_prelim <- plot_microshades(mdf_suscep, cdf_suscep) + 
  scale_y_continuous(labels = scales::percent, expand = expansion(0)) +
  facet_grid(cols = vars(time), scales = "free", space = "free") +
  theme_bw() +
  theme(legend.position = "none", plot.margin = margin(6,20,6,6)) +
  labs(title = "Order Genus - Susceptibility")

plot_grid(plot_suscep_prelim, legend_suscep,  rel_widths = c(1, .25))


#### Fantaxtic Relative Abundance ####

top_nested <- nested_top_taxa(updated_microbiome_data,
                              top_tax_level = "Order",
                              nested_tax_level = "Genus",
                              n_top_taxa = 5, 
                              n_nested_taxa = 5)
plot_nested_bar(ps_obj = top_nested$ps_obj,
                top_level = "Order",
                nested_level = "Genus")



#ordering
look <- melted_ps %>%
  group_by(Order) %>%
  reframe(tot = sum(Abundance)) %>%
  arrange(desc(tot)) %>%
  filter(!is.na(Order))

melted_ps %>% filter(Genus == "Cysteiniphilum") %>% select(Order, Family, Genus, Species)


top_asv$top_taxa %>% as_tibble() %>% select(Order) %>% distinct()
#custom

top_level <- "Order"
nested_level <- "Genus"
sample_order <- NULL

up_melt <- updated_microbiome_data %>% psmelt()

up_melt %>% as_tibble() %>% filter(Genus %in% c("MD3-55", "Rickettsia")) %>% select(Domain:Species) %>% distinct()

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
                                            Puniceicoccales = "#380135", ##111787
                                            Saprospirales = "#111787", #F7B10F
                                            Oceanospirillales = "#078090",
                                            Vibrionales = "#7C2C00",
                                            Thiotrichales = "#590404",
                                            Other = "gray75"))

# Generate the plot


mod_ggnested(psdf,
              aes_string(main_group = top_level,
                         sub_group = nested_level,
                         x = "Sample",
                         y = "Abundance"),
              main_palette = custom_palette,
         gradient_type = "tints") +
  scale_y_continuous(expand = c(0, 0)) +
  #theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90)) +
  theme_nested(theme_bw) + 
  geom_col(position = position_fill()) + 
  facet_grid(col = vars(facet_level), space = "free", scales = "free") + 
  scale_x_discrete(name = element_blank(), labels = c("Field" = "Field", "D Dose" = "Diseased", "H Dose" =
                                                        "Healthy", "T3 Diseased" = str_wrap("Disease- Exposed", 10), "T3_H_H" = str_wrap("Healthy- Exposed", 10),
                                                      "T7_H_H" = str_wrap("Healthy", 10),
                                                      "T7_D_H" = str_wrap("Disease- Exposed Healthy", 10), "T7_D_D" = "Diseased")) +
  guides(fill=guide_legend(title=substitute(bold(bd)~nb, list(bd = "Order", nb = "/ Genus")), nrow = 32),
         col=guide_legend(title=substitute(bold(bd)~nb, list(bd = "Order", nb = "/ Genus"))), nrow = 32) +
  ylab("Relative Abundance") +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=13),
        legend.title = element_text(size=13),
        strip.text = element_text(size=12))

#export 


fantaxtic_palette <- get_mod_ggnested_palette(psdf,
             aes_string(main_group = top_level,
                        sub_group = nested_level,
                        x = "Sample",
                        y = "Abundance"),
             main_palette = custom_palette,
             gradient_type = "tints")

#aquarickettsias 
psdf %>%
  filter(Sample %in% c("T7_D_D", "T7_D_H"), Genus == "MD3-55") %>%
  as_tibble() %>%
  #group_by(Sample) %>%
  summarise_at(vars(Abundance), list(t_test = ~ list(t.test(. ~ Sample)))) %>%
  rowwise() %>%
  mutate(t_stderr = t_test$stderr) %>%
  mutate(t_test = broom::tidy(t_test)) %>%
  unnest(t_test) %>%
  rename("T7_D_H" = "estimate1", "T7_D_D" = "estimate2", "DH_minus_DD" = "estimate") %>%
  select(-c(parameter, conf.low, conf.high, method, alternative))


#relative abundance
psdf %>%
  select(OTU, Genus, Sample, Abundance, old_sample) %>%
  group_by(Sample) %>%
  mutate(tot_abun = sum(Abundance)) %>%
  ungroup() %>%
  mutate(rel_abun = Abundance/tot_abun) %>%
  group_by(old_sample, Genus) %>%
  reframe(Sample, genus_rel_abun = sum(rel_abun)) %>%
  filter(Sample %in% c("T7_D_D", "T7_D_H"), Genus == "MD3-55") %>%
  summarise_at(vars(genus_rel_abun), list(t_test = ~ list(t.test(. ~ Sample)))) %>%
  rowwise() %>%
  mutate(t_stderr = t_test$stderr) %>%
  mutate(t_test = broom::tidy(t_test)) %>%
  unnest(t_test) %>%
  rename("T7_D_H" = "estimate1", "T7_D_D" = "estimate2", "DH_minus_DD" = "estimate") %>%
  select(-c(parameter, conf.low, conf.high, method, alternative))


psdf %>%
  select(OTU, Genus, Sample, Abundance, old_sample) %>%
  group_by(Sample) %>%
  mutate(tot_abun = sum(Abundance)) %>%
  ungroup() %>%
  mutate(rel_abun = Abundance/tot_abun) %>%
  group_by(old_sample, Genus) %>%
  reframe(Sample, genus_rel_abun = sum(rel_abun)) %>%
  filter(Sample %in% c("T3_H_H", "T3 Diseased"), Genus == "MD3-55") %>%
  summarise_at(vars(genus_rel_abun), list(t_test = ~ list(t.test(. ~ Sample)))) %>%
  rowwise() %>%
  mutate(t_stderr = t_test$stderr) %>%
  mutate(t_test = broom::tidy(t_test)) %>%
  unnest(t_test) %>%
  rename("T3_H_H" = "estimate1", "T3 Diseased" = "estimate2", "H_minus_D" = "estimate") %>%
  select(-c(parameter, conf.low, conf.high, method, alternative))


#### Make Venn showing ASVs to keep ####
otu_timepoint_presence <- melted_ps %>%
  mutate(across(c(exposure, final_disease_state), factor)) %>%
  filter(time %in% c('T3', 'T7') | (time == "T0" & tank == "HOMO")) %>%
  filter(Abundance > 0) %>%
  mutate(time = if_else(time == 'T0', exposure, time)) %>%
  #filter(OTU %in% bacterial_signature_asv$asv_id) %>%
  group_by(time, OTU) %>%
  summarise(n = sum(Abundance),
            .groups = 'drop') %>%
  pivot_wider(names_from = time, values_from = n, values_fill = 0L) %>%
  mutate(across(-OTU, ~. > 0))

ggvenn(otu_timepoint_presence, c('D', 'H', 'T3', 'T7')) + ggtitle("ASV Presence")

venn_all_times_and_doses <- melted_ps %>%
  mutate(across(c(exposure, final_disease_state), factor)) %>%
  filter(tank != "homogenate_fragment") %>%
  filter(Abundance > 0) %>%
  mutate(time = if_else(time == 'T0' & tank == "HOMO", exposure, time)) %>%
  #filter(OTU %in% venn_group) %>%
  group_by(time, OTU) %>%
  summarise(n = sum(Abundance),
            .groups = 'drop') %>%
  pivot_wider(names_from = time, values_from = n, values_fill = 0L) %>%
  mutate(across(-OTU, ~. > 0)) %>% 
  mutate(T0_H = ifelse(H | T0, TRUE, FALSE))

ggvenn(venn_all_times_and_doses, c('D', 'T0_H', 'T3', 'T7')) + ggtitle("Minimal Filtering")

upset(venn_all_times_and_doses %>% left_join(taxonomy_tibble, by = c("OTU" = "asv_names")), 
      c("T7", "T3", "T0", "D", "H"), 
      
      base_annotations=list(
        'Intersection size'=intersection_size(counts=T, text = aes(size = 6), fill = "slategray4")
      ),
      
      queries=list(upset_query(set='D', color="#DF0000", fill = "#DF0000"),
                   upset_query(set='T0', color="#D98EFF", fill = "#D98EFF"),
                   upset_query(set='T3', color="#B21BFF", fill = "#B21BFF"),
                   upset_query(set='T7', color="#650197", fill = "#650197"),
                   upset_query(set='H', color="#0FAB02", fill = "#0FAB02")),
      
      name='asv_names', width_ratio=0.1, min_size = 1, sort_sets = FALSE) + 
  ggtitle("Minimal Filtering")

otus_to_analyze <- filter(otu_timepoint_presence, 
                          (D & T3 & T7)) %>%
  pull(OTU)

not_t3t7_only <- filter(venn_all_times_and_doses, 
                    !c(!D & !H & !T0 & T3 & T7 & !T0_H)) %>%
  filter(!c(!D & !H & !T0 & !T3 & T7 & !T0_H)) %>%
  filter(!c(!D & !H & !T0 & T3 & !T7 & !T0_H)) %>%
  pull(OTU)
  
#### Normalize based on all samples & any other filtering ####
otu_tmm <- updated_microbiome_data %>%
  phyloseq_filter_prevalence(prev.trh = 0.2) %>%
  otu_table() %>% 
  t %>% #NOTE: *genus and family do not need the t but ASVs need the t*
  as.data.frame %>%
  as.matrix %>% 
  DGEList(remove.zeros = TRUE) %>%
  
  #Add any other filtering here
  filter_missingness(model_samples, 0.9) %>%
  filter_missing_groups(metadata, 1) %>%
 
  edgeR::calcNormFactors(method = 'TMMwsp') %>%
  filter_venn(not_t3t7_only) %>% #remove things that are only in T3 and T7
  #filter_venn(otus_to_analyze) %>% #only things in D and T3 and T7, aka potential pathogens
  filter_samples(model_samples) #%>% #remove samples not to be analyzed
  #filter_asv_meanCount(metadata, 100) #Remove ASVs with an average of less than N CPM per sample


venn_group <- otu_tmm %>%
  cpm(log = TRUE, prior.count = 0.5,
      normalized.lib.sizes = TRUE) %>%
  as_tibble(rownames = 'asv_id') %>% pull(asv_id)

#### Plot Group by number of ASVs = 0 ####
otu_tmm$counts %>%
  as_tibble(rownames = 'asv_id') %>%
  pivot_longer(cols = -asv_id,
               names_to = 'sample_id',
               values_to = 'n') %>%
  left_join(metadata, by = 'sample_id') %>%
  group_by(time, exposure, susceptability, asv_id) %>%
  summarise(pct_missing = sum(n == 0) / n(),
            .groups = 'drop_last') %>%
  summarise(mean_missing = mean(pct_missing),
            sd_missing = sd(pct_missing),
            .groups = 'drop') %>%
  ggplot(aes(x = time, y = mean_missing, ymin = mean_missing - sd_missing, 
             ymax = mean_missing + sd_missing, 
             colour = interaction(exposure, susceptability))) +
  geom_pointrange(position = position_dodge(0.5)) +
  labs(x = NULL,
       y = 'Average Percent of Samples ASVs are 0 counts from')


#### Variance weighting ####
param <- SnowParam(parallel::detectCores() - 1, "SOCK", progressbar = TRUE)
dream_weights_fullInteraction <- voomWithDreamWeights(counts = otu_tmm,
                                  formula = ~ model_comp + (1 | genotype) + (1 | tank),

                                  data = filter(metadata, !str_detect(tank, 'homo|HOMO')) %>%
                                    filter(genotype %in% longitudinal_genos) %>% #genos present in T3 and T7
                                    filter(!(exposure == "H" & final_disease_state == "D")) %>%
                                    arrange(sample_id) %>%
                                    mutate(model_comp = str_c(time, exposure, susceptability)) %>%
                                    column_to_rownames('sample_id'),
                                  BPPARAM = param,
                                  plot = TRUE)

#### ASV Modelling ####
full_data <- otu_tmm %>%
  cpm(log = TRUE, prior.count = 0.5,
      normalized.lib.sizes = TRUE) %>%
  as_tibble(rownames = 'asv_id') %>%
  pivot_longer(cols = -asv_id,
               names_to = 'sample_id',
               values_to = 'log2_cpm') %>%
  left_join(dream_weights_fullInteraction$weights %>%
              set_colnames(colnames(dream_weights_fullInteraction$E)) %>%
              set_rownames(rownames(dream_weights_fullInteraction$E)) %>%
              as_tibble(rownames = 'asv_id') %>%
              pivot_longer(cols = -asv_id,
                           names_to = 'sample_id',
                           values_to = 'weight'),
            by = c('asv_id', 'sample_id')) %>%
  left_join(as_tibble(otu_tmm$counts, rownames = 'asv_id') %>%
              pivot_longer(cols = -asv_id,
                           names_to = 'sample_id',
                           values_to = 'read_count'),
            by = c('asv_id', 'sample_id')) %>%
  
  left_join(metadata, 
            by = 'sample_id') %>%
  left_join(as_tibble(otu_tmm$samples, rownames = 'sample_id') %>%
              select(-group),
            by = 'sample_id') %>%
  left_join(tax_table(updated_microbiome_data) %>%
              as.data.frame() %>%
              as_tibble(rownames = 'asv_id'),
            by = c('asv_id'))


full_data %>% group_by(Order) %>% summarize(counts = sum(log2_cpm)) %>% arrange(desc(counts))
full_data %>% filter(Phylum == "Cyanobacteria")

#filtered for the 
write_csv(full_data, '../intermediate_files/fully_preprocessed_samples.csv.gz')


### Output ASV CPMs for doses ####
homogenate_data <- updated_microbiome_data %>%
  subset_samples(str_detect(tank, 'HOMO')) %>%
  otu_table %>%
  t() %>%
  as.data.frame %>%
  as.matrix %>% 
  DGEList(remove.zeros = FALSE) %>%
  cpm(log = TRUE, prior.count = 0.5,
      normalized.lib.sizes = TRUE) %>%
  as_tibble(rownames = 'asv_id') %>%
  pivot_longer(cols = -asv_id,
               names_to = 'sample_id',
               values_to = 'log2_cpm') %>%
  filter(asv_id %in% unique(full_data$asv_id)) %>%
  mutate(exposure = str_extract(sample_id, '[DH]'))
write_csv(homogenate_data, '../intermediate_files/homogenate_cpm.csv')

#### What's in the Homogenate Dose ####
homogenate_data %>% 
  left_join(tax_table(updated_microbiome_data) %>% 
              as.data.frame() %>%
              as_tibble(rownames = "asv_id"), by = join_by(asv_id)) %>%
  filter(log2_cpm > 5.274105) %>% #normalized zero
  select(-c(sample_id, log2_cpm, exposure)) %>%
  distinct() %>%
  group_by(Class) %>% #change taxa level here
  reframe(n = n()) %>%
  nrow()

#### Plot ASV number vs mean cpm ####
otu_tmm %>%
  cpm(log = TRUE, prior.count = 0.5,
      normalized.lib.sizes = TRUE) %>%
  as_tibble(rownames = 'asv_id') %>%
  pivot_longer(cols = -asv_id,
               names_to = 'sample_id',
               values_to = 'log2_cpm') %>%
  mutate(asv_number = str_extract(asv_id, '[0-9]+') %>% as.integer()) %>%
  group_by(asv_number) %>%
  summarise(log2_cpm_mean = mean(log2_cpm),
            log2_cpm_se = sd(log2_cpm) / sqrt(n())) %>%
  ggplot(aes(x = asv_number, y = log2_cpm_mean, 
             ymin = log2_cpm_mean - log2_cpm_se, ymax = log2_cpm_mean + log2_cpm_se)) +
  geom_pointrange() +
  geom_vline(xintercept = 200, colour = 'red')


otu_tmm %>%
  cpm(log = TRUE, prior.count = 0.5,
      normalized.lib.sizes = TRUE) %>%
  as_tibble(rownames = 'asv_id') %>%
  pivot_longer(cols = -asv_id,
               names_to = 'sample_id',
               values_to = 'log2_cpm') %>%
  mutate(asv_number = str_extract(asv_id, '[0-9]+') %>% as.integer()) %>%
  group_by(asv_number) %>%
  summarise(log2_cpm_mean = mean(log2_cpm),
            log2_cpm_se = sd(log2_cpm) / sqrt(n()),
            sum_log2_cpm = sum(log2_cpm)) %>%
  ggplot(aes(x = log2_cpm_mean, y = sum_log2_cpm)) +
  geom_point() 

otu_tmm %>%
  cpm(log = TRUE, prior.count = 0.5,
      normalized.lib.sizes = TRUE) %>%
  as_tibble(rownames = 'asv_id') %>%
  pivot_longer(cols = -asv_id,
               names_to = 'sample_id',
               values_to = 'log2_cpm') %>%
  mutate(asv_number = str_extract(asv_id, '[0-9]+') %>% as.integer()) %>%
  group_by(asv_number) %>%
  summarise(log2_cpm_mean = mean(log2_cpm),
            log2_cpm_se = sd(log2_cpm) / sqrt(n())) %>%
  filter(asv_number <= 200)

otu_tmm %>%
  cpm(log = FALSE, prior.count = 0,
      normalized.lib.sizes = FALSE) %>%
  as_tibble(rownames = 'asv_id') %>%
  pivot_longer(cols = -asv_id,
               names_to = 'sample_id',
               values_to = 'log2_cpm') %>%
  mutate(asv_number = str_extract(asv_id, '[0-9]+') %>% as.integer()) %>%
  group_by(asv_number) %>%
  summarise(mean_count = mean(log2_cpm),
            sqrt_sd = sqrt(sd(log2_cpm)),
            sum_count = sum(log2_cpm)) %>%
  ggplot(aes(x = mean_count, y = sqrt_sd, colour = asv_number < 200)) +
  geom_point() +
  geom_vline(xintercept = 100, colour = 'darkgreen') +
  scale_x_continuous(limits = c(0, 1000))



otu_tmm %>%
  cpm(log = FALSE, prior.count = 0.5,
      normalized.lib.sizes = FALSE) %>%
  as_tibble(rownames = 'asv_id') %>%
  pivot_longer(cols = -asv_id,
               names_to = 'sample_id',
               values_to = 'log2_cpm') %>%
  mutate(asv_number = str_extract(asv_id, '[0-9]+') %>% as.integer()) %>%
  group_by(asv_number) %>%
  summarise(log2_cpm_mean = mean(log2_cpm),
            sqrt_sd = sqrt(sd(log2_cpm))) %>%
  filter(log2_cpm_mean >= 100)
