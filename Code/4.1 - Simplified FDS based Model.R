# Making the model based on FDS instead of Resistance

setwd("/Users/emilytrytten/Desktop/Screenshots/Career/Vollmer Lab/GitHub/16S_Florida_Tank_Analysis/Code")

#### Libraries ####
library(vegan)
library(phyloseq)
library(EcolUtils)
library(formattable)
library(RColorBrewer)
library(magrittr)
library(lmerTest)
library(emmeans)
library(multidplyr)
library(ggupset)
library(qvalue)
library(patchwork)
library(relayer) #devtools::install_github("clauswilke/relayer")
library(ComplexUpset)
library(corrplot)
library(Hmisc)
library(broom.mixed)
library(tidyverse)

refit_models <- TRUE

#### Functions ####
cluster <- new_cluster(parallel::detectCores() - 1)
cluster_library(cluster, c('dplyr', 'lmerTest', 'emmeans', 'stringr', 'tidyr'))

fit_model <- function(formula, data, use_weights = TRUE){
  if(!use_weights){
    data$weights <- 1
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

# posthoc <- asv_models$posthoc[[1]]
process_postHoc <- function(posthoc){
  post_row <- as_tibble(posthoc) %>%
    dplyr::rename(tvalue = t.ratio,
                  pvalue = p.value) %>%
    pivot_wider(names_from = c('contrast'),
                values_from = c('estimate', 'SE', 'df', 'tvalue', 'pvalue'),
                names_vary = 'slowest')
  post_row
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

cluster_copy(cluster, c('fit_model', 'process_model', 'process_postHoc'))

# complicated functions for setting axis limits of faceted graphs
# https://stackoverflow.com/questions/63550588/ggplot2coord-cartesian-on-facets
UniquePanelCoords <- ggplot2::ggproto(
  "UniquePanelCoords", ggplot2::CoordCartesian,
  
  num_of_panels = 1,
  panel_counter = 1,
  panel_ranges = NULL,
  
  setup_layout = function(self, layout, params) {
    self$num_of_panels <- length(unique(layout$PANEL))
    self$panel_counter <- 1
    layout
  },
  
  setup_panel_params =  function(self, scale_x, scale_y, params = list()) {
    if (!is.null(self$panel_ranges) & length(self$panel_ranges) != self$num_of_panels)
      stop("Number of panel ranges does not equal the number supplied")
    
    train_cartesian <- function(scale, limits, name, given_range = NULL) {
      if (is.null(given_range)) {
        expansion <- ggplot2:::default_expansion(scale, expand = self$expand)
        range <- ggplot2:::expand_limits_scale(scale, expansion,
                                               coord_limits = self$limits[[name]])
      } else {
        range <- given_range
      }
      
      out <- list(
        ggplot2:::view_scale_primary(scale, limits, range),
        sec = ggplot2:::view_scale_secondary(scale, limits, range),
        arrange = scale$axis_order(),
        range = range
      )
      names(out) <- c(name, paste0(name, ".", names(out)[-1]))
      out
    }
    
    cur_panel_ranges <- self$panel_ranges[[self$panel_counter]]
    if (self$panel_counter < self$num_of_panels)
      self$panel_counter <- self$panel_counter + 1
    else
      self$panel_counter <- 1
    
    c(train_cartesian(scale_x, self$limits$x, "x", cur_panel_ranges$x),
      train_cartesian(scale_y, self$limits$y, "y", cur_panel_ranges$y))
  }
)

coord_panel_ranges <- function(panel_ranges, expand = TRUE, default = FALSE, clip = "on") {
  ggplot2::ggproto(NULL, UniquePanelCoords, panel_ranges = panel_ranges, 
                   expand = expand, default = default, clip = clip)
}

#
#### Data ####

#read in and format data
old_taxonomy_normalized_asv_counts <- full_data %>% #read_csv('../intermediate_files/fully_preprocessed_samples.csv.gz', show_col_types = FALSE) %>%
  mutate(time = factor(time, ordered = TRUE)) %>%
  mutate(final_disease_state = ifelse(time == "T0", "F", final_disease_state)) %>%
  filter(!c(exposure == "H" & final_disease_state == "D")) %>%
  mutate(treatment = str_c(time, exposure, final_disease_state, sep = '_'),
         time_exposure = str_c(time, exposure, sep = '_'),
         timeC = str_extract(time, '[0-9]+') %>% as.numeric,
         across(Domain:Species, str_replace_na)) %>%
  mutate(asv_number = str_extract(asv_id, '[0-9]+') %>% as.integer)

#read in unprocessed data to get old taxonomy info
microbiome_data <- read_rds("../intermediate_files/preprocess_microbiome.rds")
aggregation_level <- 'none' #or none
if(aggregation_level != 'none'){
  microbiome_data <- aggregate_taxa(microbiome_data, aggregation_level)
  taxa_names(microbiome_data) <- str_replace_all(taxa_names(microbiome_data), ' |-', '_')
} else {
  taxa_names(microbiome_data) <- str_c('ASV', 1:length(taxa_names(microbiome_data)), sep = '_')
}

#old taxonomy info
taxonomy_tibble <- tax_table(microbiome_data) %>% 
  as.data.frame %>%
  as_tibble(rownames = "asv_names")

#read in data
homogenate_data <- read_csv('../intermediate_files/homogenate_cpm.csv', show_col_types = FALSE)
updated_taxonomy <- read_csv('../intermediate_files/updated_taxonomy.csv')

#get old taxonomy info for ASVs that are all NA in updated taxonomy
all_na_in_new_tax <- updated_taxonomy %>% 
  filter(if_all(Domain:Species, ~is.na(.))) %>% 
  select(-c(Domain:Species)) %>%
  left_join(taxonomy_tibble, by = join_by("asv_id" == "asv_names"))

#combine the old taxonomy with the updated version
combined_taxonomy <- updated_taxonomy %>% 
  filter(!if_all(Domain:Species, ~is.na(.))) %>%
  full_join(all_na_in_new_tax) %>%
  mutate(Family = case_when(Genus == "Enhygromyxa" ~ "[[Nannocystaceae]]", #[[...]] to indicate manually added classifications
                            Genus == "Sedimenticola" ~ "[[incertae sedis]]",
                            TRUE ~ Family))

#update taxonomy info in data; is now fully prepared data for further analysis
normalized_asv_counts <- old_taxonomy_normalized_asv_counts %>%
  select(-c(colnames(taxonomy_tibble %>% select(-asv_names)))) %>%
  left_join(combined_taxonomy, by = join_by("asv_id"))

# posthoc_order <- c('T0.F.F', 'T3.D.D', 'T3.D.H', 'T3.H.H', 'T7.D.D', 'T7.D.H', 'T7.H.H')
simple_planned_posthocs <- list('aquarium' = c(-1, 1/6, 1/6, 1/6, 1/6, 1/6, 1/6),
                                'exposure' = c(0, 1/4, 1/4, -1/2, 1/4, 1/4, -1/2),
                                'outcome' = c(0, 1/2, -1/2, 0, 1/2, -1/2, 0))

emmeans(asv_models$model[[1]], ~treatment) %>%
  contrast(simple_planned_posthocs)


#### Make ASV Models ####

if(file.exists('../intermediate_files/mixed_model_results.rds.gz') & !refit_models){
  asv_models <- read_rds('../intermediate_files/mixed_model_results.rds.gz')
} else {
  # cluster_copy(cluster, c('posthoc_categories', 'run_posthoc'))
  cluster_copy(cluster, 'simple_planned_posthocs')
  
  asv_models <- normalized_asv_counts %>% 
    # filter(asv_id == 'ASV_65') %>%
    nest_by(across(c('asv_id', Domain:Species, Family_confidence:Species_confidence))) %>%
    partition(cluster) %>%
    mutate(fit_model(log2_cpm ~ treatment + (1 | genotype) + (1 | tank),
                     data, 
                     use_weights = FALSE),
           random_anova = list(rand(model)),
           process_model(model, re_model, random_anova),
           posthoc = list(emmeans(model, ~treatment) %>%
                            contrast(simple_planned_posthocs, 
                                     adjust = 'none') %>%
                            as_tibble)) %>%
    collect() %>%
    select(-re_model, -ends_with('global')) %>%
    ungroup %>% 
    p_adjust %>%
    relocate(anova_table, .after = model) %>% 
    relocate(posthoc, .after = random_anova)
  write_rds(asv_models, '../intermediate_files/mixed_model_results.rds.gz')
  write_csv(select(asv_models, -where(is.list)), '../intermediate_files/mixed_model_results.csv.gz')
}


#### Main Effects Upset ####
asv_models %>%
  select(asv_id, starts_with('fdr')) %>% 
  mutate(across(starts_with('fdr'), ~. < 0.05)) %>%
  
  pivot_longer(cols = -asv_id,
               names_to = c('term'),
               values_to = 'significance',
               names_prefix = 'fdr_') %>%
  filter(significance) %>%
  group_by(asv_id) %>%
  summarise(terms = list(term),
            .groups = 'drop') %>%
  
  ggplot(aes(x = terms)) +
  geom_bar() +
  scale_x_upset() +
  theme_classic() +
  theme_combmatrix(combmatrix.label.make_space = TRUE)


#### Significant Treatment Effect Upset ####
significant_models <- asv_models %>%
  filter(fdr_treatment < 0.05) %>%
  # slice(26) %>%
  rowwise %>%
  mutate(process_postHoc(posthoc)) %>%
  ungroup() %>%
  p_adjust(exclude_cols = c('treatment', 'tank', 'genotype'))


classified_asvs <- significant_models %>%
  select(asv_id, starts_with('estimate'), starts_with('fdr')) %>% # or qvalue
  select(-contains(c('treatment', 'tank', 'genotype'))) %>%
  mutate(across(starts_with('fdr'), ~. < 0.05),
         across(starts_with('estimate'), ~if_else(. < 0, -1L, 1L))) %>%
  
  pivot_longer(cols = -asv_id,
               names_to = c('.value', 'term'),
               names_pattern = '(.*)_(.*)') %>%
  rename(significance = fdr) %>%
  mutate(term = case_when(term == 'outcome' & estimate < 0 ~ 'HealthyOutcome',
                          term == 'outcome' & estimate > 0 ~ 'DiseaseOutcome',
                          
                          term == 'aquarium' & estimate < 0 ~ 'Field',
                          term == 'aquarium' & estimate > 0 ~ 'Aquaria',
                          
                          term == 'exposure' & estimate < 0 ~ 'HealthyExposed',
                          term == 'exposure' & estimate > 0 ~ 'DiseaseExposed'),
         .keep = 'unused') 


classified_asvs %>%
  filter(significance) %>%
  group_by(asv_id) %>%
  summarise(terms = list(term),
            .groups = 'drop') %>%
  
  ggplot(aes(x = terms)) +
  geom_bar() +
  scale_x_upset() +
  theme_classic() +
  theme_combmatrix(combmatrix.label.make_space = TRUE)



classified_asv_taxonomy <- select(significant_models, asv_id, Domain:Species) %>%
  left_join(classified_asvs %>%
              pivot_wider(names_from = term,
                          values_from = significance,
                          values_fill = FALSE),
            by = 'asv_id')
write_csv(classified_asv_taxonomy, '../intermediate_files/classified_significant_asvs.csv.gz')


# merged_bacterial_signature_asv <- bacterial_signature_asv %>%
#   group_by(asv_id) %>%
#   reframe(signatures = str_c(signatures, collapse = ', '), significance) %>%
#   mutate(signatures = case_when(signatures == "early_pathogen, continuous_pathogen, late_pathogen" ~ "continuous_pathogen",
#                                 signatures == "early_pathogen, early_opportunist" ~ "early_pathogen",
#                                 signatures == "late_pathogen, late_opportunist" ~ "late_pathogen",
#                                 signatures == "probiotic_t7_strict, late_opportunist" ~ "late_probiotic",
#                                 signatures == "probiotic_t7_strict" ~ "late_probiotic",
#                                 signatures == "early_pathogen, late_opportunist" ~ "early_pathogen",
#                                 signatures == "probiotic_t7_strict, early_opportunist, continuous_opportunist, late_opportunist" ~ "late_probiotic",
#                                 TRUE ~ signatures)) %>%
#   distinct()

grouped_bacterial_signature_asv <- bacterial_signature_asv %>%
  group_by(asv_id) %>%
  #reframe(signatures = str_c(signatures, collapse = ', '), significance) %>%
  mutate(signatures = case_when(str_detect(signatures, "pathogen") ~ "pathogen",
                                str_detect(signatures, "opportunist") ~ "opportunist",
                                str_detect(signatures, "crasher") ~ "crasher",
                                str_detect(signatures, "commensalist") ~ "commensalist",
                                str_detect(signatures, "probiotic") ~ "probiotic")) %>%
  distinct()

grouped_bacterial_signature_asv %>%
  group_by(asv_id) %>%
  summarise(terms = list(unique(signatures)),
            .groups = 'drop') %>%
  
  ggplot(aes(x = terms)) +
  geom_bar() +
  scale_x_upset() +
  theme_classic() +
  theme_combmatrix(combmatrix.label.make_space = TRUE)

grouped_comp_upset_prep <- grouped_bacterial_signature_asv %>%
  group_by(asv_id) %>%
  pivot_wider(names_from = signatures, values_from = significance) %>%
  summarise(commensalist = paste0(na.omit(commensalist), collapse = ","),
            crasher = paste0(na.omit(crasher), collapse = ","),
            pathogen = paste0(na.omit(pathogen), collapse = ","),
            opportunist = paste0(na.omit(opportunist), collapse = ","),
            probiotic = paste0(na.omit(probiotic), collapse = ",")) %>%
  mutate(across(commensalist:probiotic, ~ifelse(. == "", FALSE, .))) %>%
  mutate(across(commensalist:probiotic, ~as.logical(.))) %>%
  ungroup()


upset(grouped_comp_upset_prep,
  colnames(grouped_comp_upset_prep %>% select(-asv_id)), 
  matrix=(
    intersection_matrix(geom=geom_point(shape = "circle filled", size=3, stroke = 0.35, color = "gray20"))
    + scale_color_manual(
      values=c(
        "pathogen" = "#FF1E1E",
        "opportunist" = "deepskyblue3",
        "crasher" = "#FF9A00",
        "commensalist"= "#9ab010",
        "probiotic"= "#28de9e"
      )
    )
  ),
  queries=list(
    upset_query(set = "pathogen", fill = "#FF1E1E"),
    upset_query(set = "opportunist", fill = "deepskyblue3"),
    upset_query(set = "crasher", fill = "#FF9A00"),
    upset_query(set = "commensalist", fill = "#9ab010"),
    upset_query(set = "probiotic", fill = "#28de9e")
  ),
  name='ASVs', width_ratio=0.1, min_size = 0, sort_sets = FALSE,
  stripes = c(rep(c("gray91", "gray97"), 3))) +
  #stripes = c(rep(c("gray78", "gray87"), 2), "gray78", rep(c("gray97", "gray91"), 4))) +
  ggtitle("Grouped Bacterial Strategies") 


bacterial_signature_asv %>%
  ungroup() %>%
  group_by(signatures) %>%
  summarize(count = n()) %>%
  arrange(desc(count))


#venn diagram of overlaps in bacterial strats

bac_strat_venn_prep <- bacterial_signature_asv %>%
  select(-direction) %>%
  pivot_wider(names_from = signatures, values_from = significance) %>%
  mutate(across(-asv_id, ~ifelse(is.na(.), FALSE, .)))

#crashers
bac_strat_venn_prep %>%
  select(asv_id, contains("crasher")) %>%
  select(-rev_crasher_tank_assoc) %>%
  ggvenn()

#probiotics
bac_strat_venn_prep %>%
  select(asv_id, contains("probiotic")) %>%
  ggvenn()

#commensalists
bac_strat_venn_prep %>%
  select(asv_id, contains("commensalist"), contains("rev")) %>%
  ggvenn()

#### Bacterial Strategies Upset ####

#prep for comp upset
comp_upset_bac_strat <- bacterial_signature_asv %>% 
  select(-c(direction)) %>%
  pivot_wider(names_from = signatures, values_from = significance) %>%
  mutate(across(-asv_id, ~ifelse(is.na(.), FALSE, .))) %>% 
  left_join(combined_taxonomy, by = c('asv_id'))

### Upset - bacterial strategies and ASV Presence Data

upset(
  comp_upset_bac_strat %>% left_join(venn_all_times_and_doses, by = c("asv_id" = "OTU")) %>%
    mutate(across(D:T0_H, ~ifelse(is.na(.), FALSE, .))),
  colnames(comp_upset_bac_strat %>% left_join(venn_all_times_and_doses, by = c("asv_id" = "OTU")) %>% #from script 2 - asv_sample_filtering
             select(T0, D, H, contains(c('commensalist', 'crasher', 'probiotic', 'opportunist', 'pathogen')))), 
  
  base_annotations=list(
    'Intersection size'=intersection_size(counts=T, text = aes(size = 6.5),
                                          bar_number_threshold = 25,
                                          mapping=aes(fill=Family, col = Family, label = parse_number(asv_id)), col = "gray10"
    ) +
      geom_text(size = 4, position = position_stack(vjust = 0.5), col = "gray10")
  ),
  matrix=(
    intersection_matrix(geom=geom_point(shape = "circle filled", size=3, stroke = 0.35, color = "gray20"))
    + scale_color_manual(
      values=c(
        #"T7" = "#650197",
        #"T3" = "#B21BFF",
        "T0" = "#B21BFF",
        "D" = "#AE0404",
        "H" = "#0FAB02",
        "late_pathogen" = "#FF1E1E",
        #"early_pathogen" = "#FF9797",
        #"continuous_pathogen" = "maroon",
        #"continuous_opportunist" = "royalblue4",
        "late_opportunist" = "deepskyblue3",
        "early_opportunist" = "lightskyblue",
        "crasher_BMO_dysbiosis" = "#FF9A00",
        "crasher_tank_effect" = "#BD5E03",
        "crasher_dysbiosis" = "gold2",
        "probiotic_t7_general" = "#9effdd",
        "commensalist_t7"= "#9ab010",
        "probiotic_t7_responder"= "#0a8f60",
        "probiotic_t7_always"= "#28de9e"
      )
    )
  ),
  queries=list(
    #upset_query(set = "T7", fill = "#650197"),
    #upset_query(set = "T3", fill = "#B21BFF"), "#D98EFF"
    upset_query(set = "T0", fill = "#B21BFF"),
    upset_query(set = "D", fill = "#AE0404"),
    upset_query(set = "H", fill = "#0FAB02"),
    upset_query(set = "late_pathogen", fill = "#FF1E1E"),
    #upset_query(set = "early_pathogen", fill = "#FF9797"),
    #upset_query(set = "continuous_pathogen", fill = "maroon"),
    #upset_query(set = "continuous_opportunist", fill = "royalblue4"),
    upset_query(set = "late_opportunist", fill = "deepskyblue3"),
    upset_query(set = "early_opportunist", fill = "lightskyblue"),
    upset_query(set = "crasher_BMO_dysbiosis", fill = "#FF9A00"),
    upset_query(set = "crasher_tank_effect", fill = "#BD5E03"),
    upset_query(set = "crasher_dysbiosis", fill = "gold2"),
    upset_query(set = "probiotic_t7_general", fill = "#9effdd"),
    upset_query(set = "commensalist_t7", fill = "#9ab010"),
    upset_query(set = "probiotic_t7_responder", fill = "#0a8f60"),
    upset_query(set = "probiotic_t7_always", fill = "#28de9e")
  ),
  name='ASVs', width_ratio=0.1, min_size = 0, sort_sets = FALSE,
  stripes = c("#F8EAFF", "#FFE4E4", "#E7FFE5", rep(c("gray91", "gray97"), 6))) +
  #stripes = c(rep(c("gray78", "gray87"), 2), "gray78", rep(c("gray97", "gray91"), 4))) +
  ggtitle("Bacterial Strategies/Presence Data") 


### Upset - just bacterial strategies

upset(
  comp_upset_bac_strat %>% left_join(venn_all_times_and_doses, by = c("asv_id" = "OTU")),
  colnames(comp_upset_bac_strat %>% left_join(venn_all_times_and_doses, by = c("asv_id" = "OTU")) %>%
             select(contains(c('commensalist', 'crasher', 'probiotic', 'opportunist', 'pathogen')))), 
  
  base_annotations=list(
    'Intersection size'=intersection_size(counts=T, text = aes(size = 6.5),
                                          bar_number_threshold = 25,
                                          mapping=aes(fill=Family, col = Family, label = parse_number(asv_id)), col = "gray10"
    ) +
      geom_text(size = 4, position = position_stack(vjust = 0.5), col = "gray10")
  ),
  matrix=(
    intersection_matrix(geom=geom_point(shape = "circle filled", size=3, stroke = 0.35, color = "gray20"))
    + scale_color_manual(
      values=c(
        "late_pathogen" = "#FF1E1E",
        "late_opportunist" = "deepskyblue3",
        "early_opportunist" = "lightskyblue",
        "crasher_BMO_dysbiosis" = "#FF9A00",
        "crasher_tank_effect" = "#BD5E03",
        "crasher_dysbiosis" = "gold2",
        "probiotic_t7_general" = "#9effdd",
        "commensalist_t7"= "#9ab010",
        "probiotic_t7_responder"= "#0a8f60",
        "probiotic_t7_always"= "#28de9e"
      )
    )
  ),
  queries=list(
    upset_query(set = "late_pathogen", fill = "#FF1E1E"),
    upset_query(set = "late_opportunist", fill = "deepskyblue3"),
    upset_query(set = "early_opportunist", fill = "lightskyblue"),
    upset_query(set = "crasher_BMO_dysbiosis", fill = "#FF9A00"),
    upset_query(set = "crasher_tank_effect", fill = "#BD5E03"),
    upset_query(set = "crasher_dysbiosis", fill = "gold2"),
    upset_query(set = "probiotic_t7_general", fill = "#9effdd"),
    upset_query(set = "commensalist_t7", fill = "#9ab010"),
    upset_query(set = "probiotic_t7_responder", fill = "#0a8f60"),
    upset_query(set = "probiotic_t7_always", fill = "#28de9e")
  ),
  name='ASVs', width_ratio=0.1, min_size = 0, sort_sets = FALSE,
  stripes = c(rep(c("gray91", "gray97"), 4))) +
  #stripes = c(rep(c("gray78", "gray87"), 2), "gray78", rep(c("gray97", "gray91"), 4))) +
  ggtitle("Bacterial Strategies") 



#### Fancy Plots ####

## Prep Data for Homogenate Dose Panel

homogenate_models <- homogenate_data %>%
  nest_by(asv_id) %>%
  mutate(em_model = list(emmeans(lm(log2_cpm ~ exposure, data = data),
                                 ~exposure))) %>%
  mutate(dose_pairwise = em_model %>%
           contrast('pairwise', adjust = 'none') %>%
           broom::tidy(conf.int = TRUE) %>% pull(p.value)) %>%
  ungroup %>%
  mutate(dose_pairwise_adj = p.adjust(dose_pairwise, method = 'fdr')) %>%
  rowwise(asv_id) %>%
  mutate(homogenate_pred = list(em_model %>%
                                  broom::tidy(conf.int = TRUE) %>% mutate(time = "Dose", 
                                                                          graph_cat = "dose", c_time = ifelse(exposure == "D", -1.2, -1.5), final_disease_state = "na",
                                                                          facet_lab = "Doses") %>% 
                                  {. ->> set_one } %>% #save data to this var
                                  mutate(std.error = NA, df = NA, conf.low = NA, conf.high = NA, statistic = NA, p.value = dose_pairwise,
                                         c_time = c(-1.8, -0.9), final_disease_state = NA, exposure = NA) %>% #dummy set to change size of Dose Facet
                                  rbind(set_one) %>%
                                  {. ->> set_two } %>%
                                  arrange(desc(estimate)) %>%
                                  dplyr::slice(1) %>%
                                  mutate(exposure = "p_val", c_time = -1.35, estimate = estimate + 3) %>%
                                  rbind(set_two)
                                
  )) #recombine them

#set up zero line

add_zero_lines <- homogenate_models %>% ungroup() %>% select(homogenate_pred) %>% unnest(homogenate_pred) %>% 
  dplyr::slice(1:4) %>% mutate(across(everything(), ~NA)) %>% 
  mutate(facet_lab = c("Doses", "Experimental", "Doses", "Experimental"), c_time = c(-5, -5, 10, 10), 
         estimate = c(5.28, 5.294201, 5.28, 5.294201)) %>%
  mutate(sig_p = NA, time = "zero_lines")


## Make Fancy Plots

#problems:

# why is ASV 84 have a value below 5.29 for T3_d_d even though none of the log2cpm are less than 5.29

#Caused by warning in `.qf.non0()`:
#! Negative variance estimate obtained!


the_plots <- bacterial_signature_asv %>%  #merged_bacterial_signature_asv %>%
  group_by(asv_id) %>%
  mutate(grouped_signatures = case_when(str_detect(signatures, "pathogen") ~ "pathogen",
                                str_detect(signatures, "opportunist") ~ "opportunist",
                                str_detect(signatures, "crasher") ~ "crasher",
                                str_detect(signatures, "commensalist") ~ "commensalist",
                                str_detect(signatures, "probiotic") ~ "probiotic")) %>%
  group_by(asv_id, grouped_signatures) %>%
  reframe(signatures = str_c(signatures, collapse = ', '), significance) %>%
  distinct() %>%
  inner_join(significant_models,
             by = 'asv_id') %>%
  rowwise %>%
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
  left_join(homogenate_models, by = join_by(asv_id)) %>%
  mutate(plot_info = list(rbind(plot_info, homogenate_pred) %>% mutate(sig_p = ifelse(p.value < 0.05, "sig", "nonsig")))) %>%
  mutate(plot_info = list(rbind(plot_info, add_zero_lines))) %>%
  rowwise() %>%
  mutate(plot = list(
    ggplot(data = plot_info, aes(x = c_time, y = estimate, ymin = conf.low, ymax = conf.high)) +
      geom_line(data = plot_info %>% filter(time == "zero_lines"), col = "gray45") +
      
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
      
      geom_point(data = (plot_info %>% filter(exposure == "p_val")), size = 5, col = "gray20", pch = "*", aes(alpha = sig_p)) +
      
      (geom_point(data = (plot_info %>% filter(graph_cat == "dose" & final_disease_state == "na")), size = 3.7, aes(colour3 = exposure), shape = "diamond") %>%
         rename_geom_aes(new_aes = c("colour" = "colour3"))) +
      (geom_errorbar(data = (plot_info %>% filter(graph_cat == "dose" & final_disease_state == "na")), width = 0, aes(colour3 = exposure)) %>%
         rename_geom_aes(new_aes = c("colour" = "colour3"))) + 
      (geom_point(data = (plot_info %>% filter(graph_cat == "dose" & is.na(exposure))), size = 3, shape = 1, col = "black", alpha = 0)) +
      
      scale_color_manual(aesthetics = "colour1", values = c("D_H" = "#F75D5D", "D_D" = "#A70000"), guide = "legend", 
                         name = "Disease Exposed", breaks = c("D_H", "D_D"), labels = c("Healthy", "Diseased")) +
      scale_color_manual(aesthetics = "colour3", breaks = c("H", "D"), values = c("H" = "#16d9f0", "D" = "#e30e0e"), guide = "legend", 
                         name = "Doses", labels = c("H" = "Healthy", "D" = "Diseased")) +
      scale_shape_manual(values = c("D_H" = 16, "D_D" = 17, "H_H" = 16), guide = "none") +
      scale_color_manual(aesthetics = "colour2", values = c("H_H" = "#29A5B2"), guide = "legend", 
                         name = "Healthy Exposed", labels = c("Healthy")) +
      guides(colour1 = guide_legend(
        override.aes=list(linetype = c(1, 1), shape = c(16, 17))),
        colour2 = guide_legend(
          override.aes=list(linetype = c(6), shape = c(16))),
        colour3 = guide_legend(
          override.aes=list(linetype = c(0, 0)))) +
      scale_x_continuous(breaks=c(0, 3, 7)) +
      scale_alpha_manual(values = c("sig" = 1, "nonsig" = 0), guide = "none") +
      scale_linetype_manual(values = c("D_H" = 1, "D_D" = 1, "H_H" = 6), guide = "none") +
      theme_bw() +
      #theme(plot.title = element_text(face = "italic")) +
      xlab("Time") +
      ylab(expression("Normalized log"[2]*" (cpm)")) +
      labs(title = str_c(ifelse(is.na(Family), "NA", Family), " (", ifelse(is.na(Family_confidence), "NA", round(Family_confidence, digits = 0)), "%) ", ifelse(is.na(Genus), "NA", Genus)," (", ifelse(is.na(Genus_confidence), "NA", round(Genus_confidence, digits = 0)), "%) ", sep = ""),
           subtitle = str_c(ifelse(is.na(Species), "NA", Species), " (", ifelse(is.na(Species_confidence), "NA", round(Species_confidence, digits = 0)), "%); ", asv_id, "\n", signatures, sep = "")) +
      facet_grid(cols = vars(facet_lab), space = "free") + #scales = "free_x"
      coord_panel_ranges(panel_ranges = list(
        list(x=c(-1.8, -0.9)), # Dose Panel
        list(x=c(-0.75, 8)) # Experimental Panel
      ))
  )) %>%
  group_by(grouped_signatures) %>%
  summarise(combo_plots = list(wrap_plots(plot) + plot_layout(guides = 'collect') & plot_annotation(title = grouped_signatures)))


#view the plots
the_plots$combo_plots[[1]]

#### Bac Strat NMDS and PCOA ####

clean_bacstrats <- merged_bacterial_signature_asv %>%
  select(-significance) %>%
  rename("bacstrat" = "signatures")



asv_nmds <- full_data %>% 
  select(asv_id, sample_id, log2_cpm) %>%
  pivot_wider(names_from = sample_id, values_from = log2_cpm) %>%
  left_join(clean_bacstrats, by = join_by("asv_id")) %>%
  filter(!is.na(bacstrat)) %>%
  filter(!bacstrat %in% c("continuous_crasher", "late_crasher")) %>%
  select(-bacstrat) %>%
  column_to_rownames('asv_id') %>%
  t() %>%
  as.matrix()

test_nmds <- metaMDS(t(asv_nmds), distance = 'bray', k = 2, trymax = 100, autotransform = FALSE, verbose = TRUE)
sppscores(test_nmds) <- t(asv_nmds)

#shepard plot
plot(test_nmds$diss, test_nmds$dist)

# permanova

nmds_sample_metadat <- full_data %>% 
  left_join(clean_bacstrats, by = join_by("asv_id")) %>%
  filter(!is.na(bacstrat)) %>%
  filter(!bacstrat %in% c("continuous_crasher", "late_crasher")) %>%
  select(-c(bacstrat, weight, read_count, reads,lib.size, norm.factors)) %>%
  select(-(colnames(taxonomy_tibble %>% select(-asv_names)))) %>% #TAXONOMY TIBBLE USED HERE
  pivot_wider(names_from = asv_id, values_from = log2_cpm)

asv_nmds_no_colnames <- nmds_sample_metadat %>% 
  select(starts_with("ASV_")) %>%
  as.matrix()

asvs_dist <- vegdist(asv_nmds_no_colnames, method='bray')

adonis2(asvs_dist ~ time * (exposure + susceptability) + genotype + tank, data = nmds_sample_metadat, perm=999)

# genotype marginally sig, everything else sig

set.seed(0830)

adonis2(asvs_dist ~ genotype + time + tank, data = nmds_sample_metadat, perm=999)

#time and tank sig, genotype not

trial_nmds_sample_metadat <- nmds_sample_metadat %>%
  mutate(ex_res = paste(exposure, susceptability, sep = "_")) %>%
  mutate(treatment = str_c(time, exposure, susceptability, sep = '_'))

dispersion<-betadisper(asvs_dist, group=trial_nmds_sample_metadat$treatment)
permutest(dispersion)
anova(dispersion)
plot(dispersion, hull=FALSE, ellipse=TRUE, conf = 0.95)


test <- nmds_sample_metadat %>%
  select(sample_id, starts_with("ASV_")) %>%
  column_to_rownames("sample_id") %>%
  as.matrix()

test_dist <- vegdist(test, method='bray')

asvs_pcoa <- cmdscale (test_dist, eig = TRUE)
ordiplot (asvs_pcoa, display = 'sites', type = 'text')

testadon <- adonis2(test_dist ~ exposure*susceptability, data = nmds_sample_metadat, perm=999)

adonis2(test_dist ~ time*exposure*susceptability, data = nmds_sample_metadat, perm=999)


nmds_sample_metadat <- nmds_sample_metadat %>% mutate(across(time:susceptability, ~as.factor(.x))) %>% ungroup() 

nmds_sample_metadat$susceptability = as.factor(nmds_sample_metadat$susceptability)


ttti <- nmds_sample_metadat %>% mutate(treatment = paste(time, exposure, susceptability, sep = "_"), .after = susceptability) %>% as.data.frame()

ttti$treatment = as.factor(ttti$treatment)

tst <- adonis_pairwise(x = ttti, dd = test_dist, group.var = "treatment")

tst$Adonis.tab
tst$Betadisper.tab


bd_pcoa_data <- asvs_pcoa$points %>%
  as_tibble(rownames = 'sample_id') %>%
  left_join(full_data %>% select(sample_id, time, exposure, susceptability, tank, genotype) %>%
              mutate(treatment = paste(time, exposure, susceptability, sep = "_")) %>%
              mutate(shape_determ = paste(exposure, susceptability, sep = "_")) %>%
              distinct(), by = join_by("sample_id"))

bd_pcoa_centroid <- bd_pcoa_data %>% 
  group_by(treatment) %>%
  reframe(V1 = mean(V1), V2 = mean(V2))

bd_pcoa_arms <- bd_pcoa_data %>% 
  group_by(treatment) %>%
  mutate(V1_cent = mean(V1), V2_cent = mean(V2))

# better version of betadisper plot
bd_pcoa_arms %>%
  ggplot(aes(x = V1, y = V2, xend = V1_cent, yend = V2_cent)) +
  stat_ellipse(level = 0.95, aes(col = treatment)) +
  geom_point(aes(col = treatment, shape = susceptability), size = 2) + 
  geom_segment(aes(col = treatment)) +
  #scale_color_brewer(palette = "Paired", direction = -1) +
  scale_color_manual(values = c(
    "T3_H_S" = "#A6CEE3",
    "T3_H_R" = "#1F78B4",
    "T7_H_S" = "#B2DF8A",
    "T7_H_R" = "#33A02C",
    "T7_D_S" = "#FB9A99",
    "T7_D_R" = "#E31A1C",
    "T3_D_S" = "#FDBF6F",
    "T3_D_R" = "#FF7F00",
    "T0_F_S" = "#CAB2D6",
    "T0_F_R" = "#6A3D9A"
  )) +
  scale_shape_manual(values = c(
    "R" = 17,
    "S" = 16
  )) + 
  geom_text(data = bd_pcoa_centroid, 
            aes(label = treatment, x = V1, y = V2, fontface = "bold"), col = "black", inherit.aes = FALSE) +
  theme_bw() +
  labs(pch = "Susceptibility", col = "Treatment")

#plot centroid locations
bd_pcoa_data %>%
  ggplot(aes(x = V1, y = V2)) +
  stat_ellipse(geom = "polygon", alpha = 0.1, level = 0.95, aes(col = treatment, fill = treatment)) +
  geom_point(aes(col = treatment, shape = shape_determ), size = 4) + 
  geom_point(data = bd_pcoa_centroid, pch = "*")

#plot species or site alone
plot(test_nmds, "species")
orditorp(test_nmds, "species")


## NMDS with 95% CI around TIMEPOINTS

resist_cats <- full_data %>% select(sample_id, susceptability) %>% distinct()

pan3 <- scores(test_nmds)$site %>%
  as_tibble(rownames = 'asv_names') %>%
  left_join(clean_bacstrats, by = c('asv_names' = 'asv_id')) %>%
  filter(!is.na(bacstrat)) %>%
  mutate(shape_determination = case_when(
    str_detect(bacstrat, "path")  ~ "pathogen",
    str_detect(bacstrat, "opp")  ~ "opportunist",
    str_detect(bacstrat, "crash")  ~ "crasher",
    str_detect(bacstrat, "pro")  ~ "probiotic"
  )) %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  stat_ellipse(data = as_tibble(scores(test_nmds)$species, rownames = aggregation_level) %>% 
                 mutate(exp_conditions = str_split(none, "_")) %>%
                 rowwise() %>%
                 mutate(time = exp_conditions[1], exposure = exp_conditions[2], tank = exp_conditions[3]) %>%
                 left_join(resist_cats, by = join_by("none" == "sample_id")), 
               geom = "polygon", alpha = 0.1, level = 0.95, aes(col = exposure, fill = exposure)) +
  geom_point(data = as_tibble(scores(test_nmds)$species, rownames = aggregation_level) %>% 
               mutate(exp_conditions = str_split(none, "_")) %>%
               rowwise() %>%
               mutate(time = exp_conditions[1], exposure = exp_conditions[2], resist = exp_conditions[3]),
             size = 1, aes(color = time, shape = exposure)) +
  geom_point(aes(col = bacstrat, shape = shape_determination), size = 4) + #alpha = shape_determination
  geom_text(aes(label = parse_number(asv_names)), size = 3) +
  scale_shape_manual(values = c(
    "opportunist" = 16,
    "pathogen" = 17,
    "crasher" = 15,
    "probiotic" = 8,
    "H" = 6,
    "D" = 4,
    "REP1" = 3
  )) +
  scale_fill_manual(values = c(
    "T0" = "gray80",
    "T3" = "gray40",
    "T7" = "gray10",
    "D" = "red",
    "H" = "lightskyblue",
    "REP1" = "lightgreen",
    "S" = "#3DD8EA",
    "R" = "#048291"
  ), guide = "none") +
  scale_color_manual(values = c(
    "continuous_crasher" = "sandybrown",
    "continuous_pathogen" = "firebrick4",
    "late_crasher" = "darkorange2",
    "late_probiotic" = "aquamarine",
    "early_pathogen" = "hotpink",
    "late_pathogen" = "firebrick1",
    "early_opportunist" = "deepskyblue3",
    "late_opportunist" = "royalblue2",
    "T0" = "gray80",
    "T3" = "gray40",
    "T7" = "gray10",
    "D" = "red",
    "H" = "lightskyblue",
    "REP1" = "lightgreen",
    "S" = "#3DD8EA",
    "R" = "#048291"
  )) +
  guides(color = guide_legend(override.aes = list(fill = NA))) +
  theme_bw() +
  ggtitle("Exposure")

pan1 | pan2 | pan3


facet_nmds_t3 <- scores(test_nmds)$site %>%
  as_tibble(rownames = 'asv_names') %>%
  left_join(clean_bacstrats, by = c('asv_names' = 'asv_id')) %>%
  filter(!is.na(bacstrat) & !bacstrat %in% c("late_probiotic", "continuous_pathogen")) %>%
  mutate(shape_determination = case_when(
    str_detect(bacstrat, "path")  ~ "pathogen",
    str_detect(bacstrat, "opp")  ~ "opportunist",
    str_detect(bacstrat, "crash")  ~ "crasher",
    str_detect(bacstrat, "pro")  ~ "probiotic"
  )) %>%
  ggplot(aes(x = -NMDS1, y = NMDS2)) +
  stat_ellipse(data = as_tibble(scores(test_nmds)$species, rownames = "sample_id") %>% 
                 left_join((full_data %>% select(sample_id, susceptability) %>% distinct()), by = join_by("sample_id")) %>%
                 mutate(exp_conditions = str_split(sample_id, "_")) %>%
                 rowwise() %>%
                 mutate(time = exp_conditions[1], exposure = exp_conditions[2], tank = exp_conditions[3]) %>%
                 mutate(exp_res = paste(exposure, susceptability, sep = "_")) %>%
                 filter(time == "T3"), 
               geom = "polygon", alpha = 0.1, level = 0.95, aes(col = exp_res, fill = exp_res, linetype = exp_res)) +
  geom_point(data = as_tibble(scores(test_nmds)$species, rownames = "sample_id") %>% 
               left_join((full_data %>% select(sample_id, susceptability) %>% distinct()), by = join_by("sample_id")) %>%
               mutate(exp_conditions = str_split(sample_id, "_")) %>%
               rowwise() %>%
               mutate(time = exp_conditions[1], exposure = exp_conditions[2], tank = exp_conditions[3]) %>%
               mutate(exp_res = paste(exposure, susceptability, sep = "_")),
             size = 1, aes(color = time, shape = exposure)) +
  geom_point(data = as_tibble(scores(test_nmds)$species, rownames = "sample_id") %>% 
               left_join((full_data %>% select(sample_id, susceptability) %>% distinct()), by = join_by("sample_id")) %>%
               mutate(exp_conditions = str_split(sample_id, "_")) %>%
               rowwise() %>%
               mutate(time = exp_conditions[1], exposure = exp_conditions[2], tank = exp_conditions[3]) %>%
               mutate(exp_res = paste(exposure, susceptability, sep = "_")) %>%
               filter(time == "T3"),
             size = 1, aes(color = exp_res, shape = exposure)) +
  geom_point(aes(col = bacstrat, shape = shape_determination), size = 4) + #alpha = shape_determination
  geom_text(aes(label = parse_number(asv_names)), size = 3) +
  scale_linetype_manual(values = c(
    "D_S" = "dashed",
    "H_S" = "dashed",
    "D_R" = "solid",
    "H_R" = "solid"
  ), guide = "none") +
  scale_shape_manual(values = c(
    "opportunist" = 16,
    "pathogen" = 17,
    "crasher" = 15,
    "probiotic" = 18,
    "H" = 6,
    "D" = 4,
    "REP1" = 7
  )) +
  scale_fill_manual(values = c(
    "T0" = "gray80",
    "T3" = "gray40",
    "T7" = "gray10",
    "D_S" = "#F75D5D",
    "H_S" = "#3DD8EA",
    "D_R" = "#A70000",
    "H_R" = "#048291"
  ), guide = "none") +
  scale_color_manual(values = c(
    "continuous_crasher" = "sandybrown",
    "continuous_pathogen" = "firebrick4",
    "late_crasher" = "darkorange2",
    "late_probiotic" = "aquamarine2",
    "early_pathogen" = "hotpink",
    "late_pathogen" = "firebrick1",
    "early_opportunist" = "deepskyblue3",
    "late_opportunist" = "royalblue2",
    "T0" = "gray80",
    "T3" = "gray55",
    "T7" = "gray55",
    "D_S" = "#F75D5D",
    "H_S" = "#3DD8EA",
    "D_R" = "#A70000",
    "H_R" = "#048291"
  )) +
  guides(color = guide_legend(override.aes = list(fill = NA))) +
  theme_bw() +
  ylim(-0.26, 0.285) +
  xlim(-0.33, 0.25) +
  xlab("NMDS1") +
  ggtitle("T3")

facet_nmds_t7 <- scores(test_nmds)$site %>%
  as_tibble(rownames = 'asv_names') %>%
  left_join(clean_bacstrats, by = c('asv_names' = 'asv_id')) %>%
  filter(!is.na(bacstrat) & !bacstrat %in% c("late_probiotic", "continuous_pathogen")) %>%
  mutate(shape_determination = case_when(
    str_detect(bacstrat, "path")  ~ "pathogen",
    str_detect(bacstrat, "opp")  ~ "opportunist",
    str_detect(bacstrat, "crash")  ~ "crasher",
    str_detect(bacstrat, "pro")  ~ "probiotic"
  )) %>%
  ggplot(aes(x = -NMDS1, y = NMDS2)) +
  stat_ellipse(data = as_tibble(scores(test_nmds)$species, rownames = "sample_id") %>% 
                 left_join((full_data %>% select(sample_id, susceptability) %>% distinct()), by = join_by("sample_id")) %>%
                 mutate(exp_conditions = str_split(sample_id, "_")) %>%
                 rowwise() %>%
                 mutate(time = exp_conditions[1], exposure = exp_conditions[2], tank = exp_conditions[3]) %>%
                 mutate(exp_res = paste(exposure, susceptability, sep = "_")) %>%
                 filter(time == "T7"), 
               geom = "polygon", alpha = 0.1, level = 0.95, aes(col = exp_res, fill = exp_res, linetype = exp_res)) +
  geom_point(data = as_tibble(scores(test_nmds)$species, rownames = "sample_id") %>% 
               left_join((full_data %>% select(sample_id, susceptability) %>% distinct()), by = join_by("sample_id")) %>%
               mutate(exp_conditions = str_split(sample_id, "_")) %>%
               rowwise() %>%
               mutate(time = exp_conditions[1], exposure = exp_conditions[2], tank = exp_conditions[3]) %>%
               mutate(exp_res = paste(exposure, susceptability, sep = "_")),
             size = 1, aes(color = time, shape = exposure)) +
  geom_point(data = as_tibble(scores(test_nmds)$species, rownames = "sample_id") %>% 
               left_join((full_data %>% select(sample_id, susceptability) %>% distinct()), by = join_by("sample_id")) %>%
               mutate(exp_conditions = str_split(sample_id, "_")) %>%
               rowwise() %>%
               mutate(time = exp_conditions[1], exposure = exp_conditions[2], tank = exp_conditions[3]) %>%
               mutate(exp_res = paste(exposure, susceptability, sep = "_")) %>%
               filter(time == "T7"),
             size = 1, aes(color = exp_res, shape = exposure)) +
  geom_point(aes(col = bacstrat, shape = shape_determination), size = 4) + #alpha = shape_determination
  geom_text(aes(label = parse_number(asv_names)), size = 3) +
  scale_linetype_manual(values = c(
    "D_S" = "dashed",
    "H_S" = "dashed",
    "D_R" = "solid",
    "H_R" = "solid"
  ), guide = "none") +
  scale_shape_manual(values = c(
    "opportunist" = 16,
    "pathogen" = 17,
    "crasher" = 15,
    "probiotic" = 18,
    "H" = 6,
    "D" = 4,
    "REP1" = 7
  )) +
  scale_fill_manual(values = c(
    "T0" = "gray80",
    "T3" = "gray40",
    "T7" = "gray10",
    "D_S" = "#F75D5D",
    "H_S" = "#3DD8EA",
    "D_R" = "#A70000",
    "H_R" = "#048291"
  ), guide = "none") +
  scale_color_manual(values = c(
    "continuous_crasher" = "sandybrown",
    "continuous_pathogen" = "firebrick4",
    "late_crasher" = "darkorange2",
    "late_probiotic" = "aquamarine2",
    "early_pathogen" = "hotpink",
    "late_pathogen" = "firebrick1",
    "early_opportunist" = "deepskyblue3",
    "late_opportunist" = "royalblue2",
    "T0" = "gray80",
    "T3" = "gray55",
    "T7" = "gray55",
    "D_S" = "#F75D5D",
    "H_S" = "#3DD8EA",
    "D_R" = "#A70000",
    "H_R" = "#048291"
  )) +
  guides(color = guide_legend(override.aes = list(fill = NA))) +
  theme_bw() +
  ylim(-0.26, 0.285) +
  xlim(-0.33, 0.25) +
  xlab("NMDS1") +
  ggtitle("T7")

wrap_plots(facet_nmds_t3, facet_nmds_t7) + plot_layout(guides = 'collect')

## PCOA with 95% CI around BACTERIAL STRATEGIES

bacstrat_dist_mat <- vegdist(t(asv_nmds))
bacstrat_pcoa <- cmdscale (bacstrat_dist_mat, eig = TRUE)
ordiplot (bacstrat_pcoa, display = 'sites', type = 'text')

bacstrat_pcoa$points %>%
  as_tibble(rownames = 'asv_names') %>%
  left_join(clean_bacstrats, by = c('asv_names' = 'asv_id')) %>%
  filter(!is.na(bacstrat)) %>%
  mutate(shape_determination = case_when(
    str_detect(bacstrat, "path")  ~ "pathogen",
    str_detect(bacstrat, "opp")  ~ "opportunist",
    str_detect(bacstrat, "crash")  ~ "crasher",
    str_detect(bacstrat, "pro")  ~ "probiotic"
  )) %>%
  ggplot(aes(x = V1, y = V2)) +
  stat_ellipse(geom = "polygon", alpha = 0.1, level = 0.95, aes(col = bacstrat, fill = bacstrat)) +
  geom_point(aes(col = bacstrat, shape = shape_determination), size = 4) + 
  geom_text(aes(label = parse_number(asv_names)), size = 3) +
  scale_shape_manual(values = c(
    "opportunist" = 16,
    "pathogen" = 17,
    "crasher" = 15,
    "probiotic" = 8,
    "H" = 6,
    "D" = 4,
    "REP1" = 3
  )) +
  scale_fill_manual(values = c(
    "continuous_crasher" = "sandybrown",
    "continuous_pathogen" = "firebrick4",
    "late_crasher" = "darkorange2",
    "late_probiotic" = "aquamarine",
    "early_pathogen" = "hotpink",
    "late_pathogen" = "firebrick1",
    "early_opportunist" = "deepskyblue3",
    "late_opportunist" = "royalblue2"
  ), guide = "none") +
  scale_color_manual(values = c(
    "continuous_crasher" = "sandybrown",
    "continuous_pathogen" = "firebrick4",
    "late_crasher" = "darkorange2",
    "late_probiotic" = "aquamarine",
    "early_pathogen" = "hotpink",
    "late_pathogen" = "firebrick1",
    "early_opportunist" = "deepskyblue3",
    "late_opportunist" = "royalblue2",
    "T0" = "gray80",
    "T3" = "gray40",
    "T7" = "gray10"
  )) +
  guides(color = guide_legend(override.aes = list(fill = NA))) +
  theme_bw() #+
ylim(-0.3, 0.3) +
  xlim(-0.3, 0.25)


#### Emily Miscellaneous ####

## Access all ASVs in a given bacterial strategy:

asvs_by_signature <- clean_bacstrats %>%
  filter(bacstrat == "late_pathogen") %>%
  #filter(bacstrat %in% c("continuous_crasher", "late_crasher")) %>%
  pull(asv_id)


## heritability (random effect of genotype)

heritable_asvs <- asv_models %>%
  rowwise() %>%
  mutate(heritability_percent = model %>% broom.mixed::tidy("ran_pars") %>%
           mutate(variance = estimate^2, tot_var = sum(variance)) %>% 
           filter(group == "genotype") %>% mutate(herit = 100*variance/tot_var) %>%
           pull(herit), .after = asv_id) %>%
  filter(heritability_percent > 0)


#make formattable chart of heritability scores:
heritable_asvs %>% arrange(desc(heritability_percent)) %>% filter(asv_id %in% asvs_by_signature) %>%
  #filter(heritability_percent > 0.1) %>%
  filter(!asv_id %in% pathogens_not_in_dose) %>%
  select(Family, Genus, asv_id, heritability_percent) %>% 
  mutate(heritability_percent = round(heritability_percent, 3)) %>%
  #head(10) %>%
  formattable(align = c("c", "c", "c", "c"), list(
    area(col = 4) ~ color_tile("#FFEBF9", "hotpink")
  )) #%>%
export_formattable("../Figures/top10_heritability.png")


##logfold change between D and H

#make formattable chart showing logfold change of D compared to H at Dose and T7
significant_models %>%
  filter(asv_id %in% asvs_by_signature) %>%
  select(asv_id, data) %>%
  unnest(data) %>%
  select(asv_id, sample_id, log2_cpm, exposure, time) %>%
  rbind(homogenate_data %>%
          mutate(time = "dose")) %>%
  filter(time %in% c("dose", "T7")) %>%
  group_by(asv_id, time, exposure) %>%
  reframe(ave_val = mean(log2_cpm)) %>%
  ungroup() %>%
  group_by(asv_id, time) %>%
  reframe(logfold_change = ave_val[exposure == "D"] - ave_val[exposure == "H"]) %>%
  mutate(logfold_change = round(logfold_change, 2)) %>%
  pivot_wider(names_from = time, values_from = logfold_change) %>%
  select(asv_id, dose, T7) %>%
  filter(!is.na(dose) & !is.na(T7)) %>%
  rename("asv_id" = "ASV", "dose" = "Dose Logfold Change", "T7" = "T7 Logfold Change") %>%
  formattable(align = c("l", "c", "c"), list(
    area(col = c(2,3)) ~ color_tile("white", "firebrick1")
  )) %>%
  export_formattable("../Figures/put_pathogens_logfold_table.png")


## lm models predicting w T0
#TAKEAWAY: T0 abundance cannot predict disease resistance

t0_prediction_lms <- full_data %>%
  filter(time == "T0") %>%
  group_by(asv_id) %>%
  reframe(model = list(lm(log2_cpm ~ susceptability) %>% broom::tidy() %>% 
                         filter(term == "susceptabilityS"))) %>%
  unnest(model) %>%
  mutate(fdr_p = p.adjust(p.value, method = "fdr")) %>%
  filter(fdr_p < 0.05)

bacterial_signature_asv %>% filter(asv_id %in% t0_prediction_lms$asv_id)


## Correlation Matrices

#corr matrix for 8 putative pathogen candidates only

asv_corr <- full_data %>% 
  select(asv_id, sample_id, Family, log2_cpm) %>%
  pivot_wider(names_from = sample_id, values_from = log2_cpm) %>%
  mutate(combo_name = paste(Family, asv_id, sep = "_"), .after = asv_id) %>%
  left_join(clean_bacstrats, by = join_by("asv_id")) %>%
  filter(!is.na(bacstrat)) %>%
  filter(bacstrat %in% c('c("late_pathogen", "late_opportunist")', "late_pathogen")) %>%
  filter(!asv_id %in% pathogens_not_in_dose) %>%
  select(-c(bacstrat, asv_id, Family)) %>%
  column_to_rownames('combo_name') %>%
  t() %>%
  as.matrix()

asv_corr_mat <- rcorr(asv_corr, type = c("pearson","spearman"))



corrplot(asv_corr_mat$r, type="upper", order="hclust", 
         p.mat = asv_corr_mat$P, sig.level = 0.05, insig = "blank",
         col.lim = c(0,1), col = rep(brewer.pal(11,"Spectral"), 2), tl.col = "gray20")

brewer.pal(11,"Spectral")
# corr matrix for All Colwelliaceae Thalassotaleas

colwell_corr <- full_data %>% 
  filter(Genus == "Thalassotalea") %>%
  select(asv_id, sample_id, log2_cpm) %>%
  pivot_wider(names_from = sample_id, values_from = log2_cpm) %>%
  column_to_rownames('asv_id') %>%
  t() %>%
  as.matrix()

colwell_test <- rcorr(colwell_corr, type = c("pearson","spearman"))

corrplot(colwell_test$r, type="upper", order="hclust", 
         p.mat = colwell_test$P, sig.level = 0.05, insig = "blank",
         #col.lim = c(-0.2,1), col = c(rep("gray30", 10), rainbow(20)), tl.col = "gray20")
         col.lim = c(-0.2,1), col = c(rep("gray30", 6), brewer.pal(11,"Spectral")), tl.col = "gray20")


#late pathogens and opportunists


late_corr <- full_data %>% 
  select(asv_id, sample_id, Family, log2_cpm) %>%
  pivot_wider(names_from = sample_id, values_from = log2_cpm) %>%
  mutate(combo_name = paste(Family, parse_number(asv_id), sep = "_"), .after = asv_id) %>%
  left_join(clean_bacstrats, by = join_by("asv_id")) %>%
  filter(!is.na(bacstrat)) %>%
  filter(bacstrat %in% c('c("late_pathogen", "late_opportunist")', "late_pathogen", "late_opportunist")) %>%
  mutate(combo_name = ifelse(bacstrat == "late_opportunist", paste(combo_name, "_O"), paste(combo_name, "_P"))) %>%
  select(-c(bacstrat, asv_id, Family)) %>%
  column_to_rownames('combo_name') %>%
  t() %>%
  as.matrix()

late_test <- rcorr(late_corr, type = c("pearson","spearman"))

corrplot(late_test$r, type="upper", order="hclust", 
         p.mat = late_test$P, sig.level = 0.05, insig = "blank",
         #col.lim = c(-0.2,1), col = c(rep("gray30", 10), rainbow(20)), tl.col = "gray20")
         col.lim = c(-0.3,1), col = c(rep("gray30", 5), brewer.pal(11,"Spectral")), tl.col = "gray20")


#### Alpha Diversity ####

alpha_table <- microbiome::alpha(microbiome_data, index = "all") %>%
  as_tibble(rownames = 'sample_id') %>%
  inner_join(metadata, by = 'sample_id') %>%
  mutate(fragment_id = str_c(str_replace_na(exposure, 'NA'), tank, genotype, sep = '_'))

mod_alpha_tab <- alpha_table %>%
  filter(!tank %in% c("HOMO", "homogenate_fragment")) %>%
  mutate(treatment = str_c(time, exposure, susceptability, sep = '_')) %>%
  pivot_longer(cols = !c(colnames(metadata), "fragment_id", "treatment"),
               names_to = 'metric',
               values_to = 'alpha_div_value') %>%
  select(-c(final_disease_state, clone_group)) %>%
  nest_by(metric) %>%
  summarise(alpha_model = list(lmer(alpha_div_value ~ treatment + 
                                      (1 | genotype) + (1 | tank), data = data))) %>% 
  rowwise() %>% 
  mutate(sig_terms = list(anova(alpha_model) %>% 
                            rownames_to_column(var = "sig_term") %>% 
                            as_tibble() %>% 
                            rename("p_val" = `Pr(>F)`) %>%
                            mutate(fdr_p_val = p.adjust(p_val, method = 'fdr')) %>%
                            filter(fdr_p_val < 0.05) %>%
                            filter(sig_term != "time") %>%
                            pull(sig_term))) %>%
  filter(length(sig_terms) > 0) %>%
  select(-sig_terms) %>%
  mutate(alpha_type = ifelse(metric %in% c("chao1", "observed"), "richness", str_extract(metric, "[^_]+")))

## figuring out the fig

alpha_graphs <- mod_alpha_tab %>%
  rowwise() %>%
  mutate(plot_info = list(emmeans(alpha_model, ~treatment) %>%
                            broom::tidy(conf.int = TRUE) %>%
                            separate(treatment, into = c('time', 'exposure', 'susceptability')) %>%
                            mutate(graph_cat = ifelse(time == "T0", NA, 
                                                      paste(exposure, susceptability, sep = "_"))) %>%
                            {. ->> intermed } %>%
                            mutate(graph_cat = ifelse(time == "T0", paste("H", susceptability, sep = "_"), 
                                                      graph_cat)) %>%
                            dplyr::slice(rep(1:2, 1)) %>%
                            rbind(intermed) %>%
                            mutate(graph_cat = ifelse(is.na(graph_cat), paste("D", susceptability, sep = "_"), 
                                                      graph_cat)) %>%
                            mutate(c_time = parse_number(time)) %>%
                            mutate(facet_lab = "Experimental") %>%
                            mutate(c_time = ifelse(time == "T0", ifelse(susceptability == "S", c_time - 0.1, c_time + 0.1),
                                                   case_when(graph_cat == "D_S" ~ c_time - 0.35,
                                                             graph_cat == "H_S" ~ c_time - 0.15,
                                                             graph_cat == "D_R" ~ c_time + 0.15,
                                                             graph_cat == "H_R" ~ c_time + 0.35))) %>%
                            mutate(graph_cat = factor(graph_cat, levels = c("D_S", "D_R", "H_S", "H_R"), labels = c("D_S", "D_R", "H_S", "H_R"))))) %>%
  rowwise() %>%
  mutate(plot = list(
    ggplot(data = plot_info, aes(x = c_time, y = estimate, ymin = conf.low, ymax = conf.high)) +
      (geom_line(data = (plot_info %>% filter(graph_cat %in% c("D_S", "D_R"))), aes(colour1 = graph_cat, linetype = graph_cat)) %>%
         rename_geom_aes(new_aes = c("colour" = "colour1"))) + 
      (geom_line(data = (plot_info %>% filter(graph_cat %in% c("H_S", "H_R"))), aes(colour2 = graph_cat, linetype = graph_cat)) %>%
         rename_geom_aes(new_aes = c("colour" = "colour2"))) +
      (geom_errorbar(data = (plot_info %>% filter(graph_cat %in% c("D_S", "D_R"))), width = 0, aes(colour1 = graph_cat)) %>%
         rename_geom_aes(new_aes = c("colour" = "colour1"))) + 
      (geom_errorbar(data = (plot_info %>% filter(graph_cat %in% c("H_S", "H_R"))), width = 0, aes(colour2 = graph_cat)) %>%
         rename_geom_aes(new_aes = c("colour" = "colour2"))) +
      (geom_point(data = (plot_info %>% filter(graph_cat %in% c("D_S", "D_R"))), size = 3, aes(colour1 = graph_cat, pch = graph_cat)) %>%
         rename_geom_aes(new_aes = c("colour" = "colour1"))) + 
      (geom_point(data = (plot_info %>% filter(graph_cat %in% c("H_S", "H_R"))), size = 3, aes(colour2 = graph_cat, pch = graph_cat)) %>%
         rename_geom_aes(new_aes = c("colour" = "colour2"))) +
      
      
      (geom_point(data = (plot_info %>% filter(graph_cat == "dose" & susceptability == "na")), size = 3.7, aes(colour3 = exposure), shape = "diamond") %>%
         rename_geom_aes(new_aes = c("colour" = "colour3"))) +
      (geom_errorbar(data = (plot_info %>% filter(graph_cat == "dose" & susceptability == "na")), width = 0, aes(colour3 = exposure)) %>%
         rename_geom_aes(new_aes = c("colour" = "colour3"))) + 
      (geom_point(data = (plot_info %>% filter(graph_cat == "dose" & is.na(exposure))), size = 3, shape = 1, col = "black", alpha = 0)) +
      
      scale_color_manual(aesthetics = "colour1", values = c("#F75D5D", "#A70000"), guide = "legend", 
                         name = "Disease Exposed", labels = c("Susceptible", "Resistant")) +
      scale_color_manual(aesthetics = "colour3", values = c("#F10A0A", "#29A5B2"), guide = "legend", 
                         name = "Doses", labels = c("Diseased", "Healthy")) +
      scale_shape_manual(values = c(17, 16, 17, 16), guide = "none") +
      scale_color_manual(aesthetics = "colour2", values = c("#3DD8EA", "#048291"), guide = "legend", 
                         name = "Healthy Exposed", labels = c("Susceptible", "Resistant")) +
      guides(colour1 = guide_legend(
        override.aes=list(linetype = c(6, 1), shape = c(16, 17))),
        colour2 = guide_legend(
          override.aes=list(linetype = c(6, 1), shape = c(16, 17))),
        colour3 = guide_legend(
          override.aes=list(linetype = c(0, 0)))) +
      scale_x_continuous(breaks=c(0, 3, 7)) +
      
      scale_linetype_manual(values = c(1, 6, 1, 6), guide = "none") +
      theme_bw() +
      xlab("Time") +
      ylab(expression("Normalized log"[2]*" (cpm)")) +
      labs(title = metric)
  )) %>%
  group_by(alpha_type) %>%
  summarise(combo_plots = list(wrap_plots(plot) + plot_layout(guides = 'collect') & plot_annotation(title = alpha_type)))


alpha_graphs$combo_plots[[2]]


### other model type for alpha div

metrics_to_use <- c("diversity_gini_simpson", "dominance_absolute", 
                    "dominance_core_abundance", "evenness_camargo", "rarity_log_modulo_skewness")

t3t7_alpha_models <- alpha_table %>%
  filter(!tank %in% c("HOMO", "homogenate_fragment")) %>%
  mutate(treatment = str_c(time, exposure, susceptability, sep = '_')) %>%
  pivot_longer(cols = !c(colnames(metadata), "fragment_id", "treatment", "retain_sample"),
               names_to = 'metric',
               values_to = 'alpha_div_value') %>%
  select(-c(final_disease_state, clone_group, retain_sample)) %>%
  filter(time %in% c("T3", "T7")) %>%
  nest_by(metric) %>%
  summarise(t3t7_model = list(lmer(alpha_div_value ~ time*exposure*susceptability + 
                                     (1 | genotype) + (1 | tank), data = data)))

t0_alpha_models <- alpha_table %>%
  filter(!tank %in% c("HOMO", "homogenate_fragment")) %>%
  mutate(treatment = str_c(time, exposure, susceptability, sep = '_')) %>%
  pivot_longer(cols = !c(colnames(metadata), "fragment_id", "treatment", "retain_sample"),
               names_to = 'metric',
               values_to = 'alpha_div_value') %>%
  select(-c(final_disease_state, clone_group, retain_sample)) %>%
  filter(time %in% c("T0")) %>%
  nest_by(metric) %>%
  summarise(t0_model = list(lm(alpha_div_value ~ susceptability, data = data))) 

t3t7_sig_terms <- t3t7_alpha_models %>%
  rowwise() %>% 
  mutate(sig_terms = list(anova(t3t7_model) %>% 
                            rownames_to_column(var = "sig_term") %>% 
                            as_tibble() %>% 
                            rename("p_val" = `Pr(>F)`) %>%
                            mutate(fdr_p_val = p.adjust(p_val, method = 'fdr')) %>%
                            filter(fdr_p_val < 0.05) %>%
                            filter(sig_term != "time") %>%
                            pull(sig_term))) %>%
  filter(length(sig_terms) > 0) %>%
  mutate(int_term = ifelse(length(sig_terms) == 2, sig_terms[[2]], sig_terms)) %>%
  mutate(other_term = ifelse(length(sig_terms) == 2, sig_terms[[1]], NA))

t0_plot_prep <- t0_alpha_models %>%
  filter(metric %in% metrics_to_use) %>%
  filter(metric != "rarity_log_modulo_skewness") %>%
  rowwise() %>%
  mutate(t0_plot_info = list(emmeans(t0_model, ~susceptability) %>%
                               cld(Letters = letters) %>% # , adjust = 'fdr' ??
                               broom::tidy(conf.int = TRUE) %>%
                               mutate(.group = str_trim(.group)) %>%
                               mutate(c_time = 0) %>%
                               mutate(time = NA, exposure = NA, facet_lab = "Field") %>%
                               mutate(facet_lab = factor(facet_lab, levels = c("Field", "Experimental")))
  )) 

t3t7_plot_prep <- t3t7_sig_terms %>%
  filter(metric %in% metrics_to_use) %>%
  filter(metric != "rarity_log_modulo_skewness") %>%
  left_join(t0_plot_prep, by = join_by("metric")) %>%
  mutate(t3t7_plot_info = list(emmeans(t3t7_model, ~time*exposure) %>%
                                 cld(Letters = LETTERS) %>% # , adjust = 'fdr' ??
                                 broom::tidy(conf.int = TRUE) %>%
                                 mutate(.group = str_trim(.group)) %>%
                                 mutate(c_time = parse_number(time)) %>%
                                 mutate(susceptability = NA, facet_lab = "Experimental") %>%
                                 mutate(facet_lab = factor(facet_lab, levels = c("Field", "Experimental")))
  )) %>%
  mutate(both_plot_infos = list(rbind(t0_plot_info, t3t7_plot_info)))


slrp <- t3t7_plot_prep %>%
  rowwise() %>%
  mutate(plot = list(
    ggplot(data = both_plot_infos, aes(x = c_time, y = estimate, ymin = conf.low, ymax = conf.high)) +
      geom_point(data = both_plot_infos %>% filter(facet_lab == "Experimental"), 
                 aes(col = exposure, pch = exposure, size = exposure), position = position_dodge(0.5)) +
      geom_point(data = both_plot_infos %>% filter(facet_lab == "Field"), 
                 aes(col = susceptability, pch = susceptability, size = susceptability), position = position_dodge(0.5)) +
      geom_errorbar(data = both_plot_infos %>% filter(facet_lab == "Experimental"), 
                    aes(col = exposure), width = 0, position = position_dodge(0.5)) +
      geom_errorbar(data = both_plot_infos %>% filter(facet_lab == "Field"), 
                    aes(col = susceptability), width = 0, position = position_dodge(0.5)) +
      
      geom_line(data = both_plot_infos %>% filter(facet_lab == "Experimental"), 
                aes(col = exposure, linetype = exposure), position = position_dodge(0.5)) +
      geom_text(data = both_plot_infos %>% filter(facet_lab == "Experimental"), 
                aes(y = conf.high, label = .group, col = exposure),
                position = position_dodge(0.5), vjust = -1, show.legend = FALSE) +
      geom_text(data = both_plot_infos %>% filter(facet_lab == "Field"), 
                aes(y = conf.high, label = .group, col = susceptability),
                position = position_dodge(0.5), vjust = -1, show.legend = FALSE) +
      theme_bw() +
      scale_x_continuous(breaks = seq(0,7,1)) +
      scale_color_manual(values = c(
        "D" = "firebrick1",
        "H" = "seagreen2",
        "R" = "purple",
        "S" = "orange2"
      )) +
      scale_shape_manual(values = c(
        "D" = "diamond",
        "H" = "square",
        "R" = "triangle",
        "S" = "circle"
      )) +
      scale_size_manual(values = c(
        "D" = 3.5,
        "H" = 2.5,
        "R" = 2.5,
        "S" = 2.5
      )) +
      scale_linetype_manual(values = c(6,1), guide = "none") +
      facet_grid(cols = vars(facet_lab), scales = "free", space = "free") +
      ggtitle(paste(metric, " (exposure:time)")) +
      guides(color = guide_legend(
        override.aes=list(linetype = c(6, 1, 0, 0), shape = c("diamond", "square", "triangle", "circle"))))
  )) %>%
  group_by(int_term) %>%
  summarise(time_exp_plots = list(wrap_plots(plot) + plot_layout(guides = 'collect')))


slrp$time_exp_plots





#### Simplified Complex Upset ####
simple_comp_upset <- clean_bacstrats %>% 
  filter(!bacstrat %in% c("continuous_crasher", "late_crasher")) %>%
  filter(!asv_id %in% c("ASV_112", "ASV_288", "ASV_131", "ASV_15", "ASV_127")) %>% #pathogens not in D dose
  mutate(sig = TRUE) %>%
  pivot_wider(names_from = bacstrat, values_from = sig) %>%
  mutate(across(-asv_id, ~ifelse(is.na(.), FALSE, .))) %>% 
  left_join(taxonomy_tibble, by = c('asv_id' = 'asv_names')) # TAX TIBBLE USED HERE


### Upset - bacterial strategies

upset(
  simple_comp_upset,
  colnames(simple_comp_upset %>%
             select(starts_with("cont"), starts_with("late"), starts_with("early"))), 
  
  base_annotations=list(
    ' '=intersection_size(counts=T, text = aes(size = 6.5),
                          bar_number_threshold = 25,
                          mapping=aes(fill=Family, col = Family, label = parse_number(asv_id)), col = "gray10"
    ) +
      geom_text(size = 4, position = position_stack(vjust = 0.5), col = "gray10") +
      ylim(0, 19) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.text.y=element_blank(),
            axis.ticks.y=element_blank()) +
      scale_fill_manual(values = c(
        "Arenicellaceae" = "#BF4728",
        "Bdellovibrionaceae" = "#B6908B",
        "Cellvibrionaceae" = "#B5CAF0",
        "Colwelliaceae" = "#CE2220",
        "Cryomorphaceae" = "#EBDEA4",
        "Crocinitomicaceae" = "#306A26",
        "Flavobacteriaceae" = "#E67F33",
        "Fokiniaceae" = "#7EB875",
        "Francisellaceae" = "#D0B541",
        "Oligoflexaceae" = "#57A2AC",
        "P13-46" = "#7C4942",
        "Hyphomonadaceae" = "#4E78C4",
        "Puniceicoccaceae" = "#B997C7", 
        "Rhodobacteraceae" = "#824D99",
        "Saprospiraceae" = "#992572",
        "Sphingomonadaceae" = "#F3B4DE"
      ))
  ),
  matrix=(
    intersection_matrix(geom=geom_point(shape = "circle filled", size=3, stroke = 0.35, color = "gray20"))
    + scale_color_manual(
      values=c(
        "late_pathogen" = "#FF1E1E",
        "early_pathogen" = "#FF9797",
        #"continuous_pathogen" = "firebrick4",
        "late_opportunist" = "deepskyblue3",
        "early_opportunist" = "lightskyblue",
        #"continuous_crasher" = "#FF9A00",
        #"late_crasher" = "#BD5E03",
        "late_probiotic" = "#49EACC"
      )
    )
  ),
  queries=list(
    upset_query(set = "late_pathogen", fill = "#FF1E1E"),
    upset_query(set = "early_pathogen", fill = "#FF9797"),
    #upset_query(set = "continuous_pathogen", fill = "firebrick4"),
    upset_query(set = "late_opportunist", fill = "deepskyblue3"),
    upset_query(set = "early_opportunist", fill = "lightskyblue"),
    #upset_query(set = "continuous_crasher", fill = "#FF9A00"),
    #upset_query(set = "late_crasher", fill = "#BD5E03"),
    upset_query(set = "late_probiotic", fill = "#49EACC")
  ),
  name='ASVs', width_ratio=0.1, min_size = 0, wrap=TRUE, stripes = "white",
  set_sizes=(
    upset_set_size(filter_intersections=TRUE)
    + theme(axis.ticks.x=element_line())
  )) + #, sort_sets = FALSE
  ggtitle("Bacterial Strategies") 

### Upset - bacterial strategies and Presence Data

upset(
  simple_comp_upset %>% left_join(venn_all_times_and_doses, by = c("asv_id" = "OTU")),
  colnames(simple_comp_upset %>% left_join(venn_all_times_and_doses, by = c("asv_id" = "OTU")) %>%
             select(T0, D, H, starts_with("cont"), starts_with("late"), starts_with("early"))),
  base_annotations=list(
    ' '=intersection_size(counts=T, text = aes(size = 6.5),
                          bar_number_threshold = 25,
                          mapping=aes(fill=Family, col = Family, label = parse_number(asv_id)), col = "gray10"
    ) +
      geom_text(size = 4, position = position_stack(vjust = 0.5), col = "gray10") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.text.y=element_blank(),
            axis.ticks.y=element_blank()) +
      scale_fill_manual(values = c(
        "Arenicellaceae" = "#BF4728",
        "Bdellovibrionaceae" = "#B6908B",
        "Cellvibrionaceae" = "#B5CAF0",
        "Colwelliaceae" = "#CE2220",
        "Cryomorphaceae" = "#EBDEA4",
        "Crocinitomicaceae" = "#306A26",
        "Flavobacteriaceae" = "#E67F33",
        "Fokiniaceae" = "#7EB875",
        "Francisellaceae" = "#D0B541",
        "Oligoflexaceae" = "#57A2AC",
        "P13-46" = "#7C4942",
        "Hyphomonadaceae" = "#4E78C4",
        "Puniceicoccaceae" = "#B997C7", 
        "Rhodobacteraceae" = "#824D99",
        "Saprospiraceae" = "#992572",
        "Sphingomonadaceae" = "#F3B4DE"
      ))
  ),
  matrix=(
    intersection_matrix(geom=geom_point(shape = "circle filled", size=3, stroke = 0.35, color = "gray20"))
    + scale_color_manual(
      values=c(
        #"T7" = "#650197",
        #"T3" = "#B21BFF",
        "T0" = "#B21BFF",
        "D" = "#AE0404",
        "H" = "#0FAB02",
        "late_pathogen" = "#FF1E1E",
        "early_pathogen" = "#FF9797",
        #"continuous_pathogen" = "firebrick4",
        "late_opportunist" = "deepskyblue3",
        "early_opportunist" = "lightskyblue",
        #"continuous_crasher" = "#FF9A00",
        #"late_crasher" = "#BD5E03",
        "late_probiotic" = "#49EACC"
      )
    )
  ),
  queries=list(
    #upset_query(set = "T7", fill = "#650197"),
    #upset_query(set = "T3", fill = "#B21BFF"), "#D98EFF"
    upset_query(set = "T0", fill = "#B21BFF"),
    upset_query(set = "D", fill = "#AE0404"),
    upset_query(set = "H", fill = "#0FAB02"),
    upset_query(set = "late_pathogen", fill = "#FF1E1E"),
    upset_query(set = "early_pathogen", fill = "#FF9797"),
    #upset_query(set = "continuous_pathogen", fill = "firebrick4"),
    upset_query(set = "late_opportunist", fill = "deepskyblue3"),
    upset_query(set = "early_opportunist", fill = "lightskyblue"),
    #upset_query(set = "continuous_crasher", fill = "#FF9A00"),
    #upset_query(set = "late_crasher", fill = "#BD5E03"),
    upset_query(set = "late_probiotic", fill = "#49EACC")
  ),
  name='ASVs', width_ratio=0.1, min_size = 0, sort_sets = FALSE, wrap=TRUE,
  set_sizes=(
    upset_set_size(filter_intersections=TRUE)
    + theme(axis.ticks.x=element_line())
  ),
  stripes = c("#F8EAFF", "#FFE4E4", "#E7FFE5", rep(c("gray91", "gray97"), 4))) +
  #stripes = c(rep(c("gray78", "gray87"), 2), "gray78", rep(c("gray97", "gray91"), 4))) +
  ggtitle("Bacterial Strategies/Presence Data") 


### Pathogens Only

pathogen_upset <- upset(
  simple_comp_upset %>% left_join(venn_all_times_and_doses, by = c("asv_id" = "OTU")) %>% 
    filter(early_pathogen | late_pathogen),
  colnames(simple_comp_upset %>% left_join(venn_all_times_and_doses, by = c("asv_id" = "OTU")) %>%
             select(T7, T3, T0, D, H, late_pathogen, early_pathogen)),
  base_annotations=list(
    ' '=intersection_size(counts=T, text = aes(size = 6.5),
                          bar_number_threshold = 25,
                          mapping=aes(fill=Family, col = Family, label = parse_number(asv_id)), col = "gray10"
    ) +
      geom_text(size = 4, position = position_stack(vjust = 0.5), col = "gray10") +
      ylim(0,6) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.text.y=element_blank(),
            axis.ticks.y=element_blank()) +
      scale_fill_manual(values = c(
        "Arenicellaceae" = "#BF4728",
        "Bdellovibrionaceae" = "#B6908B",
        "Cellvibrionaceae" = "#B5CAF0",
        "Colwelliaceae" = "#CE2220",
        "Cryomorphaceae" = "#EBDEA4",
        "Crocinitomicaceae" = "#306A26",
        "Flavobacteriaceae" = "#E67F33",
        "Fokiniaceae" = "#7EB875",
        "Francisellaceae" = "#D0B541",
        "Oligoflexaceae" = "#57A2AC",
        "P13-46" = "#7C4942",
        "Hyphomonadaceae" = "#4E78C4",
        "Puniceicoccaceae" = "#B997C7", 
        "Rhodobacteraceae" = "#824D99",
        "Saprospiraceae" = "#992572",
        "Sphingomonadaceae" = "#F3B4DE"
      ))
  ),
  matrix=(
    intersection_matrix(geom=geom_point(shape = "circle filled", size=3, stroke = 0.35, color = "gray20"))
    + scale_color_manual(
      values=c(
        "T7" = "#650197",
        "T3" = "#B21BFF",
        "T0" = "#D98EFF",
        "D" = "#AE0404",
        "H" = "#0FAB02",
        "late_pathogen" = "#FF1E1E",
        "early_pathogen" = "#FF9797"
      )
    )
  ),
  queries=list(
    upset_query(set = "T7", fill = "#650197"),
    upset_query(set = "T3", fill = "#B21BFF"),
    upset_query(set = "T0", fill = "#D98EFF"),
    upset_query(set = "D", fill = "#AE0404"),
    upset_query(set = "H", fill = "#0FAB02"),
    upset_query(set = "late_pathogen", fill = "#FF1E1E"),
    upset_query(set = "early_pathogen", fill = "#FF9797")
  ),
  name='ASVs', width_ratio=0.1, min_size = 0, sort_sets = FALSE, wrap=TRUE,
  set_sizes=(
    upset_set_size(filter_intersections=TRUE)
    + theme(axis.ticks.x=element_line()) +
      scale_y_reverse(breaks = seq(15, 0, -5))
  ),
  stripes = c("#F8EAFF", "#F8EAFF", "#F8EAFF", "#FFE4E4", "#E7FFE5", "white", "white")) +
  #stripes = c(rep(c("gray78", "gray87"), 2), "gray78", rep(c("gray97", "gray91"), 4))) +
  ggtitle("Putative Pathogen Candidates Only") 

pathogen_upset

opportunist_upset <- upset(
  simple_comp_upset %>% left_join(venn_all_times_and_doses, by = c("asv_id" = "OTU")) %>% mutate(Family = factor(Family)) %>%
    filter(early_opportunist | late_opportunist),
  colnames(simple_comp_upset %>% left_join(venn_all_times_and_doses, by = c("asv_id" = "OTU")) %>%
             select(T7, T3, T0, D, H, late_opportunist, early_opportunist)),
  base_annotations=list(
    ' '=intersection_size(counts=T, text = aes(size = 6.5),
                          bar_number_threshold = 25,
                          mapping=aes(fill=Family, col = Family, label = parse_number(asv_id)), col = "gray10"
    ) +
      geom_text(size = 4, position = position_stack(vjust = 0.5), col = "gray10") +
      ylim(0, 9) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.text.y=element_blank(),
            axis.ticks.y=element_blank()) +
      scale_fill_manual(values = c(
        "Arenicellaceae" = "#BF4728",
        "Bdellovibrionaceae" = "#B6908B",
        "Cellvibrionaceae" = "#B5CAF0",
        "Colwelliaceae" = "#CE2220",
        "Cryomorphaceae" = "#EBDEA4",
        "Crocinitomicaceae" = "#306A26",
        "Flavobacteriaceae" = "#E67F33",
        "Fokiniaceae" = "#7EB875",
        "Francisellaceae" = "#D0B541",
        "Oligoflexaceae" = "#57A2AC",
        "P13-46" = "#7C4942",
        "Hyphomonadaceae" = "#4E78C4",
        "Puniceicoccaceae" = "#B997C7", 
        "Rhodobacteraceae" = "#824D99",
        "Saprospiraceae" = "#992572",
        "Sphingomonadaceae" = "#F3B4DE"
      ))
  ),
  
  
  
  
  matrix=(
    intersection_matrix(geom=geom_point(shape = "circle filled", size=3, stroke = 0.35, color = "gray20"))
    + scale_color_manual(
      values=c(
        "T7" = "#650197",
        "T3" = "#B21BFF",
        "T0" = "#D98EFF",
        "D" = "#AE0404",
        "H" = "#0FAB02",
        "late_opportunist" = "deepskyblue3",
        "early_opportunist" = "lightskyblue"
      )
    )
  ),
  queries=list(
    upset_query(set = "T7", fill = "#650197"),
    upset_query(set = "T3", fill = "#B21BFF"), 
    upset_query(set = "T0", fill = "#D98EFF"),
    upset_query(set = "D", fill = "#AE0404"),
    upset_query(set = "H", fill = "#0FAB02"),
    upset_query(set = "late_opportunist", fill = "deepskyblue3"),
    upset_query(set = "early_opportunist", fill = "lightskyblue")
  ),
  name='ASVs', width_ratio=0.1, min_size = 0, sort_sets = FALSE, wrap=TRUE,
  set_sizes=(
    upset_set_size(filter_intersections=TRUE)
    + theme(axis.ticks.x=element_line())
  ),
  stripes = c("#F8EAFF", "#F8EAFF", "#F8EAFF", "#FFE4E4", "#E7FFE5", "white", "white")) +
  #stripes = c("#F8EAFF", "#FFE4E4", "#E7FFE5", rep(c("gray91", "gray97"), 4))) +
  #stripes = c(rep(c("gray78", "gray87"), 2), "gray78", rep(c("gray97", "gray91"), 4))) +
  ggtitle("Putative Opportunist Candidates Only") 

opportunist_upset

wrap_plots(pathogen_upset, opportunist_upset) + 
  plot_layout(guides = 'collect')


# probiotic complex upset

probiotic_upset <- upset(
  simple_comp_upset %>% left_join(venn_all_times_and_doses, by = c("asv_id" = "OTU")) %>% mutate(Family = factor(Family)) %>%
    filter(late_probiotic),
  colnames(simple_comp_upset %>% left_join(venn_all_times_and_doses, by = c("asv_id" = "OTU")) %>%
             select(T7, T3, T0, D, H, late_probiotic)),
  base_annotations=list(
    ' '=intersection_size(counts=T, text = aes(size = 6.5),
                          bar_number_threshold = 25,
                          mapping=aes(fill=Family, col = Family, label = parse_number(asv_id)), col = "gray10"
    ) +
      geom_text(size = 4, position = position_stack(vjust = 0.5), col = "gray10") +
      ylim(0, 8) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.text.y=element_blank(),
            axis.ticks.y=element_blank()) +
      scale_fill_manual(values = c(
        "Arenicellaceae" = "#BF4728",
        "Bdellovibrionaceae" = "#B6908B",
        "Cellvibrionaceae" = "#B5CAF0",
        "Colwelliaceae" = "#CE2220",
        "Cryomorphaceae" = "#EBDEA4",
        "Crocinitomicaceae" = "#306A26",
        "Flavobacteriaceae" = "#E67F33",
        "Fokiniaceae" = "#7EB875",
        "Francisellaceae" = "#D0B541",
        "Oligoflexaceae" = "#57A2AC",
        "P13-46" = "#7C4942",
        "Hyphomonadaceae" = "#4E78C4",
        "Puniceicoccaceae" = "#B997C7", 
        "Rhodobacteraceae" = "#824D99",
        "Saprospiraceae" = "#992572",
        "Sphingomonadaceae" = "#F3B4DE"
      ))
  ),
  
  
  
  
  matrix=(
    intersection_matrix(geom=geom_point(shape = "circle filled", size=3, stroke = 0.35, color = "gray20"))
    + scale_color_manual(
      values=c(
        "T7" = "#650197",
        "T3" = "#B21BFF",
        "T0" = "#D98EFF",
        "D" = "#AE0404",
        "H" = "#0FAB02",
        "late_probiotic" = "#49EACC"
      )
    )
  ),
  queries=list(
    upset_query(set = "T7", fill = "#650197"),
    upset_query(set = "T3", fill = "#B21BFF"), 
    upset_query(set = "T0", fill = "#D98EFF"),
    upset_query(set = "D", fill = "#AE0404"),
    upset_query(set = "H", fill = "#0FAB02"),
    upset_query(set = "late_probiotic", fill = "#49EACC")
  ),
  name='ASVs', width_ratio=0.1, min_size = 0, sort_sets = FALSE, wrap=TRUE,
  set_sizes=(
    upset_set_size(filter_intersections=TRUE)
    + theme(axis.ticks.x=element_line())
  ),
  stripes = c("#F8EAFF", "#F8EAFF", "#F8EAFF", "#FFE4E4", "#E7FFE5", "white")) +
  #stripes = c("#F8EAFF", "#FFE4E4", "#E7FFE5", rep(c("gray91", "gray97"), 4))) +
  #stripes = c(rep(c("gray78", "gray87"), 2), "gray78", rep(c("gray97", "gray91"), 4))) +
  ggtitle("Putative Probiotic Candidates Only")

probiotic_upset


#& plot_annotation(title = signatures)

#### Family Groupings Table ####

simple_comp_upset %>%
  rowwise() %>%
  mutate(across(early_pathogen:late_pathogen, ~ifelse(.x == TRUE, 1, 0))) %>%
  group_by(Family) %>%
  mutate(across(early_pathogen:late_pathogen, ~sum(.x))) %>%
  select(early_pathogen:late_pathogen, Family) %>%
  distinct() %>%
  select(Family, early_pathogen, late_pathogen, early_opportunist, late_opportunist) %>%
  ungroup() %>%
  mutate(total = select(., early_pathogen:late_opportunist) %>% rowSums(na.rm = TRUE)) %>%
  arrange(desc(total)) %>%
  formattable(align = c("l", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c"), list(
    Family = formatter("span", style = ~ style(color = "gray",font.weight = "bold")),
    early_pathogen = formatter("span", style = x ~ style(color = ifelse(x < 1, "white", "black"))),
    late_pathogen = formatter("span", style = x ~ style(color = ifelse(x < 1, "white", "black"))),
    early_opportunist = formatter("span", style = x ~ style(color = ifelse(x < 1, "white", "black"))),
    late_opportunist = formatter("span", style = x ~ style(color = ifelse(x < 1, "white", "black"))),
    #late_probiotic = formatter("span", style = x ~ style(color = ifelse(x < 1, "white", "black"))),
    total = color_tile("gray", "gray")
  )) #%>%
export_formattable("../Figures/fds_exp_diffs.png")




#### Comp Upset w Sig Row ####
homog_unc_sig <- homogenate_models %>% 
  filter(dose_pairwise <= 0.05) %>% 
  select(asv_id) %>%
  mutate(unc_sig = TRUE)

homog_fdr_sig <- homogenate_models %>% 
  filter(dose_pairwise_adj <= 0.05) %>% 
  select(asv_id) %>%
  mutate(fdr_sig = TRUE)


comp_upset_w_sigs <- comp_upset_bac_strat %>%
  left_join(venn_all_times_and_doses, by = c("asv_id" = "OTU")) %>%
  left_join(homog_unc_sig, by = join_by(asv_id)) %>%
  left_join(homog_fdr_sig, by = join_by(asv_id)) %>%
  mutate(across(D:fdr_sig, ~ifelse(is.na(.), FALSE, .)))

upset(
  comp_upset_w_sigs,
  colnames(comp_upset_w_sigs %>%
             select(T0, fdr_sig, unc_sig, D, H, crasher_t7_strict, crasher_t7, crasher_t3,late_opportunist, early_opportunist, late_pathogen, continuous_pathogen, early_pathogen)), 
  
  base_annotations=list(
    'Intersection size'=intersection_size(counts=T, text = aes(size = 6.5),
                                          bar_number_threshold = 25,
                                          mapping=aes(fill=Family, col = Family, label = parse_number(asv_id)), col = "gray10"
    ) +
      geom_text(size = 4, position = position_stack(vjust = 0.5), col = "gray10")
  ),
  matrix=(
    intersection_matrix(geom=geom_point(shape = "circle filled", size=3, stroke = 0.35, color = "gray20"))
    + scale_color_manual(
      values=c(
        #"T7" = "#650197",
        #"T3" = "#B21BFF",
        "T0" = "#B21BFF",
        "unc_sig" = "#ff87cd",
        "fdr_sig" = "#ff4ab4",
        "D" = "#AE0404",
        "H" = "#0FAB02",
        "late_pathogen" = "#FF1E1E",
        "early_pathogen" = "#FF9797",
        "continuous_pathogen" = "maroon",
        "continuous_opportunist" = "royalblue4",
        "late_opportunist" = "deepskyblue3",
        "early_opportunist" = "lightskyblue",
        "crasher_t3" = "#FF9A00",
        "crasher_t7" = "#BD5E03",
        "crasher_t7_strict" = "gold2"
        #"probiotic_t7_strict" = "#49EACC"
      )
    )
  ),
  queries=list(
    #upset_query(set = "T7", fill = "#650197"),
    #upset_query(set = "T3", fill = "#B21BFF"), "#D98EFF"
    upset_query(set = "T0", fill = "#B21BFF"),
    upset_query(set = "unc_sig", fill = "#ff87cd"),
    upset_query(set = "fdr_sig", fill = "#ff4ab4"),
    upset_query(set = "D", fill = "#AE0404"),
    upset_query(set = "H", fill = "#0FAB02"),
    upset_query(set = "late_pathogen", fill = "#FF1E1E"),
    upset_query(set = "early_pathogen", fill = "#FF9797"),
    upset_query(set = "continuous_pathogen", fill = "maroon"),
    upset_query(set = "continuous_opportunist", fill = "royalblue4"),
    upset_query(set = "late_opportunist", fill = "deepskyblue3"),
    upset_query(set = "early_opportunist", fill = "lightskyblue"),
    upset_query(set = "crasher_t3", fill = "#FF9A00"),
    upset_query(set = "crasher_t7", fill = "#BD5E03"),
    upset_query(set = "crasher_t7_strict", fill = "gold2")
    #upset_query(set = "probiotic_t7_strict", fill = "#49EACC")
  ),
  name='ASVs', width_ratio=0.1, min_size = 0, sort_sets = FALSE,
  stripes = c("#F8EAFF", "#fff2fa", "#fff2fa","#FFE4E4", "#E7FFE5", rep(c("gray91", "gray97"), 4))) +
  #stripes = c(rep(c("gray78", "gray87"), 2), "gray78", rep(c("gray97", "gray91"), 4))) +
  ggtitle("Bacterial Strategies/Presence Data") 

