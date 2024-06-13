# Making the model based on FDS instead of Resistance

setwd("/Users/emilytrytten/Desktop/GitHub/16S_Florida_Tank_Analysis/Code")

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
library(relayer) #devtools::install_github("clauswilke/relayer")
library(ComplexUpset)
library(corrplot)
library(Hmisc)
library(broom.mixed)
library(taxize)
library(parallel)
library(tidyverse)
library(patchwork)

refit_models <- TRUE

#### Functions ####
cluster <- new_cluster(parallel::detectCores() - 1)
cluster_library(cluster, c('dplyr', 'lmerTest', 'emmeans', 'stringr', 'tidyr'))

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
normalized_asv_counts <- read_csv('../intermediate_files/fully_preprocessed_samples.csv.gz', show_col_types = FALSE) %>%
  mutate(time = factor(time, ordered = TRUE)) %>%
  mutate(final_disease_state = ifelse(time == "T0", "F", final_disease_state)) %>%
  filter(!c(exposure == "H" & final_disease_state == "D")) %>%
  mutate(treatment = str_c(time, exposure, final_disease_state, sep = '_'),
         time_exposure = str_c(time, exposure, sep = '_'),
         timeC = str_extract(time, '[0-9]+') %>% as.numeric,
         across(Domain:Species, str_replace_na)) %>%
  mutate(asv_number = str_extract(asv_id, '[0-9]+') %>% as.integer)

homogenate_data <- read_csv('../intermediate_files/homogenate_cpm.csv')

#### Prevalence of Cysteiniphilum ####

normalized_asv_counts %>% 
  filter(asv_id == "ASV_65") %>%
  mutate(present = ifelse(log2_cpm > 5.255240, TRUE, FALSE), .after = log2_cpm) %>%
  filter(time %in% c("T0", "T7")) %>%
  group_by(time, final_disease_state, present) %>%
  reframe(n = n())
#71% of D
#11% of H
#6% of F

#Monte Carlo Simulations

t7_mc_samples <- normalized_asv_counts %>% 
  filter(asv_id == "ASV_65") %>%
  filter(time == "T7") %>%
  mutate(present = ifelse(log2_cpm > 5.255240, 1, 0), .after = log2_cpm)

runs <- 100000

#T7 Diseased
mc_d_t7 <- sample(t7_mc_samples %>% filter(final_disease_state == "D") %>% pull(present), runs, replace = T)
sum(mc_d_t7)/runs #prev of 0.71348

sd(mc_d_t7)/sqrt(runs) #SE of 0.001429784

#T7 Healthy
mc_h_t7 <- sample(t7_mc_samples %>% filter(final_disease_state == "H") %>% pull(present), runs, replace = T)
sum(mc_h_t7)/runs #prev of 0.11259

sd(mc_h_t7)/sqrt(runs) #SE of 0.0009995724

#T0
t0_mc_samples <- normalized_asv_counts %>% 
  filter(asv_id == "ASV_65") %>%
  filter(time == "T0") %>%
  mutate(present = ifelse(log2_cpm > 5.255240, 1, 0), .after = log2_cpm)

mc_d_t0 <- sample(t0_mc_samples$present, runs, replace = T)
sum(mc_d_t0)/runs #prev of 0.06103

sd(mc_d_t0)/sqrt(runs) #SE of 0.0007570067

#### Contrasts ####

# posthoc_order <- c('T0.F.F', 'T3.D.D', 'T3.D.H', 'T3.H.H', 'T7.D.D', 'T7.D.H', 'T7.H.H')
simple_posthoc_aquar <- list('aquarium' = c(-1, 0, 0, 1/2, 0, 0, 1/2))

simple_posthoc_exp <- list('exposure' = c(0, 1/4, 1/4, -1/2, 1/4, 1/4, -1/2))

simple_posthoc_outc <- list('outcome' = c(0, 1/2, -1/4, -1/4, 1/2, -1/4, -1/4))


DD_early <- list('DD_early' = c(-1, 1, 0, 0, 0, 0, 0))

DD_late <- list('DD_late' = c(0, -1, 0, 0, 1, 0, 0))

DD_vs_DH_early <- list('DD_vs_DH_early' = c(0, 1, -1, 0, 0, 0, 0))

DD_vs_DH_late <- list('DD_vs_DH_late' = c(0, 0, 0, 0, 1, -1, 0))

two_sided_tests <- tibble(microbial_signature = c('simple_posthoc_aquar', 'simple_posthoc_exp', 'simple_posthoc_outc'),
                          contrasts = list(simple_posthoc_aquar, simple_posthoc_exp, simple_posthoc_outc),
                          direction = '=') #=

right_tests <- tibble(microbial_signature = c('DD_early', 'DD_late', 'DD_vs_DH_early', 'DD_vs_DH_late'),
                      contrasts = list(DD_early, DD_late, DD_vs_DH_early, DD_vs_DH_late),
                      direction = '>') #>

posthoc_categories <- bind_rows(two_sided_tests, right_tests) %>%
  unnest(contrasts) %>%
  mutate(contrast_name = names(contrasts)) %>%
  group_by(contrast_name, contrasts, direction) %>%
  summarise(signatures = list(c(microbial_signature)),
            .groups = 'drop') %>%
  nest(contrast = -direction)

#emmeans(asv_models$model[[1]], ~treatment) %>%
#  contrast(simple_planned_posthocs)


#### Make ASV Models ####

if(file.exists('../intermediate_files/mixed_model_results.rds.gz') & !refit_models){
  asv_models <- read_rds('../intermediate_files/mixed_model_results.rds.gz')
} else {
  cluster_copy(cluster, c('posthoc_categories', 'run_posthoc'))
  
  asv_models <- normalized_asv_counts %>% 
    mutate(tank_field = if_else(str_detect(treatment, 'F'), 'field', 'tank')) %>%
    # filter(asv_id == 'ASV_65') %>%
    nest_by(across(c('asv_id', Domain:Species))) %>% #, Family_confidence:Species_confidence
    partition(cluster) %>%
    mutate(fit_model(log2_cpm ~ treatment + (1 | genotype) + #(1 | tank),
                       (0 + dummy(tank_field, c("tank")) | tank),
                     data, 
                     use_weights = FALSE),
           random_anova = list(rand(model)),
           process_model(model, re_model, random_anova),
           posthoc = list(run_posthoc(model, posthoc_categories))) %>%
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

#number of things significant in main model
asv_models %>%
  select(asv_id, starts_with('fdr')) %>% 
  mutate(across(starts_with('fdr'), ~. < 0.05)) %>%
  
  pivot_longer(cols = -asv_id,
               names_to = c('term'),
               values_to = 'significance',
               names_prefix = 'fdr_') %>%
  filter(significance) %>%
  group_by(asv_id) %>%
  mutate(sigs = str_c(term, collapse = ", ")) %>%
  select(asv_id, sigs) %>%
  distinct() %>%
  ungroup() %>%
  group_by(sigs) %>%
  reframe(n = n())


#### Significant Treatment Effect Upset ####
significant_models <- asv_models %>%
  filter(fdr_treatment < 0.05) %>%
  # slice(26) %>%
  rowwise %>%
  mutate(process_postHoc(posthoc)) %>%
  ungroup() %>%
  p_adjust(exclude_cols = c('treatment', 'tank', 'genotype'))


simple_posthoc_sig_asvs <- significant_models %>%
  select(asv_id, starts_with('estimate'), starts_with('fdr')) %>% # fdr or qvalue
  select(asv_id, contains(c('aquar', 'outcome', 'exposure'))) %>%
  mutate(across(starts_with('fdr'), ~. < 0.05),
         across(starts_with('estimate'), ~if_else(. < 0, -1L, 1L))) %>%
  
  pivot_longer(cols = -asv_id,
               names_to = c('.value', 'signatures', 'direction'),
               names_pattern = '(.*)_(.*)_(.*)') %>%
  dplyr::rename(significance = fdr) %>%
  mutate(signatures = case_when(signatures == 'outcome' & estimate < 0 ~ 'HealthyOutcome',
                                signatures == 'outcome' & estimate > 0 ~ 'DiseaseOutcome',
                                
                                signatures == 'aquarium' & estimate < 0 ~ 'Field',
                                signatures == 'aquarium' & estimate > 0 ~ 'Aquaria',
                                
                                signatures == 'exposure' & estimate < 0 ~ 'HealthyExposed',
                                signatures == 'exposure' & estimate > 0 ~ 'DiseaseExposed'),
         .keep = 'unused') %>%
  filter(significance) %>%
  select(-c(significance, direction))

bacterial_signature_asv <- significant_models %>%
  select(asv_id, starts_with('fdr')) %>% 
  select(-contains(c('treatment', 'tank', 'genotype'))) %>%
  mutate(across(starts_with('fdr'), ~. < 0.05)) %>%
  
  pivot_longer(cols = -asv_id,
               names_to = c('term'),
               values_to = 'significance') %>%
  mutate(term = str_remove(term, 'fdr_')) %>%
  
  #Select ASVs which fit all characteristics of any given bacterial signature
  mutate(direction = str_extract(term, '[><=]'),
         term = str_remove(term, '_[><=]')) %>%
  
  left_join(unnest(posthoc_categories, contrast), 
            by = c('term' = 'contrast_name', 'direction')) %>%
  select(-contrasts) %>%
  unnest(signatures) %>%
  group_by(asv_id, signatures) %>%
  filter(all(significance)) %>%
  ungroup %>%
  select(-term) %>%
  distinct() %>%
  #only keep ASVs that match one of the simple posthocs: aquarium, exposure, or outcome
  group_by(asv_id) %>%
  reframe(asv_id, signatures = str_c(signatures, collapse = ', ')) %>%
  distinct() %>%
  filter(str_detect(signatures, "simple")) %>%
  separate_rows(signatures, sep = ", ", convert = TRUE) %>%
  #convert the simple posthocs to which direction they are
  filter(!str_detect(signatures, "simple")) %>%
  rbind(simple_posthoc_sig_asvs) %>%
  arrange(asv_id)




bacterial_signature_asv %>%
  group_by(asv_id) %>%
  summarise(terms = list(signatures),
            .groups = 'drop') %>%
  ggplot(aes(x = terms)) +
  geom_bar() +
  scale_x_upset() +
  theme_classic() +
  theme_combmatrix(combmatrix.label.make_space = TRUE)



classified_asv_taxonomy <- select(significant_models, asv_id, Domain:Species) %>%
  left_join(bacterial_signature_asv %>%
              mutate(significance = TRUE) %>%
              pivot_wider(names_from = signatures,
                          values_from = significance,
                          values_fill = FALSE),
            by = 'asv_id') %>%
  mutate(across(Field:Aquaria, ~ifelse(is.na(.), FALSE, .))) %>%
  mutate(DD_continuous = ifelse(DD_early & DD_late, TRUE, FALSE)) %>%
  mutate(DD_early = ifelse(DD_continuous, FALSE, DD_early), DD_late = ifelse(DD_continuous, FALSE, DD_late))

#write_csv(classified_asv_taxonomy, '../intermediate_files/classified_significant_asvs.csv.gz')

#how many in each category
classified_asv_taxonomy %>% 
  mutate(across(Field:DD_continuous, ~ifelse(. == TRUE, 1, 0))) %>%
  select(-colnames(taxonomy_tibble %>% select(-asv_names))) %>%
  select(-asv_id) %>%
  mutate(across(Field:DD_continuous, ~ifelse(is.na(.), 0, .))) %>%
  colSums(.)

#how many in each category without tank effect
classified_asv_taxonomy %>% 
  filter(!Aquaria & !Field) %>%
  mutate(across(Field:DD_continuous, ~ifelse(. == TRUE, 1, 0))) %>%
  select(-colnames(taxonomy_tibble %>% select(-asv_names))) %>%
  select(-asv_id) %>%
  mutate(across(Field:DD_continuous, ~ifelse(is.na(.), 0, .))) %>%
  colSums(.)

#### Bacterial Strategies Upset ####

sig_classified_asvs <- classified_asv_taxonomy %>%
  filter(!if_all(colnames(classified_asv_taxonomy %>% select(-c(asv_id, colnames(taxonomy_tibble %>%
                                                                                   select(-asv_names))))), ~. == FALSE)) %>%
  mutate(DiseaseExposed = FALSE, HealthyExposed = FALSE) %>%
  mutate(TankEffect = ifelse(Field | Aquaria, TRUE, FALSE)) %>%
  select(-c(Field, Aquaria))

#how many genera in which category
sig_classified_asvs %>% 
  filter(!TankEffect) %>% 
  pivot_longer(cols = -c(colnames(taxonomy_tibble %>% select(-asv_names)), asv_id), 
               names_to = "signature", values_to = "significance") %>%
  filter(significance) %>%
  select(-significance) %>%
  group_by(asv_id) %>%
  reframe(Family, Genus, all_sigs = str_c(signature, collapse = ", ")) %>%
  distinct() %>%
  ungroup() %>%
  group_by(Family, Genus, all_sigs) %>%
  reframe(n = n()) %>%
  arrange(all_sigs)

#how much tank effect
classified_asv_taxonomy %>%
  filter(Aquaria | Field) %>%
  select(colnames(taxonomy_tibble %>% select(-asv_names)), asv_id, Aquaria, Field) %>%
  pivot_longer(cols = c(Aquaria, Field), names_to = "TankEffect", values_to = "sig") %>%
  filter(sig) %>%
  select(Family, TankEffect) %>%
  distinct %>%
  group_by(TankEffect) %>%
  reframe(n = n())

#just overarching bac strats
upset(sig_classified_asvs %>% select(-c(contains("DD"))),
      colnames(sig_classified_asvs %>% select(-c(asv_id, colnames(taxonomy_tibble %>% select(-asv_names)), contains("DD"))) %>%
                 relocate(TankEffect, contains("Exposed"), contains("Outcome"))), 
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
            #"Field" = "#70d134",
            #"HealthyExposed" = "#78acff",
            "HealthyOutcome" = "#0d50ba",
            "TankEffect"= "#c389e0",
            #"DiseaseExposed"= "#ff7878",
            "DiseaseOutcome"= "#ba0d0d"
          )
        )
      ),
      queries=list(
        #upset_query(set = "Field", fill = "#70d134"),
        #upset_query(set = "HealthyExposed", fill = "#78acff"),
        upset_query(set = "HealthyOutcome", fill = "#0d50ba"),
        upset_query(set = "TankEffect", fill = "#c389e0"),
        #upset_query(set = "DiseaseExposed", fill = "#ff7878"),
        upset_query(set = "DiseaseOutcome", fill = "#ba0d0d")
      ),
      name='ASVs', width_ratio=0.1, sort_sets = FALSE, sort_intersections = FALSE,
      stripes = c("#fcf5ff", "#fcd5d4", "#d5d6eb", "#fcd5d4", "#d5d6eb", rep ("gray95", 3))) +
  ggtitle("Overarching Strategies")

#overarching and specific bac strats
upset(sig_classified_asvs,
      colnames(sig_classified_asvs %>% select(-c( asv_id, colnames(taxonomy_tibble %>% select(-asv_names)))) %>%
                 relocate(TankEffect, contains("Exposed"), contains("Outcome"))), 
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
            #"HealthyExposed" = "#78acff",
            "HealthyOutcome" = "#0d50ba",
            "TankEffect"= "#c389e0",
            #"DiseaseExposed"= "#ff7878",
            "DiseaseOutcome"= "#ba0d0d",
            "DD_early"= "gray80",
            "DD_late"= "gray30",
            "DD_continuous"= "gray50"
            #"time_linear" = "gray30",
            #"time_quadratic" = "gray20",
            #"time_late" = "gray60",
            #"time_early" = "gray70"
          )
        )
      ),
      queries=list(
        #upset_query(set = "HealthyExposed", fill = "#78acff"),
        upset_query(set = "HealthyOutcome", fill = "#0d50ba"),
        upset_query(set = "TankEffect", fill = "#c389e0"),
        #upset_query(set = "DiseaseExposed", fill = "#ff7878"),
        upset_query(set = "DiseaseOutcome", fill = "#ba0d0d"),
        upset_query(set = "DD_early", fill = "gray80"),
        upset_query(set = "DD_late", fill = "gray30"),
        upset_query(set = "DD_continuous", fill = "gray50")
        #upset_query(set = "time_linear", fill = "gray30"),
        #upset_query(set = "time_quadratic", fill = "gray20"),
        #upset_query(set = "time_late", fill = "gray60"),
        #upset_query(set = "time_early", fill = "gray70")
      ),
      name='ASVs', width_ratio=0.1, sort_sets = FALSE, sort_intersections=FALSE,
      stripes = c("#fcf5ff","#fcd5d4", "#d5d6eb", "#fcd5d4", "#d5d6eb", rep(c("gray95"), 3))) +
  ggtitle("All Bacterial Strategies")

# Non Tank Related bact strats  
upset(sig_classified_asvs %>% filter(!TankEffect) %>% select(-TankEffect),
      colnames(sig_classified_asvs %>% filter(!TankEffect) %>% select(-c(asv_id, colnames(taxonomy_tibble %>% select(-asv_names)), TankEffect)) %>%
                 relocate(contains("Exposed"), contains("Outcome"))), 
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
            "HealthyOutcome" = "#0d50ba",
            #"TankEffect"= "#c389e0",
            #"DiseaseExposed"= "#ff7878",
            "DiseaseOutcome"= "#ba0d0d",
            "DD_early"= "gray80",
            "DD_late"= "gray30"
            #"DD_continuous"= "gray50"
          )
        )
      ),
      queries=list(
        upset_query(set = "HealthyOutcome", fill = "#0d50ba"),
        #upset_query(set = "TankEffect", fill = "#c389e0"),
        #upset_query(set = "DiseaseExposed", fill = "#ff7878"),
        upset_query(set = "DiseaseOutcome", fill = "#ba0d0d"),
        upset_query(set = "DD_early", fill = "gray80"),
        upset_query(set = "DD_late", fill = "gray30")
        #upset_query(set = "DD_continuous", fill = "gray50")
      ),
      name='ASVs', width_ratio=0.1, sort_sets = FALSE, sort_intersections = FALSE,
      stripes = c("#fcd5d4", "#d5d6eb", "#fcd5d4", "#d5d6eb", rep(c("gray91", "gray97"), 3))) +
  ggtitle("Non Tank-Related Bacterial Strategies")


#pathogens only upset
upset(sig_classified_asvs %>% select(c(asv_id, contains("pathogen"), colnames(taxonomy_tibble %>% select(-asv_names)))) %>%
        filter(early_pathogen | late_pathogen) %>% 
        left_join(venn_all_times_and_doses %>% select(OTU, T0_H, D), by = c("asv_id" = "OTU")) %>%
        mutate(across(D:T0_H, ~ifelse(is.na(.), FALSE, .))),
      colnames(sig_classified_asvs %>% select(asv_id, contains("pathogen")) %>%
                 left_join(venn_all_times_and_doses %>% select(OTU, T0_H, D), by = c("asv_id" = "OTU")) %>% select(-asv_id)), 
      
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
            "D" = "#cf5002",
            "T0_H" = "#318f48",
            "Field" = "#70d134",
            "late_pathogen" = "#FF1E1E",
            "early_pathogen" = "#FF9797"
          )
        )
      ),
      queries=list(
        upset_query(set = "D", fill = "#cf5002"),
        upset_query(set = "T0_H", fill = "#318f48"),
        upset_query(set = "late_pathogen", fill = "#FF1E1E"),
        upset_query(set = "early_pathogen", fill = "#FF9797")
      ),
      name='ASVs', width_ratio=0.1, min_size = 1, sort_sets = FALSE, min_degree = 1,
      stripes = c(rep(c("gray91", "gray97"), 2))) +
  #stripes = c(rep(c("gray78", "gray87"), 2), "gray78", rep(c("gray97", "gray91"), 4))) +
  ggtitle("Pathogen Candidates") 


#### More Upsets ####

non_tank_sig_classified_asvs %>% select(Family) %>% distinct()

cu_family_colors <- c("Fastidiosibacteraceae" = '#e6194B', 
                      "Fokiniaceae" = '#3cb44b', 
                      "Francisellaceae" = '#42d4f4', 
                      "Oceanospirillaceae" = '#ffe119', 
                      "Colwelliaceae" = '#f58231', 
                      "Hyphomonadaceae" = '#fffac8',
                      "Thiotrichaceae" = '#469990', 
                      "Puniceicoccaceae" = '#dcbeff', 
                      "Erythrobacteraceae" = '#9A6324', 
                      "Roseobacteraceae" = '#fabed4',
                      "Cellvibrionaceae" = "purple",
                      "Vibrionaceae" = "darkred",
                      "Arenicellaceae" = "lightgreen",
                      "Planctomycetaceae" = "tan",
                      "NA" = "gray60",
                      "Unclassified" = "gray60")

non_tank_sig_classified_asvs <- sig_classified_asvs %>%
  filter(!TankEffect) %>%
  left_join(sig_homog_dose_data, by = join_by(asv_id)) %>%
  select(-TankEffect)

sig_effects_cu <- upset(sig_classified_asvs %>% select(-c(contains("DD"))),
                        colnames(sig_classified_asvs %>% select(-c(asv_id, colnames(taxonomy_tibble %>% select(-asv_names)), contains("DD"))) %>%
                                   relocate(TankEffect, contains("Exposed"), contains("Outcome"))), 
                        base_annotations=list(
                          'Intersection size'=intersection_size(counts=T, text = aes(size = 6.5),
                                                                bar_number_threshold = 25,
                                                                fill= c("#ba0d0d", "#0d50ba", "#c389e0", "#c389e0", "#c389e0") #"slategray"
                          )
                        ),
                        matrix=(
                          intersection_matrix(geom=geom_point(shape = "circle filled", size=3, stroke = 0.35, color = "gray20"))
                          + scale_color_manual(
                            values=c(
                              #"Field" = "#70d134",
                              #"HealthyExposed" = "#78acff",
                              "HealthyOutcome" = "#0d50ba",
                              "TankEffect"= "#c389e0",
                              #"DiseaseExposed"= "#ff7878",
                              "DiseaseOutcome"= "#ba0d0d"
                            )
                          )
                        ),
                        queries=list(
                          #upset_query(set = "Field", fill = "#70d134"),
                          #upset_query(set = "HealthyExposed", fill = "#78acff"),
                          upset_query(set = "HealthyOutcome", fill = "#0d50ba"),
                          upset_query(set = "TankEffect", fill = "#c389e0"),
                          #upset_query(set = "DiseaseExposed", fill = "#ff7878"),
                          upset_query(set = "DiseaseOutcome", fill = "#ba0d0d")
                        ),
                        name='ASVs', width_ratio=0.1, sort_sets = FALSE, sort_intersections = FALSE,
                        stripes = c("#fcf5ff", "#FFF1EF", "#F0F5FF", "#FFF1EF", "#F0F5FF", rep ("gray95", 3))) +
  ggtitle("Significant Effects")







# EARLY PATHOGENS
early_path_cu <- upset(non_tank_sig_classified_asvs %>% filter(DiseaseOutcome) %>% filter(DD_early) %>% mutate(DD_vs_DH_early = FALSE) %>% mutate(Family = ifelse(is.na(Family), "Unclassified", Family)),
                       colnames(non_tank_sig_classified_asvs %>% filter(DiseaseOutcome) %>% select(-c(asv_id, colnames(taxonomy_tibble %>% select(-asv_names)))) %>% mutate(DD_vs_DH_early = FALSE) %>%
                                  relocate(contains("Outcome"), contains("DD")) %>% select(-c(HealthyExposed, DiseaseExposed, DD_continuous, HealthyOutcome, sig_dose_H, DD_vs_DH_late, DD_late))), 
                       base_annotations=list(
                         'Intersection size'=intersection_size(counts=T, text = aes(size = 6.5),
                                                               bar_number_threshold = 25,
                                                               mapping=aes(fill=Family, col = Family, label = str_c("ASV", parse_number(asv_id), sep = " ")), col = "gray10"
                         ) +
                           geom_text(size = 4, position = position_stack(vjust = 0.5), col = "gray10") +
                           scale_fill_manual(values = cu_family_colors)
                       ),
                       matrix=(
                         intersection_matrix(geom=geom_point(shape = "circle filled", size=3, stroke = 0.35, color = "gray20"))
                         + scale_color_manual(
                           values=c(
                             #"HealthyOutcome" = "#0d50ba",
                             #"TankEffect"= "#c389e0",
                             #"DiseaseExposed"= "#ff7878",
                             "sig_dose_D" = "#a3088c",
                             "DiseaseOutcome"= "#ba0d0d",
                             "DD_vs_DH_early"= "darkorange3",
                             "DD_early"= "gray80"
                             #"DD_late"= "gray30"
                             #"DD_continuous"= "gray50"
                           )
                         )
                       ),
                       queries=list(
                         #upset_query(set = "HealthyOutcome", fill = "#0d50ba"),
                         #upset_query(set = "TankEffect", fill = "#c389e0"),
                         #upset_query(set = "DiseaseExposed", fill = "#ff7878"),
                         upset_query(set = "sig_dose_D", fill = "#a3088c"),
                         upset_query(set = "DiseaseOutcome", fill = "#ba0d0d"),
                         upset_query(set = "DD_vs_DH_early", fill = "darkorange3"),
                         upset_query(set = "DD_early", fill = "gray80")
                         #upset_query(set = "DD_late", fill = "gray30")
                         #upset_query(set = "DD_continuous", fill = "gray50")
                       ),
                       name='ASVs', width_ratio=0.1, sort_sets = FALSE, sort_intersections = FALSE,
                       stripes = c(rep(c("gray91", "gray97"), 3))) +
  ggtitle("Early Pathogen Candidates")

# LATE PATHOGENS
late_path_cu <- upset(non_tank_sig_classified_asvs %>% filter(DiseaseOutcome) %>% filter(DD_late) %>% mutate(Family = ifelse(is.na(Family), "Unclassified", Family)),
                      colnames(non_tank_sig_classified_asvs %>% filter(DiseaseOutcome) %>% select(-c(asv_id, colnames(taxonomy_tibble %>% select(-asv_names)))) %>%
                                 relocate(contains("Outcome")) %>% select(-c(HealthyExposed, DiseaseExposed, DD_continuous, HealthyOutcome, sig_dose_H, DD_early))), 
                      base_annotations=list(
                        'Intersection size'=intersection_size(counts=T, text = aes(size = 6.5),
                                                              bar_number_threshold = 25,
                                                              mapping=aes(fill=Family, col = Family, label = str_c("ASV", parse_number(asv_id), sep = " ")), col = "gray10"
                        ) +
                          geom_text(size = 4, position = position_stack(vjust = 0.5), col = "gray10") +
                          scale_fill_manual(values = cu_family_colors)
                      ),
                      matrix=(
                        intersection_matrix(geom=geom_point(shape = "circle filled", size=3, stroke = 0.35, color = "gray20"))
                        + scale_color_manual(
                          values=c(
                            #"HealthyOutcome" = "#0d50ba",
                            #"TankEffect"= "#c389e0",
                            #"DiseaseExposed"= "#ff7878",
                            "sig_dose_D" = "#a3088c",
                            "DiseaseOutcome"= "#ba0d0d",
                            "DD_vs_DH_late"= "darkorange3",
                            #"DD_early"= "gray80",
                            "DD_late"= "gray30"
                            #"DD_continuous"= "gray50"
                          )
                        )
                      ),
                      queries=list(
                        #upset_query(set = "HealthyOutcome", fill = "#0d50ba"),
                        #upset_query(set = "TankEffect", fill = "#c389e0"),
                        #upset_query(set = "DiseaseExposed", fill = "#ff7878"),
                        upset_query(set = "sig_dose_D", fill = "#a3088c"),
                        upset_query(set = "DiseaseOutcome", fill = "#ba0d0d"),
                        upset_query(set = "DD_vs_DH_late", fill = "darkorange3"),
                        #upset_query(set = "DD_early", fill = "gray80"),
                        upset_query(set = "DD_late", fill = "gray30")
                        #upset_query(set = "DD_continuous", fill = "gray50")
                      ),
                      name='ASVs', width_ratio=0.1, sort_sets = FALSE, sort_intersections = FALSE,
                      stripes = c(rep(c("gray91", "gray97"), 3))) +
  ggtitle("Late Pathogen Candidates")

#sig_effects_cu /
(early_path_cu | late_path_cu)

#### Fancy Plots ####

## Prep Data for Homogenate Dose Panel

homogenate_models <- homogenate_data %>%
  nest_by(asv_id) %>%
  mutate(model = list(lm(log2_cpm ~ exposure, data = data))) %>%
  mutate(homogenate_pred = list(model %>%
                                  broom::tidy(conf.int = TRUE) %>%
                                  mutate(term = case_when(term == "(Intercept)" ~ "D",
                                                          term == "exposureH" ~ "H")) %>%
                                  mutate(estimate = ifelse(term == "H", estimate[term == "H"] + estimate[term == "D"], estimate),
                                         std.error = NA,
                                         statistic = ifelse(term == "H", statistic[term == "H"] + statistic[term == "D"], statistic),
                                         conf.low = ifelse(term == "H", conf.low[term == "H"] + conf.low[term == "D"], conf.low),
                                         conf.high = ifelse(term == "H", conf.high[term == "H"] + conf.high[term == "D"], conf.high),
                                         p.value = ifelse(term == "D", p.value[term == "H"], p.value),
                                         df = NA) %>%
                                  rename("exposure" = "term") %>%
                                  mutate(time = "Dose", 
                                         graph_cat = "dose", 
                                         c_time = ifelse(exposure == "D", -1.2, -1.5), 
                                         final_disease_state = "na",
                                         facet_lab = "Doses") %>%
                                  {. ->> set_one } %>% #save data to this var
                                  mutate(std.error = NA, df = NA, conf.low = NA, conf.high = NA, statistic = NA,
                                         c_time = c(-1.8, -0.9), final_disease_state = NA, exposure = NA) %>% #dummy set to change size of Dose Facet
                                  rbind(set_one) %>%
                                  {. ->> set_two } %>% #save data to this var
                                  arrange(desc(estimate)) %>%
                                  dplyr::slice(1) %>%
                                  mutate(exposure = "p_val", c_time = -1.35, estimate = estimate + 3) %>%
                                  rbind(set_two) #recombine them
  )) %>%
  mutate(homog_p_val = homogenate_pred %>% filter(exposure == "H") %>% pull(p.value)) %>%
  ungroup %>%
  mutate(adj_homog_p_val = p.adjust(homog_p_val, method = 'fdr')) %>%
  unnest(homogenate_pred) %>%
  mutate(p.value = adj_homog_p_val) %>%
  nest(homogenate_pred = -c(asv_id, data, model, homog_p_val, adj_homog_p_val))

sig_homog_asv_list <- homogenate_models %>%
  filter(adj_homog_p_val < 0.05) %>%
  pull(asv_id)

#set up zero line

add_zero_lines <- homogenate_models %>% ungroup() %>% select(homogenate_pred) %>% unnest(homogenate_pred) %>% 
  dplyr::slice(1:4) %>% mutate(across(everything(), ~NA)) %>% 
  mutate(facet_lab = c("Doses", "Experimental", "Doses", "Experimental"), c_time = c(-5, -5, 10, 10), 
         estimate = c(5.27, 5.294201, 5.27, 5.294201)) %>%
  mutate(sig_p = NA, time = "zero_lines")


## Make Fancy Plots

the_plots <- bacterial_signature_asv %>%
  group_by(asv_id) %>%
  reframe(signatures = str_c(signatures, collapse = ', ')) %>%
  mutate(grouped_signatures = case_when(str_detect(signatures, "Aquaria") ~ "Tank Effect",
                                        str_detect(signatures, "Field") ~ "Tank Effect",
                                        str_detect(signatures, "DiseaseOutcome") & str_detect(signatures, "DD_early") & str_detect(signatures, "DD_vs_DH_early") & asv_id %in% sig_homog_asv_list ~ "Putative Early Pathogens",
                                        str_detect(signatures, "DiseaseOutcome") & str_detect(signatures, "DD_late") & str_detect(signatures, "DD_vs_DH_late") & asv_id %in% sig_homog_asv_list ~ "Putative Late Pathogens",
                                        str_detect(signatures, "DiseaseOutcome") & str_detect(signatures, "DD_late") & str_detect(signatures, "DD_vs_DH_late") & !asv_id %in% sig_homog_asv_list~ "Specialized Opportunists",
                                        str_detect(signatures, "DiseaseOutcome") ~ "Unlikely Pathogens",
                                        str_detect(signatures, "HealthyOutcome") ~ "Healthy Associated",
                                        TRUE ~ signatures)) %>%
  group_by(asv_id, grouped_signatures) %>%
  distinct() %>%
  inner_join(significant_models,
             by = 'asv_id') %>%
  mutate(taxonomy_for_graph = case_when(is.na(Order) ~ Class,
                                        is.na(Family) ~ Order,
                                        is.na(Genus) ~ Family,
                                        is.na(Species) ~ str_c(Family, Genus, sep = " "),
                                        TRUE ~ str_c(Family, Species, sep = " "))) %>%
  mutate(formatted_ASV = str_c("(ASV ", parse_number(asv_id), ")", sep = "")) %>%
  mutate(prep = ifelse(is.na(Family),
                       list(substitute(taxonomy_for_graph~italic(species_abbrev)~formatted_ASV, list(taxonomy_for_graph = taxonomy_for_graph, species_abbrev = "sp.", formatted_ASV = formatted_ASV))),
                       ifelse(is.na(Species),
                              list(substitute(italic(taxonomy_for_graph)~italic(species_abbrev)~formatted_ASV, list(taxonomy_for_graph = taxonomy_for_graph, species_abbrev = "sp.", formatted_ASV = formatted_ASV))),
                              list(substitute(italic(taxonomy_for_graph)~formatted_ASV, list(taxonomy_for_graph = taxonomy_for_graph, formatted_ASV = formatted_ASV))))
                       
                       )
           ) %>%
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
      scale_alpha_manual(values = c("sig" = 1, "nonsig" = 0), guide = "none") +
      scale_linetype_manual(values = c("D_H" = 1, "D_D" = 1, "H_H" = 6), guide = "none") +
      theme_bw() +
      #theme(plot.title = element_text(face = "italic")) +
      xlab("Time") +
      ylab(expression("Normalized log"[2]*" (cpm)")) +
      labs(title = prep) +
      #labs(title = case_when(is.na(Family) ~ bquote('hmm'~italic(.(taxonomy_for_graph))), 
      #                       TRUE ~ bquote(italic(.(taxonomy_for_graph))))) +
      # labs(title = ifelse(is.na(Family),
      #                     str_c(taxonomy_for_graph, expression(italic("sp."))), #not all italics
      #                     substitute(italic(taxa_description), list(taxa_description = taxonomy_for_graph))),
      #      subtitle = str_c("ASV ", parse_number(asv_id), sep = "")) + # all italics
      #labs(title = substitute(italic(taxa_description), list(taxa_description = str_c(Family, ifelse(str_detect(Species, " "), Species, str_c(Genus, Species, sep = " ")), sep = " ")))) +
      #labs(title = substitute(italic(taxa_description), list(taxa_description = str_c(Family, ifelse(str_detect(Species, " "), Species, str_c(Genus, Species, sep = " ")), sep = " ")))) +
      #labs(title = str_c(ifelse(is.na(Family), "NA", Family), " (", ifelse(is.na(Family_confidence), "NA", round(Family_confidence, digits = 0)), "%) ", ifelse(is.na(Genus), "NA", Genus)," (", ifelse(is.na(Genus_confidence), "NA", round(Genus_confidence, digits = 0)), "%) ", sep = ""),
      #     subtitle = str_c(ifelse(is.na(Species), "NA", Species), " (", ifelse(is.na(Species_confidence), "NA", round(Species_confidence, digits = 0)), "%); ", asv_id, "\n", signatures, sep = "")) +
      facet_grid(cols = vars(facet_lab), space = "free") + #scales = "free_x"
      scale_y_continuous(limits = c(1, 15), breaks = c(2.5, 5, 7.5, 10, 12.5, 15)) +
      coord_panel_ranges(panel_ranges = list(
        list(x=c(-1.8, -0.9)), # Dose Panel
        list(x=c(-0.75, 8)) # Experimental Panel
      ))
  )) %>%
  group_by(grouped_signatures) %>%
  summarise(combo_plots = list(wrap_plots(plot) + plot_layout(guides = 'auto'))) # & plot_annotation(title = grouped_signatures)))


#view the plots
the_plots$combo_plots[[3]]


#### Complex Upset with Doses ####

sig_homog_dose_data <- homogenate_models %>%
  unnest(homogenate_pred) %>%
  filter(exposure %in% c("D", "H")) %>%
  group_by(asv_id) %>%
  select(asv_id, adj_homog_p_val, exposure, estimate) %>%
  mutate(which_dose_higher = ifelse(estimate[exposure == "D"] > estimate[exposure == "H"], "sig_dose_D", "sig_dose_H")) %>%
  select(-c(exposure, estimate)) %>%
  distinct() %>%
  mutate(signif = ifelse(adj_homog_p_val < 0.05, TRUE, FALSE)) %>%
  pivot_wider(names_from = which_dose_higher, values_from = signif) %>%
  mutate(across(contains("sig_dose"), ~ifelse(is.na(.), FALSE, .))) %>%
  select(-adj_homog_p_val)





sig_classified_asvs %>% left_join(sig_homog_dose_data, by = join_by(asv_id)) %>% filter(sig_dose_H)
  
  
  # Non Tank Related bact strats  
  upset(non_tank_sig_classified_asvs %>% select(-sig_dose_H) %>% mutate(sig_dose_H = FALSE),
        colnames(non_tank_sig_classified_asvs %>% select(-c(asv_id, colnames(taxonomy_tibble %>% select(-asv_names)))) %>%
                   relocate(contains("Exposed"), contains("Outcome"))), 
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
              "sig_dose_H" = "#08a38e",
              "sig_dose_D" = "#a3088c",
              "HealthyOutcome" = "#0d50ba",
              #"TankEffect"= "#c389e0",
              #"DiseaseExposed"= "#ff7878",
              "DiseaseOutcome"= "#ba0d0d",
              "DD_early"= "gray60",
              "DD_late"= "gray20"
              #"DD_continuous"= "gray50"
            )
          )
        ),
        queries=list(
          upset_query(set = "sig_dose_H", fill = "#08a38e"),
          upset_query(set = "sig_dose_D", fill = "#a3088c"),
          upset_query(set = "HealthyOutcome", fill = "#0d50ba"),
          #upset_query(set = "TankEffect", fill = "#c389e0"),
          #upset_query(set = "DiseaseExposed", fill = "#ff7878"),
          upset_query(set = "DiseaseOutcome", fill = "#ba0d0d"),
          upset_query(set = "DD_early", fill = "gray60"),
          upset_query(set = "DD_late", fill = "gray20")
          #upset_query(set = "DD_continuous", fill = "gray50")
        ),
        name='ASVs', width_ratio=0.1, sort_sets = FALSE, sort_intersections = FALSE,
        stripes = c("#fcd5d4", "#d5d6eb", "#fcd5d4", "#d5d6eb", rep(c("gray91", "gray97"), 3))) +
  ggtitle("Non Tank-Related Bacterial Strategies")
  
  
#### NMDS ####

nmds_matrix <- normalized_asv_counts %>% 
  select(asv_id, sample_id, log2_cpm) %>%
  pivot_wider(names_from = sample_id, values_from = log2_cpm) %>%
  column_to_rownames('asv_id') %>%
  t() %>%
  as.matrix()

asv_nmds <- metaMDS(nmds_matrix, distance = 'bray', k = 3, trymax = 100, autotransform = FALSE, verbose = TRUE)

#shepard plot
#doesn't follow y = x line well, so MDS might be misleading -> use nmds
plot(asv_nmds$diss, asv_nmds$dist)
abline(a = 0, b = 1, col = "red")

stressplot(asv_nmds)

plot(asv_nmds)

#ugly but functional I guess
ordiplot(asv_nmds,type="n")
orditorp(asv_nmds,display="species",col="red",air=0.01)
orditorp(asv_nmds,display="sites",cex=1.25,air=0.01)

nmds_scores = as.data.frame(scores(asv_nmds)$sites)

nmds_metadata <- normalized_asv_counts %>% 
  select(asv_id, time, exposure, final_disease_state, susceptability, genotype, 
         sample_id, log2_cpm, tank, treatment) %>%
  pivot_wider(names_from = asv_id, values_from = log2_cpm) %>%
  as.data.frame()

nmds_scores$time = nmds_metadata$time
nmds_scores$exposure = nmds_metadata$exposure
nmds_scores$final_disease_state = nmds_metadata$final_disease_state
nmds_scores$susceptability = nmds_metadata$susceptability
nmds_scores$genotype = nmds_metadata$genotype
nmds_scores$treatment = str_c(nmds_scores$time, nmds_scores$exposure, nmds_scores$final_disease_state, sep = "_")

#nice nmds
ggplot(nmds_scores) +
  geom_point(aes(x = NMDS1, y = NMDS2, col = treatment, pch = time), size = 2.5) +
  stat_ellipse(aes(x = NMDS1, y = NMDS2, col = treatment)) +
  theme_bw() +
  scale_shape_manual(values = c("T0" = 16, "T3" = 15, "T7" = 17), guide = "none") +
  scale_color_manual(name = "Treatment", values = c("T0_F_F" = "#c389e0",
                     "T3_D_D" = "#E79B9B",
                     "T3_D_H" = "#97D9E1",
                     "T3_H_H" = "#95AC85",
                     "T7_D_D" = "#A70000",
                     "T7_D_H" = "#22A7B6",
                     "T7_H_H" = "#406F23")) +
  guides(color = guide_legend(
    override.aes=list(shape = c("T0_F_F" = 16,
                                "T3_D_D" = 15,
                                "T3_D_H" = 15,
                                "T3_H_H" = 15,
                                "T7_D_D" = 17,
                                "T7_D_H" = 17,
                                "T7_H_H" = 17),
                      size = 3)))


  scale_shape_manual(values = c("T0_F_F" = 1))
  
adonis2(nmds_matrix ~time*final_disease_state*exposure + genotype + tank, method = "bray", perm = 1000, data = nmds_metadata)  


#library(pairwiseAdonis)
#pairwise.adonis2(nmds_matrix ~treatment, data = nmds_metadata) %>% broom::tidy()

#### PCoA ####

dist_mat <- vegdist(nmds_matrix)
multivar_disp <- betadisper(dist_mat, nmds_metadata$treatment)

asv_pcoa <- cmdscale (dist_mat, eig = TRUE)



pcoa_data <- asv_pcoa$points %>%
  as_tibble(rownames = 'sample_id', .name_repair = "unique") %>%
  dplyr::rename("PCoA1" = "...1", "PCoA2" = "...2") %>%
  left_join(normalized_asv_counts %>% select(sample_id, treatment, time, exposure, final_disease_state, tank, genotype) %>%
              mutate(shape_determ = paste(exposure, final_disease_state, sep = "_")) %>%
              distinct(), by = join_by("sample_id")) %>% 
  mutate(PCoA1 = -PCoA1) %>%
  group_by(treatment) %>%
  mutate(PCoA1_centroid = mean(PCoA1), PCoA2_centroid = mean(PCoA2))

# better version of betadisper plot
ggplot(pcoa_data, aes(x = PCoA1, y = PCoA2, xend = PCoA1_centroid, yend = PCoA2_centroid)) +
  stat_ellipse(level = 0.95, aes(col = treatment)) +
  geom_point(aes(col = treatment, shape = treatment), size = 2.5) + 
  geom_segment(aes(col = treatment)) +
  scale_color_manual(values = c(
    "T0_F_F" = "#c389e0",
    "T3_D_D" = "#E79B9B",
    "T3_D_H" = "#97D9E1",
    "T3_H_H" = "#95AC85",
    "T7_D_D" = "#A70000",
    "T7_D_H" = "#22A7B6",
    "T7_H_H" = "#406F23"
  )) +
  scale_shape_manual(values = c(
    "T0_F_F" = 16,
    "T3_D_D" = 15,
    "T3_D_H" = 15,
    "T3_H_H" = 15,
    "T7_D_D" = 17,
    "T7_D_H" = 17,
    "T7_H_H" = 17
  ), guide = "none") + 
  geom_text(data = pcoa_data %>% select(treatment, PCoA1_centroid, PCoA2_centroid) %>% distinct(), 
            aes(label = treatment, x = PCoA1_centroid, y = PCoA2_centroid, fontface = "bold"), col = "black", inherit.aes = FALSE) +
  theme_bw() +
  labs(col = "Treatment") +
  guides(color = guide_legend(
    override.aes=list(shape = c("T0_F_F" = 16,
                                "T3_D_D" = 15,
                                "T3_D_H" = 15,
                                "T3_H_H" = 15,
                                "T7_D_D" = 17,
                                "T7_D_H" = 17,
                                "T7_H_H" = 17),
                      size = 3)))

###### OLD v

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


#### Time 0 Microbiome ####

normalized_asv_counts %>%
  filter(time == "T0")

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


#### Microshades Plot ####

mdf_prep_test1 <- altered_microbiome_data %>%
  tax_glom("Genus") %>%
  psmelt() 

#check the most abundant Orders for microshades
# mdf_prep_test1 %>%
# select(-c(Phylum:Family)) %>%
# left_join(nonoverlapping_taxonomy, by = join_by(Genus)) %>%
# group_by(Order) %>% reframe(tot_sum = sum(Abundance)) %>% arrange(desc(tot_sum))

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
  reframe(Domain, time, total, rel_abun = sum(Abundance)/total) %>%
  distinct() %>%
  left_join(nonoverlapping_taxonomy, by = join_by(Genus)) %>%
  dplyr::rename(Sample = category, Abundance = rel_abun) %>%
  mutate(Sample = factor(Sample, levels = c("Healthy", "Diseased", "T0", "T3_H", "T3_D", "T7_H_H", "T7_D_H", "T7_D_D"))) %>%
  as.data.frame()


## ORDER GENUS

color_objs_ordergenus <- create_color_dfs(mdf_processed_data, group_level = "Order", 
                                          selected_groups = c("Rickettsiales", "Alteromonadales", "Flavobacteriales", 
                                                              "Spirochaetales",  "Rhodobacterales"), cvd = TRUE)
mdf_ordergenus <- color_objs_ordergenus$mdf
cdf_ordergenus <- color_objs_ordergenus$cdf

legend_ordergenus <-custom_legend(mdf_ordergenus, cdf_ordergenus, group_level = "Order")

plot_ordergenus_prelim <- plot_microshades(mdf_ordergenus, cdf = cdf_ordergenus) + 
  scale_y_continuous(labels = scales::percent, expand = expansion(0)) +
  facet_grid(cols = vars(time), scales = "free", space = "free") +
  theme_bw() +
  theme(legend.position = "none", plot.margin = margin(6,20,6,6)) +
  #labs(title = "Order Genus") +
  scale_x_discrete(name = element_blank(), labels = c("T0" = "Field", "T3_H" = "Healthy", "T3_D" =
                                                        "Disease", "T7_H_H" = "Healthy (Healthy)",
                                                      "T7_D_H" = "Disease (Healthy)", "T7_D_D" = "Disease (Diseased)"))

plot_grid(plot_ordergenus_prelim, legend_ordergenus,  rel_widths = c(1, .25))

#### Alpha Diversity ####

alpha_table <- microbiome::alpha(altered_microbiome_data, index = "all") %>%
  as_tibble(rownames = 'sample_id') %>%
  inner_join(metadata, by = 'sample_id') %>%
  mutate(fragment_id = str_c(str_replace_na(exposure, 'NA'), tank, genotype, sep = '_'))

mod_alpha_tab <- alpha_table %>%
  filter(!tank %in% c("HOMO", "homogenate_fragment")) %>%
  mutate(final_disease_state = ifelse(exposure == "F", "F", final_disease_state)) %>%
  mutate(treatment = str_c(time, exposure, final_disease_state, sep = '_')) %>%
  pivot_longer(cols = !c(colnames(metadata), "fragment_id", "treatment"),
               names_to = 'metric',
               values_to = 'alpha_div_value') %>%
  select(-c(susceptability, resistance, clone_group)) %>%
  mutate(tank_field = if_else(str_detect(treatment, 'F'), 'field', 'tank'), .after = final_disease_state) %>%
  nest_by(metric) %>%
  summarise(alpha_model = list(lmer(alpha_div_value ~ treatment + 
                                      (1 | genotype) + (0 + dummy(tank_field, c("tank")) | tank),
                                    data = data))) %>% 
  rowwise() %>% 
  mutate(p_value = anova(alpha_model) %>% 
           rownames_to_column(var = "sig_term") %>% 
           as_tibble() %>% 
           dplyr::rename("p_val" = `Pr(>F)`) %>%
           pull(p_val)) %>%
  ungroup() %>%
  mutate(fdr_p_val = p.adjust(p_value, method = 'fdr')) %>%
  filter(fdr_p_val < 0.05) %>%
  mutate(alpha_type = ifelse(metric %in% c("chao1", "observed"), "richness", str_extract(metric, "[^_]+")))

## figuring out the fig

alpha_graphs_manuscript <- mod_alpha_tab %>%
  mutate(metric = ifelse(str_detect(metric, "chao1"), "richness_chao1", metric)) %>%
  mutate(for_manuscript = ifelse(metric %in% c("diversity_shannon", "dominance_core_abundance", "evenness_camargo", "richness_chao1"), TRUE, FALSE)) %>%
  rowwise() %>%
  mutate(plot_info = list(emmeans(alpha_model, ~treatment) %>%
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
  )) %>%
  group_by(for_manuscript) %>%
  summarise(combo_plots = list(wrap_plots(plot) + plot_layout(guides = 'collect')))

  
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


#### Comparing BH and q-value FDR Methods ####

## Q-value vs. FDR Notes

#The q-value of a test measures the proportion of false positives incurred (called the false discovery rate) 
#when that particular test is called significant.

#fdr controls the false discovery rate, the expected proportion of false discoveries 
#amongst the rejected hypotheses.

#guarantees that the true FDR rate will be less than the specified rate on average if you do an exactly similar
#experiment over and over again. So the BH approach is slightly more conservative than qvalue. 
#The BH properties hold regardless of the number of p-values, while qvalue is asymptotic, so the BH approach 
#is more robust than qvalue when the number of hypotheses being tested isn't very large.


#Using q-values allows us to decide how many false positives we are willing to accept among all the features
#that we call significant. This is particularly useful when we wish to make a large number of discoveries for
#further confirmation later on (i.e. pilot study or exploratory analyses, for example if we did a gene 
#expression microarray to pick differentially expressed genes for confirmation with real-time PCR). This is 
#also useful in genomewide studies where we expect a sizeable portion of features to be truly alternative, and
#we do not want to restrict our discovery capacity



fdr_asvs_venn <- classified_asvs %>%
  filter(significance) %>%
  rename("signatures" = "term") %>%
  group_by(asv_id) %>%
  reframe(signatures = str_c(signatures, collapse = ', '), significance) %>%
  mutate(grouped_signatures = case_when(str_detect(signatures, "Field") ~ "TankAverse",
                                        str_detect(signatures, "Aquaria") ~ "TankAssociated",
                                        TRUE ~ signatures)) %>%
  select(-signatures) %>%
  pivot_wider(names_from = grouped_signatures, values_from = significance, names_prefix = "fdr_")


qvalue_asvs_venn <- classified_asvs %>%
  filter(significance) %>%
  rename("signatures" = "term") %>%
  group_by(asv_id) %>%
  reframe(signatures = str_c(signatures, collapse = ', '), significance) %>%
  mutate(grouped_signatures = case_when(str_detect(signatures, "Field") ~ "TankAverse",
                                        str_detect(signatures, "Aquaria") ~ "TankAssociated",
                                        TRUE ~ signatures)) %>%
  select(-signatures) %>%
  distinct() %>%
  pivot_wider(names_from = grouped_signatures, values_from = significance, names_prefix = "qvalue_")


correction_comparison_cu <- fdr_asvs_venn %>%
  full_join(qvalue_asvs_venn, by = join_by(asv_id)) %>%
  mutate(across(everything(), ~ifelse(is.na(.), FALSE, .))) %>%
  ungroup()


upset(correction_comparison_cu,
      colnames(correction_comparison_cu %>% select(-asv_id)), 
      matrix=(
        intersection_matrix(geom=geom_point(shape = "circle filled", size=3, stroke = 0.35, color = "gray20"))
        + scale_color_manual(
          values=c(
            "fdr_HealthyOutcome" = "#318f48",
            "fdr_DiseaseOutcome" = "#ba0d0d",
            "fdr_TankAssociated" = "#c389e0",
            "qvalue_HealthyOutcome" = "#318f48",
            "qvalue_TankAverse" = "#70d134",
            "qvalue_TankAssociated"= "#c389e0",
            "qvalue_DiseaseOutcome"= "#ba0d0d",
            "qvalue_DiseaseExposed, DiseaseOutcome"= "#ff7878"
          )
        )
      ),
      queries=list(
        upset_query(set = "fdr_HealthyOutcome", fill = "#318f48"),
        upset_query(set = "fdr_DiseaseOutcome", fill = "#ba0d0d"),
        upset_query(set = "fdr_TankAssociated", fill = "#c389e0"),
        upset_query(set = "qvalue_HealthyOutcome", fill = "#318f48"),
        upset_query(set = "qvalue_TankAverse", fill = "#70d134"),
        upset_query(set = "qvalue_TankAssociated", fill = "#c389e0"),
        upset_query(set = "qvalue_DiseaseOutcome", fill = "#ba0d0d"),
        upset_query(set = "qvalue_DiseaseExposed, DiseaseOutcome", fill = "#ff7878")
      ),
      name='ASVs', width_ratio=0.1, min_size = 1, sort_sets = FALSE, min_degree = 1,
      stripes = c(rep("#ffe3cc", 3), rep("#ccd2ff", 6))) +
  ggtitle("Comparing FDR and q-value") 


#### Genotypes chart ####

#by exposure + FDS
full_data %>%
  select(genotype, time, exposure, final_disease_state, tank) %>%
  distinct() %>%
  group_by(time, exposure, final_disease_state) %>%
  reframe(counts = n()) %>%
  filter(time != "T0") %>% #18 genotypes T0
  pivot_wider(names_from = c(exposure, final_disease_state), values_from = counts)

#by susceptibility
full_data %>% #normalized_asv_counts %>%
  select(genotype, time, susceptability) %>%
  distinct() %>%
  group_by(time, susceptability) %>%
  reframe(counts = n()) %>%
  pivot_wider(names_from = susceptability, values_from = counts)

full_data %>%
  filter(exposure == "H" & final_disease_state == "D") %>%
  select(time, tank, genotype, clone_group, susceptability, resistance) %>%
  distinct() %>%
  arrange(time, tank, resistance)


genotype_comp_upset <- full_data %>%
  select(time, exposure, final_disease_state, genotype) %>%
  mutate(final_disease_state = ifelse(exposure == "F", "F", final_disease_state)) %>%
  mutate(combo_descriptor = str_c(exposure, final_disease_state, sep = "_")) %>%
  select(genotype, combo_descriptor) %>%
  distinct() %>%
  mutate(sig = TRUE) %>%
  pivot_wider(names_from = combo_descriptor, values_from = sig) %>%
  mutate(across(!genotype, ~ifelse(is.na(.), FALSE, .)))


upset(genotype_comp_upset,
      colnames(genotype_comp_upset %>% select(-genotype)), 
      matrix=(
        intersection_matrix(geom=geom_point(shape = "circle filled", size=3, stroke = 0.35, color = "gray20"))
        + scale_color_manual(
          values=c(
            "H_D" = "purple",
            "H_H" = "royalblue2",
            "D_H" = "#3DD8EA",
            "D_D" = "#F75D5D",
            "F_F" = "#70d134"
            
          )
        )
      ),
      queries=list(
        upset_query(set = "H_D", fill = "purple"),
        upset_query(set = "H_H", fill = "royalblue2"),
        upset_query(set = "D_H", fill = "#3DD8EA"),
        upset_query(set = "D_D", fill = "#F75D5D"),
        upset_query(set = "F_F", fill = "#70d134")
      ),
      name='ASVs', width_ratio=0.1, min_size = 1, min_degree = 1,
      stripes = c("gray90", "gray95")) +
  ggtitle("Genotypes")

# by timepoint - doesn't include H_D

genotype_timepoint_cu <- normalized_asv_counts %>%
  select(time, genotype) %>%
  distinct() %>%
  mutate(present = TRUE) %>%
  pivot_wider(names_from = time, values_from = present) %>%
  mutate(across(contains("T"), ~ifelse(is.na(.), FALSE, .)))

upset(genotype_timepoint_cu,
      colnames(genotype_timepoint_cu %>% select(-genotype) %>% select("T7", "T3", "T0")), 
      base_annotations=list(
        'Intersection size'=intersection_size(counts=T, text = aes(size = 6.5),
                                              bar_number_threshold = 25,
                                              mapping=aes(label = genotype), fill= "slategray1", col = "gray10"
        ) +
          geom_text(size = 4, position = position_stack(vjust = 0.5), col = "gray10")
      ),
      matrix=(
        intersection_matrix(geom=geom_point(shape = "circle filled", size=3, stroke = 0.35, color = "gray20"))
        + scale_color_manual(
          values=c(
            "T0" = "pink",
            "T3" = "#70d134",
            "T7" = "orange"
            
          )
        )
      ),
      queries=list(
        upset_query(set = "T0", fill = "pink"),
        upset_query(set = "T3", fill = "#70d134"),
        upset_query(set = "T7", fill = "orange")
      ),
      name='ASVs', width_ratio=0.1, min_size = 1, min_degree = 0, sort_sets = FALSE,
      stripes = c("gray90", "gray95", "gray98")) +
  ggtitle("Genotypes by Timepoint")

upset(genotype_timepoint_cu,
      colnames(genotype_timepoint_cu %>% select(-genotype) %>% select("T7", "T3", "T0")), 
      matrix=(
        intersection_matrix(geom=geom_point(shape = "circle filled", size=3, stroke = 0.35, color = "gray20"))
        + scale_color_manual(
          values=c(
            "T0" = "pink",
            "T3" = "#70d134",
            "T7" = "orange"
            
          )
        )
      ),
      queries=list(
        upset_query(set = "T0", fill = "pink"),
        upset_query(set = "T3", fill = "#70d134"),
        upset_query(set = "T7", fill = "orange")
      ),
      name='ASVs', width_ratio=0.1, min_size = 1, min_degree = 0, sort_sets = FALSE,
      stripes = c("gray90", "gray95", "gray98")) +
  ggtitle("Genotypes by Timepoint")

## by time and exposure
genotype_timepoint_exposure_cu <- normalized_asv_counts %>%
  select(genotype, time, exposure) %>%
  distinct() %>%
  mutate(present = TRUE, time_exp = str_c(time, exposure, sep = "_")) %>%
  select(genotype, time_exp, present) %>%
  pivot_wider(names_from = time_exp, values_from = present) %>%
  mutate(across(contains("T"), ~ifelse(is.na(.), FALSE, .))) %>%
  rename("T0" = "T0_F")


upset(genotype_timepoint_exposure_cu,
      colnames(genotype_timepoint_exposure_cu %>% select(-genotype) %>% select("T7_D", "T7_H", "T3_D", "T3_H", "T0")), 
      base_annotations=list(
        'Intersection size'=intersection_size(counts=T, text = aes(size = 6.5),
                                              bar_number_threshold = 25,
                                              mapping=aes(label = genotype), fill= "slategray1", col = "gray10"
        ) +
          geom_text(size = 4, position = position_stack(vjust = 0.5), col = "gray10")
      ),
      matrix=(
        intersection_matrix(geom=geom_point(shape = "circle filled", size=3, stroke = 0.35, color = "gray20"))
        + scale_color_manual(
          values=c(
            "T0" = "pink",
            "T3_H" = "#70d134",
            "T3_D" = "#317009",
            "T7_H" = "#ffa526",
            "T7_D" = "#a36308"
            
          )
        )
      ),
      queries=list(
        upset_query(set = "T0", fill = "pink"),
        upset_query(set = "T3_H", fill = "#70d134"),
        upset_query(set = "T3_D", fill = "#317009"),
        upset_query(set = "T7_D", fill = "#a36308"),
        upset_query(set = "T7_H", fill = "#ffa526")
      ),
      name='ASVs', width_ratio=0.1, min_size = 1, min_degree = 0, sort_sets = FALSE,
      stripes = c("gray90", "gray90", "gray95", "gray95", "gray98")) +
  ggtitle("Genotypes by Timepoint and Exposure")
