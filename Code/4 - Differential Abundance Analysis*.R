setwd("/Users/Emily/Desktop/GitHub/16S_Florida_Tank_Analysis/Code")

#### Libraries ####
library(tidyverse)
library(parallel)
library(multidplyr)
library(lmerTest)
library(phyloseq)
library(forcats)
library(ComplexUpset)
library(ggupset)
library(ggnested)
library(relayer) #devtools::install_github("clauswilke/relayer")
library(patchwork)
library(ggpubr)
library(vegan)
library(emmeans)
library(rempsyc)
library(glue)
library(ggtext)

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

#### Read In Data ####

#read in and format data
normalized_asv_counts <- read_csv('../intermediate_files/fully_preprocessed_samples.csv.gz', show_col_types = FALSE) %>%
  mutate(time = factor(time, ordered = TRUE)) %>%
  mutate(final_disease_state = ifelse(time == "T0", "F", final_disease_state)) %>%
  mutate(treatment = str_c(time, exposure, final_disease_state, sep = '_'),
         time_exposure = str_c(time, exposure, sep = '_'),
         timeC = str_extract(time, '[0-9]+') %>% as.numeric,
         across(Domain:Species, str_replace_na)) %>%
  mutate(asv_number = str_extract(asv_id, '[0-9]+') %>% as.integer)

#homogenate dose data
homogenate_data <- read_csv('../intermediate_files/homogenate_cpm.csv')

#taxonomy information
taxonomy_tibble <- tax_table(read_rds("../intermediate_files/updated_microbiome_data.rds")) %>% 
  as.data.frame() %>% 
  as_tibble()

#detection limit of 5.202306
detection_limit <- min(normalized_asv_counts$log2_cpm)

#detection limit in homogenate doses of 5.274105
homog_detection_limit <- min(homogenate_data$log2_cpm)

#### Miscellaneous Manuscript Metrics ####

#Supplementary Data #3 - Read Counts of each ASV per Sample (tank environment)
# normalized_asv_counts %>% 
#   select(sample_id, asv_id, read_count) %>%
#   pivot_wider(names_from = asv_id, values_from = read_count) %>%
#   as.data.frame() %>% column_to_rownames(var = "sample_id") %>%
#   as.matrix() %>%
#   write.csv("../Supplementary Data/ECT2025_supplementary_data_3.csv")

#Supplementary Data #4 - Log2 CPM of each ASV per Sample
# normalized_asv_counts %>%
#   select(sample_id, asv_id, log2_cpm) %>%
#   rbind(homogenate_data %>% select(sample_id, asv_id, log2_cpm)) %>%
#   pivot_wider(names_from = asv_id, values_from = log2_cpm) %>%
#   as.data.frame() %>% column_to_rownames(var = "sample_id") %>%
#   as.matrix() %>%
#   write.csv("../Supplementary Data/ECT2025_supplementary_data_4.csv")

#range of disease resistance scores for T0 corals
normalized_asv_counts %>%
  filter(time == "T0") %>%
  select(genotype, resistance) %>%
  distinct() %>%
  reframe(num_genotypes = length(genotype), minimum_DR = min(resistance), maximum_DR = max(resistance),
          average_DR = mean(resistance), SE_DR = sd(resistance)/sqrt(length((resistance))))

#total number of each taxa present in the homogenate doses
homogenate_data %>% 
  left_join(normalized_asv_counts %>% select(asv_id, colnames(taxonomy_tibble)) %>% distinct(), by = join_by("asv_id")) %>%
  filter(log2_cpm > homog_detection_limit + 0.000000001) %>%
  select(Genus) %>% #change taxa level here
  filter(!is.na(.)) %>%
  distinct() %>%
  nrow()

# how many healthy and diseased fragments at T7
normalized_asv_counts %>%
  filter(time == "T7") %>%
  select(fragment_id, final_disease_state) %>%
  distinct() %>%
  group_by(final_disease_state) %>%
  reframe(num_of_fragments = n())

#how many fragments at T3
normalized_asv_counts %>%
  filter(time == "T3") %>%
  select(fragment_id) %>%
  distinct() %>%
  reframe(num_of_fragments = n())

#Monte Carlo Simulations for Prevalence of Cysteiniphilum litorale

runs <- 100000

#how many times the 30 most frequently observed log2 CPM values were seen
#shows a clear detection limit
normalized_asv_counts %>% group_by(log2_cpm) %>% summarize(n = n()) %>%
  arrange(desc(n)) %>% slice(1:30) %>%
  ggplot() + geom_col(aes(x = fct_reorder(factor(log2_cpm), n), y = n), fill = "slategray3") +
  coord_flip() + theme_bw() + xlab("Log2 CPM") + ylab("Number of times observed")

#T7
t7_mc_samples <- normalized_asv_counts %>% 
  filter(asv_id == "ASV_65") %>%
  filter(time == "T7") %>%
  mutate(present = ifelse(log2_cpm > detection_limit + 0.000000001, 1, 0), .after = log2_cpm)

#T7 Diseased
mc_d_t7 <- sample(t7_mc_samples %>% filter(final_disease_state == "D") %>% pull(present), runs, replace = T)
sum(mc_d_t7)/runs #prev of 0.71063

sd(mc_d_t7)/sqrt(runs) #SE of 0.001434005

#T7 Healthy
mc_h_t7 <- sample(t7_mc_samples %>% filter(final_disease_state == "H") %>% pull(present), runs, replace = T)
sum(mc_h_t7)/runs #prev of 0.11413

sd(mc_h_t7)/sqrt(runs) #SE of 0.001005512


#T7 Healthy, Healthy Exposed
mc_h_h_t7 <- sample(t7_mc_samples %>% filter(final_disease_state == "H" & exposure == "H") %>% pull(present), runs, replace = T)
sum(mc_h_h_t7)/runs #prev of 0.11539

sd(mc_h_h_t7)/sqrt(runs) #SE of 0.001010328

#T0
t0_mc_samples <- normalized_asv_counts %>% 
  filter(asv_id == "ASV_65") %>%
  filter(time == "T0") %>%
  mutate(present = ifelse(log2_cpm > detection_limit + 0.000000001, 1, 0), .after = log2_cpm)

mc_d_t0 <- sample(t0_mc_samples$present, runs, replace = T)
sum(mc_d_t0)/runs #prev of 0.06387

sd(mc_d_t0)/sqrt(runs) #SE of 0.0007732478

#### Post Hoc Contrasts ####

# posthoc_order <- c('T0.F.F', 'T3.D.D', 'T3.D.H', 'T3.H.H', 'T7.D.D', 'T7.D.H', 'T7.H.H')
# order given by emmeans(model, ~treatment)

#compares the average of T3/T7 healthy-exposed healthy corals to T0 field-collected corals
simple_posthoc_aquar <- list('aquarium' = c(-1, 0, 0, 1/2, 0, 0, 1/2))

#compares the average of T3/T7 corals exposed to diseased homogenate to the average of T3/T7 corals exposed to healthy homogenate
simple_posthoc_exp <- list('exposure' = c(0, 1/4, 1/4, -1/2, 1/4, 1/4, -1/2))

#compares the average of T3/T7 corals that contract disease to the average of T3/T7 corals that remain asymptomatic (putatively healthy)
simple_posthoc_outc <- list('outcome' = c(0, 1/2, -1/4, -1/4, 1/2, -1/4, -1/4))

#compares fragments that eventually contract disease at T3 to field-collected corals at T0
DD_early <- list('DD_early' = c(-1, 1, 0, 0, 0, 0, 0))

#compares diseased fragments at T7 to fragments that eventually contract disease at T3
DD_late <- list('DD_late' = c(0, -1, 0, 0, 1, 0, 0))

#compares fragments that eventually contract disease to disease-exposed healthy corals at T3
DD_vs_DH_early <- list('DD_vs_DH_early' = c(0, 1, -1, 0, 0, 0, 0))

#compares diseased fragments to disease-exposed healthy corals at T7
DD_vs_DH_late <- list('DD_vs_DH_late' = c(0, 0, 0, 0, 1, -1, 0))

#set up two sided tests
two_sided_tests <- tibble(microbial_signature = c('simple_posthoc_aquar', 'simple_posthoc_exp', 'simple_posthoc_outc'),
                          contrasts = list(simple_posthoc_aquar, simple_posthoc_exp, simple_posthoc_outc),
                          direction = '=') #=
#set up one sided tests
right_tests <- tibble(microbial_signature = c('DD_early', 'DD_late', 'DD_vs_DH_early', 'DD_vs_DH_late'),
                      contrasts = list(DD_early, DD_late, DD_vs_DH_early, DD_vs_DH_late),
                      direction = '>') #>

#combine lists of contrasts
posthoc_categories <- bind_rows(two_sided_tests, right_tests) %>%
  unnest(contrasts) %>%
  mutate(contrast_name = names(contrasts)) %>%
  group_by(contrast_name, contrasts, direction) %>%
  summarise(signatures = list(c(microbial_signature)),
            .groups = 'drop') %>%
  nest(contrast = -direction)

#### Run LMER Models ####

if(file.exists('../intermediate_files/mixed_model_results.rds.gz') & !refit_models){
  asv_models <- read_rds('../intermediate_files/mixed_model_results.rds.gz')
} else {
  cluster_copy(cluster, c('posthoc_categories', 'run_posthoc'))
  
  asv_models <- normalized_asv_counts %>% 
    mutate(tank_field = if_else(str_detect(treatment, 'F'), 'field', 'tank')) %>% # categorize as tank vs. field to use in the model's dummy variable
    nest_by(across(c('asv_id', Domain:Species))) %>%
    partition(cluster) %>%
    mutate(fit_model(log2_cpm ~ treatment + (1 | genotype) +
                       (0 + dummy(tank_field, c("tank")) | tank), #dummy variable
                     data, 
                     use_weights = TRUE),
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

#### Process Post Hoc Comparisons ####

#229 ASVs are significant for treatment in the main model
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

#only ASVs significant for treatment
significant_models <- asv_models %>%
  filter(fdr_treatment < 0.05) %>% #must be significant for treatment
  rowwise() %>%
  mutate(process_postHoc(posthoc)) %>% #tidy up post hoc results
  ungroup() %>%
  p_adjust(exclude_cols = c('treatment', 'tank', 'genotype')) #FDR correct the posthoc p values

#Supplementary Data #5 - Model Outputs (T0, T3, T7)
#saveRDS(significant_models, "../Supplementary Data/ECT2025_supplementary_data_5.rds", compress = TRUE)

#interpret significance and directionality of main effects of tank, exposure, and outcome
simple_posthoc_sig_asvs <- significant_models %>%
  select(asv_id, starts_with('estimate'), starts_with('fdr')) %>%
  select(asv_id, contains(c('aquar', 'outcome', 'exposure'))) %>%
  mutate(across(starts_with('fdr'), ~. < 0.05), #significant or not
         across(starts_with('estimate'), ~if_else(. < 0, -1L, 1L))) %>% #set estimate to be either negative or positive, showing directionality
  pivot_longer(cols = -asv_id,
               names_to = c('.value', 'signatures', 'direction'),
               names_pattern = '(.*)_(.*)_(.*)') %>%
  dplyr::rename(significance = fdr) %>%
  #assign associations based on directionality
  mutate(signatures = case_when(signatures == 'outcome' & estimate < 0 ~ 'HealthyOutcome',
                                signatures == 'outcome' & estimate > 0 ~ 'DiseaseOutcome',
                                
                                signatures == 'aquarium' & estimate < 0 ~ 'Field',
                                signatures == 'aquarium' & estimate > 0 ~ 'Aquaria',
                                
                                signatures == 'exposure' & estimate < 0 ~ 'HealthyExposed',
                                signatures == 'exposure' & estimate > 0 ~ 'DiseaseExposed'),
         .keep = 'unused') %>%
  filter(significance) %>%
  select(-c(significance, direction))

#how many main model associations
simple_posthoc_sig_asvs %>%
  group_by(signatures) %>%
  reframe(n = n())


#format ASVs by what main effects and post hocs they're significant for
bacterial_signature_asv <- significant_models %>%
  select(asv_id, starts_with('fdr')) %>% 
  select(-contains(c('treatment', 'tank', 'genotype'))) %>% #remove main effects, look only at post hocs
  mutate(across(starts_with('fdr'), ~. < 0.05)) %>%
  pivot_longer(cols = -asv_id, names_to = c('term'), values_to = 'significance') %>%
  mutate(term = str_remove(term, 'fdr_')) %>%
  mutate(direction = str_extract(term, '[><=]'),
         term = str_remove(term, '_[><=]')) %>% #tidy names
  left_join(unnest(posthoc_categories, contrast), 
            by = c('term' = 'contrast_name', 'direction')) %>% #add in post hoc list
  select(-contrasts) %>%
  unnest(signatures) %>%
  group_by(asv_id, signatures) %>%
  filter(all(significance)) %>% #only keep significant post hocs
  ungroup %>%
  select(-term) %>%
  distinct() %>%
  #only keep ASVs that match one of the main effects: aquarium (tank), exposure, or outcome
  group_by(asv_id) %>%
  reframe(asv_id, signatures = str_c(signatures, collapse = ', ')) %>%
  distinct() %>%
  filter(str_detect(signatures, "simple")) %>% #main effects have simple in the name
  separate_rows(signatures, sep = ", ", convert = TRUE) %>% #separate back out
  #convert the simple posthocs to which direction they are
  filter(!str_detect(signatures, "simple")) %>% #remove the main effects
  rbind(simple_posthoc_sig_asvs) %>% #add the main effects back in with directionality indicated
  arrange(asv_id) #re-order by ASV number

#complex upset of significant signutures
bacterial_signature_asv %>%
  group_by(asv_id) %>%
  summarise(terms = list(signatures),
            .groups = 'drop') %>%
  ggplot(aes(x = terms)) +
  geom_bar() +
  scale_x_upset() +
  theme_classic() +
  theme_combmatrix(combmatrix.label.make_space = TRUE)

#add taxonomy and pivot significant ASVs wider
sig_classified_asvs <- select(significant_models, asv_id, Domain:Species) %>%
  left_join(bacterial_signature_asv %>%
              mutate(significance = TRUE) %>%
              pivot_wider(names_from = signatures,
                          values_from = significance,
                          values_fill = FALSE),
            by = 'asv_id') %>%
  filter_at(vars(DD_early:DD_vs_DH_late),all_vars(!is.na(.))) %>% #remove rows that are all NA, i.e. no significant post hocs
  mutate(DD_continuous = ifelse(DD_early & DD_late, TRUE, FALSE)) %>% #continuous is early + late
  mutate(DD_early = ifelse(DD_continuous, FALSE, DD_early), DD_late = ifelse(DD_continuous, FALSE, DD_late)) %>% #if DD_continuous 
#is true, set DD_early and DD_late to false to make it easier to see continuous ASVs
  mutate(DiseaseExposed = FALSE, HealthyExposed = FALSE) %>% #no exposure effect, manually made columns to include in complex upset
  mutate(TankEffect = ifelse(Field | Aquaria, TRUE, FALSE)) %>% #make column for tank effect, regardless of its directionality
  select(-c(Field, Aquaria)) #remove directional tank effect columns

write_csv(sig_classified_asvs, '../intermediate_files/classified_significant_asvs.csv.gz')

#how many ASVs significant for each contrast
sig_classified_asvs %>% 
  mutate(across(colnames(sig_classified_asvs %>% select(-c(asv_id, colnames(taxonomy_tibble)))), ~ifelse(. == TRUE, 1, 0))) %>%
  select(-colnames(taxonomy_tibble), -asv_id) %>%
  colSums(.)

#how many ASVs significant for each contrast after removing ASVs showing a significant effect of tank
sig_classified_asvs %>% 
  filter(!TankEffect) %>% #remove ASVs significant for tank effect
  mutate(across(colnames(sig_classified_asvs %>% select(-c(asv_id, colnames(taxonomy_tibble)))), ~ifelse(. == TRUE, 1, 0))) %>%
  select(-colnames(taxonomy_tibble), -asv_id) %>%
  colSums(.)

#view how many genera are significant for each combination of post hocs
sig_classified_asvs %>% 
  filter(!TankEffect) %>% 
  pivot_longer(cols = -c(colnames(taxonomy_tibble), asv_id), 
               names_to = "signature", values_to = "significance") %>%
  filter(significance) %>%
  select(-significance) %>%
  group_by(asv_id) %>%
  reframe(Family, Genus, all_sigs = str_c(signature, collapse = ", ")) %>%
  distinct() %>%
  ungroup() %>%
  group_by(Family, Genus, all_sigs) %>%
  reframe(n = n()) %>%
  arrange(all_sigs) %>% 
  mutate(taxa = str_c(n, Genus, sep = " ")) %>% 
  select(all_sigs, taxa)

#asv_ids for ASVs with significant post hocs
sig_classified_asvs %>% 
  filter(!TankEffect) %>% 
  pivot_longer(cols = -c(colnames(taxonomy_tibble), asv_id), 
               names_to = "signature", values_to = "significance") %>%
  filter(significance) %>%
  select(-significance) %>%
  group_by(asv_id) %>%
  reframe(asv_id, Family, Genus, all_sigs = str_c(signature, collapse = ", ")) %>%
  distinct() %>%
  ungroup() %>%
  group_by(Family, Genus, all_sigs)


#### Complex Upsets ####

#use same color palette as in relative abundance plot made using the package Fantaxtic
fantaxtic_palette <- read_csv("../intermediate_files/fantaxtic_color_palette.csv")

#add in placeholder rows for the bolded Order titles
ordered_fantaxtic_factoring <- fantaxtic_palette %>%
  group_by(Order) %>% 
  group_modify(~add_row(.x, .before = 0)) %>%
  ungroup() %>%
  mutate(subgroup_colour = ifelse(is.na(subgroup_colour), "#FFFFFF", subgroup_colour),
         Genus = ifelse(is.na(Genus), sprintf("**%s**", as.character(Order)), as.character(Genus)),
         group_subgroup = ifelse(is.na(group_subgroup), Genus, group_subgroup))

#add the bolded Order titles to the palette and format for a complex upset
full_palette <- ordered_fantaxtic_factoring$subgroup_colour
names(full_palette) <- ordered_fantaxtic_factoring$Genus

#prep data for complex upset
upset_sig_classified_asvs <- sig_classified_asvs %>% 
  left_join(fantaxtic_palette %>% select(-group_colour), by = join_by(Order, Genus)) %>%
  mutate(plot_genus = Genus) %>%
  mutate(plot_genus = case_when(is.na(subgroup_colour) & Order == "Alteromonadales" ~ "Other Alteromonadales",
                                is.na(subgroup_colour) & Order == "Flavobacteriales" ~ "Other Flavobacteriales",
                                is.na(subgroup_colour) & Order == "Oceanospirillales" ~ "Other Oceanospirillales",
                                is.na(subgroup_colour) & Order == "Rhodobacterales" ~ "Other Rhodobacterales",
                                is.na(subgroup_colour) & Order == "Saprospirales" ~ "Other Saprospirales",
                                is.na(subgroup_colour) & Order == "Spirochaetales" ~ "Other Spirochaetales",
                                is.na(subgroup_colour) & Order == "Thiotrichales" ~ "Other Thiotrichales",
                                is.na(subgroup_colour) & Order == "Verrucomicrobiales" ~ "Other Verrucomicrobiales",
                                is.na(subgroup_colour) & Order == "Vibrionales" ~ "Other Vibrionales",
                                is.na(subgroup_colour) & Order == "Rickettsiales" ~ "Other Rickettsiales",
                                TRUE ~ plot_genus)) %>% #set all genera that aren't the 4 most abundant to the "other" category for each genus
  mutate(plot_genus = ifelse(plot_genus %in% fantaxtic_palette$Genus, plot_genus, "Other")) %>%
  select(-subgroup_colour) %>%
  left_join(ordered_fantaxtic_factoring %>% select(Genus, subgroup_colour) %>% rename("plot_genus" = Genus), by = join_by(plot_genus)) %>% #add in colors
  relocate(c(plot_genus, subgroup_colour), .after = asv_id) %>%
  mutate(Order = ifelse(Order %in% fantaxtic_palette$Order, Order, "Other")) %>% #set Orders that aren't the top ones displayed here to "Other"
  group_by(Order) %>% 
  group_modify(~add_row(.x, .before = 0)) %>% #make empty row to set up bolded Order titles
  ungroup() %>%
  mutate(subgroup_colour = ifelse(is.na(subgroup_colour), "#FFFFFF", subgroup_colour),
         plot_genus = ifelse(is.na(plot_genus), sprintf("**%s**", as.character(Order)), as.character(plot_genus)),
         group_subgroup = ifelse(is.na(group_subgroup), plot_genus, group_subgroup)) %>% #format the bolded Order titles
  mutate(Order = ifelse(Order == "Puniceicoccales" & plot_genus == "Other" & subgroup_colour == "gray75", "Other", Order)) %>%
  mutate(outline = ifelse(str_detect(plot_genus, "\\*"), "no", "yes")) %>%
  mutate(across(DD_early:TankEffect, ~ifelse(is.na(.), FALSE, .))) %>% #set the newly added blank rows to FALSE so it doesn't interfere with the complex upset
  mutate(group_subgroup = ifelse(str_detect(group_subgroup, "Other") & !str_detect(group_subgroup, "\\*"), str_c(Order, " - ", plot_genus), group_subgroup)) %>%
  mutate(group_subgroup = factor(group_subgroup, ordered = T, levels = ordered_fantaxtic_factoring$group_subgroup),
         plot_genus = factor(plot_genus, ordered = T, levels = ordered_fantaxtic_factoring$Genus),
         Order = factor(Order, ordered = T)) %>%
  arrange(group_subgroup) %>%
  rename("Healthy Outcome" = HealthyOutcome, "Diseased Outcome" = DiseaseOutcome, 
         "Tank Effect" = TankEffect,
         "Healthy Exposure" = HealthyExposed, "Diseased Exposure" = DiseaseExposed)

#complex upset of main effects
upset(upset_sig_classified_asvs %>% select(-c(contains("DD"))),
      colnames(upset_sig_classified_asvs %>% select(-c(asv_id, subgroup_colour, plot_genus, outline, group_subgroup, colnames(taxonomy_tibble), contains("DD"))) %>%
                 relocate(`Tank Effect`, contains("Healthy"), contains("Outcome"))),
      base_annotations=list(
        'Intersection size'=intersection_size(counts=T, text = aes(size = 6.5),
                                              bar_number_threshold = 25,
                                              mapping=aes(fill=plot_genus, col = outline), col = "gray10"
        ) +
          scale_fill_manual(values = full_palette, limits = upset_sig_classified_asvs$plot_genus %>% unique()) +
          scale_color_manual(values = c("yes" = "gray10", "no" = "white")) +
          #theme_nested() +
          theme_nested(legend.position=c(1.26,0.2), 
                       axis.title = element_markdown(size = 14),
                       axis.text = element_markdown(size = 13),
                       legend.text = element_markdown(size = 12),
                       legend.title=element_text(size=rel(1.3))) +
          guides(fill=guide_legend(title=substitute(bold(bd)~nb, list(bd = "Order", nb = "/ Genus")),
                                   nrow = 25))
      ),
      
      matrix=(
        intersection_matrix(geom=geom_point(shape = "circle filled", size=4, stroke = 0.5, color = "gray10"))
        + scale_color_manual(
          values=c(
            "Healthy Outcome" = "#0d50ba",
            "Tank Effect"= "#c389e0",
            "Diseased Outcome"= "#ba0d0d"
          )
        )
      ),
      queries=list(
        upset_query(set = "Healthy Outcome", fill = "#0d50ba"),
        upset_query(set = "Tank Effect", fill = "#c389e0"),
        upset_query(set = "Diseased Outcome", fill = "#ba0d0d")
      ),
      name='ASVs', width_ratio=0.1, sort_sets = FALSE, sort_intersections=FALSE, 
      set_sizes=FALSE,
      themes=upset_modify_themes(
        list('intersections_matrix'=theme(axis.text=element_text(size=rel(1.3)),
                                          axis.title=element_text(size=rel(1.3))))
      ),
      intersections=list(
        'Diseased Outcome',
        'Healthy Outcome',
        c('Diseased Outcome', 'Tank Effect'),
        c('Healthy Outcome', 'Tank Effect'),
        'Tank Effect'
      )) +
  theme(plot.margin = margin(1, 11, 1, 1, "cm"))
#export 1500x900

#complex upset of main effects and post hocs with tank effect removed
upset(upset_sig_classified_asvs %>% filter(!`Tank Effect`) %>% group_by(Order) %>% filter(n() > 1), #removes headers from Orders that were filtered out
      colnames(upset_sig_classified_asvs %>% select(-c(asv_id, subgroup_colour, plot_genus, outline, group_subgroup, `Tank Effect`, colnames(taxonomy_tibble))) %>%
                 relocate(contains("Healthy"), contains("Outcome"))),
      base_annotations=list(
        'Intersection size'=intersection_size(counts=T, text = aes(size = 6.5),
                                              bar_number_threshold = 25,
                                              mapping=aes(fill=plot_genus, col = outline), col = "gray10"
        ) +
          scale_fill_manual(values = full_palette, limits = upset_sig_classified_asvs %>% filter(!`Tank Effect`) %>% group_by(Order) %>% filter(n() > 1) %>% pull(plot_genus) %>% unique()) +
          scale_color_manual(values = c("yes" = "gray10", "no" = "white")) +
          theme_nested(legend.position=c(1.06,0.2)) +
          guides(fill=guide_legend(title=substitute(bold(bd)~nb, list(bd = "Order", nb = "/ Genus")),
                                   ncol = 1))
      ),
      
      matrix=(
        intersection_matrix(geom=geom_point(shape = "circle filled", size=4, stroke = 0.5, color = "gray10"))
        + scale_color_manual(
          values=c(
            "Healthy Outcome" = "#0d50ba",
            "Diseased Outcome"= "#FF7D1A",
            "DD_early" = "gray80",
            "DD_late" = "gray30",
            "DD_vs_DH_late" = "#BA0D0D"
          )
        )
      ),
      queries=list(
        upset_query(set = "Healthy Outcome", fill = "#0d50ba"),
        upset_query(set = "Diseased Outcome", fill = "#FF7D1A"),
        upset_query(set = "DD_early", fill = "gray80"),
        upset_query(set = "DD_late", fill = "gray30"),
        upset_query(set = "DD_vs_DH_late", fill = "#BA0D0D")
      ),
      name='ASVs', width_ratio=0.1, sort_sets = FALSE, sort_intersections=FALSE, 
      set_sizes=FALSE,
      intersections=list(
        c('Healthy Outcome', 'DD_late'),
        'Healthy Outcome',
        c('Diseased Outcome', 'DD_late', 'DD_vs_DH_late'),
        c('Diseased Outcome', 'DD_late'),
        c('Diseased Outcome', 'DD_early')
      )) +
  ggtitle("Significant Main Effects") +
  theme(plot.margin = margin(1, 5.5, 1, 1, "cm"))

#### Abundance Over Time Plots ####

## Prep Data for Homogenate Dose Panel

#run homogenate linear models and prep data for plotting
homogenate_models <- homogenate_data %>%
  nest_by(asv_id) %>%
  mutate(model = list(lm(log2_cpm ~ exposure, data = data))) %>% #run lm model for homogenate doses
  mutate(homogenate_pred = list(model %>%
                                  broom::tidy(conf.int = TRUE) %>%
                                  mutate(term = case_when(term == "(Intercept)" ~ "D",
                                                          term == "exposureH" ~ "H")) %>%
                                  #D is the intercept, so you get the true values for H by adding the intercept value to the exposureH value:
                                  mutate(estimate = ifelse(term == "H", estimate[term == "H"] + estimate[term == "D"], estimate),
                                         std.error = NA,
                                         statistic = ifelse(term == "H", statistic[term == "H"] + statistic[term == "D"], statistic),
                                         conf.low = ifelse(term == "H", conf.low[term == "H"] + conf.low[term == "D"], conf.low),
                                         conf.high = ifelse(term == "H", conf.high[term == "H"] + conf.high[term == "D"], conf.high),
                                         p.value = ifelse(term == "D", p.value[term == "H"], p.value), #use the exposureH p-value for both
                                         df = NA) %>%
                                  rename("exposure" = "term") %>%
                                  mutate(time = "Dose", 
                                         graph_cat = "dose", 
                                         c_time = ifelse(exposure == "D", -1.2, -1.5), #offset the x value for H and D slightly so they don't overlap
                                         final_disease_state = "na",
                                         facet_lab = "Doses") %>%
                                  {. ->> set_one } %>% #save the data to this variable
                                  mutate(std.error = NA, df = NA, conf.low = NA, conf.high = NA, statistic = NA,
                                         c_time = c(-1.8, -0.9), final_disease_state = NA, exposure = NA) %>% #make 2 rows of NA values to increase the size of the Dose Facet so the real data points are centered in the facet
                                  rbind(set_one) %>% #add the real dose data back in
                                  {. ->> set_two } %>% #save the updated data to this variable
                                  arrange(desc(estimate)) %>% #arrange by decreasing estimate (y-values)
                                  dplyr::slice(1) %>% #keep the row with the highest estimate (y-value)
                                  mutate(exposure = "p_val", c_time = -1.35, estimate = estimate + 3) %>% #use this row with the highest y-value to set what height the stars for significance should be
                                  rbind(set_two) #add the dose data and NA rows back in
  )) %>%
  mutate(homog_p_val = homogenate_pred %>% filter(exposure == "H") %>% pull(p.value)) %>% #pull p-value
  ungroup() %>%
  mutate(adj_homog_p_val = p.adjust(homog_p_val, method = 'fdr')) %>% # FDR correct p-value
  unnest(homogenate_pred) %>%
  mutate(p.value = adj_homog_p_val) %>%
  nest(homogenate_pred = -c(asv_id, data, model, homog_p_val, adj_homog_p_val))

#Supplementary Data #6 - Model Outputs (homogenate doses)
# homogenate_data %>%
#   nest_by(asv_id) %>%
#   mutate(model = list(lm(log2_cpm ~ exposure, data = data))) %>% #run lm model for homogenate doses
#   mutate(model_stats = list(model %>%
#                                   broom::tidy(conf.int = TRUE))) %>%
#   mutate(p_val = model_stats %>% filter(term == "exposureH") %>% pull(p.value)) %>% #pull p-value
#   ungroup %>%
#   mutate(fdr_p_val = p.adjust(p_val, method = 'fdr')) %>% # FDR correct p-value
#   unnest(model_stats) %>%
#   nest(model_stats = -c(asv_id, data, model, p_val, fdr_p_val)) %>%
#   saveRDS("../Supplementary Data/ECT2025_supplementary_data_6.rds", compress = TRUE)

#34 ASVs are significantly different between the healthy and diseased doses
sig_homog_asv_list <- homogenate_models %>%
  filter(adj_homog_p_val < 0.05) %>%
  pull(asv_id)

## set up line showing normalized detection limit
add_zero_lines <- homogenate_models %>% ungroup() %>% select(homogenate_pred) %>% unnest(homogenate_pred) %>% 
  dplyr::slice(1:4) %>% mutate(across(everything(), ~NA)) %>% 
  mutate(facet_lab = c("Doses", "Experimental", "Doses", "Experimental"), c_time = c(-5, -5, 10, 10), 
         estimate = rep(detection_limit, 4)) %>%
  mutate(sig_p = NA, time = "zero_lines")


## Make Fancy Plots

abun_over_time_plots <- bacterial_signature_asv %>%
  group_by(asv_id) %>%
  reframe(signatures = str_c(signatures, collapse = ', ')) %>%
  #categorize ASVs based on post hoc signatures
  mutate(grouped_signatures = case_when(str_detect(signatures, "Aquaria") ~ "Tank Effect",
                                        str_detect(signatures, "Field") ~ "Tank Effect",
                                        str_detect(signatures, "DiseaseOutcome") & str_detect(signatures, "DD_early") & str_detect(signatures, "DD_vs_DH_early") & asv_id %in% sig_homog_asv_list ~ "Putative Early Pathogens",
                                        str_detect(signatures, "DiseaseOutcome") & str_detect(signatures, "DD_late") & str_detect(signatures, "DD_vs_DH_late") & asv_id %in% sig_homog_asv_list ~ "Putative Late Pathogens",
                                        str_detect(signatures, "DiseaseOutcome") & str_detect(signatures, "DD_late") & str_detect(signatures, "DD_vs_DH_late") & !asv_id %in% sig_homog_asv_list~ "WBD-Associated Opportunists",
                                        str_detect(signatures, "DiseaseOutcome") ~ "Unlikely Pathogens",
                                        str_detect(signatures, "HealthyOutcome") ~ "Healthy Associated",
                                        TRUE ~ signatures)) %>%
  #only keep ASVs of interest (remove tank-associated ASVs)
  filter(grouped_signatures %in% c("Putative Early Pathogens", "Putative Late Pathogens", "WBD-Associated Opportunists", "Healthy Associated")) %>%
  group_by(asv_id, grouped_signatures) %>%
  distinct() %>%
  inner_join(significant_models, by = 'asv_id') %>%
  mutate(across(Family:Species, ~ifelse(str_detect(., "NA"), NA, .))) %>% #set "NA" to NA
  #set the taxonomy to be shown on the graph as the lowest classified taxonomic level
  mutate(taxonomy_for_graph = case_when(is.na(Order) ~ Class,
                                        is.na(Family) ~ Order,
                                        is.na(Genus) ~ Family,
                                        is.na(Species) ~ str_c(Family, Genus, sep = " "),
                                        TRUE ~ str_c(Family, Species, sep = " "))) %>% #Species contains both genus and species
  mutate(formatted_ASV = str_c("(ASV ", parse_number(asv_id), ")", sep = "")) %>% #format the asv ID
  #format the taxonomy and asv ID to have the proper italics
  #if family is NA, use a higher rank and add sp. and the ASV number (no italics)
  #if species is NA, use family genus sp. (family and genus italicized)
  #if there's species-level classifications, use family genus species (all italicized)
  mutate(asv_taxa_labels = ifelse(is.na(Family),
                       list(substitute(taxonomy_for_graph~species_abbrev~formatted_ASV, 
                                       list(taxonomy_for_graph = taxonomy_for_graph, species_abbrev = "sp.", 
                                            formatted_ASV = formatted_ASV))),
                       ifelse(is.na(Species),
                              list(substitute(italic(taxonomy_for_graph)~species_abbrev~formatted_ASV, 
                                              list(taxonomy_for_graph = taxonomy_for_graph, species_abbrev = "sp.", 
                                                   formatted_ASV = formatted_ASV))),
                              list(substitute(italic(taxonomy_for_graph)~formatted_ASV, 
                                              list(taxonomy_for_graph = taxonomy_for_graph, 
                                                   formatted_ASV = formatted_ASV))))
  )) %>%
  rowwise() %>%
  mutate(plot_info = list(emmeans(model, ~treatment) %>%
                            broom::tidy(conf.int = TRUE) %>%
                            separate(treatment, into = c('time', 'exposure', 'final_disease_state')) %>%
                            #the next set of lines is adding in identical rows of the T0 points for each of the 3 combinations
                            #of exposure and outcome.  This is to allow identical T0 points to be plotted so that different
                            #colored lines for the different graph categories can be drawn from T0 to T3
                            mutate(graph_cat = ifelse(time == "T0", NA, 
                                                      paste(exposure, final_disease_state, sep = "_"))) %>%
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
  left_join(homogenate_models, by = join_by(asv_id)) %>%
  #add in homogenate dose model data
  mutate(plot_info = list(rbind(plot_info, homogenate_pred) %>% mutate(sig_p = ifelse(p.value < 0.05, "sig", "nonsig")))) %>%
  #add in data for the detection limit lines
  mutate(plot_info = list(rbind(plot_info, add_zero_lines))) %>%
  rowwise() %>%
  #make the plots
  mutate(plot = list(
    ggplot(data = plot_info, aes(x = c_time, y = estimate, ymin = conf.low, ymax = conf.high)) +
      geom_line(data = plot_info %>% filter(time == "zero_lines"), col = "gray45") + #detection limit lines
      
      # I used the relayer package to set up multiple legends
      # colour1 is disease exposed, colour2 is healthy exposed, colour3 is the homogenate doses
      
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
      
      #stars showing whether there is a significant difference between H and D doses
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
      #manually set up legend
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
      xlab(NULL) +
      ylab(NULL) +
      labs(title = asv_taxa_labels) +
      facet_grid(cols = vars(facet_lab), space = "free") + 
      scale_y_continuous(limits = c(1, 15), breaks = c(2.5, 5, 7.5, 10, 12.5, 15)) +
      coord_panel_ranges(panel_ranges = list(
        list(x=c(-1.8, -0.9)), # Dose Panel
        list(x=c(-0.75, 8)) # Experimental Panel
      ))
  )) %>%
  group_by(grouped_signatures) %>%
  #combine plots based on opportunist/pathogen signatures
  reframe(combined_plots = list(wrap_plots(plot, ncol = 2)), pathogen_plot = list(plot[asv_id == "ASV_65"]))

#extract legend
abun_over_time_legend <- ggpubr::get_legend(abun_over_time_plots$combined_plots[[1]])

#set up x and y axes labels
abun_x_axis <- cowplot::get_plot_component(ggplot() + labs(x = "Time"), "xlab-b")
abun_y_axis <- cowplot::get_plot_component(ggplot() + labs(y = expression("Normalized log"[2]*" (cpm)")), "ylab-l")

#empty space to go between plots
spacer <- ggplot() + theme_void()

#combine opportunists and putative pathogen into one plot
opp_path_design = "
EAAPL
ESSPL
EBBPL
#DDP#
"

list(
  abun_over_time_plots$pathogen_plot[[2]][[1]] + theme(legend.position="none"), # A
  abun_over_time_plots$combined_plots[[3]] + plot_layout(guides = "collect") & theme(legend.position = "none"), # B
  abun_x_axis, # D
  abun_y_axis, # E
  abun_over_time_legend, #L
  spacer, #P
  spacer #S
) %>% 
  wrap_plots() + 
  plot_layout(heights = c(1, 0.0325, 8, 0.0325), widths = c(0.25, 200, 200, 10, 50), design = opp_path_design) +
  plot_annotation(tag_levels = list(c("A", "B")))
#export 1000 x 950


#### Logfold Change Plots ####

#format ASVs to be included in logfold change plot
asvs_for_logfold <- bacterial_signature_asv %>%
  group_by(asv_id) %>%
  reframe(signatures = str_c(signatures, collapse = ', ')) %>%
  #categorize ASVs based on post hoc signatures
  mutate(grouped_signatures = case_when(str_detect(signatures, "Aquaria") ~ "Tank Effect",
                                        str_detect(signatures, "Field") ~ "Tank Effect",
                                        str_detect(signatures, "DiseaseOutcome") & str_detect(signatures, "DD_early") & str_detect(signatures, "DD_vs_DH_early") & asv_id %in% sig_homog_asv_list ~ "Putative Early Pathogens",
                                        str_detect(signatures, "DiseaseOutcome") & str_detect(signatures, "DD_late") & str_detect(signatures, "DD_vs_DH_late") & asv_id %in% sig_homog_asv_list ~ "Pathogens",
                                        str_detect(signatures, "DiseaseOutcome") & str_detect(signatures, "DD_late") & str_detect(signatures, "DD_vs_DH_late") & !asv_id %in% sig_homog_asv_list~ "White Band Disease- Associated Opportunists",
                                        str_detect(signatures, "DiseaseOutcome") ~ "Unlikely Pathogens",
                                        str_detect(signatures, "HealthyOutcome") ~ "Healthy-Associated",
                                        TRUE ~ signatures)) %>%
  #only keep ASVs of interest (remove tank-associated ASVs)
  filter(grouped_signatures %in% c("Putative Early Pathogens", "Pathogens", "White Band Disease- Associated Opportunists", "Unlikely Pathogens")) %>%
  group_by(asv_id, grouped_signatures) %>%
  inner_join(significant_models, by = 'asv_id') %>%
  mutate(across(Family:Species, ~ifelse(str_detect(., "NA"), NA, .))) #set "NA" to NA

#make list of ASVs that are late responders
late_responder_list <- significant_models %>%
  select(asv_id, `fdr_DD_late_>`) %>% #fdr_DD_late_> tests whether they are late responders
  rename("is_late_responder" = "fdr_DD_late_>") %>%
  mutate(is_late_responder = ifelse(is_late_responder < 0.05, TRUE, FALSE)) #set to true/false

#make list of stylized names with proper italics and ASV numbers
list_of_stylized_names <- normalized_asv_counts %>%
  mutate(across(Domain:Species, ~ifelse(. == "NA", NA, .))) %>%
  mutate(stylized_name = case_when(is.na(Genus) ~ glue("_{Family}_ sp. (ASV {parse_number(asv_id)})"), #if genus is NA, use family sp.
                                   !is.na(Genus) & is.na(Species) ~ glue("_{Family} {Genus}_ sp. (ASV {parse_number(asv_id)})"), #if species is NA, use family genus sp.
                                   TRUE ~ glue("_{Family} {Species}_ (ASV {parse_number(asv_id)})")), #otherwise, use family genus species
         .after = asv_id) %>%
  mutate(stylized_name_no_family = case_when(is.na(Genus) ~ glue("_{Family}_ sp. (ASV {parse_number(asv_id)})"), #if genus is NA, use family sp.
                                             !is.na(Genus) & is.na(Species) ~ glue("_{Genus}_ sp. (ASV {parse_number(asv_id)})"), #if species is NA, use genus sp.
                                             TRUE ~ glue("_{Species}_ (ASV {parse_number(asv_id)})")), #otherwise, use genus species
         .after = asv_id) %>%
  left_join(late_responder_list, by = join_by("asv_id")) %>%
  #set it so late responders will be bolded and early responders will not
  mutate(stylized_name = ifelse(is_late_responder, glue("**{stylized_name}**"), stylized_name),
         stylized_name_no_family = ifelse(is_late_responder, glue("**{stylized_name_no_family}**"), stylized_name_no_family)) %>%
  select(asv_id, stylized_name, stylized_name_no_family) %>%
  distinct()

#dose model
logfold_dose_model <- homogenate_data %>%
  filter(asv_id %in% asvs_for_logfold$asv_id) %>% #filter to ASVs of interest
  group_by(asv_id) %>%
  filter(!asv_id %in% c("ASV_372", "ASV_117", "ASV_244", "ASV_285")) %>% #data are essentially constant, can't run t-test
  summarise_at(vars(log2_cpm), list(dose_t_test = ~ list(t.test(. ~ exposure)))) %>% #run t-test
  rowwise() %>%
  mutate(t_stderr = dose_t_test$stderr, dose_t_test = broom::tidy(dose_t_test)) %>%
  unnest(dose_t_test) %>%
  rename("D" = "estimate1", "H" = "estimate2", "dose_D_minus_H" = "estimate") %>%
  mutate(fdr_dose = p.adjust(p.value)) %>% #FDR correct
  select(asv_id, dose_D_minus_H, fdr_dose, t_stderr) %>%
  pivot_longer(dose_D_minus_H, names_to = "comparison", values_to = "estimate") %>%
  rename("fdr_pval" = "fdr_dose")

#add back in ASVs that were too constant to run a t-test
logfold_dose_full_data <- logfold_dose_model %>% 
  dplyr::slice(1:4) %>% #keep first 4 rows
  mutate(across(everything(), ~NA)) %>% #make everything NA
  #set up these rows for the ASVs that were too constant for a t-test
  #set p-value to 1 and estimate for diseased dose minus healthy dose to 0
  mutate(asv_id = c("ASV_372", "ASV_117", "ASV_244", "ASV_285"), fdr_pval = c(1, 1, 1, 1), 
         comparison = c("dose_D_minus_H"), estimate = c(0, 0, 0, 0), t_stderr = c(0, 0, 0, 0)) %>%
  #combine with model results
  rbind(logfold_dose_model)


#get results from main model
logfold_model_metrics <- significant_models %>% select(asv_id, (contains("estimate") | contains("SE")) & #get estimates and standard errors
                                                         (contains("DD_vs_DH_late") | contains("exposure") | contains("outcome"))) %>% #for the expposure and outcome main effects and the DD vs DH late post hoc
  pivot_longer(cols = contains("estimate"), names_prefix = "estimate_", names_to = "comparison", values_to = "estimate") %>% #pivot estimates longer
  mutate(comparison = str_remove(comparison, '_[><=]')) %>% #remove directionality
  pivot_longer(cols = contains("SE"), names_prefix = "SE_", names_to = "comparison2", values_to = "t_stderr") %>% #pivot SE longer
  mutate(comparison2 = str_remove(comparison2, '_[><=]')) %>% #remove directionality
  filter(comparison == comparison2) %>% #pivoting twice made everything in triplicate, so only keep rows where the estimate and standard error are for the same comparison
  select(-comparison2)


model_for_plot <- significant_models %>%
  #pull whether the exposure and outcome main effects and the DD vs DH late post hoc were significant for the ASVs of interest
  select(asv_id, starts_with('fdr')) %>% 
  select(-contains(c('treatment', 'tank', 'genotype'))) %>%
  mutate(across(starts_with('fdr'), ~ifelse(. < 0.05, "sig", "nonsig"))) %>%
  pivot_longer(cols = -asv_id,
               names_to = c('comparison'),
               values_to = 'significance') %>%
  mutate(comparison = str_remove(comparison, 'fdr_')) %>%
  mutate(comparison = str_remove(comparison, '_[><=]')) %>%
  filter(asv_id %in% asvs_for_logfold$asv_id) %>% #filter for ASVs of interest
  filter(comparison %in% c("DD_vs_DH_late", "exposure", "outcome")) %>% #comparisons of interest
  left_join(logfold_model_metrics, by = join_by(asv_id, comparison)) %>% #add in estimates and standard errors
  rbind(logfold_dose_full_data %>%
          mutate(fdr_pval = ifelse(fdr_pval < 0.05, "sig", "nonsig"),
                 comparison = "dose") %>%
          rename("significance" = "fdr_pval")) %>% #add in dose information
  mutate(graph_pval = str_c(comparison, significance, sep = "_")) %>%
  filter(comparison != "exposure") %>%
  left_join(asvs_for_logfold %>% select(asv_id, grouped_signatures), by = join_by(asv_id)) %>% #add in grouped signature (i.e. pathogen or opportunist)
  mutate(grouped_signatures = ifelse(grouped_signatures == "Pathogens", "Putative Pathogen", grouped_signatures)) %>%
  mutate(grouped_signatures = factor(grouped_signatures, 
                                     levels = c("Putative Pathogen", "White Band Disease- Associated Opportunists", 
                                                "Unlikely Pathogens"))) %>%
  arrange(grouped_signatures) %>%
  left_join(list_of_stylized_names, by = join_by(asv_id)) #add in italicized names

#order ASVs by their logfold change value within their signature type to use to reorder y-axis of plot in decreasing order of logfold changes
asvs_ordered_by_outcome <- model_for_plot %>% 
  filter(comparison == "outcome") %>%
  group_by(grouped_signatures) %>%
  arrange(estimate, .by_group = TRUE) %>%
  ungroup()

#make asv_id an ordered factor
model_for_plot <- model_for_plot %>%
  mutate(asv_id = factor(asv_id, ordered = TRUE, levels = asvs_ordered_by_outcome$asv_id))


#all of the ASVs being plotted are significant for the outcome main effect
#for completeness sake of the legend, make a table where each of the 3 comparisons being examined are significant and 
#non-significant so they will show up in the legend
#pull this legend and use it for the actual logfold change plot
make_logfold_legend <- model_for_plot %>%
  filter(asv_id %in% c("ASV_65", "ASV_134")) %>%
  mutate(graph_pval = ifelse((asv_id == "ASV_134" & comparison != "dose_D_minus_H"), str_replace(graph_pval, "sig", "nonsig"), graph_pval)) %>% #set it so each comparison is both significant and non-significant
  mutate(asv_id = c(rep("test", 3), rep("test1", 3)), significance = NA, estimate = 0, t_stderr = 0, grouped_signatures = NA, stylized_name = NA, stylized_name_no_family = NA) #set all data fields to NA or 0

#make plot to pull legend from
relayer_logfold_change_for_legend <- ggplot(make_logfold_legend, aes(x = asv_id, y = estimate, ymin = estimate - t_stderr, ymax = estimate + t_stderr)) +
  geom_hline(yintercept = 0) +
  #outcome main effect
  (geom_point(data = (make_logfold_legend %>% filter(comparison == "outcome")), aes(colour1 = graph_pval, fill1 = graph_pval), pch = 23, size = 3, stroke = 1) %>%
     rename_geom_aes(new_aes = c("colour" = "colour1", "fill" = "fill1"))) + 
  scale_color_manual(aesthetics = "colour1", values = c("#FF7D1A", "#FF7D1A"), guide = "legend", 
                     name = "Outcome", breaks = c("outcome_sig", "outcome_nonsig"), labels = c("Significant", "Nonsignificant")) +
  scale_fill_manual(aesthetics = "fill1", values = c("#FF7D1A", "transparent"), guide = "legend", 
                    name = "Outcome", breaks = c("outcome_sig", "outcome_nonsig"), labels = c("Significant", "Nonsignificant")) +
  #homogenate dose
  (geom_point(data = (make_logfold_legend %>% filter(comparison == "dose")), aes(colour3 = graph_pval, fill3 = graph_pval), pch = 22, size = 4, stroke = 1) %>%
     rename_geom_aes(new_aes = c("colour" = "colour3", "fill" = "fill3"))) + 
  scale_color_manual(aesthetics = "colour3", values = c("#732dcf", "#732dcf"), guide = "legend", 
                     name = "Homogenate Doses", breaks = c("dose_sig", "dose_nonnonsig"), labels = c("Significant", "Nonsignificant")) +
  scale_fill_manual(aesthetics = "fill3", values = c("#732dcf", "transparent"), guide = "legend", 
                    name = "Homogenate Doses", breaks = c("dose_sig", "dose_nonnonsig"), labels = c("Significant", "Nonsignificant")) +
  #exposure within disease-exposed tanks post hoc
  (geom_point(data = (make_logfold_legend %>% filter(comparison == "DD_vs_DH_late")), aes(colour4 = graph_pval, fill4 = graph_pval), pch = 24, size = 3, stroke = 1) %>%
     rename_geom_aes(new_aes = c("colour" = "colour4", "fill" = "fill4"))) + 
  scale_color_manual(aesthetics = "colour4", values = c("#BA0D0D", "#BA0D0D"), guide = "legend", 
                     name = str_wrap("Outcome within Disease-Exposed Tanks", 20), breaks = c("DD_vs_DH_late_sig", "DD_vs_DH_late_nonsig"), labels = c("Significant", "Nonsignificant")) +
  scale_fill_manual(aesthetics = "fill4", values = c("#BA0D0D", "transparent"), guide = "legend", 
                    name = str_wrap("Outcome within Disease-Exposed Tanks", 20), breaks = c("DD_vs_DH_late_sig", "DD_vs_DH_late_nonsig"), labels = c("Significant", "Nonsignificant")) +
  coord_flip() +
  theme_bw() +
  guides(#outcome
    color1 = guide_legend(order = 2),
    fill1 = guide_legend(order = 2),
    #dose
    color3 = guide_legend(order = 4),
    fill3 = guide_legend(order = 4),
    #DD_DH
    color4 = guide_legend(order = 3),
    fill4 = guide_legend(order = 3))

#use this legend for logfold change plot
relayer_logfold_legend <- ggpubr::get_legend(relayer_logfold_change_for_legend) 

#make logfold change plot
logfold_change_plot <- ggplot(model_for_plot) +
  geom_hline(yintercept = 0) +
  geom_pointrange(aes(x = asv_id, y = estimate, fill = graph_pval, col = comparison, pch = comparison, ymin = estimate - t_stderr, ymax = estimate + t_stderr), 
                  size = 1, stroke = 1, position = position_dodge(width = 0.65), linewidth = 0.8) +
  scale_shape_manual(values = c(24, 22, 23)) +
  scale_color_manual(values = c("outcome" = "#FF7D1A",
                                "dose" = "#732dcf", 
                                "DD_vs_DH_late" = "#BA0D0D")) +
  scale_fill_manual(values = c("outcome_nonsig" = "white", 
                               "outcome_sig" = "#FF7D1A",
                               "exposure_nonsig" = "white", 
                               "dose_nonsig" = "white", 
                               "dose_sig" = "#732dcf", #DF369B
                               "DD_vs_DH_late_nonsig" = "white", 
                               "DD_vs_DH_late_sig" = "#BA0D0D"), guide = "none") +
  coord_flip() +
  theme_bw() +
  theme(axis.text = element_markdown(size = 10),
        axis.title = element_markdown(size = 12),
        legend.position = "none",
        strip.text = element_text(size = 12),
        strip.text.y = element_text(angle = 0)
  ) +
  scale_x_discrete(breaks = model_for_plot$asv_id, labels = model_for_plot$stylized_name_no_family) +
  facet_grid(rows = vars(grouped_signatures), space = "free", scales = "free",
             labeller = label_wrap_gen(width = 15)) +
  xlab(NULL) +
  ylab("Log Fold Change (Diseased vs. Healthy)")


#combine plot and legend
logfold_design = "
A##
A#B
A##
"

list(
  logfold_change_plot, # A
  relayer_logfold_legend # B
) %>% 
  wrap_plots() + 
  plot_layout(heights = c(0.125, 1, 0.125), widths = c(200, 5, 50), design = logfold_design)
#export 1050x700

#### NMDS ####

set.seed(68748)

#make matrix for NMDS, ASVs on columns and samples on rows
nmds_matrix <- normalized_asv_counts %>% 
  select(asv_id, sample_id, log2_cpm) %>%
  pivot_wider(names_from = sample_id, values_from = log2_cpm) %>%
  column_to_rownames('asv_id') %>%
  t() %>%
  as.matrix()

#run NMDS
asv_nmds <- metaMDS(nmds_matrix, distance = 'bray', k = 3, trymax = 100, autotransform = FALSE, verbose = TRUE)

#get NMDS scores
nmds_scores <- as.data.frame(scores(asv_nmds)$sites)

#format metadata to match NMDS matrix
nmds_metadata <- normalized_asv_counts %>% 
  select(asv_id, time, exposure, final_disease_state, susceptability, genotype, 
         sample_id, log2_cpm, tank, treatment) %>%
  pivot_wider(names_from = asv_id, values_from = log2_cpm) %>%
  as.data.frame()

#add columns for each relevant piece of metadata to the scores data frame
nmds_scores$time = nmds_metadata$time
nmds_scores$exposure = nmds_metadata$exposure
nmds_scores$final_disease_state = nmds_metadata$final_disease_state
nmds_scores$susceptability = nmds_metadata$susceptability
nmds_scores$genotype = nmds_metadata$genotype
nmds_scores$treatment = str_c(nmds_scores$time, nmds_scores$exposure, nmds_scores$final_disease_state, sep = "_")

#plot NMDS
ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2)) +
  #Day 0
  (geom_point(data = (nmds_scores %>% filter(time == "T0")), aes(colour0 = treatment, pch = time)) %>%
     rename_geom_aes(new_aes = c("colour" = "colour0"))) +
  (stat_ellipse(data = (nmds_scores %>% filter(time == "T0")), aes(colour0 = treatment), show.legend=FALSE) %>%
     rename_geom_aes(new_aes = c("colour" = "colour0"))) +
  scale_color_manual(aesthetics = "colour0", values = c("T0_F_F" = "#c389e0"), guide = "legend", 
                     name = "Day 0", breaks = c("T0_F_F"), labels = c("Field-Collected")) +
  #Day 3
  (geom_point(data = (nmds_scores %>% filter(time == "T3")), aes(colour3 = treatment, pch = time)) %>%
     rename_geom_aes(new_aes = c("colour" = "colour3"))) +
  (stat_ellipse(data = (nmds_scores %>% filter(time == "T3")), aes(colour3 = treatment), show.legend=FALSE) %>%
     rename_geom_aes(new_aes = c("colour" = "colour3"))) +
  scale_color_manual(aesthetics = "colour3", values = c("T3_D_D" = "#E79B9B", "T3_D_H" = "#97D9E1",
                                                        "T3_H_H" = "#95AC85"), guide = "legend", 
                     name = "Day 3", breaks = c("T3_H_H", "T3_D_H", "T3_D_D"), 
                     labels = c("Healthy", "Disease-Exposed, Healthy", "Eventually Develops Disease")) +
  #Day 7
  (geom_point(data = (nmds_scores %>% filter(time == "T7")), aes(colour7 = treatment, pch = time)) %>%
     rename_geom_aes(new_aes = c("colour" = "colour7"))) +
  (stat_ellipse(data = (nmds_scores %>% filter(time == "T7")), aes(colour7 = treatment), show.legend=FALSE) %>%
     rename_geom_aes(new_aes = c("colour" = "colour7"))) +
  scale_color_manual(aesthetics = "colour7", values = c("T7_D_D" = "#A70000", "T7_D_H" = "#22A7B6",
                                                        "T7_H_H" = "#406F23"), guide = "legend", 
                     name = "Day 7", breaks = c("T7_H_H", "T7_D_H", "T7_D_D"), 
                     labels = c("Healthy", "Disease-Exposed, Healthy", "Diseased")) +
  scale_shape_manual(values = c("T0" = 16, "T3" = 15, "T7" = 17), guide = "none") +
  guides(color0 = guide_legend(
    override.aes=list(shape = c("T0_F_F" = 16),
                      size = 3), order = 1),
    color3 = guide_legend(
      override.aes=list(shape = c("T3_D_D" = 15,
                                  "T3_D_H" = 15,
                                  "T3_H_H" = 15),
                        size = 3), order = 3),
    color7 = guide_legend(
      override.aes=list(shape = c(
        "T7_D_D" = 17,
        "T7_D_H" = 17,
        "T7_H_H" = 17),
        size = 3), order = 7)) +
  theme_bw()
#export 1000x550

#run permanova
permanova_results <- adonis2(nmds_matrix ~time*final_disease_state*exposure + genotype + tank, 
                             method = "bray", perm = 10000, data = nmds_metadata)

#tidy statistics
tidy_permanova_results <- broom::tidy(permanova_results) %>%
  rename("Sums of Squares" = "SumOfSqs") %>%
  mutate(term = case_when(term == "time" ~ "Time",
                          term == "final_disease_state" ~ "Disease Outcome",
                          term == "exposure" ~ "Exposure",
                          term == "genotype" ~ "Genotype",
                          term == "time:final_disease_state" ~ "Time:Disease Outcome",
                          term == "time:exposure" ~ "Time:Exposure",
                          term == "tank" ~ "Tank",
                          TRUE ~ term)) %>%
  rename("F" = "statistic")

#output statistics to a table
permanova_table <- nice_table(tidy_permanova_results, broom = "lm")
print(permanova_table, preview = "docx")

