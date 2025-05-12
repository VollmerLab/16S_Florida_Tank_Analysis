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



library(vegan)

library(EcolUtils)
library(formattable)
library(RColorBrewer)
library(magrittr)
library(emmeans)
library(ggupset)
library(qvalue)
library(relayer) #devtools::install_github("clauswilke/relayer")

library(corrplot)
library(Hmisc)
library(broom.mixed)
library(taxize)

library(rempsyc)

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

#### Miscellaneous Manuscript Metrics ####

#range of disease resistance scores for T0 corals
normalized_asv_counts %>%
  filter(time == "T0") %>%
  select(genotype, resistance) %>%
  distinct() %>%
  reframe(num_genotypes = length(genotype), minimum_DR = min(resistance), maximum_DR = max(resistance),
          average_DR = mean(resistance), SE_DR = sd(resistance)/sqrt(length((resistance))))

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

#detection limit of 5.202306
detection_limit <- normalized_asv_counts %>% select(log2_cpm) %>% min()

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

#### Main Effects Upset ####

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
          theme_nested(legend.position=c(1.06,0.2)) +
          guides(fill=guide_legend(title=substitute(bold(bd)~nb, list(bd = "Order", nb = "/ Genus")),
                                   ncol = 1))
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
      intersections=list(
        'Diseased Outcome',
        'Healthy Outcome',
        c('Diseased Outcome', 'Tank Effect'),
        c('Healthy Outcome', 'Tank Effect'),
        'Tank Effect'
      )) +
  ggtitle("Significant Main Effects") +
  theme(plot.margin = margin(1, 5.5, 1, 1, "cm"))


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

#export 2000x1100









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
  #filter(asv_id == "ASV_65") %>%
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
  #don't bother plotting tank effect, etc...
  filter(grouped_signatures %in% c("Putative Early Pathogens", "Putative Late Pathogens", "Specialized Opportunists", "Healthy Associated")) %>%
  group_by(asv_id, grouped_signatures) %>%
  distinct() %>%
  inner_join(significant_models,
             by = 'asv_id') %>%
  mutate(across(Family:Species, ~ifelse(str_detect(., "NA"), NA, .))) %>%
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
      #theme(legend.position="none") + #remove legend for combining plots
      #theme(plot.title = element_text(face = "italic")) +
      #remove x and y labels for combo plots
      xlab(NULL) +
      ylab(NULL) +
      #xlab("Time") +
      #ylab(expression("Normalized log"[2]*" (cpm)")) +
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
  reframe(combo_plots = list(wrap_plots(plot, ncol = 2)), path_plot = list(plot[asv_id == "ASV_65"]))

the_plots$path_plot[[1]]

#old way  
#summarise(combo_plots = list(wrap_plots(plot))) %>%
#pull(combo_plots) + plot_layout(guides = "none") & plot_annotation(title = grouped_signatures)))


#view the plots
the_plots$combo_plots[[1]]

plot_for_legend <- the_plots$combo_plots[[1]]

#combine plots
library(ggpubr)

the_plots$path_plot[[2]][[1]] # 700 x 450


x_axis <- cowplot::get_plot_component(ggplot() + labs(x = "Time"), "xlab-b")
y_axis <- cowplot::get_plot_component(ggplot() + labs(y = expression("Normalized log"[2]*" (cpm)")), "ylab-l")

fancy_legend <- ggpubr::get_legend(plot_for_legend)

design = "
EAAL
#HHL
FBBL
#IIL
GCCL
#DDL
"

list(
  the_plots$path_plot[[2]][[1]], # A
  the_plots$combo_plots[[1]], # B
  the_plots$combo_plots[[3]], # C
  x_axis, # D
  y_axis, # E
  y_axis, # F
  y_axis, # G
  x_axis, # H
  x_axis, # I
  fancy_legend #L
) %>% 
  wrap_plots() + 
  plot_layout(heights = c(1, 0.125, 1, 0.125, 10, 0.125), widths = c(0.25, 200, 200, 50), design = design) +
  plot_annotation(tag_levels = list(c("A", "B", "", "C")))

design = "
EAAL
EBBL
ECCL
#DD#
"

list(
  the_plots$path_plot[[2]][[1]], # A
  the_plots$combo_plots[[1]], # B
  the_plots$combo_plots[[3]], # C
  x_axis, # D
  y_axis, # E
  fancy_legend #L
) %>% 
  wrap_plots() + 
  plot_layout(heights = c(1, 1, 8, 0.125), widths = c(0.25, 200, 200, 50), design = design) +
  plot_annotation(tag_levels = list(c("A", "B", "", "C")))
#export to 1300x1700

spacer <- ggplot() + theme_void()

ggplot() + fancy_legend

list(
  the_plots$path_plot[[2]][[1]], # A
  the_plots$combo_plots[[1]], # B
  the_plots$combo_plots[[3]], # C
  x_axis, # D
  y_axis, # E
  y_axis, # F
  y_axis, # G
  x_axis, # H
  x_axis, # I
  fancy_legend #L
) %>% 
  wrap_plots() + 
  plot_layout(heights = c(1, 0.125, 1, 0.125, 10, 0.125), widths = c(0.25, 200, 200, 50), design = design) +
  plot_annotation(tag_levels = list(c("A", "B", "", "C")))

design = "
EAAL
EPPL
ECCL
ESSL
EBBL
#DD#
"

list(
  the_plots1$path_plot[[2]][[1]], # A
  the_plots1$combo_plots[[1]], # B
  the_plots1$combo_plots[[3]], # C
  x_axis, # D
  y_axis, # E
  fancy_legend, #L
  spacer, #P
  spacer #S
) %>% 
  wrap_plots() + 
  plot_layout(heights = c(1, 0.0325, 8, 0.0325, 1, 0.0125), widths = c(0.25, 200, 200, 50), design = design) +
  plot_annotation(tag_levels = list(c("A", "C", "", "B")))

# no healthies
design = "
EAAPL
ESSPL
EBBPL
#DDP#
"

list(
  the_plots$path_plot[[2]][[1]] + theme(legend.position="none"), # A
  the_plots$combo_plots[[3]] + plot_layout(guides = "collect") & theme(legend.position = "none"), # C
  x_axis, # D
  y_axis, # E
  fancy_legend, #L
  spacer, #P
  spacer #S
) %>% 
  wrap_plots() + 
  plot_layout(heights = c(1, 0.0325, 8, 0.0325), widths = c(0.25, 200, 200, 10, 50), design = design) +
  plot_annotation(tag_levels = list(c("A", "B")))

#export 1000 x 950


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

nmds_scores <-  as.data.frame(scores(asv_nmds)$sites)

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

#Nice legend separated out by day

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
                     labels = c("Healthy Control", "Exposed, Healthy", "Infected")) +
  #Day 7
  (geom_point(data = (nmds_scores %>% filter(time == "T7")), aes(colour7 = treatment, pch = time)) %>%
     rename_geom_aes(new_aes = c("colour" = "colour7"))) +
  (stat_ellipse(data = (nmds_scores %>% filter(time == "T7")), aes(colour7 = treatment), show.legend=FALSE) %>%
     rename_geom_aes(new_aes = c("colour" = "colour7"))) +
  scale_color_manual(aesthetics = "colour7", values = c("T7_D_D" = "#A70000", "T7_D_H" = "#22A7B6",
                                                        "T7_H_H" = "#406F23"), guide = "legend", 
                     name = "Day 7", breaks = c("T7_H_H", "T7_D_H", "T7_D_D"), 
                     labels = c("Healthy Control", "Exposed, Healthy", "Infected")) +
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



permanova_results <- adonis2(nmds_matrix ~time*final_disease_state*exposure + genotype + tank, 
                             method = "bray", perm = 10000, data = nmds_metadata)  

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

permanova_table <- nice_table(tidy_permanova_results, broom = "lm")

print(permanova_table, preview = "docx")
#flextable::save_as_docx(permanova_table, path = "../")

#library(pairwiseAdonis)
#pairwise.adonis2(nmds_matrix ~treatment, data = nmds_metadata) %>% broom::tidy()



