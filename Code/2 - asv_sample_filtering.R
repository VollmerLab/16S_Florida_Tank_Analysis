setwd("/Users/emilytrytten/Desktop/GitHub/16S_Florida_Tank_Analysis/Code")

#### Libraries ####
library(magrittr)
library(phyloseq)
library(microbiome)
library(vegan)
library(metagMisc)
library(edgeR)
library(variancePartition)
library(ggvenn)
library(cowplot)
library(ComplexUpset)
library(microshades)
library(lme4)
library(tidyverse)
library(emmeans)
library(relayer)
library(fantaxtic)
library(ggnested)

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

#### Data ####
aggregation_level <- 'none' #or none

microbiome_data <- read_rds("../intermediate_files/preprocess_microbiome.rds")
metadata <- sample_data(microbiome_data) %>%
  as_tibble(rownames = 'sample_id') %>%
  select(-retain_sample)

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

#### Aggregate Samples ####
if(aggregation_level != 'none'){
  microbiome_data <- aggregate_taxa(microbiome_data, aggregation_level)
  taxa_names(microbiome_data) <- str_replace_all(taxa_names(microbiome_data), ' |-', '_')
} else {
  sequences <- taxa_names(microbiome_data)
  taxa_names(microbiome_data) <- str_c('ASV', 1:length(taxa_names(microbiome_data)), sep = '_')
  names(sequences) <- taxa_names(microbiome_data)
}

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

#get old taxonomy info for ASVs that are all NA in updated taxonomy
all_na_in_new_tax <- updated_taxonomy_above80 %>% 
  filter(if_all(Domain:Species, ~is.na(.))) %>% 
  select(-c(Domain:Species)) %>%
  left_join(taxonomy_tibble, by = join_by("asv_id" == "asv_names"))

#combine the old taxonomy with the updated version
combined_taxonomy <- updated_taxonomy_above80 %>% 
  filter(!if_all(Domain:Species, ~is.na(.))) %>%
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
ncbi_classifications_by_genus <- read_csv("../intermediate_files/ncbi_classifications_for_overlaps.csv")

#most updated taxonomy, just phylum:genus
#each genus has only one described classification, is combined with old and new taxonomies
#doesnt contain NA for genus rows
nonoverlapping_taxonomy <- combined_taxonomy %>% 
  select(Phylum:Genus) %>% 
  distinct() %>% 
  filter(!is.na(Genus)) %>%
  filter(!Genus %in% multiple_classifications_list) %>%
  rbind(ncbi_classifications_by_genus)

ncbi_classifications_by_genus %>% filter(Genus == "Nannocystis")

#our ASVs with the most up to date taxonomy, use for downstream purposes
full_taxonomy <- combined_taxonomy %>%
  select(-c(Phylum:Family)) %>%
  left_join(nonoverlapping_taxonomy, by = join_by(Genus)) %>%
  relocate(Phylum:Family, .after = Domain) %>%
  filter(!asv_id %in% c(combined_taxonomy %>% filter(is.na(Genus)) %>% pull(asv_id))) %>%
  rbind(combined_taxonomy %>% filter(is.na(Genus))) %>%
  arrange(parse_number(asv_id)) %>%
  filter(!(asv_id %in% c("ASV_121", "ASV_1535", "ASV_1909", "ASV_4075", "ASV_5598") & 
             Class %in% c("Gammaproteobacteria", "Alphaproteobacteria") & Phylum != "Pseudomonadota")) %>% #remove duplicates of these w diff phyla
  filter(!(asv_id == "ASV_121" & is.na(Class))) %>%
  distinct()

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

#### Alpha Diversity ####
# pre-filtering for abundance

alpha_table <- microbiome::alpha(updated_microbiome_data, index = "all") %>%
  as_tibble(rownames = 'sample_id') %>%
  inner_join(metadata, by = 'sample_id') %>%
  mutate(fragment_id = str_c(str_replace_na(exposure, 'NA'), tank, genotype, sep = '_'))

# mod_alpha_tab <- alpha_table %>%
#   filter(!tank %in% c("HOMO", "homogenate_fragment")) %>%
#   mutate(final_disease_state = ifelse(exposure == "F", "F", final_disease_state)) %>%
#   mutate(treatment = str_c(time, exposure, final_disease_state, sep = '_')) %>%
#   pivot_longer(cols = !c(colnames(metadata), "fragment_id", "treatment"),
#                names_to = 'metric',
#                values_to = 'alpha_div_value') %>%
#   select(-c(susceptability, resistance, clone_group)) %>%
#   mutate(tank_field = if_else(str_detect(treatment, 'F'), 'field', 'tank'), .after = final_disease_state) %>%
#   nest_by(metric) %>%
#   rowwise() %>%
#   summarise(alpha_model = list(lmer(alpha_div_value ~ treatment + 
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

mod_alpha_tab <- read_rds("../intermediate_files/mod_alpha_tab.rds")

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

#### Beta Diversity ####

nmds_matrix <- melted_ps %>% 
  filter(Sample %in% model_samples) %>%
  select(OTU, Sample, Abundance) %>%
  pivot_wider(names_from = Sample, values_from = Abundance) %>%
  column_to_rownames('OTU') %>%
  t() %>%
  as.matrix()

asv_nmds <- metaMDS(nmds_matrix, distance = 'mountford', k = 2, trymax = 100, autotransform = FALSE, verbose = TRUE)

#shepard plot
#doesn't follow y = x line well, so MDS might be misleading -> use nmds
plot(asv_nmds$diss, asv_nmds$dist)
abline(a = 0, b = 1, col = "red")

stressplot(asv_nmds)

plot(asv_nmds)

nmds_scores = as.data.frame(scores(asv_nmds)$sites)

nmds_metadata <- melted_ps %>% 
  filter(Sample %in% model_samples) %>%
  select(OTU, time, exposure, final_disease_state, genotype, 
         Sample, Abundance, tank) %>%
  mutate(final_disease_state = ifelse(exposure == "F", "F", final_disease_state)) %>%
  mutate(tank_field = if_else(exposure == "F", 'field', 'tank')) %>%
  pivot_wider(names_from = OTU, values_from = Abundance) %>%
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


adonis2(nmds_matrix ~time*final_disease_state*exposure + genotype + tank, method = "mountford", perm = 999, data = nmds_metadata)

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

#custom

top_level <- "Order"
nested_level <- "Genus"
sample_order <- NULL
top_asv <- top_taxa(updated_microbiome_data, n_taxa = 10, nested_level, by_proportion = TRUE)

top_nested <- nested_top_taxa(updated_microbiome_data,
                              top_tax_level = "Order",
                              nested_tax_level = "Genus",
                              n_top_taxa = 10, 
                              n_nested_taxa = 4)

top_asv <- top_nested

# Generate a palette based  on the phyloseq object
pal <- taxon_colours(top_asv$ps_obj,
                     tax_level = top_level)

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

# Generate a base plot
ggnested(psdf,
              aes_string(main_group = top_level,
                         sub_group = nested_level,
                         x = "Sample",
                         y = "Abundance"),
              main_palette = pal) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90)) +
  theme_nested(theme_bw) + 
  geom_col(position = position_fill()) + 
  facet_grid(col = vars(facet_level), space = "free", scales = "free") + 
  scale_x_discrete(name = element_blank(), labels = c("Field" = "Field", "D Dose" = "Diseased", "H Dose" =
                                                        "Healthy", "T3 Diseased" = str_wrap("Disease- Exposed", 10), "T3_H_H" = str_wrap("Healthy- Exposed", 10),
                                                      "T7_H_H" = str_wrap("Healthy Control", 10),
                                                      "T7_D_H" = str_wrap("Disease- Exposed Healthy", 10), "T7_D_D" = "Diseased")) +
  guides(fill=guide_legend(title=substitute(bold(bd)~nb, list(bd = "Order", nb = "/ Genus"))),
         col=guide_legend(title=substitute(bold(bd)~nb, list(bd = "Order", nb = "/ Genus")))) +
  ylab("Relative Abundance")


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
  filter_samples(model_samples) %>% #remove samples not to be analyzed
  filter_asv_meanCount(metadata, 100) #Remove ASVs with an average of less than N CPM per sample


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
