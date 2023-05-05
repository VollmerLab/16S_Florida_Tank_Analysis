#Code that streamlines the process from Script #5: Disease-Associated Bacteria

setwd("~/Desktop/Screenshots/Career/Vollmer Lab/GitHub/16S_Florida_Tank_Analysis/Code")


#### Packages ####
library(ggvenn)
library(multcomp)
library(phyloseq)
library(microbiome)
library(vegan)
library(lme4)
library(afex)
library(emmeans)
library(car)
library(edgeR)
library(metagMisc)
library(ape)
library(ggdist)
library(gghalves)
library(patchwork)
library(magrittr)
library(ComplexUpset)
library(strex)
library(forcats)
library(wesanderson)
library(msa)
library(Biostrings)
library(phangorn)
library(formattable)
library(htmltools)
library(webshot)
library(tidyverse)

select <- dplyr::select

#### Functions ####

make_emmean_model <- function(model, form, alpha){
  emmeans(model, specs = form, type = 'response') %>%
    cld(Letters = LETTERS, reversed = TRUE, alpha = alpha) 
}

export_formattable <- function(f, file, width = "100%", height = NULL, 
                               background = "white", delay = 0.2){
  w <- as.htmlwidget(f, width = width, height = height)
  path <- html_print(w, background = background, viewer = NULL)
  url <- paste0("file:///", gsub("\\\\", "/", normalizePath(path)))
  webshot(url,
          file = file,
          selector = ".formattable_widget",
          delay = delay)
}

#### Read in Data ####

aggregation_level <- 'none' #or none

microbiome_data <- read_rds("../intermediate_files/preprocess_microbiome.rds")
metadata <- sample_data(microbiome_data) %>%
  as_tibble(rownames = 'sample_id') %>%
  select(-retain_sample)

if(aggregation_level != 'none'){
  microbiome_data <- aggregate_taxa(microbiome_data, aggregation_level)
  taxa_names(microbiome_data) <- str_replace_all(taxa_names(microbiome_data), ' |-', '_')
} else {
  taxa_names(microbiome_data) <- str_c('ASV', 1:length(taxa_names(microbiome_data)), sep = '_')
}

#otu_tmm
otu_tmm <- microbiome_data %>%
  phyloseq_filter_prevalence(prev.trh = 0.1) %>%
  otu_table() %>% 
  t %>% #NOTE: *genus and family do not need the t but ASVs need the t*
  as.data.frame %>%
  as.matrix %>% 
  DGEList(remove.zeros = TRUE) %>%
  edgeR::calcNormFactors(method = 'TMMwsp')

#raw model data - contains all 0 values
model_data <- cpm(otu_tmm, log = TRUE, prior.count = 2) %>%
  t %>%
  as_tibble(rownames = "sample_id") %>%
  full_join(metadata, by = "sample_id") %>%
  pivot_longer(cols = -any_of(colnames(metadata)), 
               names_to = "asv_names", values_to = "value") %>%
  mutate(across(c(exposure, final_disease_state), factor)) %>%
  filter(time %in% c('T3', 'T7') | (time == "T0" & tank == "HOMO")) %>%
  mutate(fragment_id = str_c(exposure, tank, genotype, final_disease_state)) %>%
  mutate(fragment_id = if_else(time == 'T0', 'homogenate', fragment_id))

#target microbiome data - zeroes considered and each ASV must be in 10+% of individuals
target_data <- cpm(otu_tmm, log = TRUE, prior.count = 2) %>%
  t %>%
  as_tibble(rownames = "sample_id") %>%
  full_join(metadata, by = "sample_id") %>%
  pivot_longer(cols = -any_of(colnames(metadata)), 
               names_to = "asv_names", values_to = "value") %>%
  mutate(across(c(exposure, final_disease_state), factor)) %>%
  filter(time %in% c('T3', 'T7') | (time == "T0" & tank == "HOMO")) %>%
  mutate(fragment_id = str_c(exposure, tank, genotype, final_disease_state)) %>%
  group_by(asv_names) %>%
  filter(n_distinct(sample_id[value > 7.04]) > 21) %>%
  ungroup()

#likely suspects upset prep
target_upset_data <- target_data %>%
  mutate(time = if_else(time == 'T0', exposure, time)) %>%
  group_by(time, asv_names) %>%
  summarise(n = sum(value > 7.04)) %>%
  pivot_wider(names_from = time, values_from = n, values_fill = 0L) %>%
  mutate(across(-asv_names, ~. > 0)) %>%
  left_join(taxonomy_tibble, by = "asv_names") %>%
  select(-c(Kingdom, Phylum)) %>%
  filter(D & T7 & T3)

taxonomy_tibble <- tax_table(microbiome_data) %>% 
  as.data.frame %>%
  as_tibble(rownames = "asv_names")


#### Model ####

ls_model <- model_data %>% 
  filter(time != 'T0') %>%
  inner_join(target_upset_data) %>%
  nest_by(asv_names) %>%
  # mutate(model = list(mixed(value ~time * (final_disease_state + exposure) + 
  #                        (0 + dummy(time, c('T3', 'T7')) | fragment_id), data = data, method = 'KR',
  #                        control = variancePartition:::vpcontrol)))
  mutate(model = list(mixed(value ~time * (final_disease_state + exposure) + 
                              (1 | fragment_id), data = data, method = 'KR',
                            control = variancePartition:::vpcontrol)))

#model summary for all interaction types for our subsetted data
ls_model_full <- ls_model %>%
  ungroup %>%
  rowwise(asv_names) %>%
  reframe(as_tibble(model$anova_table, rownames = 'param'),
          d_v_h = emmeans(model, ~final_disease_state) %>%
            as_tibble %>%
            select(emmean) %>%
            pull(1) %>%
            diff) %>%
  dplyr::rename(p = `Pr(>F)`) %>%
  group_by(param) %>%
  mutate(p = p.adjust(p, 'fdr')) %>%
  ungroup

#model w summary and showing significance for FDS, FDS:time, or both
ls_model_w_terms <- ls_model_full %>%
  filter(p < 0.05) %>%
  filter(param %in% c("final_disease_state", "time:final_disease_state")) %>%
  nest_by(asv_names) %>%
  rowwise() %>%
  mutate(terms = ifelse(nrow(data) < 2, NA, "both")) %>%
  unnest(data) %>%
  mutate(terms = ifelse(is.na(terms), param, terms)) %>%
  left_join(ls_model, by = join_by(asv_names))

#### Prep for Comp Upsets ####

#process data for performing complex upset - LIKELY SUSPECTS
subset_asv_comp_upset <- ls_model %>% #reduces from 305 to 249 bc 56 are significant for nothing, 133 of these are ~time
  ungroup() %>%
  rowwise(asv_names) %>%
  reframe(as_tibble(model$anova_table, rownames = 'param'),
          d_v_h =  emmeans(model, ~final_disease_state) %>%
            as_tibble %>%
            select(emmean) %>%
            pull(1) %>%
            diff) %>%
  dplyr::rename(p = `Pr(>F)`) %>%
  group_by(param) %>%
  mutate(p = p.adjust(p, 'fdr')) %>%
  ungroup %>%
  mutate(p = p < 0.05) %>%
  select(asv_names, param, p, d_v_h) %>%
  pivot_wider(names_from = 'param', values_from = p, names_prefix = 'p_', values_fill = FALSE) %>%
  left_join(taxonomy_tibble, by = join_by(asv_names)) %>%
  select(-c(Kingdom, Phylum)) %>%
  #select(-p_time) %>% #TODO decide whether to include time in upsets here
  filter(!if_all(starts_with('p_'), ~!.))
  #filter(p_final_disease_state | `p_final_disease_state:time`)


#### Complex Upsets ####
  
#Visualizing Important Intersections in Venn Diagram
  upset(target_upset_data, 
        c('D', 'H', 'T3', 'T7'), 
        annotations = list(
          # 2nd method - using ggplot
          'Order'=(
            ggplot(mapping=aes(fill=Order)) 
            + geom_bar(stat = 'count', position = 'fill') 
            + scale_y_continuous(labels = scales::percent_format())
          ) +
            ylab('Order') +
            theme(legend.position = 'top')
        ),
        name='asv_names', width_ratio=0.1, min_size = 1)

  
#Likely Suspects, 249 total; 116 w/out time
upset(subset_asv_comp_upset,
        colnames(select(subset_asv_comp_upset, starts_with('p_'))), 
        annotations = list(
          # 2nd method - using ggplot
          'Order'=(
            ggplot(mapping=aes(fill=Order)) 
            + geom_bar(stat = 'count', position = 'fill')
            + scale_y_continuous(labels = scales::percent_format())
          ) +
            ylab('Order') +
            theme(legend.position = 'top')
        ),
        name='asv_names', width_ratio=0.1, min_size = 0) +
  ggtitle("Likely Suspects")

upset(subset_asv_comp_upset,
                    colnames(select(subset_asv_comp_upset, starts_with('p_'))),
                    name='asv_names', width_ratio=0.1, min_size = 0) +
  ggtitle("Likely Suspects")
  
# All likely suspects - more abundant in Diseased
cu_disease <- upset(filter(subset_asv_comp_upset, d_v_h < 0) %>%
        select(-d_v_h),  #negative = more in disease final state, positive = more in healthy disease state
      colnames(select(subset_asv_comp_upset, starts_with('p_'))), 
      
      annotations = list(
        # 2nd method - using ggplot
        'Family'=(
          ggplot(mapping=aes(fill=Family)) 
          + geom_bar(stat = 'count', position = 'fill') 
          + scale_y_continuous(labels = scales::percent_format())
        ) +
          ylab('Order') +
          theme(legend.position = 'top')
      ),
      
      name='asv_names', width_ratio=0.1, min_size = 1) + #, max_size = 70
  ggtitle("More in Diseased")


# All likely suspects - more abundant in Healthy
cu_healthy <- upset(filter(subset_asv_comp_upset, d_v_h > 0) %>%
        select(-d_v_h),  #negative = more in disease final state, positive = more in healthy disease state
      colnames(select(subset_asv_comp_upset, starts_with('p_'))), 
      
      annotations = list(
        # 2nd method - using ggplot
        'Family'=(
          ggplot(mapping=aes(fill=Family)) 
          + geom_bar(stat = 'count', position = 'fill') 
          + scale_y_continuous(labels = scales::percent_format())
        ) +
          ylab('Family') +
          theme(legend.position = 'top')
      ),
      
      name='asv_names', width_ratio=0.1, min_size = 1) + #, max_size = 30
  ggtitle("More in Healthy")


layout <- '
AA
BB
'
wrap_plots(A = cu_disease, B = cu_healthy, design = layout)


#### Graphing ####

#venn diagram
target_data %>%
  mutate(time = if_else(time == 'T0', exposure, time)) %>%
  # count(time, asv_names) %>%
  group_by(time, asv_names) %>%
  summarise(n = sum(value > 7.04)) %>%
  pivot_wider(names_from = time, values_from = n, values_fill = 0L) %>%
  mutate(across(-asv_names, ~. > 0)) %>%
  ggvenn(c('D', 'H', 'T3', 'T7'))


##EMMEANS for VLS

#VLS FDS only
vls_plot_1 <- ls_model_w_terms %>%
  filter(terms == "final_disease_state", d_v_h < 0) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(test = list(make_emmean_model(model, as.formula("~final_disease_state"), 0.05))) %>%
  group_by(asv_names) %>%
  reframe(test =  as_tibble(test, .name_repair = "unique")) %>%
  unnest(cols = c(test)) %>%
  unnest(cols = `...1`) %>%
  mutate(.group = str_trim(.group)) %>%
  left_join(taxonomy_tibble, by = join_by(asv_names)) %>%
  select(-c("Kingdom", "Phylum")) %>%
  left_join((select(ls_model_w_terms, c(asv_names, d_v_h))), by = join_by(asv_names)) %>%
  ggplot(aes(x = fct_reorder(paste(Family, " ", Genus, " (", asv_names, ")", sep = ""), d_v_h, .desc = TRUE), y = emmean, ymin = emmean - SE, ymax = emmean + SE,
             col = final_disease_state)) +
  geom_pointrange(position = position_dodge(0.5)) +
  geom_text(aes(y = (emmean + SE), label = .group),
            position = position_dodge(0.5), vjust = -1) +
  scale_color_manual(values = c("firebrick1", "dodgerblue3")) +
  coord_flip() +
  xlab("ASV") +
  labs(title = "Final Disease State Only")

#VLS interaction only
vls_plot_2 <- ls_model_w_terms %>%
  filter(terms == "time:final_disease_state", d_v_h < 0) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(test = list(make_emmean_model(model, as.formula("~time:final_disease_state"), 0.05))) %>%
  group_by(asv_names) %>%
  reframe(test =  as_tibble(test, .name_repair = "unique")) %>%
  unnest(cols = c(test)) %>%
  unnest(cols = `...1`) %>%
  mutate(.group = str_trim(.group)) %>%
  left_join(taxonomy_tibble, by = join_by(asv_names)) %>%
  select(-c("Kingdom", "Phylum")) %>%
  left_join((select(ls_model_w_terms, c(asv_names, d_v_h))), by = join_by(asv_names)) %>%
  ggplot(aes(x = fct_reorder(paste(Family, " ", Genus, " (", asv_names, ")", sep = ""), d_v_h, .desc = TRUE), y = emmean, ymin = emmean - SE, ymax = emmean + SE,
             col = time:final_disease_state, pch = time)) +
  geom_pointrange(position = position_dodge(0.5)) +
  geom_text(aes(y = (emmean + SE), label = .group),
            position = position_dodge(0.5), vjust = -1) +
  scale_color_manual(values = c("hotpink1", "deepskyblue", "firebrick1", "dodgerblue3")) +
  coord_flip() +
  xlab("ASV") +
  labs(title = "Final Disease State:Time Interaction Only")


#Both FDS and FDS:time
vls_plot_3 <- ls_model_w_terms %>%
  filter(terms == "both", param == "time:final_disease_state", d_v_h < 0) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(test = list(make_emmean_model(model, as.formula("~time:final_disease_state"), 0.05))) %>%
  group_by(asv_names) %>%
  reframe(test =  as_tibble(test, .name_repair = "unique")) %>%
  unnest(cols = c(test)) %>%
  unnest(cols = `...1`) %>%
  mutate(.group = str_trim(.group)) %>%
  left_join(taxonomy_tibble, by = join_by(asv_names)) %>%
  select(-c("Kingdom", "Phylum")) %>%
  left_join((select(ls_model_w_terms, c(asv_names, d_v_h)) %>% distinct(asv_names, d_v_h)), by = join_by(asv_names)) %>%
  ggplot(aes(x = fct_reorder(paste(Family, " ", Genus, " (", asv_names, ")", sep = ""), d_v_h, .desc = TRUE), y = emmean, ymin = emmean - SE, ymax = emmean + SE,
             col = time:final_disease_state, pch = time)) +
  geom_pointrange(position = position_dodge(0.5)) +
  geom_text(aes(y = (emmean + SE), label = .group),
            position = position_dodge(0.5), vjust = -1) +
  scale_color_manual(values = c("hotpink1", "deepskyblue", "firebrick1", "dodgerblue3")) +
  coord_flip() +
  xlab("ASV") +
  labs(title = "Both FDS and FDS:Time")

vls_plot_1 / vls_plot_2 / vls_plot_3 + plot_layout(heights = c(4, 7, 16)) + 
  plot_annotation(title = "Very Likely Suspects") & theme_linedraw() & 
  theme(text = element_text('sans'), plot.title = element_text(size = 16))




#### Relative Abundances ####

"green"
#TODO set the taxonomic level for the plots here:
tax_level <- "Family"

ls_asv_list <- ls_model_w_terms %>% .$asv_names %>% unique()
vls_asv_list <- ls_model_w_terms %>% filter(d_v_h < 0) %>% .$asv_names %>% unique()

#ASVs in whole dataset
relabun_tot <- taxonomy_tibble %>% nest_by(!!! rlang::syms(tax_level)) %>% mutate(total_asvs = nrow(data)) %>% select(-data)

#ASVs in LS list
relabun_ls <-taxonomy_tibble %>% filter(asv_names %in% ls_asv_list) %>% nest_by(!!! rlang::syms(tax_level)) %>% 
  mutate(ls_asvs = nrow(data)) %>% select(-data)

#ASVs in VLS list
relabun_vls <-taxonomy_tibble %>% filter(asv_names %in% vls_asv_list) %>% nest_by(!!! rlang::syms(tax_level)) %>% 
  mutate(vls_asvs = nrow(data)) %>% select(-data)

relative_abundance_taxon <- relabun_ls %>% left_join(relabun_tot) %>% full_join(relabun_vls) %>%
  reframe(ls = 100*ls_asvs/total_asvs, vls = 100*vls_asvs/total_asvs, total = total_asvs) %>%
  mutate(vls = ifelse(is.na(vls), 0, vls)) %>%
  pivot_longer(c(ls, vls), names_to = "list", values_to = "percent")

relabun_ranks <- relative_abundance_taxon %>% filter(list == "ls") %>% 
  mutate(graph_order = (rank(percent))/(nrow(.) + 1)) %>% select(-c(total, list, percent))

relabun_graph <- relative_abundance_taxon %>% left_join(relabun_ranks) %>%
  ggplot(aes(x = fct_reorder((!!! rlang::syms(tax_level)), graph_order), y = percent, fill = list, 
        label = ifelse(percent == 0, "", paste(round(percent, 2), "% (of ", total, ")", sep = "")))) +
    geom_col(position = "dodge") +
    geom_text(size = 3, position = position_dodge(width = 0.9), hjust = -0.05) +
    coord_flip() +
    xlab(tax_level) +
    ylab(paste("Percent of Total ASVs in Each ", tax_level, sep = "")) +
    scale_fill_manual(values = c("orange", "firebrick1"), name = "Suspect Type", labels = c("Likely", "Very Likely")) +
    ylim(0, 5.8) +
    theme_linedraw()
  
#complimentary comp upset
taxon_comp_upset <- upset(subset_asv_comp_upset,
      colnames(select(subset_asv_comp_upset, starts_with('p_') & contains("final"))), 
      annotations = list(
        # 2nd method - using ggplot
        'Genus'=(
          ggplot(mapping=aes(fill=(!!! rlang::syms(tax_level)))) 
          + geom_bar(stat = 'count', position = 'fill')
          + scale_y_continuous(labels = scales::percent_format())
        ) +
          ylab(tax_level) +
          theme(legend.position = 'top')
      ),
      name='asv_names', width_ratio=0.1, min_size = 1, min_degree = 1) +
  ggtitle("Likely Suspects")

layout <- '
AABB#
AABB#
'
wrap_plots(A = relabun_graph, B = taxon_comp_upset, design = layout)

  
  

#### Most Abundant ASVs ####

#what's most abundant in the homogenates
bait_high_abun <- model_data %>% 
  filter(time == "T0") %>%
  group_by(asv_names, final_disease_state) %>%
  summarize(tot_abun = sum(value)) %>%
  ungroup() %>%
  group_by(final_disease_state) %>%
  arrange(desc(tot_abun)) %>%
  dplyr::slice(1:10) %>%
  left_join(taxonomy_tibble, by = join_by(asv_names)) %>%
  mutate(asv_names = parse_number(asv_names), tot_abun = round(tot_abun, 2), 
         Species = paste(Species, " (", asv_names, ")", sep = "")) %>%
  select(Order, Family, Genus, Species, final_disease_state, tot_abun) %>%
  ungroup() %>%
  pivot_wider(names_from = final_disease_state, values_from = tot_abun)

baits_most_abun <- formattable(bait_high_abun,
            align = c("l", "l", "l", "l", "c", "c", "c"), list(
              Order = formatter("span", style = ~ style(color = "gray",font.weight = "bold")),
              Family = formatter("span", style = ~ style(color = "gray",font.weight = "bold")),
              Genus = formatter("span", style = ~ style(color = "gray")),
              Species = formatter("span", style = ~ style(color = "gray")),
              D = color_tile("#FFC5C5", "#FF6767"),
              H = color_tile("#BFFFE9", "#43FFC0")
            ))

export_formattable(baits_most_abun,"most_abun_baits.png")

"red"
## Family level - MUST AGGREGATE DATA BY FAMILY TO GET THESE RESULTS
"red"

family_chart <- model_data %>% 
  rename(asv_names = "Family") %>%
  mutate(category = factor(paste(time, final_disease_state, sep = "_"),
                          levels = c("T0_H", "T3_H", "T7_H", "T0_D", "T3_D", "T7_D")), .before = time) %>%
  group_by(category, Family) %>%
  summarize(abundance = sum(value)) %>%
  ungroup() %>%
  arrange(desc(abundance)) %>%
  group_by(category) %>%
  dplyr::slice(1:20) %>%
  mutate(fam_rank = rank(-abundance)) %>%
  select(-abundance) %>%
  pivot_wider(names_from = category, values_from = fam_rank)

repeated_fams <- family_chart %>% mutate_at(vars(-c(Family)), ~ifelse(is.na(.), 0, 1)) %>% group_by(Family) %>%
  summarize(count = sum(c(T0_H, T3_H, T7_H, T0_D, T3_D, T7_D))) %>% filter(count > 1) %>% .$Family

family_chart %>% 
  filter(Family %in% repeated_fams) %>%
  formattable(align = c("l", "c", "c", "c", "c", "c", "c"), list(
    Family = formatter("span", style = ~ style(color = "gray",font.weight = "bold")),
    T0_D = color_tile("#FF3A3A", "#FFE1E1"),
    T3_D = color_tile("#FF3A3A", "#FFE1E1"),
    T7_D = color_tile("#FF3A3A", "#FFE1E1"),
    T0_H = color_tile("#00B276", "#DAFFF2"),
    T3_H = color_tile("#00B276", "#DAFFF2"),
    T7_H = color_tile("#00B276", "#DAFFF2")
  )) %>%
  export_formattable("most_abun_timepoints.png")

#### Logfold Changes ####

#same model but run on the log of the data
log_ls_model <- model_data %>% 
  filter(time != 'T0') %>%
  inner_join(target_upset_data) %>%
  nest_by(asv_names) %>%
  # mutate(model = list(mixed(value ~time * (final_disease_state + exposure) + 
  #                        (0 + dummy(time, c('T3', 'T7')) | fragment_id), data = data, method = 'KR',
  #                        control = variancePartition:::vpcontrol)))
  mutate(model = list(mixed(log(value) ~time * (final_disease_state + exposure) + 
                              (1 | fragment_id), data = data, method = 'KR',
                            control = variancePartition:::vpcontrol)))

logfold_data <- log_ls_model %>%
  ungroup() %>%
  inner_join((ls_model_w_terms %>% select(asv_names, terms)), by = join_by(asv_names), multiple = "all") %>%
  rowwise() %>%
  mutate(test = list(emmeans(model, ~final_disease_state | time) %>%
                       contrast('pairwise', adjust = 'fdr') %>% as_tibble())) %>%
  select(-c(data, model)) %>%
  unnest(test) %>%
  distinct() %>% 
  group_by(asv_names) %>% 
  mutate(diff = estimate[time == "T7"] - estimate[time == "T3"]) %>%
  ungroup() %>%
  mutate(terms = factor(terms, levels = c("final_disease_state", "time:final_disease_state", "both"),
                           labels = c("FDS", "FDS:Time", "Both"))) %>%
  left_join(taxonomy_tibble, by = join_by(asv_names))


logfold_t7_more <- ggplot(logfold_data %>% filter(diff > 0), aes(x = fct_reorder(paste(Family, " ", Genus, " (", parse_number(asv_names), ")", sep = ""), diff), y = estimate, 
             ymin = estimate - SE, ymax = estimate + SE, colour = time, pch = time)) +
  geom_hline(yintercept = 0) +
  geom_pointrange(position = position_dodge(0.5)) +
  scale_color_manual(values = c("darkorange1", "firebrick1")) +
  coord_flip() +
  xlab("") +
  labs(title = "Logfold Changes - More in T7") +
  facet_grid(rows = vars(terms), scales = "free", space = "free") +
  theme_bw()

logfold_t3_more <- ggplot(logfold_data %>% filter(diff < 0), aes(x = fct_reorder(paste(Family, " ", Genus, " (", parse_number(asv_names), ")", sep = ""), -diff), y = estimate, 
                                              ymin = estimate - SE, ymax = estimate + SE, colour = time, pch = time)) +
  geom_hline(yintercept = 0) +
  geom_pointrange(position = position_dodge(0.5)) +
  scale_color_manual(values = c("darkorange1", "firebrick1")) +
  coord_flip() +
  xlab("ASV") +
  labs(title = "Logfold Changes - More in T3") +
  facet_grid(rows = vars(terms), scales = "free", space = "free") +
  theme_bw() +
  theme(legend.position = "none")

logfold_t3_more | logfold_t7_more


#TODO add pcoa code to this doc.  upload all to github w comments.  trees w jason


