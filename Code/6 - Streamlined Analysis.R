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

my_color_tile <- function(health_state) {
  col_choices_red <- colorRampPalette(c("#FF3A3A", "#FFE1E1"))(20)
  col_choices_green <- colorRampPalette(c("#00B276", "#DAFFF2"))(20)
  
  return_col <- function(y) 
    if(health_state == "Diseased"){
      map_chr(y,function(x) case_when(x > 20  ~ "#FFF1F1",
                                      x > 0  ~ col_choices_red[x],
      ))
    } else if(health_state == "Healthy"){
      map_chr(y,function(x) case_when(x > 20  ~ "#F3FFFB",
                                      x > 0  ~ col_choices_green[x],
      ))  
    }
  
  formatter("span", 
            style = function(y) style(
              display = "block",
              padding = "0 4px",
              "border-radius" = "4px",
              "color" = ifelse( return_col(y) %in% c("#006837","#1A9850","#66BD63"),
                                csscolor("white"), csscolor("black")),
              "background-color" = return_col(y)
            )
  )
}

my_color_tile_baits <- function(health_state) {
  col_choices_red <- c(rep("white", 62), colorRampPalette(c("#FF9292","#FA0A0A"))(24))
  col_choices_green <- c(rep("white", 67), colorRampPalette(c("#8EFFD9", "#00A76F"))(28))
  
  return_col <- function(y) 
    if(health_state == "Diseased"){
      map_chr(y,function(x) case_when(x > 62  ~ col_choices_red[x],
                                      x < 62  ~ "#FFF1F1"
      ))
    } else if(health_state == "Healthy"){
      map_chr(y,function(x) case_when(x > 66  ~ col_choices_green[x],
                                      x < 66  ~ "#DDFFF3",
      ))  
    }
  
  formatter("span", 
            style = function(y) style(
              display = "block",
              padding = "0 4px",
              "border-radius" = "4px",
              "color" = ifelse( return_col(y) %in% c("#006837","#1A9850","#66BD63"),
                                csscolor("white"), csscolor("black")),
              "background-color" = return_col(y)
            )
  )
}
#
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

#write_csv((model_data %>% inner_join(target_upset_data)), "../intermediate_files/pre_model_data.csv")

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
#compare disease v healthy samples
comp_dvh <- ls_model %>%
  rowwise() %>%
  mutate(comps = list(emmeans(model, ~final_disease_state) %>%
                        contrast('pairwise', adjust = 'fdr') %>% as_tibble())) %>% 
  unnest(comps) %>%
  select(-c(SE, contrast, df, t.ratio)) %>%
  rename("d_v_h" = estimate) %>%
  mutate(adj_dvh = ifelse(p.value < 0.05, d_v_h, 0)) %>%
  mutate(up_down = ifelse(adj_dvh > 0, "up", ifelse(adj_dvh < 0, "down", "non-significant")))

vl_suspects <- comp_dvh %>%
  filter(adj_dvh < 0) %>%
  .$asv_names


#model summary for all interaction types for our subsetted data
ls_model_full <- ls_model %>%
  ungroup %>%
  rowwise(asv_names) %>%
  reframe(as_tibble(model$anova_table, rownames = 'param')) %>%
  dplyr::rename(p = `Pr(>F)`) %>%
  group_by(param) %>%
  mutate(p = p.adjust(p, 'fdr')) %>%
  ungroup %>%
  left_join(comp_dvh %>% select(-c(p.value, data, model)), by = join_by(asv_names))

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
  reframe(as_tibble(model$anova_table, rownames = 'param')) %>%
  dplyr::rename(p = `Pr(>F)`) %>%
  group_by(param) %>%
  mutate(p = p.adjust(p, 'fdr')) %>%
  ungroup %>%
  mutate(p = p < 0.05) %>%
  left_join(comp_dvh %>% select(-c(p.value, data, model)), by = join_by(asv_names)) %>%
  select(asv_names, param, p, d_v_h, up_down) %>%
  pivot_wider(names_from = 'param', values_from = p, names_prefix = 'p_', values_fill = FALSE) %>%
  left_join(taxonomy_tibble, by = join_by(asv_names)) %>%
  select(-c(Kingdom, Phylum)) %>%
  #select(-p_time) %>% #TODO decide whether to include time in upsets here
  filter(!if_all(starts_with('p_'), ~!.)) %>%
  filter(!(p_time & !p_final_disease_state & !p_exposure & !`p_time:final_disease_state` & !`p_time:exposure`))
  #filter(p_final_disease_state | `p_time:final_disease_state`)

likely_asvs <- subset_asv_comp_upset %>%
  filter(p_final_disease_state | `p_time:final_disease_state`) %>%
  .$asv_names %>% unique()

write_rds(list(likely_asvs, vl_suspects), "../intermediate_files/important_asvs.rds")

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
      
      base_annotations=list(
        'Intersection size'=intersection_size(
          mapping=aes(fill=up_down)
        ) + scale_fill_manual(values=c(
          'up'='#E41A1C', 
          'down'='#4DAF4A', 'non-significant'='#FF7F00'
        ))
      ),
      queries=list(upset_query(set='p_time:final_disease_state', color="#8400CA", fill = "#4C0075"), 
                   upset_query(set='p_final_disease_state', color="#8400CA", fill = "#4C0075")),
               
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


upset(subset_asv_comp_upset %>% filter(up_down == "up"),
      colnames(select(subset_asv_comp_upset, starts_with('p_'))), 
      
      base_annotations=list(
        'Intersection size'=intersection_size(
          mapping=aes(fill=up_down)
        ) + scale_fill_manual(values=c(
          'up'='#E41A1C', 
          'down'='#4DAF4A', 'non-significant'='#FF7F00'
        ))
      ),
      queries=list(upset_query(set='p_time:final_disease_state', color="#8400CA", fill = "#4C0075"), 
                   upset_query(set='p_final_disease_state', color="#8400CA", fill = "#4C0075")),
      
      annotations = list(
        # 2nd method - using ggplot
        'Order'=(
          ggplot(mapping=aes(fill=Family)) 
          + geom_bar(stat = 'count', position = 'fill')
          + scale_y_continuous(labels = scales::percent_format())
        ) +
          ylab('Order') +
          theme(legend.position = 'top')
      ),
      name='asv_names', width_ratio=0.1, min_size = 0) +
  ggtitle("Very Likely Suspects")
  


likely_families <-  subset_asv_comp_upset %>%
  filter(p_final_disease_state | `p_time:final_disease_state`) %>%
  .$Family %>% unique()

very_likely_families <-  subset_asv_comp_upset %>%
  filter(p_final_disease_state | `p_time:final_disease_state`, d_v_h < 0) %>%
  .$Family %>% unique()

write_rds(list(likely_families, very_likely_families), "../intermediate_files/families_of_interest.rds")
  
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
  ggplot(aes(x = fct_reorder(paste(Family, " ", Genus, " (", asv_names, ")", sep = ""), d_v_h, .desc = TRUE), 
             y = emmean, ymin = emmean - SE, ymax = emmean + SE,
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
  ggplot(aes(x = fct_reorder(paste(Family, " ", Genus, " (", asv_names, ")", sep = ""), d_v_h, .desc = TRUE), 
             y = emmean, ymin = emmean - SE, ymax = emmean + SE,
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
  ggplot(aes(x = fct_reorder(paste(Family, " ", Genus, " (", asv_names, ")", sep = ""), d_v_h, .desc = TRUE), 
             y = emmean, ymin = emmean - SE, ymax = emmean + SE,
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
top_ten <- model_data %>% 
  filter(time == "T0") %>%
  group_by(asv_names, final_disease_state) %>%
  summarize(tot_abun = sum(value)) %>%
  ungroup() %>%
  group_by(final_disease_state) %>%
  arrange(desc(tot_abun)) %>%
  dplyr::slice(1:10) %>%
  .$asv_names %>% unique()

model_data %>% 
  filter(time == "T0") %>%
  group_by(asv_names, final_disease_state) %>%
  summarize(tot_abun = sum(value)) %>%
  ungroup() %>%
  group_by(final_disease_state) %>%
  arrange(desc(tot_abun)) %>%
  filter(asv_names %in% top_ten) %>%
  left_join(taxonomy_tibble, by = join_by(asv_names)) %>%
  mutate(asv_names = parse_number(asv_names), tot_abun = round(tot_abun, 2), 
         Genus = paste(Genus, " (", asv_names, ")", sep = "")) %>%
  select(Order, Family, Genus, final_disease_state, tot_abun) %>%
  ungroup() %>%
  pivot_wider(names_from = final_disease_state, values_from = tot_abun) %>%
  rowwise() %>%
  mutate(asv_order = as.numeric(gsub("\\(([^()]+)\\)", "\\1", str_extract_all(Genus, "\\(([^()]+)\\)")[[1]]))) %>%
  arrange(asv_order) %>%
  select(-asv_order) %>%
  formattable(align = c("l", "l", "l", "l", "c", "c", "c"), list(
              Order = formatter("span", style = ~ style(color = "gray",font.weight = "bold")),
              Family = formatter("span", style = ~ style(color = "gray",font.weight = "bold")),
              Genus = formatter("span", style = ~ style(color = "gray")),
              Species = formatter("span", style = ~ style(color = "gray")),
              area(col = 5) ~ my_color_tile_baits("Diseased"),
              area(col = 4) ~ my_color_tile_baits("Healthy")
            )) %>%
  export_formattable("../Figures/most_abun_baits.png")

"red"
## Family level - MUST AGGREGATE DATA BY FAMILY TO GET THESE RESULTS
"red"

family_chart <- model_data %>% 
  dplyr::rename("Family" = asv_names) %>%
  mutate(category = factor(paste(time, final_disease_state, sep = "_"),
                          levels = c("T0_H", "T3_H", "T7_H", "T0_D", "T3_D", "T7_D")), .before = time) %>%
  group_by(category, Family) %>%
  summarize(abundance = sum(value)) %>%
  ungroup() %>%
  arrange(desc(abundance)) %>%
  group_by(category) %>%
  mutate(fam_rank = rank(-abundance, ties.method= "max"))

kept_fams <- family_chart %>%
  dplyr::slice(1:20) %>%
  .$Family %>%
  unique()

family_chart_1 <- family_chart %>%
  filter(Family %in% kept_fams) %>%
  select(-abundance) %>%
  pivot_wider(names_from = category, values_from = fam_rank) %>%
  select(Family, T0_H, T3_H, T7_H, T0_D, T3_D, T7_D)


repeated_fams <- model_data %>% 
  rename("Family" = asv_names) %>%
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
  pivot_wider(names_from = category, values_from = fam_rank) %>% 
  mutate_at(vars(-c(Family)), ~ifelse(is.na(.), 0, 1)) %>% group_by(Family) %>%
  summarize(count = sum(c(T0_H, T3_H, T7_H, T0_D, T3_D, T7_D))) %>% 
  filter(count > 1) %>% .$Family

family_chart_1 %>% 
  filter(Family %in% repeated_fams) %>%
  formattable(align = c("l", "c", "c", "c", "c", "c", "c"), list(
    Family = formatter("span", style = ~ style(color = "gray",font.weight = "bold")),
    area(col = 2:4) ~ my_color_tile("Healthy"),
    area(col = 5:7) ~ my_color_tile("Diseased")
  )) %>%
  export_formattable("../Figures/most_abun_timepoints.png")
  
#### Logfold Changes ####

#same model but run on the log of the data

logfold_data <- ls_model %>%
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


logfold_t7_more <- ggplot(logfold_data %>% filter(diff > 0), aes(x = fct_reorder(paste(Family, " ", 
                              Genus, " (", parse_number(asv_names), ")", sep = ""), diff), y = estimate, 
             ymin = estimate - SE, ymax = estimate + SE, colour = time, pch = time)) +
  geom_hline(yintercept = 0) +
  geom_pointrange(position = position_dodge(0.5)) +
  scale_color_manual(values = c("darkorange1", "firebrick1")) +
  coord_flip() +
  xlab("") +
  labs(title = "Logfold Changes - More in T7") +
  facet_grid(rows = vars(terms), scales = "free", space = "free") +
  theme_bw()

logfold_t3_more <- ggplot(logfold_data %>% filter(diff < 0), aes(x = fct_reorder(paste(Family, " ", 
                                    Genus, " (", parse_number(asv_names), ")", sep = ""), -diff), y = estimate, 
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


#TODO add pcoa code to this doc



