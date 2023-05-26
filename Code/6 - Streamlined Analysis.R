#Code that streamlines the process from Script #5: Disease-Associated Bacteria

setwd("~/Desktop/Screenshots/Career/Vollmer Lab/GitHub/16S_Florida_Tank_Analysis/Code")

#### TODO LIST ####

#global patterns package
#lm for trees
#compare my model and jason's


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
library(microshades)
library(cowplot)
library(microViz)
library(tidyverse)

select <- dplyr::select

#### Functions ####

make_emmean_model <- function(model, form, alpha){
  emmeans(model, specs = form, type = 'response') %>%
    cld(Letters = LETTERS, reversed = TRUE, alpha = alpha) 
}

make_model_plot <- function(predOut, rawData, var_inclusion){
  
  all_vars <- str_extract_all(var_inclusion, 'time|final_disease_state|exposure') %>%
    unlist
  
  if(any(str_detect(all_vars, 'time'))){
    x_var <- str_subset(all_vars, 'time.*')
  } else {
    x_var <- all_vars[1]
  }
  
  if(length(all_vars) > 1){
    colour_var <- str_subset(all_vars, x_var, negate = TRUE)
  } else {
    colour_var <- NULL
  }
  
  y_var <- str_subset(colnames(predOut), 'emmean|response')
  y_data <- str_subset(colnames(rawData), 'value')
  
  if(is.null(colour_var)){
    as_tibble(predOut) %>%
      mutate(.group = str_trim(.group)) %>%
      ggplot(aes(x = !!sym(x_var))) +
      
      stat_halfeye(data = rawData, aes(y = !!sym(y_data)),
                   adjust = 0.5, width = 0.6, .width = 0, 
                   alpha = 0.5, show.legend = FALSE,
                   fatten_point = 0, justification = -0.25, 
                   position = position_dodge(0.5), size = 0) +
      
      geom_half_point(data = rawData, aes(y = !!sym(y_data)),
                      side = 'l', range_scale = 0.1, alpha = 1,
                      position = position_dodge(width = 0.5),
                      show.legend = FALSE,
                      transformation = position_jitter(height = 0, width = 0.05)) +
      
      geom_pointrange(aes(y = !!sym(y_var),
                          ymin = lower.CL, ymax = upper.CL),
                      position = position_dodge(0.5)) +
      geom_text(aes(y = (upper.CL), label = .group),
                position = position_dodge(0.5), vjust = -0.1, 
                show.legend = FALSE) +
      labs(x = case_when(x_var %in% c('time', 'timepoint') ~ 'Time (d)',
                         x_var == 'exposure' ~ 'Exposure',
                         x_var == 'final_disease_state' ~ 'Disease State'))
    
  } else {
    as_tibble(predOut) %>%
      mutate(.group = str_trim(.group)) %>%
      ggplot(aes(x = !!sym(x_var), colour = !!sym(colour_var))) +
      
      stat_halfeye(data = rawData, aes(y = !!sym(y_data), fill = !!sym(colour_var)),
                   adjust = 0.5, width = 0.6, .width = 0, 
                   alpha = 0.5, show.legend = FALSE,
                   fatten_point = 0, justification = -0.25, 
                   position = position_dodge(0.5), size = 0) +
      
      geom_half_point(data = rawData, aes(y = !!sym(y_data), colour = !!sym(colour_var)),
                      side = 'l', range_scale = 0.1, alpha = 1,
                      position = position_dodge(width = 0.5),
                      show.legend = FALSE,
                      transformation = position_jitter(height = 0, width = 0.05)) +
      
      geom_pointrange(aes(y = !!sym(y_var),
                          ymin = lower.CL, ymax = upper.CL),
                      position = position_dodge(0.5)) +
      geom_text(aes(y = (upper.CL), label = .group),
                position = position_dodge(0.5), vjust = -0.1,
                show.legend = FALSE) +
      labs(x = case_when(x_var %in% c('time', 'timepoint') ~ 'Time (d)',
                         x_var == 'exposure' ~ 'Exposure',
                         x_var == 'final_disease_state' ~ 'Disease State'),
           colour = case_when(colour_var %in% c('time', 'timepoint') ~ 'Time (d)',
                              colour_var == 'exposure' ~ 'Exposure',
                              colour_var == 'final_disease_state' ~ 'Disease State'))
  }
  
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

plot_agg_pcoa <- function(cpm_counts){
  filtered_pcoa <- t(cpm_counts) %>%
    vegdist(method = 'euclidean') %>%
    magrittr::divide_by(1000) %>%
    pcoa()
  
  percent_variance <- filtered_pcoa$values$Eigenvalues / sum(filtered_pcoa$values$Eigenvalue)
  
  filtered_pcoa$vectors %>%
    as_tibble(rownames = 'sample_id') %>%
    dplyr::select(sample_id, Axis.1, Axis.2) %>%
    inner_join(metadata,
               by = 'sample_id') %>%
    ggplot(aes(x = Axis.1, y = Axis.2, colour = exposure, 
               shape = time, group = fragment_id)) +
    geom_point() +
    geom_path() +
    labs(x = str_c('PCoA 1 (', scales::percent(percent_variance[1]), ')'),
         y = str_c('PCoA 2 (', scales::percent(percent_variance[2]), ')'),
         title = agg_title) +
    theme_classic()
  
}

find_unique_significant_terms_rmANOVA <- function(model, alpha){
  significant_terms <- model$anova_table %>%
    as_tibble(rownames = 'param') %>%
    janitor::clean_names() %>%
    filter(pr_f < alpha) %>%
    pull(param)
  
  unique_values <- outer(significant_terms, significant_terms, str_count) %>%
    colSums() %>%
    equals(1)
  
  str_replace(significant_terms[unique_values], ':', '*')
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

write_csv((model_data %>% inner_join(target_upset_data)), "../intermediate_files/pre_model_data.csv")

#target microbiome data - zeroes considered
target_data <- cpm(otu_tmm, log = TRUE, prior.count = 2) %>%
  t %>%
  as_tibble(rownames = "sample_id") %>%
  full_join(metadata, by = "sample_id") %>%
  pivot_longer(cols = -any_of(colnames(metadata)), 
               names_to = "asv_names", values_to = "value") %>%
  mutate(across(c(exposure, final_disease_state), factor)) %>%
  filter(time %in% c('T3', 'T7') | (time == "T0" & tank == "HOMO")) %>%
  mutate(fragment_id = str_c(exposure, tank, genotype, final_disease_state)) %>%
  group_by(asv_names)
  #filter(n_distinct(sample_id[value > 7.04]) > 19)

#likely suspects upset prep
target_upset_data <- target_data %>%
  mutate(time = if_else(time == 'T0', exposure, time)) %>%
  group_by(time, asv_names) %>%
  summarise(n = sum(value > 7.12)) %>%
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

fds_dvh_summary <- comp_dvh %>%
  select(asv_names, d_v_h, up_down)

#by exposure
exposure_dvh <- ls_model %>%
  rowwise() %>%
  mutate(comps = list(emmeans(model, ~exposure) %>%
                        contrast('pairwise', adjust = 'fdr') %>% as_tibble())) %>% 
  unnest(comps) 

exposure_dvh_summary <- exposure_dvh %>%
  select(-c(SE, contrast, df, t.ratio, data, model)) %>%
  rename("exp_d_v_h" = estimate) %>%
  #mutate(adj_dvh = ifelse(p.value < 0.05, exp_d_v_h, 0)) %>%
  mutate(exp_up_down = ifelse(p.value < 0.05, ifelse(exp_d_v_h > 0, "up", "down"), "non-significant")) %>%
  select(-c(p.value))


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

likely_suspects <- ls_model_w_terms$asv_names %>% unique()

very_likely_suspects <- ls_model_w_terms %>% filter(adj_dvh > 0) %>% .$asv_names %>% unique()


#### Prep for Comp Upsets ####

#process data for performing complex upset - LIKELY SUSPECTS
subset_asv_comp_upset <- ls_model %>% #reduces from 382 to 304 bc 78 are significant for nothing, 132 of these are ~time
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
  filter(!if_all(starts_with('p_'), ~!.)) %>%
  filter(!(p_time & !p_final_disease_state & !p_exposure & !`p_time:final_disease_state` & !`p_time:exposure`))


write_rds(list(likely_suspects, very_likely_suspects), "../intermediate_files/important_asvs.rds")

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
  ggtitle("Full Subset")


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
  

likely_families <-  likely_suspects %>% as_tibble() %>% rename("asv_names" = value) %>% 
  left_join(taxonomy_tibble) %>% .$Family %>% unique()


very_likely_families <- very_likely_suspects %>% as_tibble() %>% rename("asv_names" = value) %>% 
  left_join(taxonomy_tibble) %>% .$Family %>% unique()

write_rds(list(likely_families, very_likely_families), "../intermediate_files/families_of_interest.rds")
  
#### Tables RE: Differences in FDS or Exposure####

table_prep1 <- ls_model %>% #reduces from 305 to 249 bc 56 are significant for nothing, 133 of these are ~time
  ungroup() %>%
  rowwise(asv_names) %>%
  reframe(as_tibble(model$anova_table, rownames = 'param')) %>%
  dplyr::rename(p = `Pr(>F)`) %>%
  group_by(param) %>%
  mutate(p = p.adjust(p, 'fdr')) %>%
  ungroup %>%
  #mutate(p = p < 0.05) %>%
  mutate(p = ifelse(p < 0.05, 1, 0)) %>%
  left_join(comp_dvh %>% select(-c(p.value, data, model)), by = join_by(asv_names)) %>%
  select(asv_names, param, p, d_v_h, up_down) %>%
  pivot_wider(names_from = 'param', values_from = p, names_prefix = 'p_', values_fill = FALSE) %>%
  left_join(taxonomy_tibble, by = join_by(asv_names)) %>%
  select(-c(Kingdom, Phylum)) %>%
  filter(!if_all(starts_with('p_'), ~!.)) 

sig_for_anything <- table_prep %>% .$asv_names

remove_time_only <- table_prep %>%
  filter(!(p_time & !p_final_disease_state & !p_exposure & !`p_time:final_disease_state` & !`p_time:exposure`)) %>% 
  .$asv_names


family_table <- exposure_dvh_summary %>% 
  full_join(fds_dvh_summary) %>% 
  left_join(taxonomy_tibble) %>%
  select(-c(Kingdom, Phylum)) %>%
  mutate(any_sig_int = ifelse(asv_names %in% sig_for_anything, 1, 0), 
         non_time_sig = ifelse(asv_names %in% remove_time_only, 1, 0)) %>%
  select(-c(Class, Order, Species, Genus)) %>%
  group_by(Family) %>%
  summarize(total = n(), 
            total_sig = sum(any_sig_int), 
            total_non_time = sum(non_time_sig),
            ave_exp = round(mean(exp_d_v_h), 3), 
            exp_less_D = length(exp_up_down[exp_up_down == "down"]),
            exp_more_D = length(exp_up_down[exp_up_down == "up"]), 
            exp_neutral = length(exp_up_down[exp_up_down == "non-significant"]),
            ave_fds = round(mean(d_v_h), 3), 
            fds_less_D = length(up_down[up_down == "down"]),
            fds_more_D = length(up_down[up_down == "up"]), 
            fds_neutral = length(up_down[up_down == "non-significant"]),
            all_asv = list(asv_names)
            ) %>%
  filter(!(total_sig == 0 & total_non_time == 0)) %>%
  select(-all_asv) %>%
  filter(total_non_time > 0) %>%
  arrange(desc(total_non_time)) %>%
  rename("All ASVs" = total, "All Significant" = total_sig, "Not Time" = total_non_time)
  

#make color bar nice
#unit.scale = function(x) (x - min(x)) / (max(x) - min(x))
#total = color_bar("#FA614B", fun = unit.scale)



formattable(family_table, align = c("l", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c"), list(
  Family = formatter("span", style = ~ style("font-size:10px", color = "gray",font.weight = "bold")),
  area(col = c(6, 10)) ~ color_tile("white", "#4DAF4A"),
  area(col = c(7, 11)) ~ color_tile("white", "#E41A1C"),
  area(col = c(8, 12)) ~ color_tile("white", "#FF7F00"),
  ave_exp = formatter("span", style = x ~ style(color = ifelse(x < 0, "red", "green"))),
  ave_fds = formatter("span", style = x ~ style(color = ifelse(x < 0, "red", "green")))
)) %>%
  export_formattable("../Figures/fds_exp_diffs.png")



table_prep1 %>% 
  left_join(exposure_dvh_summary) %>% 
  select(-Species) %>%
  group_by(Family) %>%
  summarize(
    total = n(),
    `   ` = " ",
    time = sum(p_time),
    fds = sum(p_final_disease_state),
    exp = sum(p_exposure),
    t_fds = sum(`p_time:final_disease_state`),
    t_exp = sum(`p_time:exposure`),
    ` ` = " ",
    exp_more_D = length(exp_up_down[exp_up_down == "up"]), 
    fds_more_D = length(up_down[up_down == "up"]), 
  ) %>%
  rowwise() %>%
  filter(!((time > 0) & (sum(fds, exp, t_fds, t_exp) < 1))) %>%
  arrange(desc(total)) %>%
  formattable(align = c("l", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c"), list(
    Family = formatter("span", style = ~ style(color = "gray",font.weight = "bold")),
    area(col = 2) ~ color_tile("white", "#539CA0"),
    area(col = 4) ~ color_tile("white", "#01B8F7"),
    area(col = 5) ~ color_tile("white", "#F63517"),
    area(col = 6) ~ color_tile("white", "#F9E820"),
    area(col = 7) ~ color_tile("white", "#A33FFF"),
    area(col = 8) ~ color_tile("white", "#02C40A"),
    area(col = c(10, 11)) ~ color_tile("white", "#E41A1C")
    )) %>%
  export_formattable("../Figures/signifs_by_family.png")

#buffer
#### Graphing ####

#venn diagram
target_data %>%
  mutate(time = if_else(time == 'T0', exposure, time)) %>%
  # count(time, asv_names) %>%
  group_by(time, asv_names) %>%
  summarise(n = sum(value > 7.12)) %>%
  pivot_wider(names_from = time, values_from = n, values_fill = 0L) %>%
  mutate(across(-asv_names, ~. > 0)) %>%
  ggvenn(c('D', 'H', 'T3', 'T7'))


##EMMEANS for LS

#LS FDS only
ls_plot_1 <- ls_model_w_terms %>%
  filter(terms == "final_disease_state") %>%
  filter(asv_names %in% very_likely_suspects) %>% #asv_names for VLS or !asv_names for LS
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
  ggplot(aes(x = fct_reorder(paste(Family, " ", Genus, " (", asv_names, ")", sep = ""), -d_v_h, .desc = TRUE), 
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
ls_plot_2 <- ls_model_w_terms %>%
  filter(terms == "time:final_disease_state") %>%
  filter(asv_names %in% very_likely_suspects) %>% #asv_names for VLS or !asv_names for LS
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
  ggplot(aes(x = fct_reorder(paste(Family, " ", Genus, " (", asv_names, ")", sep = ""), -d_v_h, .desc = TRUE), 
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
ls_plot_3 <- ls_model_w_terms %>%
  filter(terms == "both", param == "time:final_disease_state") %>%
  filter(asv_names %in% very_likely_suspects) %>% #asv_names for VLS or !asv_names for LS
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
  ggplot(aes(x = fct_reorder(paste(Family, " ", Genus, " (", asv_names, ")", sep = ""), -d_v_h, .desc = TRUE), 
             y = emmean, ymin = emmean - SE, ymax = emmean + SE,
             col = time:final_disease_state, pch = time)) +
  geom_pointrange(position = position_dodge(0.5)) +
  geom_text(aes(y = (emmean + SE), label = .group),
            position = position_dodge(0.5), vjust = -1) +
  scale_color_manual(values = c("hotpink1", "deepskyblue", "firebrick1", "dodgerblue3")) +
  coord_flip() +
  xlab("ASV") +
  labs(title = "Both FDS and FDS:Time")

#ls - garbage ASVs
ls_plot_1 / ls_plot_2 / ls_plot_3 + plot_layout(heights = c(1, 20, 15)) + 
  plot_annotation(title = "Likely Suspects") & theme_linedraw() & 
  theme(text = element_text('sans'), plot.title = element_text(size = 16))

#vls
ls_plot_1 / ls_plot_3 + plot_layout(heights = c(21, 25)) + 
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


ggplot(logfold_data %>% filter(!asv_names %in% very_likely_suspects), 
                          aes(x = fct_reorder(paste(Family, " ", 
                              Genus, " (", parse_number(asv_names), ")", sep = ""), diff), y = estimate, 
             ymin = estimate - SE, ymax = estimate + SE, colour = time, pch = time)) +
  geom_hline(yintercept = 0) +
  geom_pointrange(position = position_dodge(0.5)) +
  scale_color_manual(values = c("darkorange1", "firebrick1")) +
  coord_flip() +
  xlab("") +
  labs(title = "LS Logfold Changes") +
  facet_grid(rows = vars(terms), scales = "free", space = "free") +
  theme_bw()


ggplot(logfold_data %>% filter(asv_names %in% very_likely_suspects), aes(x = fct_reorder(paste(Family, " ", 
                   Genus, " (", parse_number(asv_names), ")", sep = ""), diff), y = estimate, 
                   ymin = estimate - SE, ymax = estimate + SE, colour = time, pch = time)) +
  geom_hline(yintercept = 0) +
  geom_pointrange(position = position_dodge(0.5)) +
  scale_color_manual(values = c("darkorange1", "firebrick1")) +
  coord_flip() +
  xlab("") +
  labs(title = "VLS Logfold Changes") +
  facet_grid(rows = vars(terms), scales = "free", space = "free") +
  theme_bw()


#TODO add pcoa code to this doc




#### Comparing Methods ####

plot_method_asvs <- read_rds("../intermediate_files/plot_method_asvs.rds")

#check for LS/VLS and linear/quadratic (plot model iter1)
comp_models <- ls_model_w_terms %>% 
  mutate(ls = TRUE, vls = ifelse(adj_dvh > 0, TRUE, FALSE)) %>%
  select(asv_names, ls, vls) %>%
  distinct() %>%
  full_join(plot_method_asvs) %>%
  mutate(across(c(ls, vls, linear, quadratic), ~ifelse(is.na(.), FALSE, .))) %>%
  left_join(taxonomy_tibble) %>%
  mutate(Name = paste(Family, Genus), .before = asv_names) %>%
  select(-c(Kingdom, Phylum, Class, Order, Family, Genus, Species)) %>%
  mutate(ordering = parse_number(asv_names)/max(parse_number(asv_names))) %>%
  group_by(terms) %>%
  arrange(ordering, .by_group = TRUE) %>%
  ungroup() %>%
  select(-c(ordering, terms)) %>%
  rowwise() %>%
  mutate(both = ifelse(any(ls,vls) & any(linear, quadratic), TRUE, FALSE), .before = ls)

#table comparing ls/vls and plot model iter1
formattable(comp_models, align = c("l", "c", "c", "c"), list(
    Name = formatter("span", style = ~ style(color = "gray",font.weight = "bold")),
    asv_names = formatter("span", style = ~ style(color = "gray",font.weight = "bold")),
    both = formatter("span", style = ~style(color = "white", display = "block", padding = "0 4px", `border-radius` = "4px", 
                                           `background-color` = ifelse(both, "#A846FF", "white"))),
    ls = formatter("span", style = ~style(color = "white", display = "block", padding = "0 4px", `border-radius` = "4px", 
                                           `background-color` = ifelse(ls, "#29CD13", "white"))),
    vls = formatter("span", style = ~style(color = "white", display = "block", padding = "0 4px", `border-radius` = "4px", 
                             `background-color` = ifelse(vls, "#29CD13", "white"))),
    linear = formatter("span", style = ~style(color = "white", display = "block", padding = "0 4px", `border-radius` = "4px", 
                                           `background-color` = ifelse(linear, "#12CFCF", "white"))),
    quadratic = formatter("span", style = ~style(color = "white", display = "block", padding = "0 4px", `border-radius` = "4px", 
                                           `background-color` = ifelse(quadratic, "#12CFCF", "white")))
  )) #%>%
  export_formattable("../Figures/compare_asv_methods.png")

# overlap in ls/vls and plot model iter1
comp_models %>% reframe(asv_names, emily_model = ifelse(any(ls,vls), TRUE, FALSE),
                          jason_model = ifelse(any(linear,quadratic), TRUE, FALSE)) %>% ggvenn()

iter2 <- read_rds("../intermediate_files/both_jasons_models.rds")[[3]]

z_stats <- iter2 %>% select(asv_names, coef_comp) %>% 
  unnest(coef_comp) %>% select(asv_names, contrast, Z) %>%
  mutate(Z = -Z) %>% pivot_wider(names_from = contrast, values_from = Z) %>% 
  rename("linear_z" = linear, "quadratic_z" = quadratic)
  
#the above table + check all sig terms
all_sig_terms_comp <- comp_models %>%
  left_join(table_prep %>% 
  select(-c(d_v_h, up_down, Class, Order, Family, Genus, Species)), by = join_by(asv_names)) %>%
  left_join(z_stats) %>% 
  mutate(linear_z = round(linear_z, 2), quadratic_z = round(quadratic_z, 2)) %>%
  #filter(any(linear, quadratic)) %>% 
  arrange(linear_z) %>% 
  select(Name, asv_names, both, ls, vls, linear, linear_z,
                 quadratic, quadratic_z, p_time, p_final_disease_state,
                 `p_time:final_disease_state`, p_exposure, `p_time:exposure`)

formattable(all_sig_terms_comp, align = c("l", "c", "c", "c"), list(
  Name = formatter("span", style = ~ style(color = "gray",font.weight = "bold")),
  asv_names = formatter("span", style = ~ style(color = "gray",font.weight = "bold")),
  both = formatter("span", style = ~style(color = "white", display = "block", padding = "0 4px", `border-radius` = "4px", 
                                          `background-color` = ifelse(both, "#A846FF", "white"))),
  ls = formatter("span", style = ~style(color = "white", display = "block", padding = "0 4px", `border-radius` = "4px", 
                                        `background-color` = ifelse(ls, "#29CD13", "white"))),
  vls = formatter("span", style = ~style(color = "white", display = "block", padding = "0 4px", `border-radius` = "4px", 
                                         `background-color` = ifelse(vls, "#29CD13", "white"))),
  linear = formatter("span", style = ~style(color = "white", display = "block", padding = "0 4px", `border-radius` = "4px", 
                                            `background-color` = ifelse(linear, "#12CFCF", "white"))),
  quadratic = formatter("span", style = ~style(color = "white", display = "block", padding = "0 4px", `border-radius` = "4px", 
                                               `background-color` = ifelse(quadratic, "#12CFCF", "white"))),
  linear_z = formatter("span", style = ~style(color = "white", display = "block", padding = "0 4px", `border-radius` = "4px", 
                                            `background-color` = ifelse(is.na(linear_z), "white", "#12CFCF"))),
  quadratic_z = formatter("span", style = ~style(color = "white", display = "block", padding = "0 4px", `border-radius` = "4px", 
                                               `background-color` = ifelse(is.na(quadratic_z), "white", "#12CFCF"))),
  p_time = formatter("span", style = ~style(color = "white", display = "block", padding = "0 4px", `border-radius` = "4px", 
                                               `background-color` = ifelse(p_time, "orange", "white"))),
  p_exposure = formatter("span", style = ~style(color = "white", display = "block", padding = "0 4px", `border-radius` = "4px", 
                                            `background-color` = ifelse(p_exposure, "orange", "white"))),
  `p_time:exposure` = formatter("span", style = ~style(color = "white", display = "block", padding = "0 4px", `border-radius` = "4px", 
                                            `background-color` = ifelse(`p_time:exposure`, "orange", "white"))),
  `p_time:final_disease_state` = formatter("span", style = ~style(color = "white", display = "block", padding = "0 4px", `border-radius` = "4px", 
                                            `background-color` = ifelse(`p_time:final_disease_state`, "orange", "white"))),
  p_final_disease_state = formatter("span", style = ~style(color = "white", display = "block", padding = "0 4px", `border-radius` = "4px", 
                                            `background-color` = ifelse(p_final_disease_state, "orange", "white")))
))

## 4 Way Venn - First/Second Model Iteration and FDS/Exposure

venn_sig_terms <- ls_model_full %>% select(asv_names, param, p) %>% mutate(p = p < 0.05) %>%
  pivot_wider(names_from = param, values_from = p) %>% 
  rowwise() %>%
  mutate(exp = ifelse(any(exposure, `time:exposure`), TRUE, FALSE), 
         fds = ifelse(any(final_disease_state, `time:final_disease_state`), TRUE, FALSE)) %>%
  select(asv_names, exp, fds)

first_iter <- read_rds("../intermediate_files/both_jasons_models.rds")[[1]]
second_iter <- read_rds("../intermediate_files/both_jasons_models.rds")[[2]]

comp_models_indepth <- model_data %>% select(asv_names) %>% unique() %>% rowwise() %>%
  mutate(first = ifelse(asv_names %in% first_iter$asv_names, TRUE, FALSE),
         second = ifelse(asv_names %in% second_iter$asv_names, TRUE, FALSE)) %>%
  full_join(venn_sig_terms, by = join_by(asv_names)) %>%
  filter(any(first, second, exp, fds)) %>%
  rename("Exposure" = exp, "FDS" = fds) %>%
  mutate(first__1 = ifelse(first, 1, 0),
         second__1 = ifelse(second, 1, 0),
         Exposure__1 = ifelse(Exposure, 1, 0),
         FDS__1 = ifelse(FDS, 1, 0)) %>%
  rowwise() %>%
  mutate(num_of_interactions = sum(first__1, second__1, Exposure__1, FDS__1)) %>% 
  arrange(desc(num_of_interactions)) %>%
  select(-c(ends_with("__1")))

comp_models_indepth %>%
  relocate(Exposure, .after = first) %>%
  relocate(second, .after = FDS) %>%
  ggvenn()

# 4 Way Venn as a Table - Z SCORES ARE FIXED: pos is more disease
comp_models_indepth %>% full_join(all_sig_terms_comp) %>% select(-starts_with("p_")) %>% 
  select(Name, asv_names, both, ls, vls, FDS, Exposure, first, second, linear, linear_z, quadratic, quadratic_z, num_of_interactions) %>%
  arrange(desc(num_of_interactions)) %>%
  filter(num_of_interactions > 1) %>%
  select(-num_of_interactions) %>%
  mutate(linear = ifelse(linear, linear_z, 0), quadratic = ifelse(quadratic, quadratic_z, 0)) %>%
  select(-c(linear_z, quadratic_z)) %>%
  left_join(taxonomy_tibble) %>%
  mutate(Name = paste(Family, Genus), .before = asv_names) %>%
  select(-c(Kingdom, Phylum, Class, Order, Family, Genus, Species, both, ls)) %>%
  relocate(c(first, linear, quadratic), .after = Exposure) %>%
  relocate(FDS, .before = second) %>%
  formattable(align = c("l", "c", "c", "c"), list(
    Name = formatter("span", style = ~ style(color = "gray",font.weight = "bold")),
    asv_names = formatter("span", style = ~ style(color = "gray",font.weight = "bold")),
    vls = formatter("span", style = ~style(color = "white", display = "block", padding = "0 4px", `border-radius` = "4px", 
                                           `background-color` = ifelse(vls, "#29CD13", "white"))),
    linear = formatter("span", style = ~style(color = "white", display = "block", padding = "0 4px", `border-radius` = "4px", 
                                              `background-color` = ifelse(linear != 0, "#12CFCF", "white"))),
    quadratic = formatter("span", style = ~style(color = "white", display = "block", padding = "0 4px", `border-radius` = "4px", 
                                                 `background-color` = ifelse(quadratic != 0, "#12CFCF", "white"))),
    FDS = formatter("span", style = ~style(color = "white", display = "block", padding = "0 4px", `border-radius` = "4px", 
                                                 `background-color` = ifelse(FDS, "#FC5629", "white"))),
    Exposure = formatter("span", style = ~style(color = "white", display = "block", padding = "0 4px", `border-radius` = "4px", 
                                                 `background-color` = ifelse(Exposure, "#FC5629", "white"))),
    first = formatter("span", style = ~style(color = "white", display = "block", padding = "0 4px", `border-radius` = "4px", 
                                                 `background-color` = ifelse(first, "#12CFCF", "white"))),
    second = formatter("span", style = ~style(color = "white", display = "block", padding = "0 4px", `border-radius` = "4px", 
                                                 `background-color` = ifelse(second, "#12CFCF", "white")))
  ))


#gather emmeans plots for disease state or exposure significance
model_comp_emmeans_plots <- ls_model %>%
  ungroup %>%
  filter(asv_names %in% (comp_models_indepth %>% filter(num_of_interactions >= 3) %>% .$asv_names)) %>%
  rowwise %>%
  mutate(terms = list(find_unique_significant_terms_rmANOVA(model, 0.05))) %>%
  unnest(terms, keep_empty = TRUE) %>%
  rowwise %>%
  mutate(em_out = list(possibly(make_emmean_model, otherwise = NULL)(model, 
                                                                     as.formula(str_c('~', terms)), 
                                                                     0.05))) %>%
  mutate(plot = list(possibly(make_model_plot, otherwise = NULL, quiet = TRUE)
                     (em_out, data, terms))) %>%
  group_by(asv_names) %>%
  reframe(plot = ifelse(any(is.na(terms)), #if NAs in data, then no plot.  else, plot
                        list(NULL),
                        list(wrap_plots(plot) &
                               plot_annotation(title = asv_names) & 
                               theme_bw())))

comp_models_plots <- read_rds("../intermediate_files/comp_models_plots_1_2.rds")

#all plots for the ASVs with 3+ interactions
all_plots_comp_models <- model_comp_emmeans_plots %>% 
    full_join(comp_models_plots) %>%
    filter(asv_names %in% (comp_models_indepth %>% filter(num_of_interactions >= 3) %>% .$asv_names)) %>%
    pivot_longer(cols = where(is.list),
                 names_to = 'name',
                 values_to = 'plot') %>%
    rowwise %>%
    mutate(missing_plot = !is.null(plot)) %>%
    ungroup %>%
    filter(missing_plot) %>%
    group_by(asv_names) %>%
    summarise(plots = list(wrap_plots(plot) &
                             plot_annotation(title = asv_names) & 
                             theme_bw()),
              n_plot = n()) %>%
  left_join(taxonomy_tibble %>% select(asv_names, Family, Genus)) %>%
  relocate(c(Family, Genus), .before = asv_names)
  
#### Patchwork PCoA Plot for ASVs, Genus, Family ####
microbiome_data <- read_rds("../intermediate_files/preprocess_microbiome.rds") %>%
  subset_samples(time %in% c('T3', 'T7'))
metadata <- sample_data(microbiome_data) %>%
  as_tibble(rownames = 'sample_id') %>%
  dplyr::select(-retain_sample) %>%
  mutate(fragment_id = str_c(str_replace_na(exposure, 'NA'), tank, genotype, sep = '_'),
         .after = sample_id)

agg_levels <- c("none", "Genus", "Family")

for(i in 1:length(agg_levels)){
  aggregation_level = agg_levels[i]
  
  if(aggregation_level != 'none'){
    microbiome_data <- aggregate_taxa(microbiome_data, aggregation_level)
    taxa_names(microbiome_data) <- str_replace_all(taxa_names(microbiome_data), ' |-', '_')
    agg_title <-  aggregation_level
  } else {
    taxa_names(microbiome_data) <- str_c('ASV', 1:length(taxa_names(microbiome_data)), sep = '_')
    agg_title <-  "ASV"
  }
  
  if(i == 1){
    otu_tmm <- microbiome_data %>%
      phyloseq_filter_prevalence(prev.trh = 0.1) %>%
      otu_table() %>% 
      t %>% #NOTE: *genus and family do not need the t but ASVs need the t*
      as.data.frame %>%
      as.matrix %>% 
      DGEList(remove.zeros = TRUE) %>%
      edgeR::calcNormFactors(method = 'TMMwsp')
    plota <- plot_agg_pcoa(cpm(otu_tmm, log = TRUE, prior.count = 2))
    
  }else{ 
    otu_tmm <- microbiome_data %>%
      phyloseq_filter_prevalence(prev.trh = 0.1) %>%
      otu_table() %>% 
      as.data.frame %>%
      as.matrix %>% 
      DGEList(remove.zeros = TRUE) %>%
      edgeR::calcNormFactors(method = 'TMMwsp')
    if(i == 2){
      plotb <- plot_agg_pcoa(cpm(otu_tmm, log = TRUE, prior.count = 2)) +
        scale_x_reverse() #+
        #scale_y_reverse()
    }else if(i == 3){
      plotc <- plot_agg_pcoa(cpm(otu_tmm, log = TRUE, prior.count = 2))
    }
  }
}
plota + plotb + plotc



#### NMDS ####

mb_data <- microbiome_data %>%
  phyloseq_transform_css %>% #normalizing by column and then log transform
  phyloseq_filter_prevalence(prev.trh = 0.1) %>% 
  otu_table %>%
  t

otu_nmds <- metaMDS(mb_data, distance = 'mountford', k = 2, trymax = 100, autotransform = FALSE, verbose = TRUE)

nmds_plot_exp <- scores(otu_nmds)$sites %>%
  as_tibble(rownames = 'sample_id') %>%
  left_join(metadata, by = 'sample_id') %>%
  
  ggplot(aes(x = NMDS1, y = NMDS2, colour = exposure, shape = time)) +
  geom_point(data = as_tibble(scores(otu_nmds)$species, rownames = aggregation_level),
             colour = 'gray60', size = 0.1, shape = 'circle') +
  
  geom_point() +
  theme_bw() +
  labs(title = "Exposure") +
  labs(color = "Exposure")

nmds_plot_fds <- scores(otu_nmds)$sites %>%
  as_tibble(rownames = 'sample_id') %>%
  left_join(metadata, by = 'sample_id') %>%
  
  ggplot(aes(x = NMDS1, y = NMDS2, colour = final_disease_state, shape = time)) +
  geom_point(data = as_tibble(scores(otu_nmds)$species, rownames = aggregation_level),
             colour = 'gray60', size = 0.1, shape = 'circle') +
  
  geom_point() +
  theme_bw() +
  labs(title = "Final Disease State") +
  labs(color = "FDS")

nmds_plot_exp | nmds_plot_fds

#### Microshades Microbe Abundances ####


mdf_prep_test1 <- microbiome_data %>%
           tax_glom("Genus") %>%
           psmelt() 

mdf_processed_data <- mdf_prep_test1 %>%
  filter(Abundance > 0) %>%
  mutate(category = paste(time, exposure, final_disease_state, 
                          sep = "_")) %>%
  mutate(category = ifelse(str_detect(category, "Field"), "T0", category)) %>%
  mutate(time = ifelse(category %in% c("T0_D_D", "T0_H_H"), "Doses", time)) %>%
  mutate(category = ifelse(category == "T0_D_D", "Diseased", ifelse(category == "T0_H_H", "Healthy", category))) %>%
  group_by(category) %>%
  mutate(total = sum(Abundance)) %>%
  ungroup() %>%
  select(-c(retain_sample)) %>%
  group_by(category, Genus) %>%
  reframe(Kingdom, Phylum, Class, Order, Family, time, total, rel_abun = sum(Abundance)/total) %>%
  distinct() %>%
  rename(Sample = category, Abundance = rel_abun) %>%
  as.data.frame()


color_objs_microbes <- create_color_dfs(mdf_processed_data, group_level = "Order", 
                                       selected_groups = c("Rickettsiales", "Enterobacterales", "Flavobacteriales", 
                                                           "Pseudomonadales",  "Rhodobacterales"), top_n_subgroups = 4,  cvd = TRUE)
mdf_microbes <- color_objs_microbes$mdf
cdf_microbes <- color_objs_microbes$cdf

legend_microbes <-custom_legend(mdf_microbes, cdf_microbes, group_level = "Order")

plot_microbes_prelim <- plot_microshades(mdf_microbes, cdf_microbes) + 
  scale_y_continuous(labels = scales::percent, expand = expansion(0)) +
  facet_grid(cols = vars(time), scales = "free", space = "free") +
  theme_bw() +
  theme(legend.position = "none", plot.margin = margin(6,20,6,6))
  
plot_grid(plot_microbes_prelim, legend_microbes,  rel_widths = c(1, .25))







