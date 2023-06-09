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





,
d_v_h =  emmeans(model, ~final_disease_state) %>%
  as_tibble %>%
  select(emmean) %>%
  pull(1) %>%
  diff




GP_legend_new <-custom_legend(mdf_GP_test, cdf_GP_test, group_level = "Order")

plot_diff <- plot_microshades(mdf_GP_test, cdf_GP_test) + 
  scale_y_continuous(labels = scales::percent, expand = expansion(0)) +
  theme(legend.position = "none")  +
  #theme(axis.text.x = element_text(size= 6)) +
  facet_wrap(~time, scales = "free_x", nrow = 1) +
  #theme(axis.text.x = element_text(size= 6)) + 
  theme(plot.margin = margin(6,20,6,6))

plot_grid(plot_diff, GP_legend_new,  rel_widths = c(1, .25))


### our data

all_model_data <-  cpm(otu_tmm, log = TRUE, prior.count = 2) %>%
  t %>%
  as_tibble(rownames = "sample_id") %>%
  full_join(metadata, by = "sample_id") %>%
  pivot_longer(cols = -any_of(colnames(metadata)), 
               names_to = "asv_names", values_to = "value") %>%
  filter(tank != "homogenate_fragment") %>%
  filter(value > 7.13) %>%
  mutate(final_disease_state = ifelse(exposure == "Field", "Field", final_disease_state)) %>%
  mutate(across(c(exposure, final_disease_state), factor)) %>%
  left_join(taxonomy_tibble) %>%
  select(-disease_state) %>%
  group_by(time, exposure, final_disease_state) %>%
  mutate(tot = sum(value)) %>%
  ungroup() %>%
  group_by(time, exposure, final_disease_state, Kingdom, Phylum, Class, Order, Family, Genus, tot) %>%
  reframe(val = sum(value), Abundance = sum(value)/tot) %>%
  distinct() %>%
  mutate(Sample = case_when(exposure == "Field" ~ "T0",
                            time == "T0" & exposure != "Field" & final_disease_state == "H" ~ "Healthy Dose",
                            time == "T0" & exposure != "Field" & final_disease_state == "D" ~ "Diseased Dose",
                            time %in% c("T3", "T7") ~ paste(time, exposure, final_disease_state, sep = "_"))) %>%
  as.data.frame()

color_objs_GP_test <- create_color_dfs(all_model_data, group_level = "Order", 
                                       selected_groups = c("Rickettsiales", "Enterobacterales", "Flavobacteriales", 
                                                           "Pseudomonadales", "Rhodobacterales"), cvd = TRUE)
mdf_GP_test <- color_objs_GP_test$mdf
cdf_GP_test<- color_objs_GP_test$cdf

GP_legend_new <-custom_legend(mdf_GP_test, cdf_GP_test, group_level = "Order")

plot_diff <- plot_microshades(mdf_GP_test, cdf_GP_test) + 
  scale_y_continuous(labels = scales::percent, expand = expansion(0)) +
  theme(legend.position = "none")  +
  #theme(axis.text.x = element_text(size= 6)) +
  facet_wrap(~time, scales = "free_x", nrow = 1) +
  #theme(axis.text.x = element_text(size= 6)) + 
  theme(plot.margin = margin(6,20,6,6))

plot_grid(plot_diff, GP_legend_new,  rel_widths = c(1, .25))

#### ####

color_test %>% group_by(Sample) %>% summarize(tot = sum(Abundance)) %>% arrange(desc(tot))

color_test %>% filter(Sample == "T3_H_H4_U67")

remotes::install_github("david-barnett/microViz")


use_asv <- 'ASV_401'

blah <- filter(iter2, asv_names %in% use_asv) %>%
  mutate(healthy_em = list(emmeans(healthy_regression, ~time) %>% tidy(conf.int = TRUE) %>% 
                             mutate(exposure = 'H', final_disease_state = 'H')),
         final_em = list(emmeans(final_outcome_model, ~time * final_disease_state) %>% tidy(conf.int = TRUE) %>%
                           mutate(exposure = 'D'))) %>%
  select(asv_names, data, healthy_em, final_em) %>%
  ungroup


bind_rows(blah$healthy_em[[1]],
          blah$final_em[[1]]) %>%
  filter((time == 'T0' & exposure == 'H') | time != 'T0') %>%
  mutate(exposure = if_else(time == 'T0', 'Field', exposure),
         final_disease_state = if_else(time == 'T0', 'Field', final_disease_state)) %>%
  
  ggplot(aes(x = time, y = estimate, ymin = conf.low, ymax = conf.high, shape = exposure,
             colour = final_disease_state)) +
  geom_pointrange(position = position_dodge(0.5), size = 1) +
  geom_point(data = mutate(blah$data[[1]],
                           final_disease_state = if_else(time == 'T0', 'Field', final_disease_state)), 
             inherit.aes = FALSE,
             aes(x = time, y = value, shape = exposure, colour = final_disease_state),
             position = position_dodge(0.5), size = 0.5)
