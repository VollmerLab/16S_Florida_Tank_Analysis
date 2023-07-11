
sig_pairs <- asv_models %>%
#filter(asv_id == "ASV_202") %>%
  rowwise %>%
  mutate(comparisons = list(emmeans(model, ~treatment) %>%
                            contrast("pairwise") %>%
                            broom::tidy(conf.int = TRUE) %>%
                            filter(contrast %in% c("T0_F_R - T0_F_S", "T3_D_R - T3_D_S",
                                                   "T3_D_R - T3_H_R", "T3_D_S - T3_H_S", 
                                                   "T3_H_R - T3_H_S",
                                                   "T7_D_R - T7_D_S", "T7_D_R - T7_H_R",
                                                   "T7_D_S - T7_H_S", "T7_H_R - T7_H_S")) %>%
                              filter(adj.p.value < 0.05))) %>%
  select(asv_id, comparisons) %>%
  unnest(comparisons)
                         

sig_pairs %>%
  select(asv_id, contrast, estimate) %>%
  mutate(contrast = case_when(contrast == "T0_F_R - T0_F_S" ~ "T0.R-S",
                              contrast == "T3_D_R - T3_D_S" ~ "T3.DR-DS",
                              contrast == "T3_D_R - T3_H_R" ~ "T3.DR-HR",
                              contrast == "T3_D_S - T3_H_S" ~ "T3.DS-HS", 
                              contrast == "T3_H_R - T3_H_S" ~ "T3.HR-HS",
                              contrast == "T7_D_R - T7_D_S" ~ "T7.DR-DS",
                              contrast == "T7_D_R - T7_H_R" ~ "T7.DR-HR",
                              contrast == "T7_D_S - T7_H_S" ~ "T7.DS-HS", 
                              contrast == "T7_H_R - T7_H_S" ~ "T7.HR-HS")) %>%
  pivot_wider(names_from = contrast, values_from = estimate) %>%
  select(asv_id, `T0.R-S`, `T3.HR-HS`, `T3.DS-HS`, `T3.DR-HR`, `T3.DR-DS`, `T7.HR-HS`, `T7.DS-HS`,
         `T7.DR-HR`, `T7.DR-DS`) %>%
  formattable(align = c("l", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c"), list(
    #Family = formatter("span", style = ~ style("font-size:10px", color = "gray",font.weight = "bold")),
    `T0.R-S` = formatter("span", style = x ~ style(color = ifelse(is.na(x), "white", ifelse(x < 0, "red", "green")))),
    `T3.DR-DS` = formatter("span", style = x ~ style(color = ifelse(is.na(x), "white", ifelse(x < 0, "red", "green")))),
    `T3.DR-HR` = formatter("span", style = x ~ style(color = ifelse(is.na(x), "white", ifelse(x < 0, "red", "green")))),
    `T3.DS-HS` = formatter("span", style = x ~ style(color = ifelse(is.na(x), "white", ifelse(x < 0, "red", "green")))),
    `T3.HR-HS` = formatter("span", style = x ~ style(color = ifelse(is.na(x), "white", ifelse(x < 0, "red", "green")))),
    `T7.DR-DS` = formatter("span", style = x ~ style(color = ifelse(is.na(x), "white", ifelse(x < 0, "red", "green")))),
    `T7.DR-HR` = formatter("span", style = x ~ style(color = ifelse(is.na(x), "white", ifelse(x < 0, "red", "green")))),
    `T7.DS-HS` = formatter("span", style = x ~ style(color = ifelse(is.na(x), "white", ifelse(x < 0, "red", "green")))),
    `T7.HR-HS` = formatter("span", style = x ~ style(color = ifelse(is.na(x), "white", ifelse(x < 0, "red", "green"))))
  )) #%>%
  #export_formattable("../Figures/name.png")

                          
                         
  