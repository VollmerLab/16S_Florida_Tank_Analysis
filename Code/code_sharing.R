clean_afex <- function(model){
  #model is an afex::mixed model
  model$anova_table %>%
    as_tibble(rownames = 'param') %>%
    janitor::clean_names() %>%
    rename_with(~str_replace_all(., '_df', 'DF')) %>%
    rename(pvalue = pr_f) %>%
    mutate(param = str_replace(param, ':', 'X'),
           param = str_replace(param, 'final_disease_state', 'finalDisease')) %>%
    pivot_wider(names_from = 'param',
                values_from = where(is.numeric),
                names_vary = 'slowest')
}
