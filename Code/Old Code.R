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