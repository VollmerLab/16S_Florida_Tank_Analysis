# 16S Florida Tank Analysis

### Filters
- minimum 10,000 reads  
- 10% prevalence via phyloseq_filter_prevalence(prev.trh = 0.1)  
- each ASV must be present in at least 10% of individuals
- Analyzing data for Timepoints 3 and 7
- Must be present in the Diseased Homogenate, T3 and T7 samples

### Remaining ASVs:
##### Sections of interest are 220 and 85 for a total of 305 ASVs
![plot](./Figures/Venn.png)


##### 305 ASVs gets reduced to 249 ASVs because 56 are significant for nothing

![plot](./Figures/comp_upset_full_subset.png)


#### More Abundant in Disease:

![plot](./Figures/comp_upset_more_disease.png)  

#### More Abundant in Healthy:

![plot](./Figures/comp_upset_more_healthy.png)

#### If you remove ASVs that are only significant for time:
##### Only 116 remain, 57 are more abundant in Disease

![plot](./Figures/comp_upset_timeless.png)

#### Relative Representation of Families of Interest
Of the total number of ASVs present in all of the data within each Order/Family, what percent of those are found in the Likely or Very Likely Suspect Lists?  

Of the 116 remaining, 57 are significant for final disease state or final disease state:time (more abundant in either healthy or diseased)

![plot](./Figures/rel_abun_by_order.png)

![plot](./Figures/rel_abun_by_family.png)

### Highest Abundances 

![plot](./Figures/most_abun_baits.png)

![plot](./Figures/most_abun_timepoints.png)


### Likely Suspects
*present in T3, T7, and Diseased bait*
*AND differs significantly for either final disease state or the interaction of final disease state and time*  


### Very Likely Suspects
*present in T3, T7, and Diseased bait*   
*AND differs significantly for either final disease state or the interaction of final disease state and time*  
*AND more abundant in diseased than healthy*

![plot](./Figures/VLS_emmeans.png)

## Logfold Changes in Bacterial Abundance
- positive logfold change means more disease

![plot](./Figures/logfold_changes.png)



### PCoA and NMDS

![plot](./Figures/PCoA_multiple_aggregations.png)

![plot](./Figures/Mountford_NMDS.png)

