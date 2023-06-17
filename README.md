# 16S Florida Tank Analysis

### Microbe Abundances
- category titles are "timepoint_exposure_susceptibility"  

![plot](./Figures/microshades_og_susceptibility.png)

- category titles are "timepoint_exposure_final disease state
- **removed** samples that were healthy exposed and contracted WBD

![plot](./Figures/microshades_ordergenus.png)

![plot](./Figures/microshades_classfamily.png)

### Current Model

- lmer model with formula of "log2_cpm ~ treatment + (1 | genotype) + (1 | tank)"
- examine fdr corrected p-values for treatment, genotype, and tank: only keep ASVs w/ significant effect of treatment
- make planned comparisons, then group ASVs into early/late/continuous and probiotic/opportunist/pathogen based on expected profiles for those strategies

![plot](./Figures/venn_nm.png)

![plot](./Figures/cu_bacterial_strategies.png)

##### significant ASVs

![plot](./Figures/bac_strat_cp1.png)

![plot](./Figures/bac_strat_cp2.png)

![plot](./Figures/bac_strat_cp3.png)

![plot](./Figures/bac_strat_cp4.png)

![plot](./Figures/bac_strat_cp5.png)

![plot](./Figures/bac_strat_cp6.png)



###### compared to random forest model:
![plot](./Figures/nm_comparison_venn.png)
(see random_forest_analysis.md for more details and plots)

### Trees

![plot](./Figures/Tree_Colwelliaceae.png)

![plot](./Figures/Tree_flavobacteriaceae.png)

![plot](./Figures/Tree_francisellaceae.png)

![plot](./Figures/Tree_rhodobacteraceae.png)

![plot](./Figures/Tree_saccharospirillaceae.png)


# OLD MODEL FIGURES

### Remaining ASVs:
##### Sections of interest are 278 and 104 for a total of 305 ASVs
![plot](./Figures/Venn.png)


##### 382 ASVs gets reduced to 304 ASVs because 78 are significant for nothing

![plot](./Figures/comp_upset_full_subset.png)  

#### If you remove ASVs that are only significant for time, only 172 remain:
##### legend is describing diseased relative to healthy ("up" means it's significantly higher in disease, etc...)

![plot](./Figures/comp_upset_ls.png)

##### only ASVs that are significantly more abundant in D than H:
###### updates to preprocessing data have increased this subset by 35
![plot](./Figures/comp_upset_vls.png)

#### Significant ASV Counts by Family
removed any families where the only significant term in the whole family is time
![plot](./Figures/signifs_by_family.png)

#### Differences in exposure and final disease state
![plot](./Figures/fds_exp_diffs.png)

### Most Abundant Families in each Timepoint and Final Disease State intersection
**aggregated by Family**  

- selected the 20 most abundant Families in each group and ranked 1-20 (1 is most abundant)  
- removed any Family that is only present in one of the six groups

![plot](./Figures/most_abun_timepoints.png)


### Likely Suspects
*present in T3, T7, and Diseased bait*  
*AND differs significantly for either final disease state or the interaction of final disease state and time*  

### Very Likely Suspects
*present in T3, T7, and Diseased bait*   
*AND differs significantly for either final disease state or the interaction of final disease state and time*  
*AND significantly more abundant in diseased than healthy*

![plot](./Figures/emmeans_VLS.png)

## Logfold Changes in Bacterial Abundance
- positive logfold change means more disease
- LS vs VLS almost perfectly separates by whether it's more abundant at T3 (LS) or T7 (VLS)

![plot](./Figures/logfold_VLS.png)

### Potential Link between Rickettsiales and Disease

![plot](./Figures/rickettsiales_corr.png)

### PCoA

![plot](./Figures/pcoa_fds.png)

![plot](./Figures/pcoa_exp.png)

### NMDS

![plot](./Figures/dual_nmds_plots.png)

