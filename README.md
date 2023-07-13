# 16S Florida Tank Analysis

### Somewhat Unfiltered Data
only filtered for 20% prevalence & less than 90% missingness

![plot](./Figures/cu459.png)

![plot](./Figures/venn459.png)

### Everything that goes into the model:
- filter 20% prevalence, less than 90% missingness, remove ASVs with an average of less than *100* cpm per sample
- removed anything that was only T3 & T7 or only T3 or only T7

![plot](./Figures/venn_min_filter.png)

![plot](./Figures/cu_min_filter.png)

## Current Model

from Vega Thurber et al 2020:
![plot](./Figures/vegathurber_dysbiosis.jpeg)

- filter 20% prevalence, less than 90% missingness, remove ASVs with an average of less than *100* cpm per sample
- removed anything that was only T3 & T7 or only T3 or only T7  
- lmer model with formula of "log2_cpm ~ treatment + (1 | genotype) + (1 | tank)"
- treatment is time, exposure, and susceptibility all combined into one descriptor
- examine fdr corrected p-values for treatment, genotype, and tank: only keep ASVs w/ significant effect of treatment
- make planned comparisons, then group ASVs into early/late/continuous and probiotic/opportunist/pathogen based on expected profiles for those strategies

**Early Pathogen:** Disease-Exposed Susceptible is greater than Disease-Exposed Resistant at T3  
AND Disease-Exposed Susceptible is greater than the average of all other T3 treatments   
AND Disease-Exposed Susceptible must be more abundant at T3 than T0  
**Late Pathogen:** Disease-Exposed Susceptible is greater than Disease-Exposed Resistant at T7  
AND Disease-Exposed Susceptible is greater than the average of all other T7 treatments   
AND Disease-Exposed Susceptible must be more abundant at T7 than T0  
**Continuous Pathogen:** meets all criteria for both early and late pathogens

**Early Opportunist:** the average of Disease-Exposed is greater than the average of Healthy-Exposed at T3  
AND Disease-Exposed is more abundant at T3 than T0  
**Late Opportunist:** the average of Disease-Exposed is greater than the average of Healthy-Exposed at T7  
AND Disease-Exposed is more abundant at T7 than T0  
**Continuous Opportunist:** meets criteria for both early and late opportunists  

**Probiotic_T7_Strict:** Disease-Exposed Resistant is more abundant than Disease-Exposed Susceptible at T7   
AND Disease-Exposed Resistant is more abundant than all other treatments at T7  
AND Disease-Exposed Resistant is more abundant at T7 than at T3  
AND Disease-Exposed Resistant is more abundant at T7 than Resistant at T0  

**Crasher_T3:** Susceptible is more abundant at T0 than Disease-Exposed Susceptible is at T3  
**Crasher_T7:** Susceptible is more abundant at T0 than Disease-Exposed Susceptible is at T7    

##### Planned Comparisons where Nothing was Found:

-- none in T0,T3 and anything in probiotic T7 but not probiotic T7 Strict did not look like a probiotic (all treatments had very similar values in most cases) --
**Probiotic_T0:** Disease-Exposed Resistant is more abundant than Disease-Exposed Susceptible at T0  
**Probiotic_T3:** Disease-Exposed Resistant is more abundant than Disease-Exposed Susceptible at T3  
**Probiotic_T7:** Disease-Exposed Resistant is more abundant than Disease-Exposed Susceptible at T7  

**Probiotic_T3_Strict:** Disease-Exposed Resistant is more abundant than Disease-Exposed Susceptible at T3  
AND Disease-Exposed Resistant is more abundant than all other treatments at T3  
AND Disease-Exposed Resistant is more abundant at T3 than Resistant at T0

**Crasher_T3_Strict:** Susceptible is more abundant at T0 than Disease-Exposed Susceptible is at T3  
AND Disease-Exposed Susceptible is less than the average of all other treatments at T3  
**Crasher_T7_Strict:** Susceptible is more abundant at T0 than Disease-Exposed Susceptible is at T7
AND Disease-Exposed Susceptible is less than the average of all other treatments at T7  

##### Model Results:

![plot](./Figures/comp_upset_bacstrat_7_11_23.png)

removed T3 and T7 from complex upset because all ASVs are present in both
![plot](./Figures/comp_upset_bacstrat_origin_7_11_23.png)

## significant ASVs - 0.01%

![plot](./Figures/bac_strat_cp1.png)

![plot](./Figures/bac_strat_cp2.png)

![plot](./Figures/bac_strat_cp3.png)

![plot](./Figures/bac_strat_cp4.png)

![plot](./Figures/bac_strat_cp5.png)

![plot](./Figures/bac_strat_cp6.png)

![plot](./Figures/bac_strat_cp7.png)

##### Different Levels of Filtering

-- in 0.01% of samples --  
![plot](./Figures/venn0.01.png)  
![plot](./Figures/bac_strat_cu0.01.png)  

-- in 0.1% of samples --  
![plot](./Figures/venn0.1.png)  
![plot](./Figures/bac_strat_cu0.1.png)  

#### ASV NMDS by bacterial strategy

![plot](./Figures/asv_nmds1.png) 
ASV as species, sample id as site

![plot](./Figures/sample_nmds.png)  
gray dots are ASVs, colored dots are individual coral fragments

### Correlation Tests
- correlation test (cor.test) between log2cpm abundance and continuous resistance value, p < 0.05 at each time point
- INITIAL DATA FILTERS: filter 20% prevalence, less than 90% missingness, remove ASVs with an average of less than *100* cpm per sample

##### Correlations based only on T0
![plot](./Figures/disease_corr_t0.png)

![plot](./Figures/disease_corr_t3.png)  

![plot](./Figures/disease_corr_t7.png)  


### Microbe Abundances
- category titles are "timepoint_exposure_susceptibility"

![plot](./Figures/microshades_og_susceptibility.png)

# OLD MODEL FIGURES

##### Exposure:Final Disease State
- category titles are "timepoint_exposure_final disease state
- **removed** samples that were healthy exposed and contracted WBD

![plot](./Figures/microshades_ordergenus.png)

![plot](./Figures/microshades_classfamily.png)


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

