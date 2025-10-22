## Code to pre-process microbiome sample data into easier to use format

setwd("/Users/Emily/Desktop/GitHub/16S_Florida_Tank_Analysis/Code")

#### Libraries ####
library(tidyverse)
library(magrittr)
library(phyloseq)

#### Read In Data ####
preprocess_metadata <- read_csv('../intermediate_files/preprocess_metadata.csv', 
                                show_col_types = FALSE)

#### Data Filtering ####

#make a list of all samples that were exposed to healthy homogenate but contracted disease
H_D_samples <- preprocess_metadata %>% filter((exposure == "H" & final_disease_state == "D")) %>%
  pull(sample_id)

#filter out host/symbiont contamination
microbiome_raw <- read_rds("../Data/updated_8_25_23_decipher_16s_ps.rds") %>%
  subset_taxa(Domain == "Bacteria" & 
                Phylum != "Cyanobacteria" &
                !is.na(Phylum) &
                Family != "Mitochondria" &
                Class != "Chloroplast" & 
                Order != "Chloroplast")

#alternate filtering that doesn't remove putative symbiont ASVs
#microbiome_raw <- read_rds("../Data/updated_8_25_23_decipher_16s_ps.rds") %>%
#  subset_taxa(Domain == "Bacteria" & 
#                !is.na(Phylum) &
#                Family != "Mitochondria")

#remove healthy exposed samples that contract disease
microbiome_raw <- subset_samples(microbiome_raw, !sample_names(microbiome_raw) %in% H_D_samples) 
#remove fragments used to make the doses
microbiome_raw <- subset_samples(microbiome_raw, time != "Acerv")

### Check initial values

#234 samples
sample_data(microbiome_raw) %>% as_tibble() %>% nrow()

#visualize distribution of read counts
sample_sums(microbiome_raw) %>%
  as_tibble() %>%
  mutate(sample = names(sample_sums(microbiome_raw)), .before = value) %>%
  rename("read_count" = "value") %>%
  ggplot() +
  geom_histogram(aes(read_count), bins = 30, fill = "#a3a3a8") +
  geom_vline(aes(xintercept = 1000), col = "#cc0606", lwd = 1, lty = "dashed") +
  theme_bw() +
  xlab("Number of Reads") +
  ylab("Number of Samples")
#export 900x500

#how many homogenate dose samples
sample_data(microbiome_raw) %>% as_tibble() %>% filter(tank == "HOMO" | is.na(exposure)) %>% nrow()

# T0 number of samples
sample_data(microbiome_raw) %>% as_tibble() %>% filter(time == "T0" & tank != "HOMO") %>% nrow()

# T0 number of genotypes
sample_data(microbiome_raw) %>% as_tibble() %>% filter(time == "T0" & tank != "HOMO") %>% pull(genotype) %>% unique() %>% length()

# T3/T7 number of samples
sample_data(microbiome_raw) %>% as_tibble() %>% filter(time == "T3" | time == "T7") %>% nrow()

# T3/T7 number of genotypes
sample_data(microbiome_raw) %>% as_tibble() %>% filter(time == "T3" | time == "T7") %>% pull(genotype) %>% unique() %>% length()

### Filter for at least 1000 reads
microbiome_raw <- prune_samples(sample_sums(microbiome_raw)>=1000, microbiome_raw)

### Check values after 1000 read filter

#227 samples
sample_data(microbiome_raw) %>% as_tibble() %>% nrow()

#how many homogenate dose samples
sample_data(microbiome_raw) %>% as_tibble() %>% filter(tank == "HOMO" | is.na(exposure)) %>% nrow()

# T0 number of samples
sample_data(microbiome_raw) %>% as_tibble() %>% filter(time == "T0" & tank != "HOMO") %>% nrow()

# T0 number of genotypes
sample_data(microbiome_raw) %>% as_tibble() %>% filter(time == "T0" & tank != "HOMO") %>% pull(genotype) %>% unique() %>% length()

# T3/T7 number of samples
sample_data(microbiome_raw) %>% as_tibble() %>% filter(time == "T3" | time == "T7") %>% nrow()

# T3/T7 number of genotypes
sample_data(microbiome_raw) %>% as_tibble() %>% filter(time == "T3" | time == "T7") %>% pull(genotype) %>% unique() %>% length()

#### Process Data ####
#process T0, T3, T7 samples
tank_samples <- sample_data(microbiome_raw) %>%
  as_tibble(rownames = 'sample_id') %>% 
  filter(tank != 'HOMO') %>% #remove homogenate doses
  mutate(exposure = if_else(time == 'T0', 'F', exposure), #set exposure to F (Field) for field-collected T0 samples
         tank = if_else(time == 'T0', 'Field', tank), #set tank for field-collected T0 samples
         sample_id_store = sample_id, #store original sample names
         sample_id = str_replace(sample_id, 'REP1_[RS]', 'F_Field')) %>% #make T0 sample names consistent with preprocess_metadata
  full_join(preprocess_metadata,
          by = join_by(sample_id, time, exposure, tank, genotype)) %>% #add metadata
  mutate(sample_id = if_else(tank == 'Field', NA_character_, sample_id)) %>% #set T0 sample names to NA
  mutate(sample_id = coalesce(sample_id, sample_id_store), .keep = 'unused') #reset T0 sample names to original names
  
#process homogenate samples
homogenate_samples <- sample_data(microbiome_raw) %>%
  as_tibble(rownames = 'sample_id') %>%
  filter(tank == 'HOMO') %>% #keep only homogenate doses
  mutate(final_disease_state = exposure)

#Recombine and put data into sample_data slot of phyloseq object
processed_microbiome <- microbiome_raw

sample_data(processed_microbiome) <- bind_rows(homogenate_samples, tank_samples) %>% 
  column_to_rownames('sample_id')

#### Output pre-processed microbiome as RDS file ####
write_rds(processed_microbiome, '../intermediate_files/preprocess_microbiome.rds')
