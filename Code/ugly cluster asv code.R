thala <- tax_table(microbiome_raw) %>% 
  as.data.frame %>%
  as_tibble(rownames = "sequence") %>%
  filter(Genus == 'Thalassotalea') %>%
  mutate(asv_id = str_c('asv_', row_number()))


thala %>%
  count(sequence) %>%
  filter(n > 1)
  

library(bioseq)
thala_seq <- thala %$%
  set_names(sequence, asv_id) %>%
  as_dna()
length(thala_seq)
seq_nchar(thala_seq)
seq_cluster(thala_seq)


library(Biostrings)
thala_seq2 <- as_DNAbin(thala_seq) %>%
  as.character %>% 
  lapply(.,paste0,collapse="") %>%
  unlist %>% 
  DNAStringSet

Biostrings::writeXStringSet(thala_seq2, '../intermediate_files/Thalassotalea_asvs.fasta')

# module load anaconda3
# #conda create -n clustalo -c bioconda clustalo
# source activate clustalo
# cd /scratch/j.selwyn
# clustalo \
# -i Thalassotalea_asvs.fasta \
# -o Thalassotalea_asvs_msa.fasta \
# --threads=${SLURM_CPUS_PER_TASK} \
# --verbose


cluster_id <- readDNAMultipleAlignment('../intermediate_files/Thalassotalea_asvs_msa.fasta') %>%
  as.character() %>%
  as_dna() %>%
  seq_cluster(threshold = 0.03, method = 'average') 

thala %>%
  mutate(cluster = cluster_id) %>%
  count(cluster) %>%
  arrange(-n)
