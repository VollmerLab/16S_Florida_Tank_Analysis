taxonomy_tibble <- tax_table(microbiome_data) %>% 
  as.data.frame %>%
  as_tibble(rownames = "asv_names")

#read in data
updated_taxonomy <- read_csv('../intermediate_files/updated_taxonomy.csv')

#replace everything less than 80% confidence with NA
updated_taxonomy_above80 <- updated_taxonomy %>%
  mutate(Domain = ifelse(Domain_confidence > 80, Domain, NA),
         Phylum = ifelse(Phylum_confidence > 80, Phylum, NA),
         Class = ifelse(Class_confidence > 80, Class, NA),
         Order = ifelse(Order_confidence > 80, Order, NA),
         Family = ifelse(Family_confidence > 80, Family, NA),
         Genus = ifelse(Genus_confidence > 80, Genus, NA),
         Species = ifelse(Species_confidence > 80, Species, NA)
  )

unclassified_below_class <-  updated_taxonomy_above80 %>%
  filter(!if_all(Domain:Species, ~is.na(.))) %>%
  filter((is.na(Order) & is.na(Family) & is.na(Genus) & is.na(Species))) %>%
  pull(asv_id)

#get old taxonomy info for ASVs that are all NA in updated taxonomy
all_na_in_new_tax <- updated_taxonomy_above80 %>% 
  filter(if_all(Domain:Species, ~is.na(.)) | asv_id %in% unclassified_below_class) %>% 
  mutate(across(contains("confidence"), ~NA)) %>% 
  select(-c(Domain:Species)) %>%
  left_join(taxonomy_tibble, by = join_by("asv_id" == "asv_names"))

#combine the old taxonomy with the updated version
combined_taxonomy <- updated_taxonomy_above80 %>% 
  filter(!(if_all(Domain:Species, ~is.na(.)) | asv_id %in% unclassified_below_class)) %>%
  full_join(all_na_in_new_tax)

#get list of genera that have multiple described taxonomies
multiple_classifications_list <- combined_taxonomy %>%
  select(Domain:Genus) %>%
  group_by(Genus) %>% 
  distinct() %>%
  summarise(n = n()) %>%
  filter(n > 1, !is.na(Genus)) %>%
  plyr::arrange(Genus) %>%
  pull(Genus)

#get the ncbi classifications for the genera with multiple described taxonomies

# ncbi_classifications_by_genus <- tax_name(multiple_classifications_list, db = "ncbi", 
#                                           get = c("Domain", "Phylum", "Class", "Order", "Family")) %>%
#   select(-c(db, Domain)) %>%
#   dplyr::rename("Genus" = "query")
# write_csv(ncbi_classifications_by_genus, "../intermediate_files/ncbi_classifications_for_overlaps.csv")
ncbi_classifications_by_genus <- read_csv("../intermediate_files/ncbi_classifications_for_overlaps.csv")

#most updated taxonomy, just phylum:genus
#each genus has only one described classification, is combined with old and new taxonomies
#doesnt contain NA for genus rows
nonoverlapping_taxonomy <- combined_taxonomy %>% 
  select(Phylum:Genus) %>% 
  distinct() %>% 
  filter(!is.na(Genus)) %>%
  filter(!Genus %in% multiple_classifications_list) %>%
  rbind(ncbi_classifications_by_genus) %>%
  distinct() %>%
  filter(!(Genus == "Oleiphilus" & Phylum == "Proteobacteria") & !(Genus == "Nannocystis" & is.na(Class)) 
         & !(Genus == "Aliiroseovarius" & Family == "Rhodobacteraceae") & !(Genus == "Porticoccus" & Phylum == "Proteobacteria"))

#our ASVs with the most up to date taxonomy, use for downstream purposes
full_taxonomy <- combined_taxonomy %>%
  select(-c(Phylum:Family)) %>%
  left_join(nonoverlapping_taxonomy, by = join_by(Genus)) %>%
  relocate(Phylum:Family, .after = Domain) %>%
  filter(!asv_id %in% c(combined_taxonomy %>% filter(is.na(Genus)) %>% pull(asv_id))) %>%
  rbind(combined_taxonomy %>% filter(is.na(Genus))) %>%
  arrange(parse_number(asv_id))
