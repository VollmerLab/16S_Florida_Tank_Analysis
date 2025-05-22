setwd("/Users/Emily/Desktop/GitHub/16S_Florida_Tank_Analysis/Code")

#### Libraries ####

library(tidyverse)
library(phyloseq)
library(taxize)
library(microbiome) #BiocManager::install("microbiome")

#### Read In Data ####
aggregation_level <- 'none' #set to either taxonomic level or none for ASV level

microbiome_data <- read_rds("../intermediate_files/preprocess_microbiome.rds")

#set taxa names and aggregate based on aggregation level
if(aggregation_level != 'none'){
  microbiome_data <- aggregate_taxa(microbiome_data, aggregation_level)
  taxa_names(microbiome_data) <- str_replace_all(taxa_names(microbiome_data), ' |-', '_')
} else {
  sequences <- taxa_names(microbiome_data)
  taxa_names(microbiome_data) <- str_c('ASV', 1:length(taxa_names(microbiome_data)), sep = '_')
  names(sequences) <- taxa_names(microbiome_data)
}


#### Update Taxonomy where Possible ####

#Taxonomy was assigned to each ASV first using a Bayesian taxonomic classifier based on the NCBI microbial 
#16S rRNA gene database and classified to the lowest taxonomic level possible with greater than 80% classification
#confidence (Gao et al., 2017). ASVs lacking confident identifications below the rank of class were reclassified 
#against the Silva database (Quast et al., 2012).

#This dual approach to taxonomic classification, while it provided more updated taxonomic classifications, also led
#to some mismatches at various taxonomic levels, such as most ASVs of a particular species having an updated family
#but one ASV of that species being listed under a different family, with all other taxonomic levels being identical
#otherwise.  This could result in those ASVs not aggregating together when examining genus or family level associations
#despite their matching identities.  This script combines the updated taxonomy (using a Bayesian taxonomic classifier)
#with the older taxonomy (Silva database) and then aims to fix any mismatches in single taxonomic levels by pulling
#the correct taxonomy from NCBI.


#the taxonomy in the phyloseq object is currently the more outdated version from the Silva database
#taxonomy_tibble pulls these Silva classifications to substitute in where the Bayesian classifier has less
#than 80% confidence
taxonomy_tibble <- tax_table(microbiome_data) %>% 
  as.data.frame %>%
  as_tibble(rownames = "asv_names")

#read in updated taxonomy from Bayesian classifier (script run in Bash)
updated_taxonomy <- read_csv('../intermediate_files/updated_taxonomy.csv')

#replace everything in updated taxonomy with less than 80% confidence with NA
updated_taxonomy_above80 <- updated_taxonomy %>%
  mutate(Domain = ifelse(Domain_confidence > 80, Domain, NA),
         Phylum = ifelse(Phylum_confidence > 80, Phylum, NA),
         Class = ifelse(Class_confidence > 80, Class, NA),
         Order = ifelse(Order_confidence > 80, Order, NA),
         Family = ifelse(Family_confidence > 80, Family, NA),
         Genus = ifelse(Genus_confidence > 80, Genus, NA),
         Species = ifelse(Species_confidence > 80, Species, NA)
  )

#make list of ASVs that are not confidently classified beyond class
unclassified_below_class <-  updated_taxonomy_above80 %>%
  filter(!if_all(Domain:Species, ~is.na(.))) %>% #remove ASVs with no confident classifications
  filter((is.na(Order) & is.na(Family) & is.na(Genus) & is.na(Species))) %>% #find ASVs that are not 
  #confidently classified beyond class
  pull(asv_id)

#get Silva taxonomy info for ASVs that are all NA in updated taxonomy
all_na_in_new_tax <- updated_taxonomy_above80 %>% 
  filter(if_all(Domain:Species, ~is.na(.)) | asv_id %in% unclassified_below_class) %>% #filter for ASVs not adequately classified
  mutate(across(contains("confidence"), ~NA)) %>% #set their confidence scores to NA
  select(-c(Domain:Species)) %>%
  left_join(taxonomy_tibble, by = join_by("asv_id" == "asv_names")) #add in original Silva classifications for these ASVs

#combine the Silva taxonomy with the updated Bayesian taxonomy
combined_taxonomy <- updated_taxonomy_above80 %>% 
  filter(!(if_all(Domain:Species, ~is.na(.)) | asv_id %in% unclassified_below_class)) %>% #remove inadequately classified ASVs
  full_join(all_na_in_new_tax) #add in the Silva version of these ASVs

#### Fix Misclassifications at Genus ####

#get list of genera that have multiple described taxonomies
multiple_classifications_list <- combined_taxonomy %>%
  select(Domain:Genus) %>%
  group_by(Genus) %>% 
  distinct() %>%
  summarise(n = n()) %>% #how many distinct taxonomies are listed for a single genus
  filter(n > 1, !is.na(Genus)) %>% #filter for genera with more than 1 taxonomy (i.e. with misclassifications)
  plyr::arrange(Genus) %>%
  pull(Genus)

#get the ncbi classifications for the genera with multiple described taxonomies

#running this code to get the classifications from NCBI will fail sometimes if your connection gets interrupted
#I just kept re-running it until it successfully made it through the whole list

# ncbi_classifications_by_genus <- tax_name(multiple_classifications_list, db = "ncbi",
#                                          get = c("Domain", "Phylum", "Class", "Order", "Family")) %>%
#    select(-c(db, Domain)) %>%
#    dplyr::rename("Genus" = "query")
# 
# write_csv(ncbi_classifications_by_genus, "../intermediate_files/ncbi_classifications_by_genus.csv")

ncbi_classifications_by_genus <- read_csv("../intermediate_files/ncbi_classifications_by_genus.csv")


#pull only the phylum:genus to match the output from NCBI and combine it to remove multiple classifications at genus level
#creates set of taxonomy from Phylum:Genus
nonoverlapping_taxonomy <- combined_taxonomy %>% 
  select(Phylum:Genus) %>% 
  distinct() %>% 
  filter(!is.na(Genus)) %>%
  filter(!Genus %in% multiple_classifications_list) %>% #filter for only genera that don't have multiple classifications
  rbind(ncbi_classifications_by_genus) %>% #add in correct classifications pulled from NCBI for the multiple classified genera
  distinct()

#add the nonoverlapping_taxonomy (at Genus level) back to ASV list
full_taxonomy <- combined_taxonomy %>%
  select(-c(Phylum:Family)) %>%
  left_join(nonoverlapping_taxonomy, by = join_by(Genus)) %>%
  relocate(Phylum:Family, .after = Domain) %>%
  #any ASVs with NAs for Genus no longer have taxonomic descriptions, so remove all ASVs with NA for genus
  filter(!asv_id %in% c(combined_taxonomy %>% filter(is.na(Genus)) %>% pull(asv_id))) %>%
  rbind(combined_taxonomy %>% filter(is.na(Genus))) %>% #add back in the taxonomy for ASVs that were NA for genus
  arrange(parse_number(asv_id)) #put in ASV number order again

#make a version of the microbiome data phyloseq to update the taxonomy in
updated_microbiome_data <- microbiome_data

#format full_taxonomy as a matrix and set it as the taxonomy info of the phyloseq
tax_table(updated_microbiome_data) <- full_taxonomy %>% 
  select(-contains("confidence")) %>%
  arrange(parse_number(asv_id)) %>%
  column_to_rownames("asv_id") %>%
  as.matrix()

#### Fix Misclassifications at Family ####

#aggregate by family
aggregated_microbiome_data_fam <- aggregate_taxa(updated_microbiome_data, "Family")

#get list of families that have multiple classifications
nonunique_taxonomy_families <- tax_table(aggregated_microbiome_data_fam) %>%
  as.data.frame() %>%
  as_tibble() %>%
  filter(str_detect(unique, "_")) %>% #Families with multiple classifications have the full taxonomic list written out
  #with underscores in the unique column whereas unique familes do not. This gives us families with multiple classifications only
  filter(unique != "Unknown Family_4") %>% #remove this family, which has an underscore in its name but does not have multiple classifications
  pull(Family) %>%
  unique()

#get accurate classifications from Domain:Order from NCBI for theses families
# ncbi_classifications_by_family <- tax_name(nonunique_taxonomy_families, db = "ncbi",
#                                            get = c("Domain", "Phylum", "Class", "Order")) %>%
#   select(-c(db, Domain)) %>%
#   dplyr::rename("Family" = "query")
# 
# write_csv(ncbi_classifications_by_family, "../intermediate_files/ncbi_classifications_by_family.csv")

ncbi_classifications_by_family <- read_csv("../intermediate_files/ncbi_classifications_by_family.csv")

#correct classifications for families that had multiple classifications
new_tax_just_non_unique_fams <- tax_table(updated_microbiome_data) %>% 
  as.data.frame() %>%
  as_tibble(rownames = "asv_names") %>%
  filter(Family %in% nonunique_taxonomy_families) %>% #filter for families with multiple classifications
  select(-c(Phylum, Class, Order)) %>%
  left_join(ncbi_classifications_by_family, by = join_by("Family")) %>% #add in the Phylum:Order from NCBI
  relocate(c(Phylum, Class, Order), .after = Domain)

#update the ASV list with the correct family level classifications
unique_family_level_taxonomy <- tax_table(updated_microbiome_data) %>% 
  as.data.frame() %>%
  as_tibble(rownames = "asv_names") %>%
  filter(!Family %in% nonunique_taxonomy_families) %>% #remove families with multiple classifications
  rbind(new_tax_just_non_unique_fams) %>% #add in correct classifications for families that had multiple classifications
  arrange(parse_number(asv_names))
  
#format family level taxonomy as a matrix and set it as the taxonomy info of the phyloseq
tax_table(updated_microbiome_data) <- unique_family_level_taxonomy %>% 
  select(-contains("confidence")) %>%
  arrange(parse_number(asv_names)) %>%
  column_to_rownames("asv_names") %>%
  as.matrix()

#### Fix Misclassifications at Order ####

#aggregate by order
aggregated_microbiome_data_ord <- aggregate_taxa(updated_microbiome_data, "Order")

#get list of orders that have multiple classifications
nonunique_taxonomy_orders <- tax_table(aggregated_microbiome_data_ord) %>%
  as.data.frame() %>%
  as_tibble() %>%
  filter(str_detect(unique, "_")) %>% #filter for orders with multiple classifications
  pull(Order) %>%
  unique()

#get accurate classifications from Phylum:Class from NCBI for theses orders
# ncbi_classifications_by_order <- tax_name(nonunique_taxonomy_orders, db = "ncbi",
#                                            get = c("Phylum", "Class")) %>%
#   select(-c(db)) %>%
#   dplyr::rename("Order" = "query")

# write_csv(ncbi_classifications_by_order, "../intermediate_files/ncbi_classifications_by_order.csv")

ncbi_classifications_by_order <- read_csv("../intermediate_files/ncbi_classifications_by_order.csv")

#correct classifications for orders that had multiple classifications
new_tax_just_non_unique_orders <- tax_table(updated_microbiome_data) %>% 
  as.data.frame() %>%
  as_tibble(rownames = "asv_names") %>%
  filter(Order %in% nonunique_taxonomy_orders) %>% #filter for orders with multiple classifications
  select(-c(Phylum, Class)) %>%
  left_join(ncbi_classifications_by_order, by = join_by("Order")) %>% #add in the Phylum:Class from NCBI
  relocate(c(Phylum, Class), .after = Domain)

#update the ASV list with the correct order level classifications
unique_order_level_taxonomy <- tax_table(updated_microbiome_data) %>% 
  as.data.frame() %>%
  as_tibble(rownames = "asv_names") %>%
  filter(!Order %in% nonunique_taxonomy_orders) %>% #remove orders with multiple classifications
  rbind(new_tax_just_non_unique_orders) %>% #add in correct classifications for orders that had multiple classifications
  arrange(parse_number(asv_names))

#format order level taxonomy as a matrix and set it as the taxonomy info of the phyloseq
tax_table(updated_microbiome_data) <- unique_order_level_taxonomy %>% 
  select(-contains("confidence")) %>%
  arrange(parse_number(asv_names)) %>%
  column_to_rownames("asv_names") %>%
  as.matrix()

#### Fix Misclassifications at Class ####

#aggregate by class
aggregated_microbiome_data_class <- aggregate_taxa(updated_microbiome_data, "Class")

#get list of classes that have multiple classifications
nonunique_taxonomy_classes <- tax_table(aggregated_microbiome_data_class) %>%
  as.data.frame() %>%
  as_tibble() %>%
  filter(str_detect(unique, "_")) %>% #filter for classes with multiple classifications
  pull(Class) %>%
  unique()

#get accurate classifications for Phylum from NCBI for theses classes
# ncbi_classifications_by_class <- tax_name(nonunique_taxonomy_classes, db = "ncbi",
#                                            get = c("Phylum")) %>%
#   select(-c(db)) %>%
#   dplyr::rename("Class" = "query")
# 
# write_csv(ncbi_classifications_by_class, "../intermediate_files/ncbi_classifications_by_class.csv")

ncbi_classifications_by_class <- read_csv("../intermediate_files/ncbi_classifications_by_class.csv")

#correct classifications for classes that had multiple classifications
new_tax_just_non_unique_classes <- tax_table(updated_microbiome_data) %>% 
  as.data.frame() %>%
  as_tibble(rownames = "asv_names") %>%
  filter(Class %in% nonunique_taxonomy_classes) %>% #filter for classes with multiple classifications
  select(-c(Phylum)) %>%
  left_join(ncbi_classifications_by_class, by = join_by("Class")) %>% #add in the Phylum from NCBI
  relocate(c(Phylum), .after = Domain)

#update the ASV list with the correct class level classifications
unique_class_level_taxonomy <- tax_table(updated_microbiome_data) %>% 
  as.data.frame() %>%
  as_tibble(rownames = "asv_names") %>%
  filter(!Class %in% nonunique_taxonomy_classes) %>% #remove classes with multiple classifications
  rbind(new_tax_just_non_unique_classes) %>% #add in correct classifications for classes that had multiple classifications
  arrange(parse_number(asv_names)) %>%
  mutate(Class = ifelse(is.na(Class) & Order %in% c("Polyangiales", "Haliangiales", "Nannocystales"), "Polyangia", Class))
  #manually set the Class for these 3 orders where NA, as class is their only NA taxonomic rank

#format class level taxonomy as a matrix and set it as the taxonomy info of the phyloseq
tax_table(updated_microbiome_data) <- unique_class_level_taxonomy %>% 
  select(-contains("confidence")) %>%
  arrange(parse_number(asv_names)) %>%
  column_to_rownames("asv_names") %>%
  as.matrix()

#### Write Updated Taxonomy RDS ####

write_rds(updated_microbiome_data, "../intermediate_files/updated_microbiome_data.rds")

#total number of each taxa present in the dataset
tax_table(updated_microbiome_data) %>% 
  as.data.frame() %>% 
  as_tibble() %>%
  select(Genus) %>% #change taxa level here
  filter(!is.na(.)) %>%
  distinct() %>%
  nrow()
