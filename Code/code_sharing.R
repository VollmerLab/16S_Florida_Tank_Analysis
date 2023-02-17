#### Filtering and normalizing data
otu_tmm <- microbiome %>%
  phyloseq_filter_prevalence(prev.trh = 0.1) %>%
  # phyloseq_transform_css() %>%
  otu_table() %>% 
  t %>% #need this if no phyloseq_transform_css
  as.data.frame %>%
  as.matrix %>%
  DGEList %>%
  edgeR::calcNormFactors(method = 'TMM') #change method to TMM



filter_transcripts <- function(data, min_cpm, prop_samples){
  keep <- rowMeans(cpm(data, log = TRUE) > min_cpm) >= prop_samples
  message('Genes Removed by Filter: ', scales::comma(table(keep)[1]))
  message('Genes Kept by Filter: ', scales::comma(table(keep)[2]))
  # message('Percet of Genes of Interest Kept after Filter: ', 
  #         scales::percent(sum(rownames(data[keep, ]$counts) %in% genes_of_interest$gene_id) / length(genes_of_interest$gene_id)))
  # message('Percet of Genes in Linkage Region Kept after Filter: ', 
  #         scales::percent(sum(rownames(data[keep, ]$counts) %in% linkage_genes$gene_id) / length(linkage_genes$gene_id)))
  
  data[keep, keep.lib.sizes = FALSE]
}

gene_normalize_factors <- pivot_wider(raw_gene_expression, names_from = sequence_id,
                                      values_from = gene_count) %>%
  # filter(str_detect(gene_id, '^__', negate = TRUE)) %>%
  column_to_rownames('gene_id') %>%
  as.matrix() %>% 
  DGEList(remove.zeros = TRUE) %>%
  filter_transcripts(0.1, 0.1) %>%
  calcNormFactors(method = 'TMMws')