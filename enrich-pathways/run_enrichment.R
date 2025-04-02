library(optparse)
library(data.table)
library(dplyr)
library(tidyr)
library(glue)

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
}

if (!requireNamespace("ReactomePA", quietly = TRUE)) {
  BiocManager::install("ReactomePA")
}

if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  BiocManager::install("clusterProfiler")
}

library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)

# Define command-line arguments
option_list <- list(
  make_option(c("--input"), type="character", help="Directory with phenotype specific files of genes and p-values"),
  make_option(c("--output"), type="character", help="Output Directory")
)

# Parse command-line arguments
opt <- parse_args(OptionParser(option_list=option_list))

print(opt$input)
print(opt$output)


gwas_data = fread(glue(opt$input), sep='\t')
  
gene_symbols <- gwas_data$gene
  
# Convert gene symbols to Entrez Gene IDs
entrez_ids <- bitr(gene_symbols, 
  fromType = "SYMBOL", 
  toType = "ENTREZID", 
  OrgDb = org.Hs.eg.db)
  
# Print the conversion result
entrez_id <- merge(entrez_ids, gwas_data, by.x='SYMBOL', by.y='gene', how = 'inner')
gwas_data <- entrez_id
  
# Rank genes by p-value
gene_list <- gwas_data$'pvals'
names(gene_list) <- gwas_data$ENTREZID
gene_list <- sort(gene_list)
  
# Perform Reactome GSEA using the ReactomePA package
reactome_gsea <- enrichPathway(gene = names(gene_list), 
  pvalueCutoff = 0.05, 
  organism = "human", 
  readable = TRUE)
  
# Create a ranked list of genes (you can use p-value, effect size, etc.)
gene_list <- gwas_data$pvals
names(gene_list) <- entrez_ids$ENTREZID
print(names(gene_list))
gene_list <- sort(gene_list)
  
# Perform KEGG enrichment analysis
kegg_enrich <- enrichKEGG(gene = names(gene_list), 
                          organism = "hsa",
                          keyType = 'kegg',
                          pvalueCutoff = 0.05)

# Save results
write.table(as.data.frame(kegg_enrich), file=glue('{opt$output}_kegg.tsv'), sep='\t',quote=FALSE, row.names=FALSE)
write.table(as.data.frame(reactome_gsea), file=glue('{opt$output}_reactome.tsv'),sep='\t', quote=FALSE, row.names=FALSE)
write.table(entrez_id, file=glue('{opt$output}_entrez_map.tsv'), sep='\t', quote=FALSE,row.names=FALSE)
  
# Plot
png(glue('{opt$output}_kegg.png'), width = 480, height = 480) 
  print(dotplot(kegg_enrich))
dev.off() 
  
png(glue('{opt$output}_reactome.png'), width = 480, height = 480) 
  print(dotplot(reactome_gsea))
dev.off() 