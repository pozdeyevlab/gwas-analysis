# Reactome and Kegg Enrichment

To try and understand functional pathways associated with out top-genes we used a simple enrichment analysis to find significant pathways. 

## Steps
1) For each phenotype provide a tsv with two columns:gene, pvals. 
2) Run [ReactomePA::enrichPathway()](https://www.rdocumentation.org/packages/ReactomePA/versions/1.16.2/topics/enrichPathway) for each phenotype. 
3) Save tables for Reactome & Kegg results as well as dot plots.
4) Map EntrezID's back to gene names for easier interpretation.

## Required Inputs
All that is required is a tsv with two columns. The first containing gene names and the second p-values. An example has been provided in `example_inputs`. 

To run this workflow, fill out the `config.yaml` file with your data.

## Outputs
This will create the following files in your specified output directory for every phenotype in your configfile
1) {output_dir}/{phenotype}_kegg.tsv = Significant kegg pathways (nominal p-value of 0.05)
2) {output_dir}/{phenotype}_kegg_with_gene_names.tsv = 
3) {output_dir}/{phenotype}_kegg.png = Dot plot for significant kegg pathways
4) {output_dir}/{phenotype}_entrez_map.tsv = Map between entrez id's and gene names
5) {output_dir}/{phenotype}_reactome.tsv = Significant kegg pathways (nominal p-value of 0.05)
6) {output_dir}/{phenotype}_reactome.png = Dot plot for significant reactome pathways

## Environment Build & Running the Workflow
```bash
conda env create -f environment.yml
conda activate snkmk_enrich
poetry install
snakemake --cores 2 --configfile config.yaml --use-conda --conda-frontend conda --dry-run
# Assuming no errors
snakemake --cores 2 --configfile config.yaml --use-conda --conda-frontend conda
```

## Note about the example data
The provided input data is entirely made up and is not expected to yield meaningful results.