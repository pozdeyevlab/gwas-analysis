# Add significant genes from FUSION and PrediXscan to lead variants
After [combining lead variants](https://github.com/pozdeyevlab/gwas-analysis/tree/main/combine-leads) from across all ancestries we ran FUSION and PrediXscan, in order to better understand the genes which were mapped by our lead variants. 

## Required Inputs
* Final outputs from [combining lead variants](https://github.com/pozdeyevlab/gwas-analysis/tree/main/combine-leads) 
* Ancestry stratified [FUSION](TODO ADD FUSION REPO) results
* Ancestry stratified [PrediXScan](TODO ADD Fusion Repo) results
* Homo_sapiend.GRCh38.113.gtf file downloaded from [ensembl - here](https://ftp.ensembl.org/pub/current_gtf/homo_sapiens/)

**There is no example data, there is no map file required, only a lead variant file, a FUSION file, and a PrediXscan file for each phenotype. Add or remove pheotypes by adding to the nested dictionary in the config. Every key in the inputs dictionary will become an output file. It is expacted that there is at least one FUSION or PrediXscan file. Please name as `fusion_ancestry` or `predixscan_ancestry`**

## Outputs
* For every phenotype summary file provided the same file will be created, with the addition of two columns. 'Gene.FUSION', and 'GENE.PrediXscan'. These columns will contain comma separated lists of significant genes whose eQTL variant is within 1Mb of the lead variant. 

## Environment Build & Running the Workflow
```{bash}
conda env create -f environment.yml
conda activate merge_twas
poetry install

# Check dag after filling in config.yaml
snamemake --cores 1 --configfile config.yaml --dry-run
snamemake --cores 1 --configfile config.yaml
```
