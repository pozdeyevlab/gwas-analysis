# FUSION Implementation

We ran [FUSION](http://gusevlab.org/projects/fusion/) to better understand the genes mapped from our lead variants. We used precomputed [GTEx v8 All Samples model for Thyroid](https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/ALL/GTExv8.ALL.Thyroid.tar.gz) (n = 574). 

## Requiured Inputs
* The inputs to this analysis are munged sumstats which are readily available from the [run-metal](https://github.com/pozdeyevlab/gwas-analysis/tree/main/run-metal) pipeline. The necessary files will be found in `{metal_output_dir}/ldsc_munge/*_summstats.gz`
* Precomputed [GTEx model](https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/ALL/GTExv8.ALL.Thyroid.tar.gz)
* [FUSION software](http://gusevlab.org/projects/fusion/#installation)
* Homo_sapiend.GRCh38.113.gtf file downloaded from [ensembl - here](https://ftp.ensembl.org/pub/current_gtf/homo_sapiens/)
* EUR 1K Genomes LD Reference Pannel


## Environment Build & Running the Workflow
```{bash}
conda env create -f environment.yml
conda activate merge_twas
poetry install

# Check dag after filling in config.yaml
snamemake --cores 1 --configfile config.yaml --dry-run
snamemake --cores 1 --configfile config.yaml
```


