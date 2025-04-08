# PRSCS Implementation

This repository holds the code to needed to run [PRSCS](https://github.com/getian107/PRScs) on the meta analysis results. 

## Required Inputs
1) METAL results file (*_1.tbl) from leave CCPM out runs
2) LD reference file (ancestry specific) we used the EUR pannel from the [1K Genomes Project](https://www.dropbox.com/s/mt6var0z96vb6fv/ldblk_1kg_eur.tar.gz?dl=0)
3) The number of GWAS (cases+controls) per phenotype
4) RSID map file, made from dbSNP 
6) A bim file (validatin set) -- list of variants in your test data

The config.yaml requires a map file (for formatting see example_inputs/map.tsv). Please fill out config file before proceeding. Example data is not provided due to size of references. 


## Outputs
1) For ever row in the map file you will see a final 'lifted' results file in the output directory. These can be fed into plink2 scores to calculate PRS on your desired population. [plink2](https://www.cog-genomics.org/plink/2.0/) we used PLINK v2.00a6LM AVX2 Intel (9 Jun 2024). 

## Environment Build & Running the Workflow
```{bash}
conda env create -r environment.yml
conda activate run_prscs
poetry install

# Clone PRScs
git clone https://github.com/getian107/PRScs.git

# Check dag
snakemake --cores 1 --configfile config.yaml --dry-run
snakemake --cores 1 --configfile config.yaml
```
