# Run SUMMIX

## This repository holds the code used to run summix on all gnomAD aligned summary stats, and meta-analyzed mixed ancestry results for generation of custom ldsc results. 

## Required Inputs:
* map file of meta data used in `gwasqc`
* gnomad aligned files -- if using in conjunction with `gwasqc` these will be the files ending with *_aligned_results.tsv

## Running the workflow:

The commands below will build the required environment, activate it, and create the DAG. If the DAG is not output to the console check that all paths are correct in `config.yaml`. 

This repository comes with test data, use this an example for your own input data. 

```bash
conda env create -f environment.yml
conda activate summix_snkmk
snakemake --cores 1 --configfile config.yaml --snakefile Snakefile --use-conda --conda-frontend conda --dry-run

# The output TSV can be found in summix_out/summix_results/
```
