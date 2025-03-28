# Run SUMMIX

## This repository holds the code used to run [Summix2](https://github.com/hendriau/Summix) on all gnomAD aligned summary stats, and meta-analyzed mixed ancestry results for generation of custom LDSC reference files. 

## Required Inputs:
* map file of meta data used in `gwasqc`
* gnomad aligned files -- if using in conjunction with `gwasqc` these will be the files ending with *_aligned_results.tsv

## Running the workflow:

The commands below will build the required environment, activate it, and create the DAG. If the DAG is not output to the console check that all paths are correct in `config.yaml`. 

This repository comes with test data, use this an example for your own input data. 

## Snakemake workflow for all aligned results
```bash
conda env create -f environment.yml
conda activate snkmk_summix
poetry install
snakemake --cores 1 --configfile config.yaml --snakefile Snakefile --use-conda --conda-frontend conda --dry-run
```

## Simple bash script for running summix on metal results
```bash
conda env create -f r_environment.yml
conda activate summix
bash run_summix_on_metal.sh
```
