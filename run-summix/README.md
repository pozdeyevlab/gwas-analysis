# SUMMIX For Custom LDSC-References

## This repository holds the code used to run summix on all gnomAD aligned summary stats.

## Required Inputs:
* map file of meta data used in `gwasqc`
* gnomad aligned files -- if using in conjunction with `gwasqc` these will be the files ending with *_aligned_results.tsv

## Running the workflow:

The commands below will build the required environment, activate it, and create the DAG. If the DAG is not output to the console check that all paths are correct in `config.yaml`. 

This repository comes with test data, use this an example for your own input data. 

```bash
conda env create -f environment.yml
conda activate summix
snakemake --cores 1 --configfile config.yaml --snakefile Snakefile --use-conda --conda-frontend conda --dry-run
```
