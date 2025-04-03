# Combine Meta-Analysis Lead Variants
Helper workflow used to combine the lead variants from all ancestry statified meta-analysis, and make final summary files.

## Required Inputs
This code is highly specialized to take the final output files from [run-metal](https://github.com/pozdeyevlab/gwas-analysis/tree/main/run-metal). Based on the [config file in run-metal](https://github.com/pozdeyevlab/gwas-analysis/blob/main/run-metal/config.yaml) files will be found in {output_dir}/final_summary_with_genes/*_with_genes_and_regions.tsv. 

Example inputs are provided in `example_inputs`, the config file is set up to use these example inputs.

*There is a simple bash script `example_inputs/make_map.py` to help automate making the map file*

## Outputs
Counts -- {output_dir}/final_tables/{phenotype}_{sex}_count.tsv
Summary -- {output_dir}/final_tables/{phenotype}_{sex}_all_ancestries.tsv
Most Significant -- {output_dir}/final_tables/{phenotype}_{sex}_most_significant.tsv

## Environment Build & Running the Workflow
```{bash}
conda env create -f environment.yml
conda activate combine
poetry install

# After filling out config file test the DAG
snakemake --cores 1 --configfile config.yaml --dry-run

# Run code
snakemake --cores 1 --configfile config.yaml
```