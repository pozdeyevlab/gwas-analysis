# Covariate Adjusted LD Scores

## AoU Disclaimer
The input data required for this process is not publicly available and researchers who do not have "Controlled Tier Access' will need to go through the proper registration and training before attempting to complete the steps outlined below. Furthermore, this code is not meant to be run locally as it is specifically designed for the AoU workbench. 

## Required Inputs
1) The AoU Published ancestry predictions which can be found here 'gs://~/ancestry_preds.tsv' (for those with controlled tier access). Please have `ancestry_preds.tsv` in your working directory. 
2) Summix2 predictions for your cohort of interest. These can be collected using [run-summix](https://github.com/pozdeyevlab/gwas-analysis/tree/main/run-summix). These will need to be manually written into `modules/dsub_subset_and_map_bgens.py`.
3) Whole genome bgens and sample files. We used the outputs from [vds-filter](https://github.com/pozdeyevlab/gwas-analysis/tree/main/vds-filter)

## Ouptus

## Running the Code
Step 1 = extracts 5000 random genomes with population structure from Summix
Step 2 = PCA_for_LDSC.ipynb  performs QC of the 5000 genomes and calculated principal components
Step 3 = Finally, dsub_cov_ldsc_nikita.ipynb calculates covariate-adjusted LD scores


