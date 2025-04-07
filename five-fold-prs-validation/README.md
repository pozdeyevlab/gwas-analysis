# Five-Fold Cross Validation on PRS Performance
Custom workflow, used to create the final PRS summary tables. Specifically, comparing the PRS generated from lead variants vs the PRS from [PRSCS](https://github.com/getian107/PRScs).

## Required Inputs
*manual prs refers to prs generated from lead variants*
1) Per phenotype a manually calculated PRS file, with the following columns: IID, PRS
2) Per phenotype the scores files from PRSCS, with the following columns: IID, PRS
3) Phenotype definitions for each phenotype, with the following columns: IID, <disease>
4) Covariante file, workflow expects IID, age, sex, Genotyping_Batch, PCs1-10, and similarity_group

### Importand Note
Based on the phenotype dictionary and the provided weights dictionary in the config file a matrix of AUC's will be produced. All combinations will also have p-values generated. This file can be parsed to extract comparisons of interest.

*Randomly generated example data is provided in example_inputs. Please use this data to format your own inputs*

## Outputs

## Environment Build & Running the Workflow
```{bash}
conda env create -f environment.yml
conda activate five_fold
poetry install

# Check dag
snakemake --cores 1 --configfile config.yaml --use-conda --conda-frontend conda --dry-run
snakemake --cores 1 --configfile config.yaml --use-conda --conda-frontend conda
```
