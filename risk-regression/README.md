# Regression Analysis of Risk Strata and Thyroid Cancer PRS

## Overview
Patient charts from CCPM (Colorado Center for Personalized Medicine) were manually reviewed to create a spreadsheet of the risk strata listed below. To assess any association between PRS (thyroid cancer vs all & tc vs bng) we conducted a regression analysis and report all betas and p-values for age, sex, and prs. 

## Required Inputs
The input to this workflow is a tsv with all covariates [age,sex,PC's1-10,Genotyping_Batch,PRS] and risk strata. 

In order to protect patient confidentiality we cannot share the input data.

## Outputs
For each PRS, the following will be generated:
1) A table with p-values, and beta values for all covariates for all risk factors
2) Formatted results with only sex,age,and prs

## Environment Build & Running the Workflow
```{bash}
conda env create -f environment.yml
conda activate risk_regression
poetry install
Rscript run_regression.R --input prs_and_risk.tsv --output test_thyroid_cancer.tsv --prs thyroid_cancer --invnorm FALSE
Rscript run_regression.R --input prs_and_risk.tsv --output test_tc_vs_bng.tsv --prs tc_vs_bng --invnorm FALSE
python format.py --input-files test_thyroid_cancer.tsv test_tc_vs_bng.tsv --output regression_results.tsv 
```