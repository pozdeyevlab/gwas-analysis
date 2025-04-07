"""
Load data with PRS and covariates
Load phenotype definitions
Do 5-fold logistic regression cross-validation on the following models
    1) GT + PCs
    2) GT + PCs + Age + Sex
    3) GT + PCs + PRS
    4) GT + PCs + Age + Sex + PRS
Mixed ancestry and EUR
"""
import polars as pl 
import numpy as np
import sys
import pandas as pd
from pathlib import Path
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score

# get list of prs's
# Make matrix for comparing prs & predicting
list1 = ['thyroid_cancer','thyroid_cancer_GBMI_I', 'benign_nodular_goiter', 'tc_vs_bng', 'graves_disease','hypothyroidism', 'lymphocytic_thyroiditis']

list2 = ['thyroid_cancer', 	'benign_nodular_goiter', 'tc_vs_bng_all', 'tc_vs_bng_path', 'graves_disease', 'hypothyroidism', 'lymphocytic_thyroiditis']

matrix = [(x, y) for x in list1 for y in list2]

with open('final_manual_mixed_five_fold_out.tsv', mode = 'w+') as f:
    f.write(f'prs\tpredicting\tcovariates\tmean_auc\tmean_auc_se\n')
    for m in matrix:
        prs = m[0]
        pheno = m[1]
        #df = pl.read_csv(f'OneDrive_1_3-19-2025/{prs}.prs', separator=',')
        df = pl.read_csv(f'OneDrive_1_3-19-2025/{prs}.prs', separator=',')#.filter(pl.col('similarity_group').str.contains('EUR'))
        pheno = pl.read_csv(f'phenos/{pheno}.pheno', separator='\t', null_values='NA')
        if prs == 'tc_vs_bng':
            df = df.rename({"combined":"PRS"})
    
        # Merge data
        merged = df.join(pheno, on = 'IID', how = 'inner')
        hardcall = merged.columns[-1]
        merged = merged.filter(~(pl.col(hardcall).is_null())).to_pandas()
    
        covariates = [['Genotyping_Batch', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10'], ['Genotyping_Batch', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10', 'age', 'sex'], ['Genotyping_Batch', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10', 'age', 'sex', 'PRS'], ['Genotyping_Batch', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10', 'PRS'], ['age', 'sex'], ['PRS']]


        names = [['gt_pc'], ['gt_pc_age_sex'], ['gt_pc_age_sex_prs'], ['gt_pc_prs'], ['age_sex_only'], ['prs_only']]


        for covariate in list(zip(covariates,names)):
            name = covariate[1]
            print(name)
            y = merged[hardcall].to_numpy()
            X = merged[covariate[0]].to_numpy()
            # Initialize the logistic regression model
            logreg = LogisticRegression(max_iter = 1000)
            # Perform five-fold cross-validation
            print(X)
            auc_scores = cross_val_score(logreg, X, y, cv=5, scoring='roc_auc')
            # Mean nd stder
            mean_auc = auc_scores.mean()
            std_auc = auc_scores.std()
            print(mean_auc)
            f.write(f'{prs}\t{hardcall}\t{name[0]}\t{mean_auc}\t{std_auc}\n')




