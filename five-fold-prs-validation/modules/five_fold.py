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
from typing import List

import defopt
import sys
import polars as pl
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score


def five_fold(
    *,
    prs_files: List[str],
    prs_names: List[str],
    pheno_files: List[str],
    pheno_names: List[str],
    output: str,
    covar_file: str,
    ancestry: str,
    delong_out: str,
) -> None:
    """
    :param prs_files:
    :param prs_names:
    :param pheno_files:
    :param pheno_names:
    :param output:
    :param delong_out:
    """
    # get list of prs's
    # Make matrix for comparing prs & predicting
    names_matrix = [(x, y) for x in prs_names for y in pheno_names]

    files_matrix = [(x, y) for x in prs_files for y in pheno_files]

    dfs = []
    with open(output, mode="w+") as f:
        f.write(f"prs\tpredicting\tcovariates\tmean_auc\tmean_auc_se\n")
        for m in list(zip(names_matrix, files_matrix)):
            prs = m[0][0]
            pheno = m[0][1]
            prs_file = m[1][0]
            pheno_file = m[1][1]
            weight = pl.read_csv(prs_file, separator="\t")

            pheno = pl.read_csv(pheno_file, separator="\t", null_values="NA")

            covariate_df = pl.read_csv(covar_file, separator="\t")
            if 'PRS' in covariate_df.columns:
                covariate_df = covariate_df.drop("PRS")
            if ancestry != "mixed":
                covariate_df = covariate_df.filter(
                    pl.col("similarity_group").str.contains(ancestry.upper())
                )

            # Merge data
            merged = weight.join(pheno, on="IID", how="inner")

            hardcall = merged.columns[-1]

            # Add covars
            merged = merged.join(covariate_df, on="IID", how="inner")
            print(merged.columns)
            
            dfs.append(
                merged.rename({"PRS": f"{prs}_prs"}).select(
                    "IID", f"{prs}_prs", hardcall
                )
            )
            merged = merged.filter(~(pl.col(hardcall).is_null())).to_pandas()

            covariates = [
                [
                    "Genotyping_Batch",
                    "PC1",
                    "PC2",
                    "PC3",
                    "PC4",
                    "PC5",
                    "PC6",
                    "PC7",
                    "PC8",
                    "PC9",
                    "PC10",
                ],
                [
                    "Genotyping_Batch",
                    "PC1",
                    "PC2",
                    "PC3",
                    "PC4",
                    "PC5",
                    "PC6",
                    "PC7",
                    "PC8",
                    "PC9",
                    "PC10",
                    "age",
                    "sex",
                ],
                [
                    "Genotyping_Batch",
                    "PC1",
                    "PC2",
                    "PC3",
                    "PC4",
                    "PC5",
                    "PC6",
                    "PC7",
                    "PC8",
                    "PC9",
                    "PC10",
                    "age",
                    "sex",
                    "PRS",
                ],
                [
                    "Genotyping_Batch",
                    "PC1",
                    "PC2",
                    "PC3",
                    "PC4",
                    "PC5",
                    "PC6",
                    "PC7",
                    "PC8",
                    "PC9",
                    "PC10",
                    "PRS",
                ],
                ["age", "sex"],
                ["PRS"],
            ]

            names = [
                ["gt_pc"],
                ["gt_pc_age_sex"],
                ["gt_pc_age_sex_prs"],
                ["gt_pc_prs"],
                ["age_sex_only"],
                ["prs_only"],
            ]

            for covariate in list(zip(covariates, names)):
                name = covariate[1]

                y = merged[hardcall].to_numpy()
                X = merged[covariate[0]].to_numpy()
                # Initialize the logistic regression model
                logreg = LogisticRegression(max_iter=1000)
                # Perform five-fold cross-validation

                auc_scores = cross_val_score(logreg, X, y, cv=5, scoring="roc_auc")
                # Mean nd stder
                mean_auc = auc_scores.mean()
                std_auc = auc_scores.std()

                f.write(f"{prs}\t{hardcall}\t{name[0]}\t{mean_auc}\t{std_auc}\n")
    final = dfs[0]
    columns = final.columns
    for df in dfs:
        for c in df.columns:
            if c not in columns:
                columns.append(c)
                final = final.join(df.select("IID", c), on="IID", how="left")

    final.write_csv(delong_out, separator="\t")


if __name__ == "__main__":
    defopt.run(five_fold)
