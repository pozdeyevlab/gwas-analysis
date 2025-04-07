"""
Formats outputs from modules/five_fold.py into an easy to read matrix, for each covariavate combo
"""
from pathlib import Path

import defopt
import polars as pl


def format_five_fold(*, input_file: Path, output: str) -> None:
    """
    :param input_file:
    :param output:
    """
    df = pl.read_csv(input_file, separator="\t")
    print(df)
    covariates = df["covariates"].unique()
    for covar in covariates:
        tmp = df.filter(pl.col("covariates") == covar).unique()
        # Format values
        tmp = tmp.with_columns(
            (pl.col("mean_auc") + (1.96 * pl.col("mean_auc_se")))
            .round(4)
            .alias("upper_ci")
        )
        tmp = tmp.with_columns(
            (pl.col("mean_auc") - (1.96 * pl.col("mean_auc_se")))
            .round(4)
            .alias("lower_ci")
        )
        tmp = tmp.with_columns(
            (
                pl.col("mean_auc").round(4).cast(str)
                + " "
                + "("
                + pl.col("lower_ci").cast(str)
                + ","
                + pl.col("upper_ci").cast(str)
                + ")"
            ).alias("value")
        )
        pivot_df = tmp.pivot(values="value", index="prs", columns="predicting")
        pivot_df.write_csv(f"{output}_{covar}_five_fold_matrix.tsv", separator="\t")


if __name__ == "__main__":
    defopt.run(format_five_fold)
