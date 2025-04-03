"""
Format the results from regresison analysis
"""

from pathlib import Path
from typing import List

import defopt
import polars as pl

# pylint: disable=R0914, R0913, R0903, C0301


def format(*, input_files: List[str], output: Path) -> pl.DataFrame:
    """
    :param input_files: Input path to regresison outputs from different PRS tests
    :param output: Where to write results
    """
    dfs = []
    for file in input_files:
        tmp = pl.read_csv(file, null_values="NA", separator="\t").rename(
            {"Estimate": "beta", "Std..Error": "SE", "Pr...t..": "pval"}
        )
        dfs.append(tmp)
    df = pl.concat(dfs, how="vertical")

    df = df.with_columns(
        (
            pl.col("beta").map_elements("{:.4g}".format, return_dtype=pl.String).cast(str)
            + " ("
            + pl.col("SE").map_elements("{:.4g}".format, return_dtype=pl.String).cast(str)
            + ")"
            + " "
            + pl.col("pval").map_elements("{:.5g}".format, return_dtype=pl.String).cast(str)
        ).alias("test")
    )

    df = df.filter(pl.col("covariates").is_in(["PRS", "sex", "age"]))
    df = df.with_columns(
        (
            pl.when(pl.col("pval") > 0.05)
            .then(pl.col("test").cast(str) + "^")
            .otherwise(pl.col("test"))
        ).alias("test")
    )

    df = df.with_columns(
        (
            pl.when((pl.col("pval") < 0.05))
            .then(pl.col("test").cast(str) + "*")
            .otherwise(pl.col("test"))
        ).alias("test")
    )

    df = df.with_columns(
        (
            pl.when((pl.col("pval") < 0.01))
            .then(pl.col("test").cast(str) + "*")
            .otherwise(pl.col("test"))
        ).alias("test")
    )

    df = df.with_columns(
        (
            pl.when(pl.col("pval") < pl.col("bonferroni"))
            .then(pl.col("test").cast(str) + "*")
            .otherwise(pl.col("test"))
        ).alias("test")
    )

    grouped = df.select("risk", "prs", "covariates", "test").group_by(pl.col("risk"))
    dfs = []
    for group in grouped:
        dfs.append((group[1]))

    final = pl.concat(dfs, how="vertical")

    final.write_csv(output, separator="\t")


if __name__ == "__main__":
    defopt.run(format)
