"""Make QQ & Manhattan Plots"""

from pathlib import Path

import defopt
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import pyarrow
from qqman import qqman
from scipy import stats

# pylint: disable=C0301
# pylint: disable=R0914
# pylint: disable=R0913
# pylint: disable=R0915
# pylint: disable=C0121
# pylint: disable=W0611


def plot(
    *,
    file_path: str,
    pval_col: str,
    manhattan_out: Path,
) -> None:
    """
    :param file_path: Path to alignment results files
    :param pval_col: Column with p-value
    :param manhattan_out: Write manhattan plots to this file
    """
    # Columns to read from each file
    columns_to_read = [
        "Aligned_AF",
        pval_col,
        "STUDY_ID",
        "outlier_stdev",
        "CHR_gnomad",
        "POS_gnomad",
        "REF_gnomad",
        "ALT_gnomad",
        "GNOMAD_AN_Flag",
    ]

    # Read specific columns from all files into a single polars DF
    print(f"File used: {file_path}")
    combined_df = _read_specific_columns(file_path, columns_to_read, pval_col)

    # Manhattan plots for common and rare variants & qqplots
    figure, axes = plt.subplots(nrows=2, ncols=3, figsize=(50, 20))
    rows = axes.shape[0]
    columns = axes.shape[1]

    for row in list(range(0, rows - 1)):
        for column in list(range(0, columns)):
            if column == 0:
                _make_plot(
                    axes=axes,
                    col=column,
                    row=row,
                    df=combined_df,
                    pval_col=pval_col,
                    maf_cutoff="rare",
                )
            if column == 1:
                _make_plot(
                    axes=axes,
                    col=column,
                    row=row,
                    df=combined_df,
                    pval_col=pval_col,
                    maf_cutoff="common",
                )
            if column == 2:
                _make_plot(
                    axes=axes,
                    col=column,
                    row=row,
                    df=combined_df,
                    pval_col=pval_col,
                    maf_cutoff="both",
                )

    figure.tight_layout()
    plt.savefig(
        manhattan_out,
        format="png",
        dpi=600,
    )
    plt.clf()
    plt.close()


def _make_plot(
    *,
    axes: np.ndarray,
    col: int,
    row: int,
    df: pl.DataFrame,
    pval_col: str,
    maf_cutoff: str,
) -> None:
    df = df.with_columns(
        (pl.col("CHR_gnomad").str.replace("chr", "")).alias("CHR_gnomad")
    )
    df = df.with_columns((pl.col("POS_gnomad").cast(int)).alias("POS_gnomad"))
    df = df.with_columns((pl.col("STUDY_ID").str.replace("X", "23")).alias("STUDY_ID"))
    df = df.with_columns(
        (pl.col("CHR_gnomad").str.replace("X", "23")).cast(int).alias("CHR_gnomad")
    ).sort("CHR_gnomad", "POS_gnomad")

    if maf_cutoff.lower() == "common":
        df_filtered = df.filter(pl.col("Aligned_AF") > 0.05)
    if maf_cutoff.lower() == "rare":
        df_filtered = df.filter(pl.col("Aligned_AF") <= 0.05)
    if maf_cutoff.lower() == "both":
        df_filtered = df

    # Must convert to pandas for compatability with qqman
    pandas_df = df_filtered.to_pandas()

    # Calculate chisquare
    pandas_df["CHISQ"] = stats.chi2.isf(pandas_df[pval_col], df=1)

    # Calculate lambda
    lambda_value = round(np.median(pandas_df["CHISQ"]) / 0.4549364231195724, 5)
    common_threshold = 0.05

    if maf_cutoff.lower() == "common":
        title = f"{maf_cutoff.capitalize()} Variants (MAF > {common_threshold})\nƛ = {lambda_value}\nN: {pandas_df.shape[0]}"
    if maf_cutoff.lower() == "rare":
        title = f"{maf_cutoff.capitalize()} Variants (MAF ≤ {common_threshold})\nƛ = {lambda_value}\nN: {pandas_df.shape[0]}"
    if maf_cutoff.lower() == "both":
        title = (
            f"{maf_cutoff.capitalize()}\nƛ = {lambda_value}\nN: {pandas_df.shape[0]}"
        )

    # Create plots
    if df.shape[0] > 0:
        qqman.manhattan(
            pandas_df,
            ax=axes[row, col],
            col_chr="CHR_gnomad",
            col_bp="POS_gnomad",
            col_p=pval_col,
            col_snp="STUDY_ID",
            title=title,
        )
        qqman.qqplot(
            pandas_df,
            ax=axes[row + 1, col],
            title=title,
            col_p=pval_col,
        )


def _read_specific_columns(file, columns, pval_col):
    columns = list(set(columns))

    df = pl.read_csv(file, columns=columns, separator="\t", dtypes={"CHR_gnomad": str})
    print(df)

    df = df.with_columns(
        pl.col("CHR_gnomad").cast(str).str.replace("chr", "").alias("CHR_gnomad")
    )
    df = df.filter(pl.col(pval_col).is_not_null())
    df = df.filter(pl.col("GNOMAD_AN_Flag") == 0).filter(
        pl.col("outlier_stdev") == "No"
    )

    if pval_col == "LOG10P":
        df = df.with_columns((10 ** (-1 * pl.col(pval_col))).alias(pval_col))

    return df


if __name__ == "__main__":
    defopt.run(plot)
