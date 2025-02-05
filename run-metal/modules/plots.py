"""
This module makes manhattan & qq plots from META outputs. 
"""

from pathlib import Path
from typing import List, Union, Optional

import defopt
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import pyarrow
from qqman import qqman
from scipy import stats

# pylint: disable=R0914, R0913, R0915, C0121, line-too-long, unused-import


def plot(
    *,
    file_path: Path,
    manhattan_out: Path,
    phenotype: str,
    biobank_min: int,
    cases: Union[int, float],
    controls: Union[int, float],
    or_filter: float = None,
    af_filter: float = None,
) -> None:
    """
    :param file_path: Path to alignment results files
    :param manhattan_out: Write manhattan plots to this file
    :param phenotype: Disease
    :param biobank_min: Minimum number of biobanks
    :param cases: Total number of cases
    :param controls: Total number of controls
    :params or_filter: OR max
    :params af_filter: AF min
    """
    # Columns to read from each file
    columns_to_read = [
        "MarkerName",
        "Allele1",
        "Allele2",
        "P-value",
        "Freq1",
        "Direction",
        "Effect",
    ]

    # Read specific columns from all files into a single polars DF
    combined_df = _read_specific_columns(file_path, columns_to_read, biobank_min)

    # Filter for power analysis
    if (af_filter is not None) and (or_filter is not None):
        combined_df = combined_df.filter(
            pl.col("adjusted_metal_beta").exp() < or_filter
        ).filter(pl.col("adjusted_metal_af") > af_filter)

    # Manhattan & QQ plots for common, rare and mixed variants
    figure, axes = plt.subplots(nrows=2, ncols=1, figsize=(20, 20))
    rows = axes.shape[0]
    if not combined_df.is_empty():
        for row in list(range(0, rows - 1)):
            _make_plot(
                axes=axes,
                row=row,
                df=combined_df,
                pval_col="P-value",
                phenotype=phenotype,
                cases=cases,
                controls=controls,
            )

    # figure.tight_layout()
    plt.savefig(
        manhattan_out,
        format="png",
    )
    plt.clf()
    plt.close()


def _make_plot(
    *,
    axes: np.ndarray,
    row: int,
    df: pl.DataFrame,
    pval_col: str,
    phenotype: str,
    cases: int,
    controls: int,
):
    """
    Given the average minor allele frequency per variant the data frame is filtered, chisq is calculated, ambda is calculated and plots are created.

    Args:
        axes: Numpy array for matplotlib plots
        row: Row for which to place plot
        df: Polars data frame with meta anlyzed information
        pval_col: Column that contains the meta p-value
        pehenotype: Disease
    """
    df_filtered = (
        df.with_columns((pl.col("CHR").str.replace("chr", "")).alias("CHR"))
        .with_columns((pl.col("CHR").str.replace("X", "23")).alias("CHR"))
        .sort(pl.col("CHR").cast(int))
    )

    if not df_filtered.is_empty():
        # Must convert to pandas for compatability with qqman
        # df_filtered = df_filtered.sort('CHR', 'POS')
        pandas_df = df_filtered.to_pandas()
        pandas_df["CHR"] = pandas_df["CHR"].astype(int)
        pandas_df = pandas_df.sort_values("CHR")

        # Calculate chisquare
        pandas_df["CHISQ"] = stats.chi2.isf(pandas_df["P-value"], df=1)
        print(pandas_df)

        # Calculate lambda
        lambda_value = round(np.median(pandas_df["CHISQ"]) / 0.4549364231195724, 5)
        title = f"{' '.join(i.capitalize() for i in phenotype.split('_'))} Meta Analyzed GWAS Results\nƛ = {lambda_value}    N Variants: {pandas_df.shape[0]}\nCases: {cases}   Controls: {controls}"

        # Create plots
        if pandas_df.shape[0] > 0:
            qqman.manhattan(
                pandas_df,
                ax=axes[row],
                col_chr="CHR",
                col_bp="POS",
                col_p=pval_col,
                col_snp="MarkerName",
                label_size=25,
                title="",
                xtick_size=20,
                ytick_size=20,
                suggestiveline=False,
            )
            qqman.qqplot(
                pandas_df,
                ax=axes[row + 1],
                col_p=pval_col,
                label_size=25,
                title="",
                xtick_size=15,
                ytick_size=15,
            )
            axes[row].set_title(title, fontsize=30, fontweight="bold")
            axes[row].tick_params(axis="x", rotation=90)
            axes[row + 1].tick_params(axis="both", labelsize=12)
            axes[row + 1].grid(True, linestyle="--", linewidth=0.5)


def _read_specific_columns(
    file: Path, columns: List[str], biobank_min: int
) -> pl.DataFrame:
    """
    Reads a list of columns into a polars data frame, removes rows in which the p-value is null & or NaN. If the desired columns do not exist an empty df is returned. In order to split the data based on allele frequency the avg AF is calculated per variant, and named "AVG_AF"

    Args:
        file: The file to read in (assumed to be tab seperated)
        columns: List of columns to read
        biobank_min: Minimum number of biobanks
    """
    try:
        df = pl.read_csv(file, columns=columns, separator="\t")
    except pl.exceptions.ColumnNotFoundError:
        return pl.DataFrame({})
    if not df.is_empty():
        df = df.filter(pl.col("P-value").is_not_null())
        df = df.filter((pl.col("P-value") != np.NaN))

    # Extract chr and pos from I
    df = df.with_columns(
        pl.col("MarkerName")
        .str.replace("X", "23")
        .alias("MarkerName")).with_columns(
            pl.col("MarkerName")
            .str.split(":")
            .map_elements(lambda arr: arr[2].upper())
            .alias("ref_from_id"),
            pl.col("MarkerName")
            .str.split(":")
            .map_elements(lambda arr: arr[3].upper())
            .alias("alt_from_id"),
            pl.col("MarkerName")
            .str.split(":")
            .map_elements(lambda arr: arr[1])
            .alias("POS"),
            pl.col("MarkerName")
            .str.split(":")
            .map_elements(lambda arr: arr[0].replace("X", "23"))
            .alias("CHR")
            .cast(int),
        )
    

    df = (
        df.with_columns(pl.col("Allele1").str.to_uppercase().alias("Allele1"))
        .with_columns(pl.col("Allele2").str.to_uppercase().alias("Allele2"))
        .with_columns(
            pl.when(
                (pl.col("ref_from_id") == pl.col("Allele1").str.to_uppercase())
                & (pl.col("alt_from_id") == pl.col("Allele2").str.to_uppercase())
            )
            .then(True)
            .otherwise(False)
            .alias("reported_alleles_match_id")
        )
        .with_columns(
            pl.when(pl.col("reported_alleles_match_id"))
            .then(pl.col("Freq1"))
            .otherwise(1 - pl.col("Freq1"))
            .alias("adjusted_metal_eaf")
        )
        .with_columns(
            pl.when(pl.col("reported_alleles_match_id"))
            .then(pl.col("Effect"))
            .otherwise(-1 * pl.col("Effect"))
            .alias("adjusted_metal_beta")
        )
    )

    biobank_count = df.with_columns(pl.col("Direction").str.len_bytes().alias("temp"))[
        "temp"
    ].min()
    if biobank_min < biobank_count:
        df = df.filter((pl.col("Direction").str.count_matches(r"\+|\-")) >= biobank_min)
        return df
    else:
        df = df.filter((pl.col("Direction").str.count_matches(r"\+|\-")) >= 2)
        return df


if __name__ == "__main__":
    defopt.run(plot)
