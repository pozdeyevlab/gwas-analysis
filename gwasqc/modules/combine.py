"""Combine aligned results, find outliers, and add AN flags from gnomad"""

from pathlib import Path
from typing import List

import defopt
import mahalanobis
import matplotlib.pyplot as plt
import polars as pl

# pylint: disable=C0301 # line too long
# pylint: disable=R0914 # Too many local variables
# pylint: disable=R0913 # Too many arguments
# pylint: disable=R0915 # Too many statements
# pylint: disable=C0121 # Comparison


def combine(
    *,
    file_paths: List[Path],
    gnomad_flag_dir: Path,
    output_file: Path,
) -> None:
    """
    Helper for generating scatterplots for effect allele frequencies between reference and study variants

    :param file_paths: Path to alignment results files
    :param gnomad_flag_dir: Directory to blacklisted variant files
    :param output_file: Where plots will be written
    """
    # Read all files into a single polars DF and add gnomad glag
    print(f"Files used:\n{list(file_paths)}")

    # Add option to handle AoU missing frequencies
    combined_pl = _combine_files(file_paths, gnomad_flag_dir)
    missing_count = combined_pl.filter(pl.col("Aligned_AF").is_null()).shape[0]
    print(
        f"There are {missing_count}/{combined_pl.shape[0]} variants are missing AF from raw summary stat"
    )
    combined_pl = combined_pl.with_columns(
        pl.when(pl.col("Aligned_AF").is_null())
        .then(pl.col("AF_gnomad"))
        .otherwise(pl.col("Aligned_AF"))
        .alias("Aligned_AF")
    )

    print(combined_pl)
    # Calculate outliers with mahalobis distances
    outlier_pl = mahalanobis.calculate(
        aligned_pl=combined_pl.select(["STUDY_ID", "Aligned_AF", "AF_gnomad"])
    )

    outlier_pl = pl.concat([outlier_pl, combined_pl], how="align")
    print(outlier_pl.columns)
    print(
        f'\nMahalanobis Summary:\nTotal aligned varinats with Mahalanobis distance greater than three standard deviations from the mean: {outlier_pl.filter(outlier_stdev="Yes").shape[0]}/{outlier_pl.shape[0]}'
    )
    print(
        f'Total aligned varinats with Mahalanobis distance p-value less than 0.001: {outlier_pl.filter(outlier_pval="Yes").shape[0]}/{outlier_pl.shape[0]}\n'
    )

    # Write output
    outlier_pl.write_csv(output_file, include_header=True, separator="\t")


def _combine_files(files, gnomad_flag_dir):
    dfs = []
    for file in files:
        chrom = f"{file}".split("/")[-1].split("_")[0]
        blacklist_df: pl.DataFrame = _get_blacklist_variants(
            gnomad_flag_dir=gnomad_flag_dir, chrom=chrom
        )
        df = pl.read_csv(
            file, separator="\t", null_values="NA", dtypes={"CHR_gnomad": str}
        )
        # Add low AN flag to df
        df = _add_an_flag(study_df=df, gnomad_df=blacklist_df)
        dfs.append(df)
    # Align column before concat
    columns = []
    for df in dfs:
        df_cols = df.columns
        [columns.append(c) for c in df_cols]
    all_columns = set(columns)
    print(f"\nAll Columns: {all_columns}")
    fixed_dfs = []
    for df in dfs:
        df_col = df.columns
        cols_to_add = all_columns.difference(set(df_col))
        if len(cols_to_add) > 0:
            for col in list(cols_to_add):
                df = df.with_columns((pl.lit(None).cast(str)).alias(col))
        cols = df.columns
        for c in ["#chrom", "chromosome", "CHR", "CHROM"]:
            if c in cols:
                df = df.drop(c)
        fixed_dfs.append(df)
    concat_df = pl.concat(fixed_dfs, how="align")
    print(concat_df)
    return concat_df


def _get_blacklist_variants(gnomad_flag_dir: Path, chrom: str) -> pl.DataFrame:
    if chrom == 23:
        chrom = "X"
        gnomad_tsv = list(gnomad_flag_dir.glob(f"flagged_variants_*chr{chrom}.tsv"))[0]
    else:
        gnomad_tsv = list(gnomad_flag_dir.glob(f"flagged_variants_*chr{chrom}.tsv"))[0]
    df = pl.read_csv(gnomad_tsv, separator="\t")
    return df


def _add_an_flag(study_df: pl.DataFrame, gnomad_df: pl.DataFrame) -> pl.DataFrame:
    # Remove chr from chrom col if present
    gnomad_df = _make_id_column(
        new_column_name="GNOMAD_ID",
        polars_df=gnomad_df,
        CHR="CHR",
        POS="POS",
        REF="REF",
        ALT="ALT",
    )
    study_df = _make_id_column(
        new_column_name="STUDY_ID",
        polars_df=study_df,
        CHR="CHR_gnomad",
        POS="POS_gnomad",
        REF="REF_gnomad",
        ALT="ALT_gnomad",
    )

    study_df = study_df.with_columns(
        pl.when((pl.col("STUDY_ID").is_in(list(set(gnomad_df["GNOMAD_ID"])))))
        .then(1)
        .otherwise(0)
        .alias("GNOMAD_AN_Flag")
    )
    return study_df


def _make_id_column(
    *,
    new_column_name: str,
    polars_df: pl.DataFrame,
    CHR: str,
    POS: str,
    REF: str,
    ALT: str,
) -> pl.DataFrame:
    """
    Helper function to create id column
    """
    id_column = (
        polars_df[CHR].str.replace("chr", "")
        + pl.lit(":")
        + polars_df[POS].cast(str)
        + pl.lit(":")
        + polars_df[REF]
        + pl.lit(":")
        + polars_df[ALT]
    )
    polars_df = polars_df.with_columns(id_column.alias(new_column_name))
    return polars_df


if __name__ == "__main__":
    defopt.run(combine)
