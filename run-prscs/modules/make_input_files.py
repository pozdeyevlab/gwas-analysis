"""
Module to read do meta analysis
"""

from pathlib import Path
from typing import List, Optional

import defopt
import polars as pl

# pylint: disable = C0301
# pylint: disable = R0903 # Too few public methods
# pylint: disable = R1728


def make_bim(*, meta_results_file: Path, rsid_file: Path, output_sst: str) -> None:
    """
    :param meta_analysis_file: Path to meta analyzed results file
    :param output_bim: Path to write bim file
    :param output_sst: Path to write PRScs formatted summary stat to
    """
    # Read in METAL analyzed results and create chrom, pos, ref and alt from MarkerName column
    metal_results = (
        pl.read_csv(meta_results_file, separator="\t")
        .with_columns(
            pl.col("MarkerName")
            .str.split(":")
            .map_elements(lambda arr: arr[0], return_dtype=str)
            .alias("chrom")
        )
        .with_columns(
            pl.col("MarkerName")
            .str.split(":")
            .map_elements(lambda arr: arr[1], return_dtype=str)
            .alias("pos")
        )
        .with_columns(
            pl.col("MarkerName")
            .str.split(":")
            .map_elements(lambda arr: arr[2], return_dtype=str)
            .alias("ref")
        )
        .with_columns(
            pl.col("MarkerName")
            .str.split(":")
            .map_elements(lambda arr: arr[3], return_dtype=str)
            .alias("alt")
        )
    )

    metal_results = metal_results.with_columns(
        pl.when(pl.col("ref").str.to_lowercase() != pl.col("Allele1"))
        .then(1 - pl.col("Freq1").cast(float))
        .otherwise(pl.col("Freq1"))
        .alias("Freq1")
    )

    metal_results = metal_results.with_columns(
        pl.when(pl.col("ref").str.to_lowercase() != pl.col("Allele1"))
        .then(-1 * pl.col("Effect").cast(float))
        .otherwise(pl.col("Effect"))
        .alias("Effect")
    )

    # Add RSID's based on reference file extracted from gnomAD VCF
    chroms = list(set(metal_results["chrom"]))
    try:
        chroms.remove("X")
        chroms.remove("23")
    except:
        print("X not in chroms")
    dfs = []
    for chrom in chroms:
        temp = metal_results.filter(pl.col("chrom") == chrom)
        rsid = (
            pl.scan_csv(rsid_file, separator="\t", has_header=False)
            .filter(pl.col("column_1").is_in(temp["MarkerName"]))
            .rename({"column_1": "MarkerName", "column_2": "rsid"})
            .filter(pl.col("MarkerName").is_in(temp["MarkerName"]))
            .collect()
        )
        temp_merged = rsid.join(temp, on="MarkerName", how="inner")
        dfs.append(temp_merged)

    merged = pl.concat(dfs, how="vertical")
    merged = merged.with_columns(extra=0)
    # bim_file = merged.select(["chrom", "rsid", "extra", "pos", "ref", "alt"])
    sst_file = merged.select(["rsid", "ref", "alt", "Effect", "StdErr"]).rename(
        {"rsid": "SNP", "ref": "A1", "alt": "A2", "Effect": "BETA", "StdErr": "SE"}
    )

    # bim_file.write_csv(output_bim, separator="\t", include_header=False)
    sst_file.write_csv(output_sst, separator="\t", include_header=True)

    return None


if __name__ == "__main__":
    defopt.run(make_bim)
