"""
Module to format gnomad aligned summary stats for METAL
Removes variants that:
1) Have a gnomAD AN warning
2) Have significantly different AF compared to gnomAD (based on mahalanobis distance)
3) Have MAC (minor allele count less than 20)
4) gnomAD variant does not PASS QC
5) Have MAF less than 0.0005 (0.05%)
"""

from pathlib import Path
from typing import Optional, Union

import defopt
import polars as pl

# pylint: disable=R0914, R0913, R0903, C0301


def prep_qc_files(
    *,
    input_file: str,
    metal_output: str,
    output_path: Path,
) -> pl.DataFrame:
    """
    :param input_files: Input path to regenie or saige summary stats
    :param metal_output: Input path to regenie or saige summary stats
    :param output_path: Path to file to write unusable variants to
    """
    # Read and correct metal results:
    # Read in metal results
    metal_df = (
        (
            pl.read_csv(metal_output,  separator="\t")
            .with_columns(
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
            .with_columns(pl.col("Allele1").str.to_uppercase().alias("Allele1"))
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
        )
        .with_columns(
            pl.when(pl.col("reported_alleles_match_id"))
            .then(pl.col("Effect"))
            .otherwise(-1 * pl.col("Effect"))
            .alias("adjusted_metal_beta")
        )
    ).rename({"MarkerName": "STUDY_ID", 'adjusted_metal_beta':'beta', 'StdErr':'se', 'P-value':'pval'}).select(["STUDY_ID",'beta', 'se', 'pval'])
    
    # Read in old details
    sumstat = pl.read_csv(input_file, separator='\t').drop(['beta', 'se', 'pval'])
    combined_df = metal_df.join(sumstat, on = 'STUDY_ID', how = 'inner')
    combined_df = combined_df.select(
        [
            "STUDY_ID",
            "beta",
            "se",
            "Aligned_AF",
            "MAC",
            "pval",
            "REF_gnomad",
            "ALT_gnomad",
            "AF_gnomad",
            "Potential_Strand_Flip",
            "outlier_stdev",
            "per_variant_N",
            "missingness",
        ]
    )
    combined_df.write_csv(output_path, separator="\t")


if __name__ == "__main__":
    defopt.run(prep_qc_files)
