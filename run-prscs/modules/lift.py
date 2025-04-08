"""
Module to revert PRSCS results back to GRCh38
"""

from pathlib import Path
from typing import List, Optional

import defopt
import polars as pl

# pylint: disable = C0301
# pylint: disable = R0903 # Too few public methods
# pylint: disable = R1728


def lift(*, prscs_results_file: Path, bim_file: Path, output: str) -> None:
    """
    :param prscs_results_file: Path to prscs results file
    :param output_bim: Path to write lifted results file
    :param bim_file: Bim used as input to prscs
    """
    prscs_results = pl.read_csv(
        prscs_results_file,
        separator="\t",
        has_header=False,
        new_columns=["CHROM", "RSID", "POS", "REF", "ALT", "BETA"],
    )
    print(prscs_results)
    bim_file = pl.read_csv(
        bim_file,
        separator="\t",
        has_header=False,
        new_columns=["CHROM", "RSID", "EXTRA", "GRCh38_POS", "REF", "ALT"],
    )
    print(bim_file)
    # Add GRCH38 position to PRSCS output file
    prscs_results = (
        prscs_results.join(bim_file, on="RSID", how="inner")
        .with_columns(
            (
                "chr"
                + pl.col("CHROM").cast(str)
                + ":"
                + pl.col("GRCh38_POS").cast(str)
                + ":"
                + pl.col("REF")
                + ":"
                + pl.col("ALT")
            ).alias("ID")
        )
        .select(["ID", "CHROM", "RSID", "GRCh38_POS", "REF", "ALT", "BETA"])
        .unique()
        .write_csv(output, separator="\t", include_header=False)
    )
    return None


if __name__ == "__main__":
    defopt.run(lift)
