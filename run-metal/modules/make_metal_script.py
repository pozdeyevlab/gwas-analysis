"""
This module creates the script which is used to run inverse-variance weighted meta-analysis using METAL. The output can be found in {output_dir}/metal_results{phenotype}.metal
"""

from pathlib import Path
from typing import List

import defopt
import polars as pl

# pylint: disable=R0914, R0913, R0903, C0301


def write_metal(
    *, input_files: List[str], gc: str, output_path: Path, metal_out: str
) -> pl.DataFrame:
    """
    :param input_files: Input path to regenie or saige summary stats
    :param output_path: Path to file to write unusable variants to
    :param gc: Genomic Control ON or OFF
    """
    # For every file add command to metal script
    with open(output_path, mode="w+", encoding="utf-8") as out:
        out.writelines(
            f"SCHEME STDERR\nSEPARATOR TAB\nAVERAGEFREQ ON\nMINMAXFREQ ON\nGENOMICCONTROL {gc}\n CUSTOMVARIABLE per_variant_N\nLABEL per_variant_N as per_variant_N\n"
        )
        for file in input_files:
            out.writelines(
                f"MARKER STUDY_ID\nALLELE REF_gnomad ALT_gnomad\nEFFECT beta\nSTDERR se\nFREQLABEL Aligned_AF\nPROCESS {file}\n\n"
            )
        out.write(f"OUTFILE {metal_out}\nANALYZE\nQUIT")


if __name__ == "__main__":
    defopt.run(write_metal)

