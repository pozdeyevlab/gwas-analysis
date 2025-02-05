"""
Module to read in gnomad reference and flag variants that have an AN less than 50% of the max AN
"""
import csv
import re
from pathlib import Path
from typing import List, Optional, Union

import defopt
import polars as pl

# pylint: disable = C0301


def read_reference(
    *,
    gnomad_tsv: Path,
    output_tsv: Path,
) -> None:
    """"
    :param gnomad_tsv: Path to gnomAD reference file
    :param output_tsv: Path to output file

    """
    header: List[str] = _get_header(reference=gnomad_tsv)

    an_column = "[76]AN"

    # Read in only necessary columns
    names_dict = _make_names_dict(
        allele_number=an_column, header=header
    )

    # Read in only necessary columns
    dtypes_dict = _make_dtypes_dict(
        allele_number=an_column, header=header
    )

    try:
        gnomad_pl: pl.DataFrame = pl.read_csv(
            gnomad_tsv,
            separator="\t",
            columns=list(names_dict.keys()),
            dtypes=dtypes_dict,
            infer_schema_length=10000,
        )
    except pl.exceptions.ComputeError as error:
        print(
            f"There is an error reading in the gnomad reference file {gnomad_tsv}\nError: {error}"
        )

    gnomad_pl = gnomad_pl.rename(names_dict)
    gnomad_pl = gnomad_pl.filter(pl.col("AN") > 0)

    # Filter for AN < 50% of max(AN)
    max_an = gnomad_pl["AN"].max()
    gnomad_pl = gnomad_pl.filter(pl.col("AN") < max_an/2)

    # Add ID column
    # Create the 'ID' column by concatenating values from 'CHR', 'POS', 'REF', and 'ALT'
    id_column = (
        gnomad_pl["CHR"].str.replace('chr', '')
        + pl.lit(":")
        + gnomad_pl["POS"].cast(str)
        + pl.lit(":")
        + gnomad_pl["REF"]
        + pl.lit(":")
        + gnomad_pl["ALT"]
    )

    # Create a new DataFrame with the 'ID' column added
    gnomad_pl = gnomad_pl.with_columns(id_column.alias("ID"))

    # Write to tsv
    gnomad_pl.write_csv(output_tsv, separator="\t")


def _make_names_dict(
    allele_number: str, header: List[str]
) -> dict:
    dtypes = {}

    # Constants
    chromsome = _search_patterns_in_header(pattern=r"CHROM", header=header)[0]
    position = _search_patterns_in_header(pattern=r"POS", header=header)[0]
    ref = _search_patterns_in_header(pattern=r"REF", header=header)[0]
    alt = _search_patterns_in_header(pattern=r"ALT", header=header)[0]
    filter_flag = _search_patterns_in_header(pattern=r"FILTER", header=header)[0]

    # Add to empty dictionary
    dtypes[chromsome] = "CHR"
    dtypes[position] = "POS"
    dtypes[ref] = "REF"
    dtypes[alt] = "ALT"
    dtypes[filter_flag] = "FILTER"
    dtypes[allele_number] = "AN"

    return dtypes


def _make_dtypes_dict(
    allele_number: str, header: List[str]
) -> dict:
    dtypes = {}

    # Constants
    chromsome = _search_patterns_in_header(pattern=r"CHROM", header=header)[0]
    position = _search_patterns_in_header(pattern=r"POS", header=header)[0]
    ref = _search_patterns_in_header(pattern=r"REF", header=header)[0]
    alt = _search_patterns_in_header(pattern=r"ALT", header=header)[0]
    filter_flag = _search_patterns_in_header(pattern=r"FILTER", header=header)[0]

    # Add to empty dictionary
    dtypes[chromsome] = pl.Utf8
    dtypes[position] = pl.Int32
    dtypes[ref] = pl.Utf8
    dtypes[alt] = pl.Utf8
    dtypes[filter_flag] = pl.Utf8
    dtypes[allele_number] = pl.Int32

    return dtypes


def _get_header(reference: Path) -> List[str]:
    with open(reference, "r", encoding="utf-8") as file:
        tsv_reader = csv.reader(file, delimiter="\t")
        header = next(tsv_reader)
        return header


def _search_patterns_in_header(pattern: str, header: List[str]) -> List[str]:
    matching_elements = [element for element in header if re.search(pattern, element)]
    filtered_result = [
        element for element in matching_elements if "joint" not in element.lower()
    ]
    return filtered_result


if __name__ == "__main__":
    defopt.run(read_reference)
