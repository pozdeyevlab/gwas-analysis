"""
Module to read in gnomad reference
This module completes the following steps:
1) Searches for the AN (allele number) and
AF (allele frequency) columns in the gnomad header
that correspond to the sex and ancestry of the summary stat (if applicable)
2) Reads in only the required columns
3) Returns a polars df with the following column names 'CHR', 'POS', 'REF', 'ALT', 'AN', 'AF'

"""
import csv
import re
import sys
from datetime import datetime
from pathlib import Path
from typing import List, Optional, Union

import defopt
import polars as pl

# pylint: disable = C0301


def read_reference(
    *,
    gnomad_tsv: Path,
    sex: str,
    ancestry: str,
    positions: List[int],
) -> None:
    """
    :param sex: Sex if any
    :ancestry: Three letter ancestry code if any
    :param gnomad_tsv: Path to gnomAD reference file
    :param positions: List of positions from gwas summary stats
    """
    gnomad_get_start = datetime.now()
    print(f"Starting to search for AN and AF:\n{datetime.now()}")
    print(f"Provided sex = {sex}\nProvided ancestry = {ancestry}")
    header: List[str] = _get_header(reference=gnomad_tsv)

    # Search for AN & AF columns according to sex and ancestry
    an_columns = find_correct_columns(
        header, "AN", sex=_convert_sex(sex), ancestry=_convert_ancestry(ancestry)
    )
    print(f"Found AN = {an_columns}")

    af_columns = find_correct_columns(
        header, "AF", sex=_convert_sex(sex), ancestry=_convert_ancestry(ancestry)
    )
    print(f"Found AF = {af_columns}")

    # Read in only necessary columns
    names_dict = _make_names_dict(
        allele_number=an_columns, allele_frequency=af_columns, header=header
    )

    # Read in only necessary columns
    dtypes_dict = _make_dtypes_dict(
        allele_number=an_columns, allele_frequency=af_columns, header=header
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
    gnomad_pl = gnomad_pl.filter(pl.col("POS_gnomad").is_in(positions))
    gnomad_pl = gnomad_pl.filter(pl.col("AN_gnomad") > 0)

    # Add ID column
    # Create the 'ID_gnomad' column by concatenating values from 'CHR', 'POS', 'REF', and 'ALT'
    id_column = (
        gnomad_pl["CHR_gnomad"].cast(str)
        + pl.lit(":")
        + gnomad_pl["POS_gnomad"].cast(str)
        + pl.lit(":")
        + gnomad_pl["REF_gnomad"]
        + pl.lit(":")
        + gnomad_pl["ALT_gnomad"]
    )

    # Create a new DataFrame with the 'ID_gnomad' column added
    gnomad_pl = gnomad_pl.with_columns(id_column.alias("ID_gnomad"))

    gnomad_end = datetime.now()
    total = gnomad_end - gnomad_get_start
    print(
        f"\nFinished searching and filtering gnomAD reference in {total}:\n{datetime.now()}\n"
    )
    print(gnomad_pl)
    return gnomad_pl


def _make_names_dict(
    allele_number: str, allele_frequency: str, header: List[str]
) -> dict:
    dtypes = {}

    # Constants
    chromsome = _search_patterns_in_header(pattern=r"CHROM", header=header)[0]
    position = _search_patterns_in_header(pattern=r"POS", header=header)[0]
    ref = _search_patterns_in_header(pattern=r"REF", header=header)[0]
    alt = _search_patterns_in_header(pattern=r"ALT", header=header)[0]
    filter_flag = _search_patterns_in_header(pattern=r"FILTER", header=header)[0]

    # Add to empty dictionary
    dtypes[chromsome] = "CHR_gnomad"
    dtypes[position] = "POS_gnomad"
    dtypes[ref] = "REF_gnomad"
    dtypes[alt] = "ALT_gnomad"
    dtypes[filter_flag] = "FILTER"
    dtypes[allele_number] = "AN_gnomad"
    dtypes[allele_frequency] = "AF_gnomad"

    return dtypes


def _make_dtypes_dict(
    allele_number: str, allele_frequency: str, header: List[str]
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
    dtypes[allele_frequency] = pl.Float32

    return dtypes


def find_correct_columns(header: List[str], field: str, sex: str, ancestry: str) -> str:
    """
    Searches for the correct AN or AF column in the gnomad reference given the sex and ancestry of the gwas summary stat and returns the column name corresponding to the gnomad reference.

    Args:
        header: The header from the gnoamd reference file
        field: AN (allele number) or AF (allele frequency)
        sex: Female, male or nan
        ancestry: Ancestry code or nan
    """
    matches = _search_patterns_in_header(
        pattern=_generate_pattern(field, sex, ancestry), header=header
    )
    try:
        matched = matches[0]
        return matched
    except IndexError:
        print(
            f"No matches were found for field:{field} sex:{sex} ancestry:{ancestry} in this header:\n{header}\nUsing sex:None and ancestry:None instead"
        )
        try:
            matches = _search_patterns_in_header(
                pattern=_generate_pattern(field, None, None), header=header
            )
            matched = matches[0]
            return matched
        except IndexError:
            print(
                f"No matches were found for field:{field} sex:{sex} ancestry:{ancestry} in this header:\n{header}"
            )
            sys.exit()


def _get_header(reference: Path) -> List[str]:
    with open(reference, "r", encoding="utf-8") as file:
        tsv_reader = csv.reader(file, delimiter="\t")
        header = next(tsv_reader)
        return header


def _generate_pattern(
    field: str, sex: Optional[str], ancestry: Optional[str]
) -> Union[str, ValueError]:
    # Constructing the regex pattern dynamically
    if sex is not None and ancestry is not None:
        return rf"(?:{field}|{sex}|{ancestry}).*(?:{field}|{sex}|{ancestry}).*(?:{field}|{sex}|{ancestry})"
    if sex is not None and ancestry is None:
        return rf"(?:{field}).*(?:{sex})"
    if sex is None and ancestry is not None:
        return rf"(?=.*{field})(?=.*{ancestry})(?!.*_(XX|XY))"
    if sex is None and ancestry is None:
        return rf"(?=.*{field})"
    return ValueError


def _search_patterns_in_header(pattern: str, header: List[str]) -> List[str]:
    matching_elements = [element for element in header if re.search(pattern, element)]
    filtered_result = [
        element for element in matching_elements if "joint" not in element.lower()
    ]
    return filtered_result


def _convert_sex(sex: str) -> Optional[str]:
    if sex.lower() == "female":
        return "XX"
    if sex.lower() == "male":
        return "XY"
    return None


def _convert_ancestry(ancestry: str) -> Optional[str]:
    if ancestry.lower() == "nan":
        return None
    if ancestry.lower() == "eur":
        return "nfe"
    return ancestry.lower()


if __name__ == "__main__":
    defopt.run(read_reference)
