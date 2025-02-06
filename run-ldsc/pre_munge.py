"""
Prep METAL output for GSEM::LDSC.munge()
"""

from pathlib import Path
from typing import List, Optional
from math import exp
import defopt
import polars as pl
import attr
import re


@attr.s(frozen=False, auto_attribs=True, kw_only=True)
class Columns:
    """Commmon names for columns from summary stats"""

    chrom: Optional[str]
    pos: Optional[str]
    non_effect_allele: Optional[str]
    effect_allele: Optional[str]
    eaf: Optional[str]
    se: Optional[str]
    beta: Optional[str]
    pval: Optional[str]
    variant_id: Optional[str]
    total_n: Optional[str] = None
    case_n: Optional[str] = None
    control_n: Optional[str] = None
    imputation: Optional[str] = None


# pylint: disable=R0914, R0913, R0903, C0301
def prep_metal_output(
    *,
    gwas_results: Path,
    gwas_software: str,
    output_file: Path,
    rsid_map: Path,
    chrom_col: str,
    pos_col: str,
    ea: str,
    eaf: str,
    non_ea: str,
    beta: str,
    pval: str,
    se: str,
    variant_id: str,
    n_case: str,
    n_control: str,
    n_total: str,
    calculated_total: float
) -> pl.DataFrame:
    """
    :param gwas_results: Path to metal results file (*1_.tbl)
    :param gwas_software: Path to metal results file (*1_.tbl)
    :param chrom_col: Name of chrom col
    :param postion: NAme of pos col
    :param ea: effect allele
    :param eaf: str
    :param non_ea: Ref
    :param beta: str
    :param pval: str
    :param se: str
    :param variant_id: str
    :param n_case: str
    :param n_control: str
    :param n_total: str
    :param output_file: Path to file to write unusable variants to
    :param rsid_map: Path to tsv with two columns ID & RSID (chr:pos:ref:alt\trsid)
    :param calculated_total: provide number of total
    """
    found_columns: Columns = Columns(
        chrom=_search_header_for_positions(chrom_col),
        pos=_search_header_for_positions(pos_col),
        non_effect_allele=_search_header_for_positions(non_ea),
        effect_allele=_search_header_for_positions(ea),
        eaf=_search_header_for_positions(eaf),
        beta=_search_header_for_positions(beta),
        pval=_search_header_for_positions(pval),
        total_n=_search_header_for_positions(n_total),
        case_n=_search_header_for_positions(n_case),
        control_n=_search_header_for_positions(n_control),
        se=_search_header_for_positions(se),
        variant_id=_search_header_for_positions(variant_id),
    )
    # Read in gwas summary stats based on analysis software 'REGENIE' or 'SAIGE'
    try:
        if gwas_software.lower() == "saige":
            raw_df: pl.DataFrame = read_gwas(gwas_results, found_columns, sep="\t")
        elif gwas_software.lower() == "regenie":
            raw_df: pl.DataFrame = read_gwas(gwas_results, found_columns, sep=" ")
        elif gwas_software.lower() == "nan":
            print("No gwas-software was provided, assuming data is tab separated\n")
            raw_df: pl.DataFrame = read_gwas(gwas_results, found_columns, sep="\t")
        else:
            raise ValueError(
                f"The provided GWAS software: {gwas_software} is not supported -- available options are regenie or saige\n"
            )
    except FileNotFoundError:
        print(f"File '{gwas_results}' not found.")

    # Check that an ID column exists, if not add one named 'MarkerID'
    if found_columns.variant_id is None:
        id_column = (
            raw_df[found_columns.chrom].cast(str).replace("chr", "").replace("23", "X")
            + pl.lit(":")
            + raw_df[found_columns.pos].cast(str)
            + pl.lit(":")
            + raw_df[found_columns.non_effect_allele]
            + pl.lit(":")
            + raw_df[found_columns.effect_allele]
        )
        raw_df = raw_df.with_columns(id_column.alias("MarkerID"))
        found_columns.variant_id = "MarkerID"

    # Remove chr from ID column if present, and enforce all separators be ':'
    raw_df = raw_df.with_columns(
        pl.col(found_columns.variant_id)
        .str.replace_all("chr", "")
        .alias(found_columns.variant_id)
    )
    raw_df = raw_df.with_columns(
        pl.col(found_columns.variant_id)
        .str.replace_all("_", ":")
        .alias(found_columns.variant_id)
    )
    raw_df = raw_df.with_columns(
        pl.col(found_columns.variant_id)
        .str.replace_all("/", ":")
        .alias(found_columns.variant_id)
    )
    if gwas_software == "regenie":
        sep = " "
    else:
        sep = "\t"

    # Attach rsid
    metal_pl = raw_df
    metal_pl = metal_pl.with_columns(
        pl.col(found_columns.variant_id)
        .str.replace_all("chr", "")
        .alias(found_columns.variant_id)
    ).with_columns(
        "chr" + pl.col(found_columns.variant_id).alias(found_columns.variant_id)
    )

    print("Metal Output")
    print(metal_pl)
    rsids = (
        pl.scan_csv(
            rsid_map,
            separator="\t",
            has_header=False,
            new_columns=(found_columns.variant_id, "rsid"),
        )
        .filter(
            (pl.col(found_columns.variant_id).is_in(metal_pl[found_columns.variant_id]))
        )
        .collect()
    )
    print("RSID Output")
    print(rsids)

    # Check the N columns
    print(found_columns)
    if (found_columns.total_n is not None) and (found_columns.case_n is not None) and (found_columns.control_n is not None):
        metal_pl = metal_pl.with_columns((pl.col(found_columns.case_n)+pl.col(found_columns.control_n)).alias('N'))
    if (found_columns.total_n is None) and (found_columns.case_n is not None) and (found_columns.control_n is not None):
        metal_pl = metal_pl.with_columns((pl.col(found_columns.case_n)+pl.col(found_columns.control_n)).alias('N'))
    if (found_columns.total_n is None) and (found_columns.case_n is None) and (found_columns.control_n is None):
        metal_pl = metal_pl.with_columns(N=calculated_total)
    # Join tables
    #final = rsids.join(metal_pl, on=found_columns.variant_id, how="inner")
    final = metal_pl.rename({'SNP_ID':'rsid'})
    if found_columns.pval == "LOG10P":
        final = final.with_columns((10 ** (-1 * pl.col("LOG10P"))).alias("LOG10P"))
    final = (
        final.drop(found_columns.variant_id)
        .rename(
            {
                "rsid": "ID",
                found_columns.pval: "Pval",
                found_columns.eaf: "EAF",
                found_columns.non_effect_allele: "A1",
                found_columns.effect_allele: "A2",
                found_columns.pos: "POS",
                found_columns.chrom: "CHR",
                found_columns.beta: "BETA",
                found_columns.se: "SE",

            }
        )
        .select(
            [
                "ID",
                "A1",
                "A2",
                "POS",
                "CHR",
                "BETA",
                "EAF",
                "Pval",
                "SE",
                "N",
            ]
        )
    )

    # Calculate ODDS Ratio & Add counts
    final = final.with_columns(OR=pl.col("BETA").exp())
    final.write_csv(output_file, separator="\t")


def _search_header_for_positions(col_name: str) -> Optional[str]:
    if col_name is None or col_name == "nan":
        return None
    return col_name


def read_gwas(filename: Path, column_map: Columns, sep: str) -> pl.DataFrame:
    """
    Checks if the provided summary stat file is gzipped, if necessary de-compresses file. Then scans file, only reading in lines with information pertaining to the provided chromosome.

    Args:
        filename: Path to summary stat file
        column_map: Instance of class Columns for that summary stat
        chromosome: Which chromosome to read in
        sep: File deliminator
    """
    # check if file is gzipped
    # convert chromosomes to string so that X and Y can be handled in addition to int type 1-22
    if re.search(".gz$", f"{filename}"):
        gwas = pl.read_csv(
            filename,
            separator=sep,
            null_values=["NA"],
            truncate_ragged_lines=True,
            dtypes={column_map.chrom: str, column_map.pos: pl.Float64},
        )

    else:
        gwas = (
            pl.scan_csv(
                filename,
                separator=sep,
                null_values=["NA"],
                truncate_ragged_lines=True,
                dtypes={column_map.chrom: str, column_map.pos: pl.Float64},
            )
            .collect()
        )

    return gwas.with_columns(
        pl.col(column_map.pos).cast(pl.Int64()).alias(column_map.pos)
    )


if __name__ == "__main__":
    defopt.run(prep_metal_output)
