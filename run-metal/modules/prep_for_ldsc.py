"""
This module prepares metal outputs for linkage-disequilibrium score regression (LDSC). 

1)	If the ref & alt alleles in the snp id do not match the reported ref and alt allele then the beta value is inversed to reflect the true effect allele. 
2)	rsID’s are attached to snp id’s
3)	Remove variants that do not meet the required number of biobanks
4)	Calculates the odds ratio of the beta value
5)	Writes a tab separated file, compatible with the first step of ldsc, munge.py
"""

from math import exp
from pathlib import Path
from typing import List, Optional, Union

import defopt
import polars as pl

# pylint: disable=R0914, R0913, R0903, C0301


def prep_metal_output(
    *,
    metal_results: Path,
    cases: Union[int,float],
    controls: Union[int,float],
    output_path: Path,
    rsid_map: Path,
    biobank_min: int,
    or_filter: float = None,
    af_filter: float = None,
) -> pl.DataFrame:
    """
    :param metal_results: Path to metal results file (*1_.tbl)
    :param cases: Number of total cases
    :param controls: Number of total controls
    :param output_path: Path to file to write unusable variants to
    :param rsid_map: Path to tsv with two columns ID & RSID (chr:pos:ref:alt\trsid)
    :param biobank_min: The number of biobanks each variant must be found in
    :param or_filter: OR minimum
    :param af_filter: AF minimum
    """
    metal_pl = pl.read_csv(metal_results, separator="\t")
    metal_pl = (
        metal_pl.with_columns(
            pl.col("MarkerName")
            .str.split(":")
            .map_elements(lambda arr: arr[2])
            .alias("REF_FROM_ID")
        )
        .with_columns(
            pl.col("MarkerName")
            .str.split(":")
            .map_elements(lambda arr: arr[3])
            .alias("ALT_FROM_ID")
        )
        .with_columns(
            pl.when((pl.col("REF_FROM_ID") == pl.col("Allele1").str.to_uppercase()))
            .then(pl.col("Effect"))
            .otherwise(-1 * pl.col("Effect"))
            .alias("BETA")
        )
        .with_columns(
            pl.when((pl.col("REF_FROM_ID") == pl.col("Allele1").str.to_uppercase()))
            .then(pl.col("Freq1"))
            .otherwise(1 - pl.col("Freq1"))
            .alias("EAF")
        )
        .with_columns(pl.col("MarkerName").str.replace("chr", "").alias("MarkerName"))
        .with_columns("chr" + pl.col("MarkerName").alias("MarkerName"))
    )
    if (af_filter is not None) and (or_filter is not None):
        metal_pl = metal_pl.filter(pl.col('EAF') > af_filter)
        metal_pl = metal_pl.filter(pl.col('BETA').exp() < or_filter)

    biobank_count = metal_pl.with_columns(pl.col('Direction').str.len_bytes().alias('temp'))['temp'].max()
    print(biobank_min)
    if int(biobank_min) <= int(biobank_count):
        metal_pl = metal_pl.filter(
            (pl.col("Direction").str.count_matches(r"\+|\-")) >= biobank_min
        )
    else:
        metal_pl = metal_pl.filter(
            (pl.col("Direction").str.count_matches(r"\+|\-")) >= 2
        )
    print("Metal Output")
    print(metal_pl)

    rsids = (
        pl.read_csv(rsid_map, separator="\t")
        .with_columns(
            (
                pl.col("Chr").cast(str)
                + ":"
                + pl.col("Start").cast(str)
                + ":"
                + pl.col("Ref")
                + ":"
                + pl.col("Alt")
            ).alias("MarkerName")
        )
        .rename({"avsnp151": "rsid"})
    )

    print("RSID Output")
    print(rsids)

    # Join tables
    final = rsids.join(metal_pl, on="MarkerName", how="inner")
    print(final)
    final = (
        final.with_columns(
            pl.col("MarkerName")
            .str.split(":")
            .map_elements(lambda arr: arr[0])
            .alias("CHR"),
            pl.col("MarkerName")
            .str.split(":")
            .map_elements(lambda arr: arr[1])
            .alias("POS"),
        )
        .drop("MarkerName")
        .rename(
            {
                "rsid": "MarkerName",
                "StdErr": "SE",
                "ALT_FROM_ID": "A2",
                "REF_FROM_ID": "A1",
                "P-value": "Pval",
                "per_variant_N": "N"
            }
        )
        .select(
            ["MarkerName", "A1", "A2", "POS", "CHR", "BETA", "EAF", "Pval", "SE", "N"]
        )
    )

    # Calculate ODDS Ratio & Add counts
    final = final.with_columns(OR=pl.col("BETA").exp())
    final = final.with_columns(CASES=cases)
    final = final.with_columns(CONTROLS=controls)

    final.write_csv(output_path, separator="\t")


def read_gwas(
    filenames: List[Path], columns: List[str], chrom: str, pos: str, ref: str, alt: str
) -> Optional[pl.DataFrame]:
    """
    Checks if the provided summary stat file is gzipped, if necessary de-compresses file. Then scans file, only reading in lines with information pertaining to the provided chromosome.

    Args:
        filename: Path to summary stat file
        column_map: Instance of class Columns for that summary stat
        sep: File deliminator
    """
    # Check if the total n, or case & control n columns exist
    dfs = []
    for filename in filenames:
        df = pl.read_csv(
            filename,
            separator="\t",
            columns=columns,
            null_values=["NA"],
            dtypes={chrom: str, pos: pl.Float64},
        )
        dfs.append(df)

    gwas = pl.concat(dfs, how="vertical")

    return (
        gwas.with_columns(pl.col(pos).cast(pl.Int64()).alias(pos))
        .filter((pl.col(chrom).cast(str) != "23"))
        .filter((pl.col(chrom).cast(str) != "X"))
        .with_columns(
            SNPID=pl.col(chrom).str.replace("chr", "")
            + ":"
            + pl.col(pos).cast(str)
            + ":"
            + pl.col(ref)
            + ":"
            + pl.col(alt)
        )
    )


def read_rsid_reference(
    rsid_map: Path, study_pl: pl.DataFrame, snp_col: str
) -> pl.DataFrame:
    # Scan rsid map file for variants found in study data
    rsid = (
        pl.scan_csv(
            rsid_map,
            separator="\t",
            has_header=False,
            new_columns=[snp_col, "RSID"],
            low_memory=True,
        )
        .filter((pl.col(snp_col).is_in(study_pl[snp_col])))
        .collect()
    )
    print(rsid)

    before_count = study_pl.shape[0]

    # Merge to add rsid
    study_pl = study_pl.join(rsid, on=snp_col, how="inner")

    print(
        f"{before_count - study_pl.shape[0]}/{before_count} variants will be excluded as they do not have matching rsids"
    )
    return study_pl


if __name__ == "__main__":
    defopt.run(prep_metal_output)

