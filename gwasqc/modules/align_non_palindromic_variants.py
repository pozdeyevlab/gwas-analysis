"""
Align non-palindromic variants

This module is used to align non-palindromic variants. 
In order to compare variants from multiple biobanks it is crucial to align variants to a common reference. This module handles non-palindromic variants and checks for the following scenarios:

1) Exact match - study reference == gnoamd reference & study alternate == gnomad alternate
2) Inverse match  -  study reference == gnomad alternate & study alternate == gnomad reference
3) Exact transcribed match - transcribed study reference == gnoamd reference & transcribed study alternate == gnomad alternate
4) Inverse transcribed match - transcribed study reference == gnoamd alternate & transcribed study alternate == gnomad reference

This module also handles complementary variants. Meaning a study variant in whcih both exact and inverse matches are valid gnomad variants. In these cases the gnoamd variant with the allele frequency closest to the study frequency is assigned as the true match. 

Alignment Logic:
In order to avoid OOM errors this module uses set comparison to determine matches between the study and gnomad. For each possible match (exact, inverse, exact transcribed, and inverse transcribed) two sets of variant ID's are created. The first is consisted for all scenarios, the gnomad set, the second changes according to the scenario, this is from the study. The set intersection is used to create an new data frame, and the alignment method writted to the column 'Alignment Method'. After all scenarios have been checked data frames containing 
"""
from collections import namedtuple
from typing import List, Optional

import attr
import defopt
import filter_gwas
import numpy as np
import polars as pl

# pylint: disable=C0301
# pylint: disable=R0914 # too many local variables
# pylint: disable=R0915 # Too many statements
# pylint: disable=R0913 # too many arguments
# pylint: disable=R0903 # too few public methods


@attr.s(frozen=False, auto_attribs=True, kw_only=True)
class AlignmentResults:
    """Class to represent the alignment outcome"""

    aligned_and_merged: Optional[pl.DataFrame]
    unaligned: pl.DataFrame


def harmonize(
    *,
    gnomad_pl: pl.DataFrame,
    gwas_pl: pl.DataFrame,
    col_map: filter_gwas.Columns,
) -> None:
    """
    :param gnomad_pl: Reference Data Frame
    :param gwas_pl: Study Data Frame
    :param col_map: Class of columns (see inititation of class in modules/harmonize.pt)
    """
    # For book-keeping print the number of non-palindromic complimentary matches to stdout
    # Find the number of non-palindromic inverse pairs in gnomad (chr10:100:A:G & chr10:100:G:A)
    new = _make_id_column_gnomad(
        new_column_name="Flipped_Allele_ID",
        polars_df=gnomad_pl,
    )

    matches_count = new.filter(pl.col("ID_gnomad").is_in(pl.col("Flipped_Allele_ID")))

    print(
        f"\nNon-Palindromic Summary:\nTotal complementary non-palindromic variants: {matches_count.shape[0]}"
    )
    ################################ Start of Alignment ################################
    list_of_results: List[AlignmentResults] = []
    # Exact match (ref=ref & alt=alt)
    exact: AlignmentResults = align_alleles(
        gnomad_df=gnomad_pl,
        gwas_df=gwas_pl,
        id_column=col_map.variant_id,
        method="exact_match",
    )
    # Calculate the variants with an abs difference in AF between study and gnomad greater than 0.1
    #likely_exact = (
    #    exact.aligned_and_merged.with_columns(
    #        ABS_DIF_AF=abs(pl.col("AF_gnomad") - pl.col(col_map.eaf))
    #    )
    #    .filter(pl.col("ABS_DIF_AF") < 0.05)
    #    .drop("ABS_DIF_AF")
    #)

    # Inverse (ref=alt & alt=ref)
    gwas_pl = _make_id_column(
        new_column_name="Flipped_Allele_ID",
        col_map=col_map,
        first_allele=col_map.effect_allele,
        second_allele=col_map.non_effect_allele,
        polars_df=gwas_pl.join(exact.aligned_and_merged, on=col_map.variant_id, how="anti"),
    )

    inverse: AlignmentResults = align_alleles(
        gnomad_df=gnomad_pl,
        gwas_df=gwas_pl,
        id_column="Flipped_Allele_ID",
        method="inverse_match",
    )
    if (inverse.aligned_and_merged is not None) and (exact.aligned_and_merged is not None):
        results = handle_complementary_variants(
            exact_pl=exact.aligned_and_merged,
            inverse_pl=inverse.aligned_and_merged,
            col_map=col_map,
        )
        exact.aligned_and_merged = results.exact
        inverse.aligned_and_merged = results.inverse
        list_of_results.append(exact)
        list_of_results.append(inverse)
    if (inverse.aligned_and_merged is not None) and (exact.aligned_and_merged is None):
        list_of_results.append(inverse)
    if (inverse.aligned_and_merged is None) and (exact.aligned_and_merged is not None):
        list_of_results.append(exact)

    # Transcribed exact (ref = transcribed(ref) alt = transcribed(alt))
    gwas_pl = _make_id_column(
        new_column_name="Transcribed_ID",
        col_map=col_map,
        first_allele=col_map.transcribed_non_effect_allele,
        second_allele=col_map.transcribed_effect_allele,
        polars_df=inverse.unaligned,
    )

    transcribed: AlignmentResults = align_alleles(
        gnomad_df=gnomad_pl,
        gwas_df=gwas_pl,
        id_column="Transcribed_ID",
        method="transcribed_match",
    )
    list_of_results.append(transcribed)

    # Transcribed inverse (ref = transcribed(alt) & alt = transcribed(ref))
    gwas_pl = _make_id_column(
        new_column_name="Transcribed_Flipped_ID",
        col_map=col_map,
        first_allele=col_map.transcribed_effect_allele,
        second_allele=col_map.transcribed_non_effect_allele,
        polars_df=transcribed.unaligned,
    )

    transcribed_flipped: AlignmentResults = align_alleles(
        gnomad_df=gnomad_pl,
        gwas_df=gwas_pl,
        id_column="Transcribed_Flipped_ID",
        method="transcribed_fliped_match",
    )
    list_of_results.append(transcribed_flipped)

    # Concatenate results
    list_of_dfs: List[pl.DataFrame] = []
    for result in list_of_results:
        df = result.aligned_and_merged
        if df is not None:
            list_of_dfs.append(df)

    stacked_pl = pl.concat(list_of_dfs, how="diagonal")

    # If needed flip beta and allele frequency
    stacked_pl = stacked_pl.with_columns(
        pl.when(
            (
                pl.col("Alignment_Method").is_in(
                    ["inverse_match", "transcribed_flipped_match"]
                )
            )
        )
        .then(-1 * pl.col(col_map.beta) if col_map.beta is not None else np.NaN)
        .otherwise(pl.col(col_map.beta) if col_map.beta is not None else np.NaN)
        .alias("Aligned_Beta")
    )

    stacked_pl = stacked_pl.with_columns(
        pl.when(
            (
                pl.col("Alignment_Method").is_in(
                    ["inverse_match", "transcribed_flipped_match"]
                )
            )
        )
        .then(1 - pl.col(col_map.eaf) if pl.col(col_map.eaf) is not None else np.NaN)
        .otherwise(pl.col(col_map.eaf))
        .alias("Aligned_AF")
    )

    # Write Abs(AF_study - AF_ref)
    stacked_pl = stacked_pl.with_columns(
        (abs(pl.col("AF_gnomad") - pl.col("Aligned_AF"))).alias("ABS_DIF_AF")
    )
    # Print Summary
    print(
        f"Total aligned non-palindromic variants with method 'exact_match': {stacked_pl.filter(pl.col('Alignment_Method') == 'exact_match').shape[0]}"
    )

    print(
        f"Total aligned non-palindromic variants with method 'inverse_match':{stacked_pl.filter(pl.col('Alignment_Method') == 'inverse_match').shape[0]}"
    )

    print(
        f"Total aligned non-palindromic variants with method 'transcribed_match': {stacked_pl.filter(pl.col('Alignment_Method') == 'transcribed_match').shape[0]}"
    )

    print(
        f"Total aligned non-palindromic variants with method 'inverse_match':{stacked_pl.filter(pl.col('Alignment_Method') == 'transcribed_flipped_match').shape[0]}"
    )

    return stacked_pl


def _make_id_column(
    *,
    new_column_name: str,
    col_map: filter_gwas.Columns,
    first_allele: str,
    second_allele: str,
    polars_df: pl.DataFrame,
) -> pl.DataFrame:
    """
    Helper function to create id columns with different effect and non effect allele combinaitons
    """
    id_column = (
        polars_df[col_map.chrom].cast(str)
        + pl.lit(":")
        + polars_df[col_map.pos].cast(str)
        + pl.lit(":")
        + polars_df[first_allele]
        + pl.lit(":")
        + polars_df[second_allele]
    )
    polars_df = polars_df.with_columns(id_column.alias(new_column_name))
    return polars_df


def _make_id_column_gnomad(
    *,
    new_column_name: str,
    polars_df: pl.DataFrame,
) -> pl.DataFrame:
    """
    Helper function to create id columns with differenc effect and non effect allele combinaitons
    """
    id_column = (
        polars_df["CHR_gnomad"].cast(str)
        + pl.lit(":")
        + polars_df["POS_gnomad"].cast(str)
        + pl.lit(":")
        + polars_df["ALT_gnomad"]
        + pl.lit(":")
        + polars_df["REF_gnomad"]
    )
    polars_df = polars_df.with_columns(id_column.alias(new_column_name))
    return polars_df


def align_alleles(
    *, gnomad_df: pl.DataFrame, gwas_df: pl.DataFrame, id_column: str, method: str
) -> AlignmentResults:
    """
    Uses set comparison to determine allele alignment between the reference (gnomad) variant id's and the gwas summary stat id's.

    Args:
        gnomad_df: Reference dataframe
        gwas_df: Gwas summary stat dataframe
        id_column: The name of the id_column in the gwas dataframe
        method: The descriptor for what alignemtn method is being tested
    """
    # Remove chr from both ID's if present
    gnomad_df = gnomad_df.with_columns(pl.col("ID_gnomad").str.replace("chr", "").alias("ID_gnomad"))

    gwas_df = gwas_df.with_columns(
        pl.col(id_column).str.replace("chr", "").alias(id_column)
    )

    # Create sets of id columns and find overlap
    gnomad_set: set = set(gnomad_df["ID_gnomad"])
    gwas_set: set = set(gwas_df[id_column])
    aligned_set: set = gnomad_set & gwas_set

    if len(aligned_set) > 0:
        # Filter gwas data for only aligned variants
        aligned_df: pl.DataFrame = gwas_df.filter(pl.col(id_column).is_in(aligned_set))
        unaligned_df: pl.DataFrame = gwas_df.filter(
            ~(pl.col(id_column).is_in(aligned_set))
        )
        aligned_df = aligned_df.with_columns(pl.lit(method).alias("Alignment_Method"))

        # Merge aligned_df with gnomad_df on id column
        joined_df = aligned_df.join(
            gnomad_df,
            left_on=id_column,
            right_on="ID_gnomad",
            how="left",
            suffix="_gnomad",
        )
        # Create instance of results class
        return AlignmentResults(aligned_and_merged=joined_df, unaligned=unaligned_df)

    unaligned_df = gwas_df
    return AlignmentResults(aligned_and_merged=None, unaligned=unaligned_df)


def handle_complementary_variants(
    *, exact_pl: pl.DataFrame, inverse_pl: pl.DataFrame, col_map: filter_gwas.Columns
) -> namedtuple:
    """
    # TODO add docstring
    Args:
        exact_pl: Data frame with exact matches
        inverse_pl: Data frame with inverse matches
        col_map; Class of column names
    """
    print(exact_pl)
    exact_pl_subset = exact_pl.with_columns(
        ABS_DIF_AF=abs((pl.col("AF_gnomad") - pl.col(col_map.eaf)))
    ).select(col_map.variant_id, col_map.eaf, "ABS_DIF_AF", "AF_gnomad", "REF_gnomad", "ALT_gnomad")

    inverse_pl_subset = inverse_pl.with_columns(
        ABS_DIF_AF=abs((pl.col("AF_gnomad") - (1 - pl.col(col_map.eaf))))
    ).select(col_map.variant_id, col_map.eaf, "ABS_DIF_AF", "AF_gnomad", "REF_gnomad", "ALT_gnomad")

    joined = (
        exact_pl_subset.join(inverse_pl_subset, on=col_map.variant_id, how="inner")
        .filter((pl.col("ABS_DIF_AF") > pl.col("ABS_DIF_AF_right")))
        .select(col_map.variant_id)
    )

    print(
        f"Total aligned non-palindromic complementary variants with method 'inverse_match': {len(joined)}"
    )

    # Remove variants from exact match that are more likely to be an inverse match
    exact_pl = exact_pl.filter(
        ~(pl.col(col_map.variant_id).is_in(joined[col_map.variant_id]))
    ).with_columns(test=pl.col('CHR_gnomad').cast(str)+':'+pl.col('POS_gnomad').cast(str)+':'+pl.col('REF_gnomad')+':'+pl.col('ALT_gnomad'))

    # Propery format inverse.aligned_and_merged
    inverse_pl = inverse_pl.filter(
        (~(pl.col(col_map.variant_id).is_in(exact_pl[col_map.variant_id])))
        | (pl.col(col_map.variant_id).is_in(joined[col_map.variant_id]))
    ).with_columns(test=pl.col('CHR_gnomad').cast(str)+':'+pl.col('POS_gnomad').cast(str)+':'+pl.col('REF_gnomad')+':'+pl.col('ALT_gnomad'))

    # Handle Potential Duplicates
    print(len(set(inverse_pl['test']) & set(exact_pl['test'])))
    results = namedtuple("results", ["exact", "inverse"])
    return results(exact_pl, inverse_pl)


if __name__ == "__main__":
    defopt.run(harmonize)

