"""
Handles Palindromic variants
"""

from collections import namedtuple
from typing import List, Optional

import attr
import defopt
import filter_gwas
import numpy as np
import polars as pl

# pylint: disable=C0301 # line too long
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
    :param col_map: Class of columns
    """
    list_of_results: List[AlignmentResults] = []

    # Check for variants with exact alignment (ref=ref & alt=alt)
    exact: AlignmentResults = align_alleles(
        gnomad_df=gnomad_pl,
        gwas_df=gwas_pl,
        id_column=col_map.variant_id,
        method="exact_match",
    )

    # Inverse Exact
    gwas_pl = _make_id_column(
        new_column_name="Flipped_Allele_ID",
        col_map=col_map,
        first_allele=col_map.effect_allele,
        second_allele=col_map.non_effect_allele,
        polars_df=gwas_pl,
    )
    inverse: AlignmentResults = align_alleles(
        gnomad_df=gnomad_pl,
        gwas_df=gwas_pl,
        id_column="Flipped_Allele_ID",
        method="inverse_match",
    )

    if (inverse.aligned_and_merged is not None) and (
        exact.aligned_and_merged is not None
    ):
        results = handle_palindromic_variants(
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

    # Concatenate results
    list_of_dfs: List[pl.DataFrame] = []
    for result in list_of_results:
        df = result.aligned_and_merged
        if not df.is_empty():
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

    # Print the number of palindromic aligned variants
    exact_align = stacked_pl.filter(pl.col("Alignment_Method") == "exact_match")[
        col_map.palindromic_flag
    ].sum()
    inverse_align = stacked_pl.filter(pl.col("Alignment_Method") == "inverse_match")[
        col_map.palindromic_flag
    ].sum()

    print(
        f"\nPalindromic Summary:\nTotal aligned palindromic variants: {sum(stacked_pl[col_map.palindromic_flag])}/{sum(gwas_pl[col_map.palindromic_flag])}\nTotal aligned palindromic varinats with method 'exact_match': {exact_align}\nTotal aligned palindromic varinats with method 'inverse_match': {inverse_align}\nTotal aligned palindromic varinats with allele frequencies between 0.4 and 0.6: {sum(stacked_pl[col_map.palindromic_af_flag])}"
    )

    # Calculate:
    # Absolute difference in AF
    # Fold change
    # Alternate allele comparison
    # AF-GWAS < 0.4 and AF-gnomAD > 0.6, AF-GWAS > 0.6 and AF-gnomAD < 0.4

    stacked_pl = (
        stacked_pl.with_columns(
            (abs(pl.col("AF_gnomad") - pl.col("Aligned_AF"))).alias("ABS_DIF_AF")
        )
        .with_columns(
            ((pl.col("AF_gnomad") / pl.col("Aligned_AF"))).alias("FOLD_CHANGE_AF")
        )
        .with_columns(
            (
                pl.when(
                    (
                        ((abs((1 - pl.col("AF_gnomad")) - pl.col("Aligned_AF"))))
                        < pl.col("ABS_DIF_AF")
                    )
                )
                .then(1)
                .otherwise(0)
                .alias("GNOMAD_EAF_FLAG")
            )
        )
        .with_columns(
            (
                pl.when(
                    ((pl.col("AF_gnomad") < 0.4) & (pl.col("Aligned_AF") > 0.6))
                    | ((pl.col("AF_gnomad") > 0.6) & (pl.col("Aligned_AF") < 0.4))
                )
                .then(1)
                .otherwise(0)
                .alias("4_6")
            )
        )
        .with_columns(
            (
                pl.when(pl.col("4_6") == 1)
                .then(pl.lit("4_6_af"))
                .otherwise(
                    (
                        pl.when(pl.col("GNOMAD_EAF_FLAG") == 1)
                        .then(pl.lit("alt_af"))
                        .otherwise(
                            pl.when((pl.col("FOLD_CHANGE_AF") > 2))
                            .then(pl.lit("fold_change"))
                            .otherwise(pl.lit("PASS"))
                        )
                    )
                )
                .alias("Potential_Strand_Flip")
            )
        )
        .drop(["4_6", "GNOMAD_EAF_FLAG"])
    )

    fold_change_count = (
        stacked_pl.with_columns(
            pl.when(pl.col("FOLD_CHANGE_AF") > 2).then(1).otherwise(0).alias("x")
        )
        .filter(x=1)
        .shape[0]
    )
    print(
        f"Total aligned palindromic varinats with a fold change greater than 2 (gnomad_af/gwas_af): {fold_change_count}"
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
    Helper function to create id columns with differenc effect and non effect allele combinaitons
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
    gnomad_df = gnomad_df.with_columns(
        (pl.col("ID_gnomad").str.replace("chr", "")).alias("ID_gnomad")
    )

    gwas_df = gwas_df.with_columns(
        (pl.col(id_column).str.replace("chr", "")).alias(id_column)
    )

    # Create sets of id columns and find overlap
    gnomad_set: set = set(gnomad_df["ID_gnomad"])
    gwas_set: set = set(gwas_df[id_column])
    aligned_set: set = gnomad_set & gwas_set

    if len(aligned_set) > 0:
        # Filter gwas data for only aligned variants
        aligned_df: pl.DataFrame = gwas_df.filter(pl.col(id_column).is_in(aligned_set))
        unaligned_df: pl.DataFrame = gwas_df.filter(
            ~pl.col(id_column).is_in(aligned_set)
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


def handle_palindromic_variants(
    *, exact_pl: pl.DataFrame, inverse_pl: pl.DataFrame, col_map: filter_gwas.Columns
) -> namedtuple:
    """
    # TODO add docstring
    Args:
        exact_pl: Data frame with exact matches
        inverse_pl: Data frame with inverse matches
        col_map; Class of column names
    """
    exact_pl_subset = exact_pl.with_columns(
        ABS_DIF_AF=abs((pl.col("AF_gnomad") - pl.col(col_map.eaf)))
    ).select(col_map.variant_id, "ABS_DIF_AF")

    inverse_pl_subset = inverse_pl.with_columns(
        ABS_DIF_AF=abs((pl.col("AF_gnomad") - (1 - pl.col(col_map.eaf))))
    ).select(col_map.variant_id, "ABS_DIF_AF")

    joined = (
        exact_pl_subset.join(inverse_pl_subset, on=col_map.variant_id, how="inner")
        .filter((pl.col("ABS_DIF_AF") > pl.col("ABS_DIF_AF_right")))
        .select(col_map.variant_id)
    )

    print(
        f"Palindromic variants with a lower abs difference in allele frequency (compared to gnomad) when aligned via inverse match vs exact match: {len(joined)}"
    )

    # Remove variants from axact match that are more likely to be an inverse match
    exact_pl = exact_pl.filter(~(pl.col(col_map.variant_id).is_in(joined)))

    # Propery format inverse.aligned_and_merged
    inverse_pl = inverse_pl.filter(
        (~(pl.col(col_map.variant_id).is_in(exact_pl[col_map.variant_id])))
        | (pl.col(col_map.variant_id).is_in(joined))
    )

    results = namedtuple("results", ["exact", "inverse"])
    return results(exact_pl, inverse_pl)


if __name__ == "__main__":
    defopt.run(harmonize)
