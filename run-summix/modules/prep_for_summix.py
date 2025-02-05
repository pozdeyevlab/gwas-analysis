"""
Module to format gnomad aligned summary stats for SUMMIX
Removes variants that:
1) Have a gnomAD AN warning
2) Have significantly different AF compared to gnomAD (based on mahalanobis distance)
3) Have MAC (minor allele count less than 20)
4) gnomAD variant does not PASS QC
"""

from pathlib import Path
from typing import Optional

import defopt
import polars as pl

# pylint: disable=R0914, R0913, R0903, C0301


def prep_qc_files(
    *,
    input_file: str,
    ea: str,
    non_ea: str,
    se: str,
    output_path: Path,
    chrom_col: str,
    pos_col: str,
    case_count: int,
    control_count: int,
    pval_col: str,
    mac_filter: int,
    case_count_col: Optional[str],
    control_count_col: Optional[str],
    total_count_col: Optional[str]
) -> pl.DataFrame:
    """
    :param input_files: Input path to regenie or saige summary stats
    :param gwas_software: Software used to generate summary stats
    :param ea: Column name of effect allele
    :param non_ea: Column name of non-effect allele
    :param se: Column name of se column
    :param output_path: Path to file to write unusable variants to
    :param chrom_col: Column name of variant if (Marker or ID)
    :param pos_col: Path to file to write unusable variants to
    :param case_count: Total number of people
    :param control_count: Total number of people
    :param pval_col: Column of original p-value
    :param mac_filter: Filter for minor allele frequency
    :param case_count_col: Case count column
    :param control_count_col: Control count col
    :param total_count_col: Total count col
    """
    # Read in all chromosomes for each biobank, phenptype, and sex
    print(input_file)
    total_n = case_count + control_count
    combined_df = _read_specific_columns(
        file=input_file,
        chrom=chrom_col,
        pos=pos_col,
        ref=non_ea,
        alt=ea,
        total_n=total_n,
        mac_filter=mac_filter,
        pval_col=pval_col,
        case_count_col=case_count_col,
        control_count_col=control_count_col,
        total_count_col=total_count_col,
    )

    # handle log10 pvalues:
    if pval_col == "LOG10P":
        combined_df = combined_df.with_columns(
            (10 ** (-1 * pl.col(pval_col).cast(float))).alias("pval")
        )
    else:
        combined_df = combined_df.rename({pval_col: "pval"})

    # Write output
    combined_df = combined_df.drop("beta")
    combined_df = combined_df.rename({se: "se", "Aligned_Beta": "beta"})
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


def _read_specific_columns(
    *,
    file: str,
    chrom: str,
    pos: str,
    ref: str,
    alt: str,
    total_n: int,
    mac_filter: int,
    pval_col: str,
    case_count_col: Optional[str],
    control_count_col: Optional[str],
    total_count_col: Optional[str]
) -> pl.DataFrame:
    """
    For each file read in specific columns and combine
    """
    # Read & Format
    old_chrom = chrom
    chrom = "CHR_gnomad"
    df = (
        pl.read_csv(
            file,
            separator="\t",
            dtypes={old_chrom: str, pos: int, ref: str, alt: str, chrom: str},
        )
        .with_columns(CHR=pl.col(chrom).str.replace("chrom", "").replace("X", "23"))
        .filter((pl.col("FILTER") == "PASS"))
        .filter((pl.col("GNOMAD_AN_Flag") == 0))
        .filter((pl.col("outlier_stdev") == "No"))
    )
    if (case_count_col != "nan") and (control_count_col != "nan"):
        df = df.with_columns(
            (pl.col(case_count_col) + pl.col(control_count_col)).alias("per_variant_N")
        ).with_columns(pl.col("per_variant_N").fill_null(total_n))

    if total_count_col != "nan":
        df = df.rename({total_count_col: "per_variant_N"}).with_columns(
            pl.col("per_variant_N").fill_null(total_n)
        )
    if (
        (total_count_col == "nan")
        and (case_count_col == "nan")
        and (control_count_col == "nan")
    ):
        df = df.with_columns(per_variant_N=total_n)

    df = (
        df.with_columns(
            (pl.when(pl.col("Aligned_AF").cast(float) > 0.5))
            .then(
                ((1 - pl.col("Aligned_AF").cast(float)) * pl.col("per_variant_N")) * 2
            )
            .otherwise((pl.col("Aligned_AF") * pl.col("per_variant_N")) * 2)
            .alias("MAC")
        )
        .filter(pl.col("MAC") >= mac_filter)
        .drop(["FILTER", "BETA"])
    )

    # Filter for MAF > 0.0005
    df = df.filter(
        (pl.col("Aligned_AF").cast(float) > 0.0005) & (pl.col("Aligned_AF") < 0.9995)
    )

    deduplicated_df = df.sort("STUDY_ID").unique(
        keep="first", maintain_order=True, subset=["STUDY_ID"]
    )

    # Calculate missingness
    deduplicated_df = deduplicated_df.with_columns(
        ((total_n - pl.col("per_variant_N")) / total_n).alias("missingness")
    )
    print(deduplicated_df)
    print(deduplicated_df.columns)

    # Remove Pvals that are nan
    deduplicated_df = deduplicated_df.filter(pl.col(pval_col).is_not_null())

    # Remove variants with potential strand flip
    deduplicated_df = deduplicated_df.with_columns(
        pl.col("Potential_Strand_Flip").fill_null("PASS")
    ).filter(pl.col("Potential_Strand_Flip") == "PASS")

    return deduplicated_df


if __name__ == "__main__":
    defopt.run(prep_qc_files)
