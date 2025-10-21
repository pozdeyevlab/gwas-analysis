"""
This module takes output from METAL and creates a list of genome wide significant variants that do not overlap with each other in 500kb increments. This methodology was inspired by the methods used in the GBMI I flagship paper. 

For each chromosome variants are sorted by position, then starting with the most significant variant (lowest p-value) variants that are within +/-500kb of the most significant variant are eliminated. This is done until all genome wide significant variants are either marked as being within a 500kb range of another significant variant OR are in fact the most significant variant in said region.
"""

from pathlib import Path

import defopt
import numpy as np
import polars as pl

import cochrans

# pylint: disable=R0914, R0913, R0903, C0301


def sig_loci(
    *,
    input_file: str,
    output_path: Path,
    significance_threshold: float,
    phenotype_catalog_path: Path,
) -> pl.DataFrame:
    """
    :param input_file: Input path to regenie or saige summary stats
    :param significance_threshold: Desired significance_threshold for filtering significant variants
    :param output_path: Path to file to write unusable variants to
    :param phenotype_catalog_path Path to gwas catalog of published associated variants
    """
    ################################################################################################
    # PART ONE:
    # Find the +/- 500kb regions
    ################################################################################################
    id_col = "STUDY_ID"
    p_value = "P-value"
    # Read in all chromosomes for each biobank, phenptype, and sex
    df = _read_input(
        file=input_file,
        id_col=id_col,
    ).with_columns((pl.col("chr").str.replace_all("X", "23")).alias("chr"))
    sig_data = df.filter((pl.col("P-value").cast(float)) <= (significance_threshold))
    n_sig_variants = sig_data.shape[0]
    sig_data = cochrans.cochran_q(df=sig_data)
    sig_data = sig_data.with_columns(
        (pl.col("cochran_q_p_value") < (1 / n_sig_variants)).alias(
            "cochran_heterogeneity"
        )
    )

    chroms = set(list(sig_data["chr"]))

    regions_dfs = []

    for chrom in chroms:
        regions_df = _find_regions(
            sig_threshold=significance_threshold,
            raw_df=sig_data.filter(chr=chrom),
            pval_col=p_value,
            pos_col="pos",
            chr_col="chr",
            id_col=id_col,
            flank=1000000,
        )
        regions_dfs.append(regions_df)
    merged_df = pl.concat(regions_dfs, how="vertical")

    # Read in catalog data
    catalog = pl.read_csv(
        phenotype_catalog_path,
        null_values=["NA", "NR"],
        separator="\t",
    ).with_columns(
        (pl.col("CHR_ID").str.replace_all("X", "23").cast(int)).alias("CHR_ID")
    )

    # Add known or potentially novel to final data
    merged_chroms = set(merged_df["chr"])
    merged_dfs = []
    for chr in merged_chroms:
        df = merged_df.filter(pl.col("chr").cast(int) == int(chr))
        chr_catalog = catalog.filter(pl.col("CHR_ID") == int(chr))
        df = df.with_columns(
            (
                pl.when(
                    (
                        sum(
                            [
                                (
                                    (pl.col("start_region").cast(int) < i)
                                    & (pl.col("end_region").cast(int) > i)
                                )
                                for i in chr_catalog["CHR_POS"]
                            ]
                        )
                        >= 1
                    )
                )
                .then(pl.lit("known"))
                .otherwise(pl.lit("potentially_novel"))
                .alias("known_in_gwas_catalog")
            )
        )
        merged_dfs.append(df)
    merged_df = pl.concat(merged_dfs, how="vertical")

    print(
        f"\n{merged_df.shape[0]} unique regions, based on {sig_data.shape[0]} variants with a pvalue less than 5e-8"
    )

    # Add regions to significant metal results
    final_data = (
        sig_data.join(merged_df, on="chr", how="inner")
        .filter(
            (pl.col("pos") >= pl.col("start_region"))
            & (pl.col("pos") <= pl.col("end_region"))
        )
        .drop(["tmp", "counts"])
    )

    # Filter to retain most significant variant per region
    lowest_pval = final_data.group_by(
        ["chr", "start_region", "end_region"], maintain_order=True
    ).agg(pl.col('P-value').min())
    print(lowest_pval)
    final_data = lowest_pval.join(
        final_data, on=["chr", "start_region", "end_region", "P-value"]
    ).unique(subset = ["chr", "start_region", "end_region", "P-value"], keep="first", maintain_order=True)

    ################################################################################################
    # PART TWO:
    # Look for variants in the catalog that are within the boundaries of the found regions
    ################################################################################################
    catalog = (
        pl.read_csv(
            phenotype_catalog_path,
            separator="\t",
            null_values=["NA"],
            columns=[
                "CHR_ID",
                "CHR_POS",
                "MAPPED_GENE",
                "P-VALUE",
                "SNP_GENE_IDS",
                "SNPS",
                "PUBMEDID",
                "FIRST AUTHOR",
                "LINK",
            ],
        )
        .rename({"CHR_ID": "chr", "P-VALUE": "CATALOG_PVALUE"})
        .with_columns((pl.col("chr").str.replace_all("X", "23")).alias("chr"))
    )

    ### Add Significant Data ###

    annotated_catalog = (
        (merged_df.join(catalog, on="chr", how="inner"))
        .with_columns(
            pl.when(
                (pl.col("CHR_POS").cast(int) > pl.col("start_region"))
                & (pl.col("CHR_POS").cast(int) < pl.col("end_region"))
            )
            .then(pl.lit("known"))
            .otherwise(pl.lit("unknown"))
            .alias("test")
        )
        .filter(pl.col("test") == "known")
        .drop("test")
        .with_columns(
            (
                pl.col("chr").cast(str)
                + "-"
                + pl.col("start_region").cast(str)
                + "-"
                + pl.col("end_region").cast(str)
            ).alias("ID")
        )
    ).select(
        [
            "PUBMEDID",
            "FIRST AUTHOR",
            "LINK",
            "MAPPED_GENE",
            "SNP_GENE_IDS",
            "SNPS",
            "CATALOG_PVALUE",
            "ID",
        ]
    )

    annotated_catalog = annotated_catalog.group_by("ID", maintain_order=True).agg(
        **{col: pl.col(col) for col in annotated_catalog.columns if col != "ID"}
    )

    cols = annotated_catalog.drop("ID").columns
    for col in cols:
        annotated_catalog = (
            annotated_catalog.with_columns(pl.col(col).list.unique().alias(col))
            .with_columns(
                col,
                annotated_catalog[col]
                .map_elements(lambda x: _list_to_str(x), skip_nulls=False)
                .alias(f"second_{col}"),
            )
            .drop(col)
        )
        annotated_catalog = annotated_catalog.rename({f"second_{col}": col})

    final_data = final_data.with_columns(
        (
            pl.col("chr").cast(str)
            + "-"
            + pl.col("start_region").cast(str)
            + "-"
            + pl.col("end_region").cast(str)
        ).alias("ID")
    )

    final_result = pl.concat([final_data, annotated_catalog], how="align").with_columns(
        pl.when(pl.col("start_region") < 0)
        .then(1)
        .otherwise(pl.col("start_region"))
        .alias("start_region")
    )
    #final_result = final_result.drop('P-value').rename({'str_pval':'P-value'})
    final_result.write_csv(output_path, separator="\t")


def _list_to_str(lst) -> str:
    lst = [f"{i}" for i in lst]
    lst = list(set(lst))
    return ";".join(lst)


def _find_regions(
    *,
    sig_threshold: float,
    raw_df: pl.DataFrame,
    pval_col: str,
    pos_col: str,
    chr_col: str,
    id_col: str,
    flank: int = 500000,
) -> pl.DataFrame:
    # Filter data frame for significant variants only & sort by chr and pos
    raw_df = (
        raw_df.filter(pl.col(pval_col).cast(float) < sig_threshold)
        .sort([chr_col, pos_col])
        .with_columns(pl.col(chr_col).cast(int).alias(chr_col))
        .with_columns(pl.col(pos_col).cast(int).alias(pos_col))
    )

    # Calcualte differences for all chr and pos in df
    diff_chrom: np.ndarray = np.diff(raw_df[chr_col])

    diff_pos: np.ndarray = np.diff(raw_df[pos_col])
    # start = c(1, which(diff(y) != 0 | diff(x) <= (flank | diff(x) >= flank) + 1)
    # end = c(start - 1, length(x))

    # Define conditions
    conditions: np.ndarray = (
        (diff_chrom != 0) | (diff_pos <= (-1 * flank)) | (diff_pos >= flank)
    )
    conditions = np.append(True, conditions)

    start_indices = np.where(conditions)[0]

    end_indices = np.append(start_indices - 1, len(conditions) - 1)[1:]

    starts = raw_df[id_col].gather(start_indices)

    ends = raw_df[id_col].gather(end_indices)

    # Assumes name column is chr:pos:ref:alt
    extracted_start = np.array([int(x.split(":")[1]) for x in starts])
    extracted_end = np.array([int(x.split(":")[1]) for x in ends])

    # Print Summary
    chrom = raw_df[chr_col][0]
    print(f"chr{chrom} has {len(extracted_start)} distinct regions")
    return pl.DataFrame(
        {
            "start_region": extracted_start - 500000,
            "end_region": extracted_end + 500000,
            "chr": np.array([str(chrom)] * len(extracted_start)),
            "left_most_sig_var": starts,
            "right_most_sig_var": ends,
        }
    )


def _read_input(*, file: str, id_col: str) -> pl.DataFrame:
    """
    For each file read in specific columns to avoid mem load
    """
    # Read & Format
    raw_df = pl.read_csv(file, separator="\t", dtypes={'P-value':float}) #.with_columns(pl.col('P-value').cast(float).alias('float_pval')).rename({'P-value':'str_pval', 'float_pval':'P-value'})

    raw_df = raw_df.with_columns(
        pl.col(id_col).str.split(":").map_elements(lambda arr: arr[0]).alias("chr"),
        pl.col(id_col)
        .str.split(":")
        .map_elements(lambda arr: arr[1])
        .alias("pos")
        .cast(int),
        pl.col(id_col)
        .str.split(":")
        .map_elements(lambda arr: arr[2])
        .alias("ref_from_id"),
        pl.col(id_col)
        .str.split(":")
        .map_elements(lambda arr: arr[3])
        .alias("alt_from_id"),
    )

    return raw_df


if __name__ == "__main__":
    defopt.run(sig_loci)

