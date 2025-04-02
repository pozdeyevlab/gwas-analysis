"""
This module is used to summarize all lead variants into one complete non-overlapping list. For each phenotype, by chromosome we start with the MOST significant variant, then search +/- 500kb if no other variant is found great, otherwise that variant is kept but flagged as 'FOUND in other region'
"""

import sys
from pathlib import Path
from typing import List

import defopt
import numpy as np
import polars as pl

# pylint: disable=R0914, R0913, R0903, C0301


def sig_loci(
    *,
    files: str,
    ancestries: str,
    sexes: str,
    output_most_sig: str,
    output_count: str,
    output_all_ancestries: str,
) -> pl.DataFrame:
    """
    :param files: Input path to regenie or saige summary stats
    :param output_file: Path to file to write unusable variants to
    """
    ################################################################################################
    # PART ONE:
    # Load Data
    ################################################################################################
    id_col = "ID"
    p_value = "p-value"

    # Loop through phenotypes
    dfs = []
    for file in list(zip(files.split(","), sexes.split(","), ancestries.split(","))):
        ancestry = file[2]
        sex = file[1]
        df = pl.read_csv(file[0], separator="\t").with_columns(
            pl.lit(f"{ancestry}_{sex}").alias("gwas")
        )
        cols_to_keep = [col for col in df.columns if "missingness" not in col]
        df = df.select(cols_to_keep)
        dfs.append(df)
    sig_data = (
        pl.concat(dfs, how="align")
        .rename({"temp_ID": "ID"})
        .unique(subset=["ID", "p-value", "beta (ALT)"], keep="first")
    )
    chroms = set(list(sig_data["CHR"]))
    print(chroms)

    regions_dfs = []
    for chrom in chroms:
        regions_df = _find_regions(
            raw_df=sig_data.filter(CHR=chrom),
            pos_col="POS(hg38)",
            chr_col="CHR",
            id_col=id_col,
            flank=1000000,
        )
        regions_dfs.append(regions_df)
    merged_df = pl.concat(regions_dfs, how="align")
    print(merged_df)
    # Add regions to significant metal results
    merged_df = merged_df.with_columns((pl.col("chr").cast(int))).rename({"chr": "CHR"})
    final_data = sig_data.join(merged_df, on="CHR", how="inner").filter(
        (pl.col("POS(hg38)") >= pl.col("start_region"))
        & (pl.col("POS(hg38)") <= pl.col("end_region"))
    )
    print(final_data)
    dfs = []
    temp = final_data.group_by(["start_region", "end_region"])
    for t in temp:
        most_sig = t[1].select("p-value").min()
        df = t[1].with_columns(
            pl.when(pl.col("p-value") == most_sig)
            .then(pl.lit("most_significant"))
            .otherwise(pl.lit("tagged"))
            .alias("STATUS")
        )
        dfs.append(df)

    final_data = pl.concat(dfs, how="vertical").rename({"gwas": "Meta-analysis"})

    # Sort
    final_data = final_data.with_columns((pl.col("CHR").cast(int)).alias("CHR"))
    final_data = final_data.with_columns(
        (pl.col("POS(hg38)").cast(int)).alias("POS(hg38)")
    )
    final_data = final_data.sort(["CHR", "POS(hg38)"])

    # Set the column names
    final_data = final_data.select(
        [
            "Endpoint",
            "Meta-analysis",
            "Category",
            "ID",
            "CHR",
            "POS(hg38)",
            "REF",
            "ALT(effect)",
            "rsID",
            "Func.refGene",
            "Gene.refGene",
            "GeneDetail.refGene",
            "ExonicFunc.refGene",
            "allele frequency (ALT)",
            "beta (ALT)",
            "se",
            "p-value",
            "p-value for heterogeneity (Cochran’s Q)",
            "Ancestry*",
            "Number of data sets (per biobank per ancestry)",
            "direction_ALT",
            "number of biobanks",
            "number biobanks with significant association (p-value < 5E-8)",
            "minimum p-value among all biobanks",
            "QC_strand_Flip (flag to indicate whether this variant had a strand flip in any data set)",
            "QC_Freq diff than gnomAD (flag to indicate whether this variant failed the AF QC (compared to gnomAD) in any data set)",
            "PUBMEDID",
            "FIRST AUTHOR",
            "MAPPED_GENE",
            "SNP_GENE_IDS",
            "SNPS",
            "CATALOG_PVALUE",
            "STATUS",
        ]
    )

    final_data.filter(pl.col("STATUS") == "most_significant").sort(
        by=["CHR", "POS(hg38)"]
    ).write_csv(output_most_sig, separator="\t")

    final_data.write_csv(output_all_ancestries, separator="\t")

    final_data.group_by("STATUS", "Category", "Meta-analysis").len().write_csv(
        output_count, separator="\t"
    )


def _find_regions(
    *,
    raw_df: pl.DataFrame,
    pos_col: str,
    chr_col: str,
    id_col: str,
    flank: int = 500000,
) -> pl.DataFrame:
    # Filter data frame for significant variants only & sort by chr and pos
    raw_df = (
        raw_df.sort([chr_col, pos_col])
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
    # print(f"chr{chrom} has {len(extracted_start)} distinct regions")
    return pl.DataFrame(
        {
            "start_region": extracted_start - 500000,
            "end_region": extracted_end + 500000,
            "chr": np.array([str(chrom)] * len(extracted_start)),
            "left_most_sig_var": starts,
            "right_most_sig_var": ends,
        }
    )


if __name__ == "__main__":
    defopt.run(sig_loci)
