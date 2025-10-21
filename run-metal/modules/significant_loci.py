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
    ld_file: str,
    output_path: Path,
    significance_threshold: float,
    phenotype_catalog_path: Path,
) -> pl.DataFrame:
    """
    :param ld_file: LD Clump file
    :param input_file: Metal result
    :param significance_threshold: Desired significance_threshold for filtering significant variants
    :param output_path: Path to file to write unusable variants to
    :param phenotype_catalog_path Path to gwas catalog of published associated variants
    """
    ################################################################################################
    # PART ONE:
    # Find the +/- 500kb regions
    ################################################################################################
    id_col = "STUDY_ID"

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
    # Read in all chromosomes for each biobank, phenptype, and sex
    ld = pl.read_csv(ld_file, separator='\t', null_values=['.'], dtypes={'P':str}).rename({'#CHROM':'chr', 'POS':'pos'}).with_columns((pl.col("chr").str.replace_all("X", "23")).alias("chr")).with_columns(pl.col('ID').str.replace_all('_', ':').str.replace('chr', '').alias('ID'))
    ld = ld.with_columns(pl.col('SP2').fill_null(pl.col('ID').str.replace_all(':', '_')).alias('SP2'))
    
    #ld = pl.read_csv(ld_file, separator='\t', dtypes={'P':str}).rename({'#CHROM':'chr', 'POS':'pos'}).with_columns((pl.col("chr").str.replace_all("X", "23")).alias("chr")).with_columns(pl.col('ID').str.replace_all('_', ':').str.replace('chr', '').alias('ID'))
    ld = ld.with_columns(pl.col('SP2').str.split(',').alias('tmp'))
    ld = ld.with_columns(pl.col('tmp').map_elements(lambda arr: int(min([item.split('_')[1] for item in arr]))).alias('start_region'))

    ld = ld.with_columns(pl.col('tmp').map_elements(lambda arr: int(max([item.split('_')[1] for item in arr]))).alias('end_region')).drop('tmp')
    print(ld)
    # Read in catalog data
    catalog = pl.read_csv(
        phenotype_catalog_path,
        null_values=["NA", "NR"],
        separator="\t",
        dtypes={'CHR_ID':str}
    ).filter(pl.col('P-VALUE')<5e-8).with_columns(
        (pl.col("CHR_ID").str.replace_all("X", "23").cast(int)).alias("CHR_ID")
    )

    # Add known or potentially novel to final data
    merged_chroms = set(ld["chr"])
    merged_dfs = []
    for chr in merged_chroms:
        df = ld.filter(pl.col("chr").cast(int) == int(chr))

        chr_catalog = catalog.filter(pl.col("CHR_ID") == int(chr))
        df = df.with_columns(
            (
                pl.when(
                    (
                        sum(
                            [
                                (
                                    (500000>(abs(pl.col("pos").cast(int) - i)))
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
    merged_df = merged_df.rename({'ID':'STUDY_ID'})

    print(
        f"\n{merged_df.shape[0]} unique regions, based on {sig_data.shape[0]} variants with a pvalue less than 5e-8"
    )
    # Add regions to significant metal results
    final_data = sig_data.join(merged_df, on="STUDY_ID", how="inner")
    final_data = final_data.drop(["tmp", "counts", 'P-value']).rename({'P':'P-value'})
    
    #final_data = final_data.select('chr','start_region','end_region','P-value','STUDY_ID','Allele1','Allele2','Freq1','Effect','StdErr','Direction','per_variant_N','pos','ref_from_id','alt_from_id','pos_right','SP2','known_in_gwas_catalog')
    ################################################################################################
    # PART TWO:
    # Look for variants in the catalog that are within the boundaries of the found regions
    ################################################################################################
    catalog = (
        pl.read_csv(
            phenotype_catalog_path,
            separator="\t",
            null_values=["NA"],
            dtypes={'CHR_ID':str},
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
                "DISEASE/TRAIT"
            ],
        )
        .rename({"CHR_ID": "chr", "P-VALUE": "CATALOG_PVALUE"})
        .with_columns((pl.col("chr").str.replace_all("X", "23")).alias("chr"))
    ).filter(pl.col('CATALOG_PVALUE')<5e-8)

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
            "DISEASE/TRAIT"
        ]
    )

    annotated_catalog = annotated_catalog.group_by("ID", maintain_order=True).agg(
        **{col: pl.col(col) for col in annotated_catalog.columns if col != "ID"}
    )

    cols = annotated_catalog.drop("ID").columns
    for col in cols:
        print(col)
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

    final_result = final_data.join(annotated_catalog, on="ID", how = 'left').with_columns(
        pl.when(pl.col("start_region") < 0)
        .then(1)
        .otherwise(pl.col("start_region"))
        .alias("start_region")
    )

    final_result.write_csv(output_path, separator="\t")


def _list_to_str(lst) -> str:
    lst = [f"{i}" for i in lst]
    lst = list(lst)
    return ";".join(lst)


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


