"""
Combine fusion results and annotate with gene names 
"""

import defopt
import polars as pl
import pyranges as pr
from typing import List


def gene_names(*, gtf: str, fusion_files: List[str], output_file: str) -> pl.DataFrame:
    """
    :param gtf: GTF File
    :param fusion_files: FUSION results
    :param output_file: String output
    """

    # Set up Ensembl code
    gr = pr.read_gtf(gtf)
    df = gr.df
    df2 = pl.from_pandas(df)
    df3 = df2.select(
        "Chromosome", "Start", "End", "gene_name", "gene_id", "gene_version", "Feature"
    ).filter(pl.col("Feature") == "start_codon")

    df2 = (
        df2.select(
            "Chromosome",
            "Start",
            "End",
            "gene_name",
            "gene_id",
            "gene_version",
            "Feature",
        )
        .filter(~pl.col("gene_name").is_in(df3["gene_name"]))
        .filter(pl.col("Feature") == "gene")
    )
    df3 = pl.concat([df3, df2], how="vertical")

    # Filter results to exclude unsuccesful TWAS
    fusion = []
    for file in fusion_files:
        tmp = pl.read_csv(file, separator = "\t")
        fusion.append(tmp)
    
    fusion_results = pl.concat(fusion, how = 'vertical')
    fusion_results = fusion_results.filter(~pl.col('TWAS.P').is_null())

    # Separate genes and variants
    fusion_results = fusion_results.with_columns(pl.col('ID').str.split('.').map_elements(lambda arr: arr[0]).alias('gene_id'))

    fusion_results = fusion_results.with_columns(pl.col('ID').str.split('.').map_elements(lambda arr: arr[1]).alias('gene_version'))

    # Merge to get gene name
    final_results = fusion_results.merge(pl.from_pandas(df).select("gene_id","gene_name"), on = 'gene_id', how = "left")
    final_results.write_csv(output_file, separator='\t')



if __name__ == "__main__":
    defopt.run(gene_names)
