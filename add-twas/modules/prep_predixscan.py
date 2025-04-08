"""
Goal: Haoyu provided prediscan results separatred by ancestry for all 6 phenotyes. 
We need to get an estimated position for these genes using the ensembl gtf data. 
For all 6 phenotypes 
"""

import defopt
import polars as pl
import pyranges as pr


def prep(*, gtf: str, predixscan_file: str) -> pl.DataFrame:
    """
    :param gtf: GTF File
    :param predixscan_file: Predixscan output
    """

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

    predi = pl.read_csv(predixscan_file, separator=",", null_values="NA")
    predi = predi.filter(pl.col("pvalue") < (0.05 / predi.shape[0]))
    predi = predi.join(df3, on="gene_name", how="left").unique()
    return predi


if __name__ == "__main__":
    defopt.run(prep)
