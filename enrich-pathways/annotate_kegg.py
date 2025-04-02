"""
Make Kegg results easier to interpret by mapping EntrezID number back to gene names
"""

from pathlib import Path

import defopt
import polars as pl

# pylint: disable=R0914, R0913, R0903, C0301


def match_kegg(*, kegg_results: str, kegg_map: str, output: str) -> None:
    """
    :param kegg_results: Output Kegg TSV
    :param kegg_map: Map between Kegg EntrezID's and Gene Names
    """
    print(kegg_results)
    print(kegg_map)
    print(output)
    map_df = pl.read_csv(kegg_map, separator="\t")
    kegg = pl.read_csv(kegg_results, separator="\t")
    kegg = kegg.with_columns(
        pl.col("geneID")
        .map_elements(lambda arr: map_back(arr, map_df), return_dtype=pl.String)
        .alias("geneID")
    )
    kegg.write_csv(output, separator="\t")


def map_back(entrez: str, map_df: pl.DataFrame) -> str:
    """
    Helper to map IDs to gene names
    """
    entrez_list = entrez.split("/")
    entrez_mapped = []
    for entrez_id in entrez_list:
        symbol = map_df.filter(pl.col("ENTREZID").cast(str) == entrez_id)[
            "SYMBOL"
        ].to_list()[0]
        entrez_mapped.append(symbol)
    final_str = "/".join(entrez_mapped)
    return final_str


if __name__ == "__main__":
    defopt.run(match_kegg)
