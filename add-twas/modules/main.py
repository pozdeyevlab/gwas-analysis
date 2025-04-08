"""
Main script to organize combining FUSION & PrediXscan
"""

import re
from pathlib import Path
from typing import List

import defopt
import match_fusion
import match_predixscan
import polars as pl
import prep_predixscan


def main(
    *,
    gwas: Path,
    fusions: List[str],
    predixscans: List[str],
    ancestries: List[str],
    gtf: Path,
    output: str,
) -> None:
    """
    :param gwas: GWAS-Meta Summary file
    :param fusions: Comma separated list of fusion files (by ancestry)
    :param predixscans: Comma separated list of predixscan files (by ancestry)
    :param ancestries: Comma separared list of ancestries
    :param gtf: Emsembl gtf file
    """
    # For every ancestry add significant genes then combine results
    final_files = []
    gwas_df = pl.read_csv(gwas, separator="\t", null_values="NA")
    for ancestry in ancestries:
        # Filter gwas results
        tmp = gwas_df.filter(pl.col("Meta-analysis").str.contains(ancestry))

        # Prep Fusion
        try:
            pattern = rf".*?{re.escape(ancestry)}.*"
            fusion = [
                re.match(pattern, s, re.IGNORECASE).group()
                for s in fusions
                if re.match(pattern, s, re.IGNORECASE) is not None
            ][0]

            fusion = pl.read_csv(fusion, separator="\t")

            fusion = fusion.filter(pl.col("TWAS.P") < (0.05 / fusion.shape[0]))
            fusion = fusion.filter(~(pl.col("V3").str.contains("gene_source")))

            # Combine results
            gwas_with_fusion = match_fusion.match_fusion(
                summary=tmp,
                fusion=fusion,
            )
            # empty_column = pl.Series("Gene.PrediXscan", [None] * gwas_with_fusion.height)
            # gwas_with_fusion.insert_column(gwas_with_fusion.shape[1], empty_column)
            final_files.append(gwas_with_fusion)

        except IndexError:
            print(
                f"Ancestry: {ancestry} does not have a fusion input. If this is incorrect please check the config file"
            )

        # Prep predixscan
        try:
            pattern = rf".*?{re.escape(ancestry)}.*"
            pred_file = [
                re.match(pattern, s, re.IGNORECASE).group()
                for s in predixscans
                if re.match(pattern, s, re.IGNORECASE) is not None
            ][0]
            # pred_file = [s for s in predixscans if ancestry in s][0]
            predi = prep_predixscan.prep(gtf=gtf, predixscan_file=pred_file)

            gwas_with_predixscan = match_predixscan.match_predixscan(
                summary=tmp, predixscan=predi
            )

            # empty_column = pl.Series("Gene.FUSION", [None] * gwas_with_predixscan.height)
            # gwas_with_predixscan.insert_column(gwas_with_predixscan.shape[1], empty_column)
            final_files.append(gwas_with_predixscan)

        except IndexError:
            print(
                f"Ancestry: {ancestry} does not have a PrediXscan input. If this is incorrect please check the config file"
            )

    # Remove ancestries from starting file and add new genes
    for ancestry in ancestries:
        original = gwas_df.filter(~(pl.col("Meta-analysis").str.contains(ancestry)))

    # Add dummy columns to gwas_df
    empty_column = pl.Series("Gene.FUSION", [None] * original.height)
    original.insert_column(original.shape[1], empty_column)
    empty_column = pl.Series("Gene.PrediXscan", [None] * original.height)
    original.insert_column(original.shape[1], empty_column)

    order = [
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
        "Gene.FUSION",
        "Gene.PrediXscan",
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
    formatted = pl.concat(final_files, how="align").select(order)

    results = pl.concat([formatted, original.select(order)], how="vertical").sort(
        by=["CHR", "POS(hg38)"]
    )
    results.write_csv(output, separator="\t")


if __name__ == "__main__":
    defopt.run(main)
