"""
Match fusion results to gwas
"""

import defopt
import polars as pl


def match_predixscan(
    *,
    summary: pl.DataFrame,
    predixscan=pl.DataFrame,
) -> pl.DataFrame:
    """
    Args:
        --
        --
    Returns:
        --
    """
    chrs = summary["CHR"].unique()
    matches = []
    pathway_genes = []
    for chr in chrs:
        first_fusion = predixscan.filter(pl.col("Chromosome").cast(int) == chr)

        tmp_summary = summary

        positions = tmp_summary["POS(hg38)"]
        for position in positions:
            tmp_fusion = first_fusion.with_columns(
                pl.when(abs(pl.col("Start") - position) < 1000000)
                .then(1)
                .otherwise(0)
                .alias("potential_match")
            )
            tmp_fusion = tmp_fusion.filter(pl.col("potential_match") == 1)

            if not tmp_fusion.is_empty():

                # Create a ';' separated list of significant genes within 1Mb of the lead variant
                genes = tmp_fusion["gene_name"].to_list()
                pvals = tmp_fusion["pvalue"].to_list()

                lead_gene = tmp_summary.filter(pl.col("POS(hg38)") == position)[
                    "Gene.refGene"
                ].to_list()[0]

                lead_gene_pval = tmp_summary.filter(pl.col("POS(hg38)") == position)[
                    "p-value"
                ].to_list()[0]

                # Format genes for addition to summary file
                matched_df = tmp_summary.filter(pl.col("POS(hg38)") == position)
                added = f"{','.join(list(set(genes)))}".replace(" ", "")
                s = pl.Series("Gene.PrediXscan", [added])
                matches.append(matched_df.insert_column(matched_df.shape[1], s))

                # Combine lists for pathway annotation
                genes.append(lead_gene)
                pvals.append(lead_gene_pval)
                for i in list(zip(genes, pvals)):
                    pathway_genes.append({"gene": i[0], "pvals": i[1]})

            if tmp_fusion.is_empty():
                # Combine lists for pathway annotation
                lead_gene = tmp_summary.filter(pl.col("POS(hg38)") == position)[
                    "Gene.refGene"
                ].to_list()[0]
                lead_gene_pval = tmp_summary.filter(pl.col("POS(hg38)") == position)[
                    "p-value"
                ].to_list()[0]
                pathway_genes.append({"gene": lead_gene, "pvals": lead_gene_pval})

                # Format genes for addition to summary file
                matched_df = tmp_summary.filter(pl.col("POS(hg38)") == position)
                empty_column = pl.Series("Gene.PrediXscan", [None] * matched_df.height)
                matched_df.insert_column(matched_df.shape[1], empty_column)
                matched_df = matched_df.with_columns(
                    pl.col("Gene.PrediXscan").cast(str).alias("Gene.PrediXscan")
                )
                matches.append(matched_df)

    final_matched = pl.concat(matches, how="vertical")

    pathways = pl.DataFrame(pathway_genes).sort("pvals").group_by("gene").first()

    # pathways.write_csv(pathway_output, separator="\t")

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
    return final_matched.select(order)


if __name__ == "__main__":
    defopt.run(match_predixscan)
