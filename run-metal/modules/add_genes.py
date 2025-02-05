"""
This module is the final step in this pipeline. It creates a human readable summary of the potentially novel and already reported variants.

Variants are classified as potentially novel if there are no catalog variants which fall within the start and end region of the variants identified in modules/significant_locy.py

Categories in output:
Endpoint (phenotype)
Category (previously reported)
Pos(hg38)
REF
ALT(effect)
rsID
Func.refGene
Gene.refGene
GeneDetail.refGene
ExonicFunc.refGene
allele frequency (ALT)
beta (ALT)
se
p-value
p-value for heterogeneity (Cochran’s Q)
Ancestry*
Number of data sets (per biobank per ancestry)
direction_ALT
number of biobanks
number biobanks with significant association (p-value < 5E-8)
minimum p-value among all biobanks
QC_strand_Flip (flag to indicate whether this variant had a strand flip in any data set)
QC_Freq diff than gnomAD (flag to indicate whether this variant failed the AF QC (compared to gnomAD) in any data set)
PUBMEDID
FIRST AUTHOR
LINK
MAPPED_GENE
SNP_GENE_IDS
SNPS
CATALOG_PVALUE
Missingness per dataset
"""

# pylint: disable=locally-disabled, multiple-statements, fixme, line-too-long

from pathlib import Path

import defopt
import numpy as np
import polars as pl


def add_genes(
    *,
    significant_loci: Path,
    annovar_genes: Path,
    sig_loci_id_col: str,
    output: Path,
) -> pl.DataFrame:
    """
    :param significant_loci: Path to raw significant loci ({output_dir}/significant_loci/{phenotype}.tsv)
    :param annovar_genes: Path to annovar results ({output_dir}/annovar/{phenotype}.hg38_multianno.txt})
    :param sig_loci_id_col: Name of id col in sig loci tsv
    :param output: Where to write results
    """
    ################################################################################################
    # PART ONE:
    # Read in annovar results and make ID col -> chr:pos:ref:alt
    ################################################################################################
    annovar_pl = (
        pl.read_csv(
            annovar_genes,
            separator="\t",
            columns=[
                "Chr",
                "Start",
                "Ref",
                "Alt",
                "Func.refGene",
                "Gene.refGene",
                "GeneDetail.refGene",
                "ExonicFunc.refGene",
                "avsnp151"
            ],
        )
        .with_columns(
            (
                pl.col("Chr").cast(str)
                + ":"
                + pl.col("Start").cast(str)
                + ":"
                + pl.col("Ref")
                + ":"
                + pl.col("Alt")
            ).alias(sig_loci_id_col)
        )
        .select(
            [
                sig_loci_id_col,
                "Func.refGene",
                "Gene.refGene",
                "GeneDetail.refGene",
                "ExonicFunc.refGene",
                "avsnp151",
            ]
        )
    )
    annovar_pl = annovar_pl.filter(annovar_pl.is_unique()).rename({"avsnp151": "rsid"})
    print(annovar_pl)
    ################################################################################################
    # PART TWO:
    # Read in significant_loci and concatenate files
    ################################################################################################
    sig_pl = pl.read_csv(significant_loci, separator="\t", null_values="None")
    print(sig_pl)

    final_pl = (
        (pl.concat([sig_pl, annovar_pl], how="align"))
        .with_columns(
            pl.when((pl.col("Allele1") != pl.col("ref_from_id")))
            .then(
                pl.col("Direction")
                .str.replace_all(r"\+", "pos")
                .str.replace_all(r"\-", "neg")
            )
            .otherwise(pl.col("Direction"))
            .alias("Direction_fixed")
        )
        .with_columns(
            (
                pl.col("Direction_fixed")
                .str.replace_all("pos", "-")
                .str.replace_all("neg", "+")
            ).alias("direction_ALT")
        )
        .with_columns(
            pl.col("direction_ALT")
            .str.count_matches(r"\+|\-")
            .alias("Number of data sets (per biobank per ancestry)")
        )
    ).filter(pl.col('chr').is_not_null())
    print(final_pl)


    ################################################################################################
    # PART THREE:
    # Format the final output
    ################################################################################################

    # Count The Number of Biobanks with P-value less than 5e-8
    col_names = final_pl.columns
    sig_pval_count = [
        _count_significant_biobanks(row, col_names) for row in final_pl.iter_rows()
    ]
    final_pl = final_pl.with_columns(
        pl.Series(
            "number biobanks with significant association (p-value < 5E-8)",
            sig_pval_count,
        )
    )

    # Find Minimum P-value
    col_names = final_pl.columns
    min_pvals = [_min_pval_biobanks(row, col_names) for row in final_pl.iter_rows()]
    final_pl = final_pl.with_columns(
        pl.Series("minimum p-value among all biobanks", min_pvals)
    )

    # Find Potential Strand Flip
    col_names = final_pl.columns
    potential_flips = [
        _find_potential_strand_flip(row, col_names) for row in final_pl.iter_rows()
    ]
    final_pl = final_pl.with_columns(
        pl.Series(
            "QC_strand_Flip (flag to indicate whether this variant had a strand flip in any data set)",
            potential_flips,
        )
    )

    # Find QC Allele Frequency Differences
    col_names = final_pl.columns
    af_flags = [_find_af_qc(row, col_names) for row in final_pl.iter_rows()]
    final_pl = final_pl.with_columns(
        pl.Series(
            "QC_Freq diff than gnomAD (flag to indicate whether this variant failed the AF QC (compared to gnomAD) in any data set)",
            af_flags,
        )
    )

    # Record ancestry per variant
    col_names = final_pl.columns
    ancestries = [_find_ancestries(row, col_names) for row in final_pl.iter_rows()]
    final_pl = final_pl.with_columns(pl.Series("Ancestry*", ancestries))

    # Add number of biobanks
    col_names = final_pl.columns
    biobanks = [
        _count_number_of_biobanks(row, col_names) for row in final_pl.iter_rows()
    ]
    final_pl = final_pl.with_columns(pl.Series("number of biobanks", biobanks))

    # Add Endpoint
    pheno = str(significant_loci).split("/")[-1].split("_significant_loci.tsv")[0]

    final_pl = final_pl.with_columns(Endpoint=pl.lit(pheno))

    # Add Category
    final_pl = final_pl.with_columns(
        pl.when(pl.col("known_in_gwas_catalog") == "known")
        .then(pl.lit("Previously reported"))
        .otherwise(pl.lit("Novel"))
        .alias("Category")
    )

    # Make previous records more readable
    columns = [
        "PUBMEDID",
        "FIRST AUTHOR",
        "LINK",
        "MAPPED_GENE",
        "SNP_GENE_IDS",
        "SNPS",
        "CATALOG_PVALUE",
    ]

    for col in columns:
        final_pl = final_pl.with_columns(
            pl.col(col).str.replace_all(";", ",").alias(col)
        )

    # Rename columns and select order
    missing_cols = [x for x in final_pl.columns if x.startswith("missing")]
    final_cols = [
        "Endpoint",
        "Category",
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
        "LINK",
        "MAPPED_GENE",
        "SNP_GENE_IDS",
        "SNPS",
        "CATALOG_PVALUE",
    ]
    for col in missing_cols:
        final_cols.append(col)

    final_pl = (
        final_pl.drop("REF")
        .rename(
            {
                "POS": "POS(hg38)",
                "alt_from_id": "ALT(effect)",
                "ref_from_id": "REF",
                "rsid": "rsID",
                "adjusted_metal_eaf": "allele frequency (ALT)",
                "adjusted_metal_beta": "beta (ALT)",
                "StdErr": "se",
                "P-value": "p-value",
                "cochran_q_p_value": "p-value for heterogeneity (Cochran’s Q)",
            }
        )
        .select(final_cols)
    )

    # Write final output
    final_pl.write_csv(output, separator="\t")


def _count_significant_biobanks(row, col_names) -> int:
    # Extract pvalues per biobank
    pvals = np.array(
        [row[i] for i in range(len(row)) if col_names[i].startswith("pval")]
    ).astype(float)
    sig_count = sum(pvals < 5e-8)

    return sig_count


def _count_number_of_biobanks(row, col_names) -> int:
    # Extract beta (arbitraty to get biobank name)
    betas = [row[i] for i in range(len(row)) if col_names[i].startswith("beta_")]
    columns = [
        col_names[i] for i in range(len(row)) if col_names[i].startswith("beta_")
    ]

    # Loop through biobank and grab ancestry
    biobanks = []
    for item in list(zip(columns, betas)):
        if item[1] is not None:
            biobank = item[0].split("_")[1]
            biobanks.append(biobank)

    return len(set(biobanks))


def _min_pval_biobanks(row, col_names) -> str:
    # Extract pvalues per biobank
    pvals = np.array(
        [row[i] for i in range(len(row)) if col_names[i].startswith("pval")]
    ).astype(float)
    pvals = pvals[np.logical_not(np.isnan(pvals))]
    # min_pval = "{:.3e}".format(min(pvals))
    min_pval = min(pvals)
    print(min_pval)

    return min_pval


def _find_potential_strand_flip(row, col_names) -> str:
    # Extract Potential Strand Flip per biobank
    strand_flips = [
        str(row[i]).replace("None", "PASS")
        for i in range(len(row))
        if col_names[i].startswith("Potential_Strand_Flip")
    ]
    # Count items in list that are 'PASS'
    pass_count = strand_flips.count("PASS")

    if pass_count == len(strand_flips):
        return "PASS"

    not_pass = len(strand_flips) - pass_count
    return f"FAIL: {not_pass} biobank(s) had a potential strand flip"


def _find_af_qc(row, col_names) -> str:
    # Extract AF Differences
    af_flags = [
        str(row[i]).replace("None", "No")
        for i in range(len(row))
        if col_names[i].startswith("outlier_stdev")
    ]

    # Count items in list that are 'No'
    pass_count = af_flags.count("No")

    if pass_count == len(af_flags):
        return "FALSE"

    not_pass = len(af_flags) - pass_count
    return f"TRUE: {not_pass} biobank(s) had a significantly differenct AF compared to gnomAD"


def _find_ancestries(row, col_names) -> str:
    # Extract beta (arbitraty to get biobank name)
    betas = [row[i] for i in range(len(row)) if col_names[i].startswith("beta_")]
    columns = [
        col_names[i] for i in range(len(row)) if col_names[i].startswith("beta_")
    ]

    # Loop through biobank and grab ancestry
    ancestries = []
    for item in list(zip(columns, betas)):
        if item[1] is not None:
            ancestry = (
                item[0].split("_")[-1].split(".tsv")[0].replace("nan", "mixed").lower()
            )
            ancestries.append(ancestry)

    return ",".join(set(ancestries))


if __name__ == "__main__":
    defopt.run(add_genes)
