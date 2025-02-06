"""
Harmonization module
This module does the following:
    1) Calls filter_gwas.filter_summary_stats() to read in the summary stat file as a polars df and apply qc filters
    2) Calls get_gnomad_ref.read_reference() to read in the gnomad reference file as a polars df searching for ancestry and sex specific data
    3) Removes variants that do not have a matching reference position, and writes variants to output
    4) Creates the following additional ID columns in the gwas summary df
        a) Flipped_Allele_ID = chr-pos-effect_allele-non_effect_allele
        b) Transcribed_ID = chr-pos-transcribed_non_effect-transcribed_effect
        c) Transcribed_Flipped_ID = chr-pos-transcribed_effect-transcribed_non_effect
    5) Call align_palindromic_variants.harmonize() and align_non_palindromic_variants.harmonize() respectively.
    6) If the gwas variant aligns via inverse match or transcribed flipped match then aligned beta becomes (-1 * original beta) and the aligned af becomes (1 - original af)
    7) Results are concatenated and written to <aligned/file_path>aligned_to_gnomad.tsv
    8) All unaligned variants are written to <aligned/file_path>un_aligned.tsv
"""

from datetime import datetime
from pathlib import Path
from typing import Optional, Union

import align_non_palindromic_variants
import align_palindromic_variants
import attr
import defopt
import filter_gwas
import get_gnomad_ref
import mahalanobis
import polars as pl

# pylint: disable=C0301 # line too long
# pylint: disable=R0914 # too many local variables
# pylint: disable=R0915 # Too many statements
# pylint: disable=R0913 # too many arguments
# pylint: disable=R0903 # too few public methods


@attr.s(frozen=False, auto_attribs=True, kw_only=True)
class AlignmentResults:
    """Class to represent the alignment outcome"""

    aligned_and_merged: Optional[pl.DataFrame]
    unaligned: pl.DataFrame


def harmonize(
    *,
    gwas_results: Path,
    gwas_software: str,
    output_file: Path,
    output_unaligned: Path,
    gnomad_ref_dir: Path,
    sex: str,
    ancestry: str,
    chromosome: Union[int, str],
    chrom_col: Optional[str],
    position: Optional[str],
    ea: Optional[str],
    non_ea: Optional[str],
    eaf: Optional[str],
    beta: Optional[str],
    pval: Optional[str],
    se: Optional[str],
    variant_id: str,
    output_no_position: str,
    impute: Optional[str] = None,
    n_case: Optional[str] = None,
    n_control: Optional[str] = None,
    n_total: Optional[str] = None,
) -> None:
    """
    :param gwas_results: Input path to regenie or saige summary stats
    :param gwas_software: Define if 'regenie' or 'saige' was used to generate summary stats
    :param output_file: Where plots will be written
    :param sex: Sex if any
    :param chromosome: The chromosome to align
    :param ancestry: Three letter ancestry code if any (EUR, FIN, AFR, etc)
    :param gnomad_ref_dir: Path to gnomAD reference files
    :param chrom: Column name of chromosome in summary stat
    :param position: Column name of genomic position in summary stat
    :param ea: Column name of effect allele in summary stat
    :param non_ea: Column name of non-effect allele in summary stat
    :param eaf: Column name of effect allele frequency in summary stat
    :param beta: Column name of beta in summary stat
    :param se: Column name of standard error in summary stat
    :param n_case: Column name of case N in summary stat
    :param n_control: Column name of control N in summary stat
    :param n_total: Column name of total N in summary stat
    :param impute: Column name of imputation value (if available) in summary stat
    :param pval: Column name of p-value in summary stat
    :param variant_id: Column name of variant id (if any, cannot be rsids) in summary stat
    :param output_unaligned: Path to file for unaligned variants

    """
    # Read in summary stats for specific chromosome
    start_original = datetime.now()
    gwas_results_class: filter_gwas.Results = filter_gwas.filter_summary_stats(
        gwas_results=gwas_results,
        gwas_software=gwas_software,
        chrom=chrom_col,
        position=position,
        ea=ea,
        non_ea=non_ea,
        eaf=eaf,
        beta=beta,
        pval=pval,
        se=se,
        variant_id=variant_id,
        chromosome=chromosome,
        impute=impute,
        n_case=n_case,
        n_control=n_control,
        n_total=n_total,
        unusable_path=Path(
            f'{output_no_position.replace("_no_position_in_gnomad", "_invalid_variants")}'
        ),
    )

    gwas_pl: pl.DataFrame = gwas_results_class.summary_stats
    col_map: filter_gwas.Columns = gwas_results_class.column_map

    print(f"\nStarting harmonization for chr{chromosome}:")

    # Based on chromosome, read in the corresponding gnomad reference and find appropriate AN and AF columns

    # EXCEED has 1-22 and 23 which is listed as 'X' when looking at the 'ID' col

    if chromosome == 23:
        chromosome = "X"
        gwas_pl = gwas_pl.with_columns(
            (pl.col(col_map.chrom).replace("23", "X")).alias(col_map.chrom)
        )
        gnomad_tsv = list(gnomad_ref_dir.glob(f"gnomad_*chr{chromosome}.tsv"))[0]
    else:
        gnomad_tsv = list(gnomad_ref_dir.glob(f"gnomad_*chr{chromosome}.tsv"))[0]

    total_variants = gwas_pl.shape[0]
    gnomad_pl: pl.DataFrame = get_gnomad_ref.read_reference(
        sex=sex,
        ancestry=ancestry,
        gnomad_tsv=gnomad_tsv,
        positions=list(set(gwas_pl[col_map.pos])),
    )

    # Find the number of non-existant positions
    dif = set(gwas_pl[col_map.pos]).difference(set(gnomad_pl["POS_gnomad"]))

    # Write file of non matching positions
    _write_or_append_to_file(
        output_no_position, gwas_pl.filter(pl.col(col_map.pos).is_in(dif))
    )

    # Remove un-matched positions from gwas_df to reduce search space
    gwas_pl = gwas_pl.filter(~pl.col(col_map.pos).is_in(dif))

    # Number of varinats in gwas with matching position in gnomad
    chr_possible = gwas_pl.shape[0]
    print(
        f"Alignment:\n{chr_possible}/{total_variants} Variants in the gwas summary file {gwas_results} have matching positions in the reference file {gnomad_tsv}"
    )

    #################################### Start Alignment Logic #####################################
    # Find the number of non-palindromic inverse pairs in gnomad (chr10:100:A:G & chr10:100:G:A)
    print(f'gnomad_pl:\n{gnomad_pl}')
    print(f'gwas_pl:\n{gwas_pl}')
    palindromic_results: pl.DataFrame = align_palindromic_variants.harmonize(
        gnomad_pl=gnomad_pl,
        gwas_pl=gwas_pl.filter(pl.col(col_map.palindromic_flag) == 1),
        col_map=col_map,
    )

    non_palindromic_results: pl.DataFrame = align_non_palindromic_variants.harmonize(
        gnomad_pl=gnomad_pl,
        gwas_pl=gwas_pl.filter(pl.col(col_map.palindromic_flag) == 0),
        col_map=col_map,
    )

    stacked_pl = pl.concat(
        [palindromic_results, non_palindromic_results], how="diagonal"
    )

    # Re-order stacked_pl so that is is easier to read
    reorder = []

    first = ['CHR_gnomad', 'POS_gnomad', 'REF_gnomad', 'ALT_gnomad', 'AF_gnomad', col_map.eaf, 'Aligned_AF', 'Aligned_Beta', 'Alignment_Method']

    reorder.extend(first)
    [reorder.append(x) for x in stacked_pl.columns if x not in first]
    stacked_pl = stacked_pl.select(reorder)

    # Make unaligned output files
    unaligned = gwas_pl.filter(
        ~(pl.col(col_map.variant_id).is_in(stacked_pl[col_map.variant_id]))
    )
    _write_or_append_to_file(output_file, stacked_pl)
    _write_or_append_to_file(output_unaligned, unaligned)
    end = datetime.now()
    total = end - start_original
    print(f"\nCompleted alignment for {gwas_results} chr{chromosome} in {total}\n")


def _write_or_append_to_file(file_path: Path, data: pl.DataFrame) -> None:
    data.write_csv(file_path, separator="\t", include_header=True)


if __name__ == "__main__":
    defopt.run(harmonize)
