"""
This module takes in raw output directly from METAL and using the list of input files provided in the info file re-attaches p-values, betas, mac, se, missingness, N, mahalanobis outlier results (from gwasqc), potential strand flip (from gwasqc), and gnomad aligned effect allele frequency (from gwasqc) to each variant from metal according to direction. This module also adjusts beta values and allele frequencies output by METAL according to the match between the allele order reported in the snp id, and the REF & ALT alleles reported in the columns Allele1 and Allele2. 

Example:
If the direction is ‘++?’ then the first two biobanks listed in the info file created by metal will be recorded and the input data for those biobanks will be read and appropriately recorded. Otherwise, the values will be null, in the case of missing biobanks (‘?’). 

Outputs:
This module creates two main outputs which are coded as ‘summary_df’ and ‘metal_df’ in the code below. The first, ‘summary_df’, contains all genome wide significant variants from METAL and is written to {output_dir}/metal_summary/{phenotype}.tsv. The latter, ‘metal_df’, contains all genome wide significant variants with the summary stats per biobanks. 
"""

import subprocess
from pathlib import Path
from typing import List, Optional

import defopt
import polars as pl

# pylint: disable=R0914, R0913, R0915, C0121, line-too-long, unused-import


def summarize(
    *,
    metal_input_dir: Path,
    metal_results: Path,
    output: Path,
    info: Path,
    output_raw: Path,
    biobank_min: Optional[int],
    or_filter: float = None,
    af_filter: float = None,
) -> None:
    """
    :param metal_input_dir: Directory with metal input files
    :param metal_results: Path to alignment results files
    :param output: Where to write results
    :param info: Info path from metal output
    :param output_raw : If a path is provided the metal output file is written to theat path
    :param biobank_min: Optionally remove variants that are not present in a user defined minimum of biobanks
    :param or_filter: OR Minimum
    :param af_filter: AF Minimum
    """
    # Get ordered list of input stats from metal info file
    names_list = extract_info_from_file(info)

    # Columns to read from each file
    columns_to_read = [
        "MarkerName",
        "Allele1",
        "Allele2",
        "P-value",
        "Direction",
        "Freq1",
        "Effect",
        "StdErr",
        "per_variant_N",
    ]

    # Read in metal results
    metal_df = (
        (
            pl.read_csv(metal_results, columns=columns_to_read, separator="\t")
            .filter(pl.col("P-value") < 5e-8)
            .with_columns(
                pl.col("MarkerName")
                .str.split(":")
                .map_elements(lambda arr: arr[2].upper())
                .alias("ref_from_id"),
                pl.col("MarkerName")
                .str.split(":")
                .map_elements(lambda arr: arr[3].upper())
                .alias("alt_from_id"),
                pl.col("MarkerName")
                .str.split(":")
                .map_elements(lambda arr: arr[1])
                .alias("POS"),
                pl.col("MarkerName")
                .str.split(":")
                .map_elements(lambda arr: arr[0].replace("X", "23"))
                .alias("CHR")
                .cast(int),
            )
            .with_columns(pl.col("Allele1").str.to_uppercase().alias("Allele1"))
            .with_columns(pl.col("Allele2").str.to_uppercase().alias("Allele2"))
            .with_columns(
                pl.when(
                    (pl.col("ref_from_id") == pl.col("Allele1").str.to_uppercase())
                    & (pl.col("alt_from_id") == pl.col("Allele2").str.to_uppercase())
                )
                .then(True)
                .otherwise(False)
                .alias("reported_alleles_match_id")
            )
        )
        .with_columns(
            pl.when(pl.col("reported_alleles_match_id"))
            .then(pl.col("Freq1"))
            .otherwise(1 - pl.col("Freq1"))
            .alias("adjusted_metal_eaf")
        )
        .with_columns(
            pl.when(pl.col("reported_alleles_match_id"))
            .then(pl.col("Effect"))
            .otherwise(-1 * pl.col("Effect"))
            .alias("adjusted_metal_beta")
        )
    ).rename({"MarkerName": "STUDY_ID"})
    if (af_filter is not None) and (or_filter is not None):
        metal_df = metal_df.filter(pl.col('adjusted_metal_df').exp()<or_filter)
        metal_df = metal_df.filter(pl.col('adjusted_metal_af')>af_filter)
    summary_df = metal_df

    # Read in the input to each name
    raw_files = []
    for name in names_list:
        try:
            raw_tsv = list(metal_input_dir.glob(f"{name}.tsv"))[0]
            raw_files.append(raw_tsv)
        except IndexError:
            adj_dir = Path('/pl/active/pozdeyevlab/GBMI/GWAS/GWAS/QC_Pipeline_Output/FINAL_METAL_INPUTS_PRE_ADJ/metal_results')
            raw_tsv = list(adj_dir.glob(f"{name}.tsv"))[0]
            print(name)
            raw_files.append(raw_tsv)
    dfs = []
    count = 0

    for file in raw_files:
        columns = [
            "STUDY_ID",
            "Aligned_AF",
            "MAC",
            "pval",
            "se",
            "beta",
            "Potential_Strand_Flip",
            "outlier_stdev",
            "per_variant_N",
            "missingness",
        ]
        df = pl.read_csv(file, separator="\t", columns=columns).filter(
            pl.col("STUDY_ID").is_in(metal_df["STUDY_ID"])
        )

        columns.remove("STUDY_ID")
        biobank = file.parts[-1].replace(".tsv", "")
        for col in columns:
            df = df.rename({col: f"{col}_{biobank}"})

        dfs.append(df)
        count = count + 1

    raw_data = pl.concat(dfs, how="align")
    metal_df = metal_df.join(raw_data, on="STUDY_ID", how="inner")

    # metal_df = metal_df.join(details_from_raw, on = 'STUDY_ID')
    cols = [col for col in metal_df.columns if col.startswith("MAC")]
    metal_df = metal_df.with_columns(
        metal_df.select(cols).sum_horizontal().alias("total_mac")
    )

    # metal_df = metal_df.join(details_from_raw, on = 'STUDY_ID')
    cols = [col for col in metal_df.columns if col.startswith("per_variant_N")]
    metal_df = metal_df.with_columns(
        metal_df.select(cols).sum_horizontal().alias("total_N")
    )
    biobank_count = metal_df.with_columns(
        pl.col("Direction").str.len_bytes().alias("temp")
    )["temp"].max()
    if biobank_min <= biobank_count:
        metal_df = metal_df.filter(
            (pl.col("Direction").str.count_matches(r"\+|\-")) >= biobank_min
        )
    else:
        metal_df = metal_df.filter(
            (pl.col("Direction").str.count_matches(r"\+|\-")) >= 2
        )

    # Filter summary accordingly
    print(summary_df)
    print(metal_df)
    summary_df = summary_df.filter(pl.col("STUDY_ID").is_in(metal_df["STUDY_ID"]))
    summary_df.write_csv(output, separator="\t")
    metal_df.write_csv(output_raw, separator="\t")


def extract_info_from_file(info_file: Path) -> List[str]:
    """
    Reads info file from METAL in order to get list of input files. Should return a list of names, if an error occurs an empty list will be returned.

    Args:
        info_file: Path to file output from metal ending in .info
    """
    # Define the command to execute
    command = ["grep", "Input File", f"{info_file}"]
    command += ["|", "cut", "-d':'", "-f", "3", "|", "rev"]
    command += ["|", "cut", "-d'/'", "-f", "1", "|", "rev"]
    command += ["|", "cut", "-d'.'", "-f", "1"]
    try:
        # Run the command and capture its output
        result = subprocess.run(
            " ".join(command), shell=True, capture_output=True, text=True, check=True
        )

        # Split the output into lines and save them into an ordered list
        output_lines = result.stdout.split("\n")

        # Remove empty strings from the list
        output_lines = [line.strip() for line in output_lines if line.strip()]

        return output_lines

    except subprocess.CalledProcessError as e:
        # Handle if the command exits with a non-zero status
        print("Error executing command:", e)
        return []


if __name__ == "__main__":
    defopt.run(summarize)

