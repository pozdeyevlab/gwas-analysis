"""Make plots for deciding qc cut-offs"""
from pathlib import Path

import defopt
import matplotlib.pyplot as plt
import polars as pl

# pylint: disable=C0301 # line too long
# pylint: disable=R0914 # Too many local variables
# pylint: disable=R0913 # Too many arguments
# pylint: disable=R0915 # Too many statements
# pylint: disable=C0121 # Comparison


def plot(
    *,
    file_path: Path,
    output_file: Path,
) -> None:
    """
    Helper for generating scatterplots for effect allele frequencies between reference and study variants

    :param file_path: Path to alignment results files
    :param output_file: Where plots will be written
    """
    # Columns to read from each file
    columns_to_read = [
        "Aligned_AF",
        "AF_gnomad",
        "palindromic_af_flag",
        "gwas_is_palindromic",
        "FILTER",
        "Alignment_Method",
        "ABS_DIF_AF",
        "CHR_gnomad",
        "POS_gnomad",
        "REF_gnomad",
        "ALT_gnomad",
        "outlier_pval",
        "outlier_stdev",
        "mahalanobis",
        "GNOMAD_AN_Flag"
    ]

    # Read specific columns from all files into a single polars DF
    print(f"File used:{file_path}")
    combined_df = _read_specific_columns(file_path, columns_to_read)

    # Filter aligned variants for all possible qc values
    combined_df = combined_df.filter((pl.col("FILTER") == "PASS"))

    # Prepare data for palindromic and filter scatterplots
    ref_eaf_non_palindromic = combined_df.filter(
        pl.col("gwas_is_palindromic") == False
    )["AF_gnomad"]
    ref_eaf_palindromic = combined_df.filter(pl.col("gwas_is_palindromic") == True)[
        "AF_gnomad"
    ]

    # Set up scatterplot figure
    figure, axes = plt.subplots(nrows=3, ncols=2, figsize=(40, 40))

    # Palindromic vs Non-Palindromic
    alignment_methods = combined_df["Alignment_Method"].unique()
    palindrome = combined_df["gwas_is_palindromic"].unique()

    for alignment_method in alignment_methods:
        # Scatterplots
        axes[0, 0].scatter(
            combined_df.filter((pl.col("Alignment_Method") == alignment_method))["AF_gnomad"],
            combined_df.filter((pl.col("Alignment_Method") == alignment_method))[
                "Aligned_AF"
            ],
            label=alignment_method.replace("_", " "),
            s=10,
        )

    axes[0, 0].set_title(
        f"Alignment Methods Palindromic & Non-Palindromic\nN:{combined_df.shape[0]}",
        fontsize=20,
    )
    axes[0, 0].set_xlabel("gnomad EAF", fontsize=15)
    axes[0, 0].set_ylabel("study EAF", fontsize=15)
    axes[0, 0].legend(loc="upper right")

    # Palindromic vs Non-Palindromic
    for pal in palindrome:
        # Scatterplots
        axes[0, 1].scatter(
            combined_df.filter((pl.col("gwas_is_palindromic") == pal))["AF_gnomad"],
            combined_df.filter((pl.col("gwas_is_palindromic") == pal))["Aligned_AF"],
            label=pal,
            s=10,
        )

    axes[0, 1].set_title(
        f"Palindromic vs Non Palindromic\nN:{combined_df.shape[0]}", fontsize=20
    )
    axes[0, 1].set_xlabel("gnomad EAF", fontsize=15)
    axes[0, 1].set_ylabel("study EAF", fontsize=15)
    axes[0, 1].legend(loc="upper right")

    # AN gnomAD Filter
    axes[1, 0].scatter(
        combined_df.filter(pl.col("GNOMAD_AN_Flag") == 0)["AF_gnomad"],
        combined_df.filter(pl.col("GNOMAD_AN_Flag") == 0)["Aligned_AF"],
        label="AN >= .5(max(AN))",
        s=10,
    )
    axes[1, 0].scatter(
        combined_df.filter(pl.col("GNOMAD_AN_Flag") == 1)["AF_gnomad"],
        combined_df.filter(pl.col("GNOMAD_AN_Flag") == 1)["Aligned_AF"],
        label="AN < .5(max(AN))",
        s=10,
    )
    axes[1, 0].set_title(
        f"AN Filter Based On GnomAD Warning\nFlagged:{len(combined_df.filter(pl.col('GNOMAD_AN_Flag') == 1)['AF_gnomad'])}",
        fontsize=20,
    )
    axes[1, 0].set_xlabel("gnomad EAF", fontsize=15)
    axes[1, 0].set_ylabel("study EAF", fontsize=15)
    axes[1, 0].legend(loc="upper right")

    # Mahalanobis Distance Outliers
    for out in ["Yes", "No"]:
        # Scatterplots
        axes[1, 1].scatter(
            combined_df.filter((pl.col("outlier_stdev") == out))["AF_gnomad"],
            combined_df.filter((pl.col("outlier_stdev") == out))["Aligned_AF"],
            label=out,
            s=10,
        )

    axes[1, 1].set_title(
        f"Mahalanobis Outliers (Stdev Based)\nOutliers:{combined_df.filter(outlier_stdev='Yes').shape[0]}",
        fontsize=20,
    )
    axes[1, 1].set_xlabel("gnomad EAF", fontsize=15)
    axes[1, 1].set_ylabel("study EAF", fontsize=15)
    axes[1, 1].legend(loc="upper right")

    # Mahalanobis & An Flags
    axes[2, 0].scatter(
        combined_df.filter((pl.col("outlier_stdev") == out) & (pl.col("GNOMAD_AN_Flag") == 0))["AF_gnomad"],
        combined_df.filter((pl.col("outlier_stdev") == out) & (pl.col("GNOMAD_AN_Flag") == 0))["Aligned_AF"],
        label=out,
        s=10,
    )

    axes[2, 0].set_title(
        f'Mahalanobis Outliers Removed\n& AN Flagged Variants Removed\nN Remaining:{combined_df.filter((pl.col("outlier_stdev") == out) & (pl.col("GNOMAD_AN_Flag") == 0)).shape[0]}',
        fontsize=20,
    )
    axes[2, 0].set_xlabel("gnomad EAF", fontsize=15)
    axes[2, 0].set_ylabel("study EAF", fontsize=15)
    axes[2, 0].legend(loc="upper right")

    # Histogram of AF differences
    axes[2, 1].hist(combined_df["mahalanobis"], bins=100, edgecolor="black")
    axes[2, 1].set_title(
        "Distribution of\nMahalanobis Distances",
        fontsize=20,
    )
    figure.tight_layout()
    plt.savefig(
        output_file,
        format="png",
        dpi=300,
    )
    plt.clf()
    plt.close()


def _read_specific_columns(file, columns):
    df = pl.read_csv(file, columns=columns, separator="\t", dtypes={"CHR_gnoamd": str})
    return df


if __name__ == "__main__":
    defopt.run(plot)
