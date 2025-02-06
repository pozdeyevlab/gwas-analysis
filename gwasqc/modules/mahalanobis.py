"""
Calculates mahalanobis distance for gwas-af vs gnomad-af
"""
import defopt
import numpy as np
import polars as pl
from scipy.spatial import distance
from scipy.stats import chi2

# pylint: disable=C0301 # line too long
# pylint: disable=R0914 # too many local variables
# pylint: disable=R0915 # Too many statements
# pylint: disable=R0913 # too many arguments
# pylint: disable=R0903 # too few public methods


def calculate(
    *,
    aligned_pl: pl.DataFrame,
) -> None:
    """
    :param aligned_pl: Data frame with aligned variants
    """
    data_np = aligned_pl.select("AF_gnomad", "Aligned_AF").to_numpy()

    # Calculate covariance matrix
    covariance_matrix = np.cov(data_np.T)

    # Calculate Mahalanobis distances for all rows
    mahalanobis_distances = []
    for row in data_np:
        #a = abs(row[0] - row[1])
        gwas = np.array([row[0], row[0]])
        gnomad = np.array([row[1], row[1]])

        mahalanobis_distance = distance.mahalanobis(
            gwas, gnomad, np.linalg.inv(covariance_matrix)
        )
        mahalanobis_distances.append(mahalanobis_distance)

    # Add values back to polars df
    aligned_pl = aligned_pl.with_columns(
        pl.Series("mahalanobis", mahalanobis_distances)
    )

    # Calculate outliers based on mean + (3 * std)
    std = aligned_pl.std()
    avg = aligned_pl.mean()
    aligned_pl = (
        aligned_pl.with_columns(
            pl.when(
                pl.col("mahalanobis") > (avg["mahalanobis"] + (3 * std["mahalanobis"]))
            )
            .then(pl.lit("Yes"))
            .otherwise(pl.lit("No"))
            .alias("outlier_stdev")
        )
        .with_columns(mahalanobis_mean=avg["mahalanobis"])
        .with_columns(mahalanobis_stdev=std["mahalanobis"])
    )

    # Calcualte outlier as pval < 0.001
    aligned_pl = aligned_pl.with_columns(
        mahalanobis_pval=(1 - aligned_pl["mahalanobis"].apply(lambda x: chi2.cdf(x, 1)))
    )
    aligned_pl = aligned_pl.with_columns(
        pl.when(pl.col("mahalanobis_pval") < 0.001)
        .then(pl.lit("Yes"))
        .otherwise(pl.lit("No"))
        .alias("outlier_pval")
    )

    return aligned_pl


if __name__ == "__main__":
    defopt.run(calculate)
