"""
Calculate cochran q p-value to test for heterogeneity of beta values across all input datasets
"""

import defopt
import numpy as np
import polars as pl
import scipy

# pylint: disable=R0914, R0913, R0903, C0301


def cochran_q(*, df: pl.DataFrame) -> pl.DataFrame:
    """
    :param df: Data Frame
    """
    ################################################################################################
    # PART ONE:
    # Read in significnat loci
    ################################################################################################
    # Calculate Cochran's Q p-value for each row
    col_names = df.columns
    cochran_q_p_values = [_cochran_q_p_value(row, col_names) for row in df.iter_rows()]

    # Add Cochran's Q p-value to the DataFrame
    df = df.with_columns(pl.Series("cochran_q_p_value", cochran_q_p_values))
    print(df["cochran_q_p_value"])
    return df


def _cochran_q_p_value(row, col_names):
    # Extract effect sizes and standard errors from the row
    betas = np.array(
        [row[i] for i in range(len(row)) if col_names[i].startswith("beta")]
    ).astype(float)

    ses = np.array(
        [row[i] for i in range(len(row)) if col_names[i].startswith("se")]
    ).astype(float)

    # Remove NA values from both beta and se
    betas = betas[~np.isnan(betas)]
    ses = ses[~np.isnan(ses)]

    # Check if we have enough valid observations for Cochran's Q test
    if len(betas) < 2:
        return np.nan

    weights = 1 / ses**2
    sum_inv_var = sum(weights)
    effs_inv_var = weights * betas
    beta_meta = sum(effs_inv_var) / sum_inv_var
    effs_size_org = betas
    p_value = _het_test(effs_size_org, weights, beta_meta)

    return p_value


def _het_test(effs_sizes, weights, effs_size_meta):
    k = len(effs_sizes)

    effs_sizes_array = np.array(effs_sizes)
    weights_array = np.array(weights)
    eff_dev = weights_array * ((effs_sizes_array - effs_size_meta) ** 2)
    sum_eff_dev = np.sum(eff_dev)

    return scipy.stats.distributions.chi2.sf(sum_eff_dev, k - 1)


if __name__ == "__main__":
    defopt.run(cochran_q)
