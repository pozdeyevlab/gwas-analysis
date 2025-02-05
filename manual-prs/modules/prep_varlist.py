"""
Calculate AUC from manually calculated PRS
"""

from pathlib import Path

import defopt
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
from sklearn import metrics

# pylint: disable = C0301
# pylint: disable = R0903 # Too few public methods
# pylint: disable = R1728


def prep_varlist(
    *,
    significant_loci: str,
    output_path: Path,
) -> None:
    """
    :param significant_loci: Path to metal results (expect raw output)
    :param output_path: Where to write tsv (no header, snp id and beta)
    """
    # Read in all required data
    metal_pl = pl.read_csv(significant_loci, separator="\t").select(['STUDY_ID','Allele1','Allele2','Effect']).with_columns((pl.col('STUDY_ID').str.replace('chr', '')).alias('STUDY_ID'))
    print(metal_pl)
    # Ensure that beta is correct direction
    metal_pl = metal_pl.with_columns(pl.col('STUDY_ID').str.split(':').map_elements(lambda arr: arr[2]).alias('REF_FROM_ID')).with_columns(pl.when((pl.col('REF_FROM_ID').str.to_uppercase() == pl.col('Allele1'))).then(pl.col('Effect')).otherwise(-1*pl.col('Effect')).alias('CORRECTED_BETA')).select(['STUDY_ID', 'CORRECTED_BETA'])

    # Print for log
    print(metal_pl)

    # Write to output
    metal_pl.write_csv(output_path, separator='\t', include_header=False)

    return None


if __name__ == "__main__":
    defopt.run(prep_varlist)
