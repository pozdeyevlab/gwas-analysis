"""
Module used to convert matrix table to bgen. Cannot be run on pre-emptible workers. 
"""

import gc
import os
import subprocess

import bokeh.io
import defopt
from bokeh import *
from bokeh.io import output_notebook, show
from bokeh.resources import INLINE

bokeh.io.output_notebook(INLINE)
import hail as hl
from bokeh.io import output_notebook, show
# Google cloud resources
from google.cloud import storage

hl.init(default_reference="GRCh38")


def convert_matrix_table(*, input_mt_dir: str, output_dir: str) -> None:
    """
    :param input_mt_dir: Google Bucket path to matrix tables created from 'vds_to_matrix_table.py'
    :param outptu_dir: Google Bucket path to output directory
    """

    chromosomes = [f"chr{i}" for i in list(range(1, 23))]
    for chrom in chromosomes:
        print(f"Starting {chrom}")
        mt = hl.read_matrix_table(f"{input_mt_dir}/{chrom}_filtered.mt")
        mt = mt.filter_rows(mt.filters.size() > 0, keep=False)
        mt = mt.filter_rows((mt.info.AC[0] >= 20) & (mt.info.AC[1] >= 20))

        # Adding dosages; all hard calls in WGS dataset
        homref_gp = hl.literal([1.0, 0.0, 0.0])
        het_gp = hl.literal([0.0, 1.0, 0.0])
        homvar_gp = hl.literal([0.0, 0.0, 1.0])

        # adding rsid
        mt = mt.annotate_rows(
            rsid=mt.locus.contig
            + "_"
            + hl.str(mt.locus.position)
            + "_"
            + mt.alleles[0]
            + "_"
            + mt.alleles[1]
        )

        # Adding dosages; all hard calls in WGS dataset
        homref_gp = hl.literal([1.0, 0.0, 0.0])
        het_gp = hl.literal([0.0, 1.0, 0.0])
        homvar_gp = hl.literal([0.0, 0.0, 1.0])

        mt = mt.annotate_entries(
            GP=hl.case()
            .when(mt.GT.is_hom_ref(), homref_gp)
            .when(mt.GT.is_het(), het_gp)
            .default(homvar_gp)
        )

        # Include GP & rsid
        hl.export_bgen(mt, f"{output_dir}/{chrom}_filtered", gp=mt.GP, rsid=mt.rsid)
        hl.index_bgen(f"{output_dir}/{chrom}_filtered.bgen", reference_genome="GRCh38")

        print(f"Successfully converted {chrom}")
        return None


if __name__ == "__main__":
    defopt.run(convert_matrix_table)
