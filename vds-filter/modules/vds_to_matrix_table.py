"""
Module used to filter whole genome variant data set to subset of user provided variants, then write data to hail matrix table. Can be run on pre-emptible workers. 
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


def filter_vds(
    *, variant_file: str, vds_storage_location: str, output_dir: str
) -> None:
    """
    :param variant_file: Google Bucket path to tab seperated file which contains the following columns;[CHR,POS,REF,ALT]
    :param vds_storage_location: Google Bucket path to whole genome variant data set (vds)
    :param outptu_dir: Google Bucket path to output directory
    """
    table = hl.import_table(variant_file, delimiter="\t", impute=True)

    # Join two columns together with ":"
    table = table.annotate(
        locus=hl.locus(table.CHR, table.POS, reference_genome="GRCh38")
    )
    table = table.annotate(alleles=hl.array([table.REF, table.ALT]))
    table = table.key_by("locus", "alleles")

    # Load VDS
    vds_storage_location = os.getenv("WGS_VDS_PATH")
    wgs_vds = hl.vds.read_vds(vds_storage_location)

    # Loop thgough and filter
    chromosomes = [f"chr{chr}" for chr in list(range(1, 23))]
    for chrom in chromosomes:
        print(f"Starting {chrom}")
        filtered_table = table.filter(table.CHR == chrom)
        vds = hl.vds.filter_chromosomes(wgs_vds, keep=chrom)
        vds = hl.vds.split_multi(vds)
        vds = hl.vds.filter_variants(vds, variants_table=filtered_table, keep=True)
        mt = vds.variant_data
        mt = mt.transmute_entries(FT=hl.if_else(mt.FT, "PASS", "FAIL"))
        mt = hl.vds.to_dense_mt(hl.vds.VariantDataset(vds.reference_data, mt))
        mt = mt.annotate_rows(info=hl.agg.call_stats(mt.GT, mt.alleles))
        mt = mt.drop(
            "as_vqsr",
            "tranche_data",
            "truth_sensitivity_snp_threshold",
            "truth_sensitivity_indel_threshold",
            "snp_vqslod_threshold",
            "indel_vqslod_threshold",
        )
        mt.write(f"{output_dir}/{chrom}_filtered.mt", overwrite=True)
        print(f"Successfully filtered {chrom}\n")

    return None

if __name__ == "__main__":
    defopt.run(filter_vds)
