"""
This module formats the output from modules/significant_loci.py as a bed file so that it is compatible with annovar software. 
"""

from pathlib import Path

import defopt
import polars as pl


# pylint: disable = C0301
def make_bed(
    *,
    input_file: Path,
    output: Path,
) -> None:
    """
    :param input_file: Path
    :param output: Path
    """
    chr_col = "chr"
    pos_col = "pos"
    id_col = "MarkerName"
    ref_col = "ref"
    alt_col = "alt"
    try:
        input_pl: pl.DataFrame = pl.read_csv(
            input_file,
            separator="\t",
            columns=[id_col],
        ).with_columns(
            pl.col(id_col)
            .str.split(":")
            .map_elements(lambda arr: arr[2].upper())
            .alias(ref_col),
            pl.col(id_col)
            .str.split(":")
            .map_elements(lambda arr: arr[3].upper())
            .alias(alt_col),
            pl.col(id_col)
            .str.split(":")
            .map_elements(lambda arr: arr[1])
            .alias(pos_col),
            pl.col(id_col)
            .str.split(":")
            .map_elements(lambda arr: arr[0].replace("X", "23"))
            .alias(chr_col)
            .cast(int),
        )
    except pl.exceptions.ComputeError as error:
        print(
            f"There is an error reading in the input file {input_file}\nError: {error}"
        )
    print(input_pl)
    input_pl = input_pl.with_columns(pl.col(pos_col).cast(int).alias(pos_col))
    print(input_pl)
    input_pl = input_pl.with_columns(
        pl.when(((pl.col(ref_col)).str.len_bytes()) > ((pl.col(alt_col)).str.len_bytes()))
        .then(pl.col(pos_col) + (((pl.col(ref_col)).str.len_bytes() - 1)))
        .otherwise(pl.col(pos_col))
        .alias("END")
    )
    input_pl = input_pl.select([chr_col, pos_col, "END", ref_col, alt_col])

    input_pl.write_csv(output, separator=" ", include_header=False)


if __name__ == "__main__":
    defopt.run(make_bed)
