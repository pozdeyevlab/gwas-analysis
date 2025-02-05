import polars as pl
metal_results = 'GBMI_G'
metal_df = (
        (
            pl.read_csv(metal_results, separator="\t")
            .with_columns(
                pl.col("direction")
                .str.count_matches(r"?", literal=True)
                .alias("missing_biobanks")
            )
            .with_columns(pl.col("direction").str.len_chars().alias("n_biobanks"))
            .filter(
                (pl.col("n_biobanks").cast(int) - pl.col("missing_biobanks").cast(int))
                >= 2
            )
            .drop(["n_biobanks", "missing_biobanks"])
        )
    ).rename({"MarkerName": "STUDY_ID", 'Effect': 'beta', 'ref_from_id': 'REF', 'alt_from_id': 'ALT'}).select(['STUDY_ID', 'beta', 'REF', 'ALT', 'P-value', 'POS', 'CHR']).write_csv('gbmi_I_tc_two_biobanks.assoc', separator=' ')
