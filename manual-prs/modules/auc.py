"""
Calculate AUC from manually calculated PRS
"""

from pathlib import Path
import sys
import defopt
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
from sklearn import metrics

# pylint: disable = C0301
# pylint: disable = R0903 # Too few public methods
# pylint: disable = R1728


def auc(
    *,
    phen_col: str,
    gt_file: Path,
    beta_file: Path,
    phenotype_file: Path,
    plot_output: str,
    summary_output: str,
    phenotype: str,
    header_file: Path
) -> None:
    """
    :param gt_file: Path to genotype data
    :param header_file: File with one column 'ArbID' followed by all sample ID's from merged psam
    :param phenotype_file: Hardcall phenotype files
    :param plot_output: Where to write plot png
    :param summary_output: Where to record performance metrics
    :param beta_file: File with variant IDs and beta values
    :param phen_col: Name of phenotype column in phenotype file
    :param phenotype: Name of phenotype
    """
    # Read in all required data
    header = pl.read_csv(header_file, separator='\t', has_header=False, new_columns = ['header'])['header'].to_list()
    genotypes = pl.read_csv(gt_file, separator="\t", dtypes={'betas':float}, has_header=False, new_columns = header, null_values = 'NA')
    betas = pl.read_csv(
        beta_file, separator="\t", has_header=False, new_columns=["ID", "beta"]
    ).with_columns(("chr" + pl.col("ID").cast(str).str.replace("chrchr", "chr")).alias("ID"))
    genotypes = genotypes.join(betas, on="ID", how="inner")
    phenotypes = pl.read_csv(phenotype_file, separator="\t", null_values='NA').select(['ArbID', phen_col]).drop_nulls()

    # Data check
    print(genotypes)

    # Calculate genotype values 0, 1, 2
    temp = genotypes.drop(['ID', 'beta'])
    print(temp)
    final = temp.with_columns(pl.all().cast(int) * genotypes['beta']).sum().transpose(include_header=True, column_names=['prs']).rename({'column':'ArbID'}).with_columns(pl.col('ArbID').cast(int).alias('ArbID')).join(phenotypes, on = 'ArbID', how='inner')
    print(final)
    
    # Calculate AUC
    fpr, tpr, threshold = metrics.roc_curve(
        final[phen_col].to_numpy(), final["prs"].to_numpy()
    )
    roc_auc = metrics.auc(fpr, tpr)
    auc = metrics.roc_auc_score(final[phen_col].to_numpy(), final["prs"].to_numpy())
    
    print(auc)
    raw_output = summary_output.replace('.tsv', '_raw.tsv')
    print(raw_output)
    final.write_csv(raw_output, separator='\t') 
    # Save data
    # Check if file exists
    new_data = {"phen_ids": phenotype, "auc": auc}
    pl.DataFrame(new_data).write_csv(summary_output, separator="\t")

    # Plot & Save png
    name = [i.capitalize() for i in phenotype.split("_")]
    name = [i.replace('Gbmi', 'GBMI') for i in name]
    name = [i.replace('Ii', f'II/n') for i in name]
    name = [i.replace('Vs', 'vs') for i in name]
    name = ' '.join(name)
    name = f'Thyroid Cancer GBMI I\nLeave CCPM Out'
    
    fig = plt.figure()
    ax = fig.add_subplot()
    plt.title(f"{name}", fontsize = 15)
    plt.plot(fpr, tpr, "b", label="AUC = %0.4f" % roc_auc)
    plt.legend(loc="lower right")
    plt.plot([0, 1], [0, 1], "r--")
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    plt.ylabel("True Positive Rate", fontsize = 15)
    plt.xlabel("False Positive Rate", fontsize = 15)
    ax.set_aspect('equal', adjustable='box')
    plt.savefig(plot_output)

if __name__ == "__main__":
    defopt.run(auc)
