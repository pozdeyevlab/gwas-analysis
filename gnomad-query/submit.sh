#!/bin/bash
#SBATCH --nodes=1
#SBATCH --partition=amilan
#SBATCH --ntasks=12
#SBATCH --job-name="plots_test_test"
#SBATCH -o "/scratch/alpine/swhite3@xsede.org/plots/plots_test.out"
#SBATCH -e "/scratch/alpine/swhite3@xsede.org/plots/plots_test.err"
#SBATCH --account=amc-general

cd /projects/swhite3@xsede.org/gnomad-query

module load anaconda
module load bcftools
conda activate gnomad_query
poetry install

snakemake --cores $nproc --configfile config.yaml

