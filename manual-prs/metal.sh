#!/bin/bash
#SBATCH --nodes=1
#SBATCH --partition=amem
#SBATCH --qos=mem
#SBATCH --ntasks=1 # requests 48 cores
#SBATCH --mem=600G # requests 500GB of RAM
#SBATCH --time=24:00:00 # requests the node for 24 hours
#SBATCH --job-name="manual_prs"
#SBATCH -o "/scratch/alpine/swhite3@xsede.org/manual_prs/metal.out"
#SBATCH -e "/scratch/alpine/swhite3@xsede.org/manual_prs/metal.err"
#SBATCH --account=amc-general

cd /projects/swhite3\@xsede.org/manual_prs

module load anaconda
conda activate metal
module load bcftools

snakemake --cores 10 --configfile metal.yaml
