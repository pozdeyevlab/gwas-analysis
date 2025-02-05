# Manual PRS

"Manual PRS" refers to the process of taking the lead variants identified from meta-analysis, extracting genotypes from a bed file, calculating PRS as the summation of (genotype * effect size). The calculated PRS is then compared to the phenotype file to generate an AUC. 

# Required inputs

1) Directory of input pgens by chromsome (pgen, pvar, psam)
2) Phenotype files (IID, FID, phenotype[0,1])
3) Significant variants output from `run-metal/module/significant_loci.py`
4) Map file (example provided `map_file`)
5) Bcftools is required (and not invcluded in the conda install)

# Run Workflow
```bash
conda env create -f environment.yml
conda activate manual_prs
poetry install

# Check that you also have bcftools installed abd available
bcftools --help

# Run pipeline
snakemake --cores 1 --configfile config.yaml
```