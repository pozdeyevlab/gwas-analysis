# gwas-qc
A custom pipeline for harmonizing REGENIE/SAIGE summary stats from a variety of bio-banks.

# How to use

### Environment & Dependency Set Up
If the steps below do not work please contact samantha.l.white@cuanschutz.edu for assistance. 
```bash
git clone https://github.com/pozdeyevlab/gwasqc.git
cd gwasqc
conda env create -f environment.yml
conda activate gwas_qc
poetry install
# To test correct install try
python modules/harmonize.py --help
```

### Necessary Input files

1) **meta analysis map file**: See details below in the _map file_ section for specific format of this file. Unfortunately, any user introduced errors are likely to orginate from improper formatting in this file.
2) **gnomAD reference files**: A directory of tab separated files named 'gnomad_ref_chr<chrom>.tsv' which contain reference loci, reference and alternate alleles, and allele frequencies according to sex and ancestry. These files can be generates using the pipeline found here.
[here](https://github.com/pozdeyevlab/gnomad-query.git)
3) **gnomAD flagged variant files**: A directory of tab separated files named 'flagged_variants_chr<chrom>.tsv'. These files can be generated following the pipeline found here.
[here](https://github.com/pozdeyevlab/gnomad-query.git)
4) **config.yaml**: See details below in _config_ section.

After setting up everything above please test that the dag looks correct with `--dry-run` before proceeding to full run.

```bash
# adjust cores to your system or use 'all' to run in series use 1 core
# test
snakemake --cores 10 --configfile config.yaml --dry-run

# run
snakemake --cores 10 --configfile config.yaml
```

## Input & Output

### Input Files

#### Map file with the following column names (case sensitive):
|Column Name    |Description     |
|---------------|----------------|
|BIOBANK        |name of biobank (must be one word, we apologize for the onconvenience but 'Some_Bank' will yield buggy results instead please use 'SomeBank'!)|
|PHENOTYPE      |disease phenotype|
|SEX|male, female or both - any other string will cause an error|
|ANCESTRY|EAS, EUR, HIS, AMR.. etc<br/>If the ancestry is unknown or the gwas is multi-ethnic leave the cell blank<br/>If the provided ancestry is not available in the reference data then mixed ancestry will be used by default|
|GWAS_SOFTWARE|saige or regenie<br/>If saige then summary stats are assumed to be tab separated</br>If regenie then summary stats are assumed to be space separated|
|CHROM|column name in summary stat containing the chromosome|
|POS|column name in summary stat containing the position of the variant|
|Effect_Allele|column name in summary stat containing the effect allele|
|Non_Effect_Allele|column name in summary stat containing the non-effect allele|
|Effect_AF|column name in summary stat containing the effect allele frequency|
|BETA|column name in summary stat containing the beta value|
|SE|column name in summary stat containing the standard error|
|PVAL|column name in summary stat containing the p-value|
|ID|column name in summary stat containing the variant ID</br>assumed to be in chr:pos:ref:alt format|
|IMPUTE|column name in summary stat containing the imputation values (Optional)|
|Total_N|column name in summary stat containing the total N|
|Case_N|column name in summary stat containing the number of case samples|
|Control_N|column name in summary stat containing the number of control samples|
|PATH|path to summary stat file|
|Case_Count|number of cases reported from biobank|
|Control_Count|number of controls reported from biobank|

#### Config file with the following key-pair values:
|Key|Value|
|---|-----|
|map_file|"path/to/map/file”|
|gnomad_ref_dir|“path/to/gnomad/reference/file/directory”<br/>pipeline for generating reference files can be found [here](https://github.com/pozdeyevlab/gnomad-query)|
|gnomad_flag_dir|"path/to/gnomad/an_flagged/file/directory" some variants in gnomad have a quality warning when the alternate allele is found in less than half of all participants, these variants have been flagged for qc|
|output_path|“where/to/write/output/files/”|

### Output Files
#### Outputs:
|Output Type|Path|
|-----------|----|
|Aligned variants|`<output_path>/<biobank>_<phenotype>_<sex>_<ancestry>_aligned_to_gnomad.tsv`|
|Un-Aligned variants|`<output_path>/<biobank>_<phenotype>_<sex>_<ancestry>_not_aligned_in_gnomad.tsv`|
|Variants without a matching reference position|`<output_path>/<biobank>_<phenotype>_<sex>_<ancestry>_no_position_in_gnomad.tsv`|
|Un-usable variants|`<output_path>/<biobank>_<phenotype>_<sex>_<ancestry>_invalid_variants.tsv`<br/>If alternate allele is invalid ie `CN0` then those variants are removed from analysis|
|Plots|`<output_path>/plots/<biobank>_<phenotype>_<sex>_<ancestry>.png`<br/>scatterplots, histograms, and bar charts to help decide qc cut-offs|

## Pipeline Overview

1. For each summary stat file per chromosome:
    1. Read the chromosome specific gnomAD reference file into memory, scan only for positions that are present in the gwas summary file. Write variants without a matching position in gnomAD to a scratch file. 
   
    2. For summary stat variants with a matching position in gnomAD:
   
        1. Enforce variant id format <chrom>:<pos>:<non_effect_allele>:<effect_allele>
            1. Variants with missing alleles or structural variants are thrown out.
        2. Create flags for the following:
            * p-value is 0
            * -1e6 < BETA < 1e6
            * -1e6 < SE < 1e6
            * If imputation data is available, imputation is less than imputation threshold (default: 0.3)
            * Variant alleles are palindromic
            * Variant alleles are palindromic and effect allele frequency is between 0.4 - 0.6
            * Effect and non effect alleles match the order found in the variant ID
    3. Run `modules/align_palindromic_variants.py` on palindromic variants and `modules/align_non_palindromic_variants.py` on non palindromic variants. The reason for this distinction is that palindromic variants are not checked for transcribed matches, as that would be redundant.
    4. For **non palindromic variants**
        1. Check for an exact match
            * Example: 1:100:A:G & 1:100:A:G
        2. Check for an inverse match 
            * Example: 1:100:G:A & 1:100:A:G
            * **IMPORTANT**: The AF & Beta become 1 - AF & -1 * BETA
        3. For variants that align BOTH by exact and inverse match, find the absolute difference in allele frequency (AF) between the gnomAD reported AF and the gwas reported AF. The method that procuces the smallest difference in AF is saved as the true match. These counts are saved to the log file. 
        4. If there is not an exact or inverse match then the varaint is tested for an exact transcribed match
            * Example: 1:100:T:C & 1:100:A:G
        5. If there is not an exact, inverse, or transcribed exact match then the variant is tested for an inverse transcribed match 
            * Example: 1:100:C:T & 1:100:A:G
            * **IMPORTANT**:The AF & Beta become 1 - AF & -1 * BETA
    5. For **palindromic variants**
        1. Check for an exact match
            * Example: 1:100:A:T & 1:100:A:T
        2. Check for an inverse match 
            * Example: 1:100:T:A & 1:100:A:T
            * **IMPORTANT**: The AF & Beta become 1 - AF & -1 * BETA
        3. For variants that align BOTH by exact and inverse match, find the absolute difference in allele frequency (AF) between the gnomAD reported AF and the gwas reported AF. The method that procuces the smallest difference in AF is saved as the true match. These counts are saved to the log file. 
   
    5. Write the aligned, unaligned, and missing position variants to three respective output files

3) Per input gwas summary stat create the following plots:

       * Scatterplot: colored by alignment method
       * Scatterplot: colored by palinfromic or not palindromic
       * Scatterplot: colored by chisq p-val > 0.05 or chisq p-val < 0.05
       * Scatterplot: Plot of aligned variants with all possible filters applied
       * Histogram: abs(Study_AF - Reference_AF)
       * Bar plot: abs(Study_AF - Reference_AF) binned into intervals of 0.01 (0 - 1.0)
