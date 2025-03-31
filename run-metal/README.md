# Inverse Variance Weighted Meta-Analysis & Downstream Processing
This pipeline is specifically designed to be used in conjunction with the [gwasqc](https://github.com/pozdeyevlab/gwasqc) pipeline. If you wish to use this workflow with your own summary stats, please ensure that your input summary stats have the required columns, discussed _Required USER Provided Input Files_. 

## Environment & Dependency Set Up
If the steps below do not work please contact samantha.l.white@cuanschutz.edu for assistance. 
```bash
git clone ''
cd metal-ldsc
conda env create -f environment.yml
conda activate metal
poetry install
```

## What is done in this pipeline?
1) Inputs are grouped by phenotype, meaning that multiple phenotypes (aka endpoints) can be run at the same time. 
2) GWAS summary stats are formatted to be analyzed by METAL: `modules/prep_for_metal.py`
    * Enforces the following naming column names
        * STUDY_ID
        * beta
        * se
        * Aligned_AF
        * MAC
        * pval
        * REF_gnomad
        * ALT_gnomad
        * AF_gnomad
        * Potential_Strand_Flip
        * outlier_stdev
        * per_variant_N
        * missingness
    * Applies the following QC
        1) Remove 'chrom' from chromosome column
        2) Rename 'X' to '23' if present in chromosome column
        3) Filter for variants that have passed the gnomAD (v1.40) QC retained from alignment in the "FILTER" column
        4) Filter for varaints that have not aligned to a flagged gnomAD (v1.40) variant
        5) Remove variants that were calculated to have a difference in allele frequency (GWAS vs gnomad) greater than 3 stdev from the mean Mahalanobis distance 
        6) Assign 'per_variant_N' based on the available data
        7) Calculate minor allele count (MAC) and filter accordingly
        8) Filter for MAF > 0.0005
        9) Ensure there are no duplicate varaints
        10) Caluclate per varaint missingness
        11) Remove variants with a missing p-value
        12) Remove variants with potenital strand flips
3) If the column meta-data column `lower_ci_95_intercept` is not Null and the value is greater than 1.0 then p-values will be adjusted for genomic correction using the LDSC intercept. 
4)  Dynamically writes the script to run METAL: `modules/write_metal_script.py`
    * METAL is default to run inv-var weigthed meta-analysis
    * If you are curious what this file looks like, it is written here: `<output_dir>/run-metal/<pehnotype>.metal`
5) Runs METAL: `~line 319 in Snakefile`
    * METAL binary is included in this repository
6) Summarize metal outputs: `modules/summary.py`
    * Takes raw METAL output and does the following:
        * Filters for genome wide significant variants (significant level must be set by the user in the config file, we used 5e-8)
        * Filters for variants that are present in the user supplied 'biobank_min' we used 3. Meaning that if a variant was significnat but only present in 2 input summary stats, then it would be removed from downstream analysis. 
    * After filtering, re-attaches summary stats from input files to meta analyzed variants
    * Sometimes one input summary stat may have a larger than expected impact on the meta-analysis. For that reason we re-attach statistics (beta, af, p-val, etc) from the input summary stats to the meta-analysis results file. 
    * Re-annotation is done using logic from the 'Direction' column and the info file output by METAL
7) Create Manhattan & QQ Plots from metal output: `modules/plots.py`
    * Generates manhattan/qq plots to easily assess results
    * Only plots variants that are found in the user provided minimum number of biobanks. 
8) Run LDSC Regression & Genomic Correlation: `modules/prep_for_ldsc.py`
    * In order to run ldsc rsid's must be attached to all variants from METAL
    * The munge & ldsc python (2) scripts are called directly in the Snakefile
        * LDSC is run for individual phenotypes and together to calculate the genetic correlation matrix
        * See instructions below on setting up LDSC & required references
9) Caculates the most significant variant per non-overlapping genomic region: `modules/significnt_loci.py`
    * Based on methods used in the GBMI flagship paper we developed an algorithm which identifies start and end coordinates surrounding the most significnat variant, which may be tagging the causal variant. This list of significant variants is then annotated with ANNOVAR and compared to a list of known variants. 
10) Run ANNOVAR to attach the nearest gene to list of significant variants: `modules/prep_for_annovar.py`
    * Significant variants are formatted for ANNOVAR and then annotation is called directly from the Snakefile
    * See below for instructions on installing ANNOVAR and required references
11) Create the final summary: `modules/add_genes.py`
    * To simplify interpretation of the results a final summary is created in `<output_directory>/final_summary_with_genes/<phenotype>_with_genes_and_regions.tsv`. This file contains information on the meta-analyzed p-value, beta, cochran q p-value, direction etc, as well as PMID's, Authors, and article links for previously published variants. 
        * Variants found in the GWAS catalog which have a position that lies between the start and end regions of out tagged variant are called 'Previously Known' if no known variants fall within the regions of the tagged varaint then that variant & gene are labeled 'Novel'. 

## Download External Software & Reference Files
_Links below were valid as of August 2024_
1) **LDSC scripts**: Download LDSC by following their instructions [here](https://github.com/bulik/ldsc)
2) **Ancestry specific linkage disequilibrium score reference files**: European references can be download from a variety of sources including free and and requester pays methods see this [github issues post](https://github.com/bulik/ldsc/issues/371) for more information. 
3) **Hapmap reference file**: A list of snps which will be used to subsed data for ldsc regression. This file can be downloaded [here](https://ibg.colorado.edu/cdrom2021/Day06-nivard/GenomicSEM_practical/eur_w_ld_chr/w_hm3.snplist)
4) **ANNOVAR**: The ANNOVAR perl scripts cannot be included in this repository, as registration is required. To download ANNOVAR please follow the instructions [here](https://www.openbioinformatics.org/annovar/annovar_download_form.php). Once the perl scripts have been obtained follow the instructions below to download reference files, see the ANNOVAR [documentation](https://annovar.openbioinformatics.org/en/latest/user-guide/startup/) for further details:
```bash
./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene humandb/

./annotate_variation.pl -buildver hg38 -downdb cytoBand humandb/

./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar exac03 humandb/ 

./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar avsnp151 humandb/ 

./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar dbnsfp30a humandb/
```
(`annovar/example/gene_xref.txt` should be included upon initial download)


## Required USER Provided Input Files
1) **Map File**: This is a tab separated file containing all necessary meta-data (see included example `map_file.tsv` and detailed layout below). Unfortunately, any user introduced errors are likely to orginate from improper formatting in this file.

It is VERY important to only include summary stats in your map file which you wish to include in your meta-analysis. This workflow automatically groups summary stats by phenotype, but not by ancestry or sex. Meaning that if you include sex and or ancestry stratified analysis in the map file everything will by meta-analyzed under the same phenotype. 

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
|population_composition_afr|afr summix output (not required)|
|population_composition_amr|amr summix output (not required)|
|population_composition_eas|eas summix output (not required)|
|population_composition_nfe|nfe summix output (not required)|
|population_composition_mid|mid summix output (not required)|
|population_composition_sas|sas summix output (not required)|
|intercept|intercept from [ldsc results](https://github.com/pozdeyevlab/gwas-analysis/tree/main/run-ldsc)|
|lower_ci_95_intercept|lower 95% CI of LDSC intercept|

2) **Known Variant Catalog(s)**: To determine if genes have been previously reported in the GWAS catalog a tsv per phenotype is requred with the following columns:
   
|Column Name    |Description     |
|---------------|----------------|
|PUBMEDID |PUBMEDID of associated publication|
|Journal |Name of journal|
|LINK|Link to publication|
|DISEASE/TRAIT|Phenotype name|
|CHR_ID|Chromosome|
|CHR_POS|GRCh38 Position|
|REPORTED GENE(S)|Associated Gene(s)|
|SNPS|Variant ID|

3) **config.yaml**: The config file is used to summply snakemake with the necessary user provided inputs.
   
|Key |Value |
|----|------|
|map_file |"map_file.tsv” Path to user provided map file|
|genomic_control |ON or OFF, for METAL|
|qc_input_dir |Path to top directory of summary stats. In order to be read by the pipeline aligned summary stats must follow this naming convention `{biobank}_{phenotype}_{sex}_{ancestry}_aligned_results.tsv`|
|prep_dir |To avoid unnecessarily re-writing metal input files provide the path where they _will be_ or _already_ have been written EX: "Pipeline_Output/METAL_INPUTS"|
|adj_dir |Where to run selective genomic correction using ldsc intercept|
|output_dir|“where/to/write/output/files/” Path to where the user wants to write output directories|
|metal_binary|"matal/metal" Path to executable [METAL](https://github.com/statgen/METAL) binary (comes pre-installed in repository)|
|rsid_map_file|Path to file with SNPID & rsid's see details on creating that file here[here]()|
|hapmap|Path to w_hm3.snplist|
|ref_ld_chr|Path to LD score reference directory|
|w_ld_chr|Path to LD score weights directory|
|sig_threshold|5e-8 The maximum p-value allowed for variants to be considered significant|
|phenotype_catalog_dict|Dictionary with phenotype:path to known variant catalog|
|disease_prevalence_dict|Dictionary of population prevalence of phenotypes (see example config.yaml)|
|gwas_catalog_disease_col |Name of disease col in catalog "DISEASE/TRAIT"|
|mac_filter|Integer for minimum minor allele count (we used 20)|
|post_metal_filter|Integer for minimum minor allele count of total across all meta-analysis (we used 20)|
|biobank_min|Minimum number of biobanks that each meta-analyzed variant must be found in to be considered in downstream analysis (we used 4)|
|annovar_perl|Path to `table_annovar.pl` (from ANNOVAR)|
|humandb_dir|Path to ANNOVAR databse reference directory|
|xref|Path to xref file, likely `annovar/example/gene_xref.txt`|
|ldsc_env|LDSC requires a python2 environment, which can be found in `ldsc_env.yml`|
|cli_munge_script|Path to `ldsc/munge_sumstats.py`|
|ldsc_script|Path to `ldsc/ldsc.py`|


# Running the Pipeline
After setting up everything above please test that the dag looks correct with `--dry-run` before proceeding to full run.

```bash
# adjust cores to your system or use 'all', to run in series use 1 core
# test
snakemake --cores 10 --configfile config.yaml --use-conda --conda-frontend conda --dry-run

# run
snakemake --cores 10 --configfile config.yaml --use-conda --conda-frontend conda
```

# Output Files
## Outputs:
|Output Type|Path|
|-----------|----|
|Metal Prep Files|`<output_path>/metal_prep/<biobank>_<phenotype>_<sex>_<ancestry>.tsv`|
|Metal Results|`<output_path>/metal_results/<phenotype>.metal & <phenotype>_1.tbl & <phenotype>_1.info`|
|Metal Plots|`<output_path>/metal_plots/<phenotype>.png`|
|Metal Significant Variant Summary|`<output_path>/metal_summary/<phenotype>.tsv`|
|LDSC Pre Munge|`<output_path>/ldsc_pre_munge/<phenotype>.tsv`|
|LDSC Munge|`<output_path>/ldsc_munge/<phenotype>.sumstats.gz`|
|LDSC Results|`<output_path>/cli_ldsc/<phenotype>.log`| 
|LDSC Matrix Results|`<output_path>/cli_ldsc_matrix/<phenotype>.log`|
|Final Summary With Genes and Known Variants|`<output_path>/final_summary_with_genes/<phenotype>.log`|

