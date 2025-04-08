# FUSION Implementation

We ran [FUSION](http://gusevlab.org/projects/fusion/) to better understand the genes mapped from our lead variants. We used precomputed [GTEx v8 All Samples model for Thyroid](https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/ALL/GTExv8.ALL.Thyroid.tar.gz) (n = 574). 

## Requiured Inputs
1) The inputs to this analysis are intermediate files realily available from the [run-metal](https://github.com/pozdeyevlab/gwas-analysis/tree/main/run-metal) pipeline. The necessary files will be found in `{metal_output_dir}/ldsc_pre_munge/*_with_rsid.tsv`
2) Precomputed [GTEx model](https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/ALL/GTExv8.ALL.Thyroid.tar.gz)
3) [FUSION software](http://gusevlab.org/projects/fusion/#installation)