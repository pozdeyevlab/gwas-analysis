#!/bin/bash

# bcftools_query.sh
#
# Usage:
#   bcftools_query.sh -v "gnomed_v4.vcf" -o "output_chr21.tsv" 
#
# Helper script to generate allele frequency table from gnomAD

set -o errexit
set -o nounset

# Read in named command line args 
while getopts ":v:o:c:t:" opt; do
  case $opt in
    v) vcf="$OPTARG"
    ;;
    o) out="$OPTARG"
    ;;
    c) chrom="$OPTARG"
    ;;
    t) temp="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    exit 1
    ;;
  esac

  case $OPTARG in
    -*) echo "Option $opt needs a valid argument"
    exit 1
    ;;
  esac
done

# Download vcf from gnomAD
gsutil cp $vcf $temp/$chrom.vcf.gz

# Query vcf
bcftools query -H -f'%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%AF\t%AF_XX\t%AF_XY\t%AF_afr_XX\t%AF_afr_XY\t%AF_afr\t%AF_ami_XX\t%AF_ami_XY\t%AF_ami\t%AF_amr_XX\t%AF_amr_XY\t%AF_amr\t%AF_asj_XX\t%AF_asj_XY\t%AF_asj\t%AF_eas_XX\t%AF_eas_XY\t%AF_eas\t%AF_fin_XX\t%AF_fin_XY\t%AF_fin\t%AF_mid_XX\t%AF_mid_XY\t%AF_mid\t%AF_nfe_XX\t%AF_nfe_XY\t%AF_nfe\t%AF_raw\t%AF_remaining_XX\t%AF_remaining_XY\t%AF_remaining\t%AF_sas_XX\t%AF_sas_XY\t%AF_sas\t%AF_joint_XX\t%AF_joint_XY\t%AF_joint\t%AF_joint_afr_XX\t%AF_joint_afr_XY\t%AF_joint_afr\t%AF_joint_ami_XX\t%AF_joint_ami_XY\t%AF_joint_ami\t%AF_joint_amr_XX\t%AF_joint_amr_XY\t%AF_joint_amr\t%AF_joint_asj_XX\t%AF_joint_asj_XY\t%AF_joint_asj\t%AF_joint_eas_XX\t%AF_joint_eas_XY\t%AF_joint_eas\t%AF_joint_fin_XX\t%AF_joint_fin_XY\t%AF_joint_fin\t%AF_joint_mid_XX\t%AF_joint_mid_XY\t%AF_joint_mid\t%AF_joint_nfe_XX\t%AF_joint_nfe_XY\t%AF_joint_nfe\t%AF_joint_raw\t%AF_joint_remaining_XX\t%AF_joint_remaining_XY\t%AF_joint_remaining\t%AF_joint_sas_XX\t%AF_joint_sas_XY\t%AF_joint_sas\t%AF_grpmax\t%AF_grpmax_joint\t%AN\t%AN_XX\t%AN_XY\t%AN_afr_XX\t%AN_afr_XY\t%AN_afr\t%AN_ami_XX\t%AN_ami_XY\t%AN_ami\t%AN_amr_XX\t%AN_amr_XY\t%AN_amr\t%AN_asj_XX\t%AN_asj_XY\t%AN_asj\t%AN_eas_XX\t%AN_eas_XY\t%AN_eas\t%AN_fin_XX\t%AN_fin_XY\t%AN_fin\t%AN_mid_XX\t%AN_mid_XY\t%AN_mid\t%AN_nfe_XX\t%AN_nfe_XY\t%AN_nfe\t%AN_raw\t%AN_remaining_XX\t%AN_remaining_XY\t%AN_remaining\t%AN_sas_XX\t%AN_sas_XY\t%AN_sas\t%AN_joint_XX\t%AN_joint_XY\t%AN_joint\t%AN_joint_afr_XX\t%AN_joint_afr_XY\t%AN_joint_afr\t%AN_joint_ami_XX\t%AN_joint_ami_XY\t%AN_joint_ami\t%AN_joint_amr_XX\t%AN_joint_amr_XY\t%AN_joint_amr\t%AN_joint_asj_XX\t%AN_joint_asj_XY\t%AN_joint_asj\t%AN_joint_eas_XX\t%AN_joint_eas_XY\t%AN_joint_eas\t%AN_joint_fin_XX\t%AN_joint_fin_XY\t%AN_joint_fin\t%AN_joint_mid_XX\t%AN_joint_mid_XY\t%AN_joint_mid\t%AN_joint_nfe_XX\t%AN_joint_nfe_XY\t%AN_joint_nfe\t%AN_joint_raw\t%AN_joint_remaining_XX\t%AN_joint_remaining_XY\t%AN_joint_remaining\t%AN_joint_sas_XX\t%AN_joint_sas_XY\t%AN_joint_sas\t%AN_grpmax\t%AN_grpmax_joint\n' $temp/$chrom.vcf.gz > $out

# Remove vcf
rm $temp/$chrom.vcf.gz