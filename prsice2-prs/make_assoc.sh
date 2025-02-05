#awk '{print "chr"$1" "$8" "$11" "$12" "$13" "$14" "$17}' /pl/active/pozdeyevlab/GBMI/GWAS/GWAS/publishable_metal_results/ldsc_adjusted_meta_analysis_three_biobanks_eur_male_leave_ccpm_out_metal_results_benign_nodular_goiter.tsv | sed 's/chrSTUDY_ID/STUDY_ID/g' > eur_male_lco_benign_nodular_goiter.assoc

#awk '{print "chr"$1" "$8" "$11" "$12" "$13" "$14" "$17}' /pl/active/pozdeyevlab/GBMI/GWAS/GWAS/publishable_metal_results/ldsc_adjusted_meta_analysis_three_biobanks_eur_male_leave_ccpm_out_metal_results_tc_vs_bng.tsv | sed 's/chrSTUDY_ID/STUDY_ID/g' > eur_male_lco_tc_vs_bng.assoc

#awk '{print "chr"$1" "$8" "$11" "$12" "$13" "$14" "$17}' /pl/active/pozdeyevlab/GBMI/GWAS/GWAS/publishable_metal_results/ldsc_adjusted_meta_analysis_three_biobanks_eur_male_leave_ccpm_out_metal_results_graves.tsv | sed 's/chrSTUDY_ID/STUDY_ID/g' > eur_male_lco_graves.assoc

#awk '{print "chr"$1" "$8" "$11" "$12" "$13" "$14" "$17}' /pl/active/pozdeyevlab/GBMI/GWAS/GWAS/publishable_metal_results/ldsc_adjusted_meta_analysis_three_biobanks_eur_male_leave_ccpm_out_metal_results_lymphocytic_thyroiditis.tsv | sed 's/chrSTUDY_ID/STUDY_ID/g' > eur_male_lco_lymphocytic_thyroiditis.assoc

#awk '{print "chr"$1" "$8" "$11" "$12" "$13" "$14" "$17}' /pl/active/pozdeyevlab/GBMI/GWAS/GWAS/publishable_metal_results/ldsc_adjusted_meta_analysis_three_biobanks_eur_male_leave_ccpm_out_metal_results_hypothyroidism.tsv | sed 's/chrSTUDY_ID/STUDY_ID/g' > eur_male_lco_hypothyroidism.assoc

awk '{print "chr"$1" "$8" "$11" "$12" "$13" "$14" "$17}' /pl/active/pozdeyevlab/GBMI/GWAS/GWAS/publishable_metal_results/ldsc_adjusted_meta_analysis_with_mvp_four_biobanks_mixed_both_leave_ccpm_out_metal_results_benign_nodular_goiter.tsv | sed 's/Effect/beta/g' |sed 's/chrSTUDY_ID/STUDY_ID/g' > mixed_both_lco_benign_nodular_goiter.assoc

awk '{print "chr"$1" "$8" "$11" "$12" "$13" "$14" "$17}' /pl/active/pozdeyevlab/GBMI/GWAS/GWAS/publishable_metal_results/ldsc_adjusted_meta_analysis_with_mvp_four_biobanks_mixed_both_leave_ccpm_out_metal_results_hypothyroidism.tsv | sed 's/Effect/beta/g'| sed 's/chrSTUDY_ID/STUDY_ID/g' > mixed_both_lco_hypothyroidism.assoc

awk '{print "chr"$1" "$8" "$11" "$12" "$13" "$14" "$17}' /pl/active/pozdeyevlab/GBMI/GWAS/GWAS/publishable_metal_results/ldsc_adjusted_meta_analysis_with_mvp_four_biobanks_mixed_both_leave_ccpm_out_metal_results_graves.tsv | sed 's/Effect/beta/g'| sed 's/chrSTUDY_ID/STUDY_ID/g' > mixed_both_lco_graves.assoc

awk '{print "chr"$1" "$8" "$11" "$12" "$13" "$14" "$17}' /pl/active/pozdeyevlab/GBMI/GWAS/GWAS/publishable_metal_results/ldsc_adjusted_meta_analysis_with_mvp_four_biobanks_mixed_both_leave_ccpm_out_metal_results_thyroid_cancer.tsv | sed 's/Effect/beta/g'| sed 's/chrSTUDY_ID/STUDY_ID/g' > mixed_both_lco_thyroid_cancer.assoc

awk '{print "chr"$1" "$8" "$11" "$12" "$13" "$14" "$17}' /pl/active/pozdeyevlab/GBMI/GWAS/GWAS/publishable_metal_results/ldsc_adjusted_meta_analysis_three_biobanks_mixed_male_leave_ccpm_out_metal_results_hypothyroidism.tsv | sed 's/Effect/beta/g'| sed 's/chrSTUDY_ID/STUDY_ID/g' > mixed_male_lco_hypothyroidism.assoc

awk '{print "chr"$1" "$8" "$11" "$12" "$13" "$14" "$17}' /pl/active/pozdeyevlab/GBMI/GWAS/GWAS/publishable_metal_results/ldsc_adjusted_meta_analysis_three_biobanks_mixed_male_leave_ccpm_out_metal_results_thyroid_cancer.tsv | sed 's/Effect/beta/g'| sed 's/chrSTUDY_ID/STUDY_ID/g' > mixed_male_lco_thyroid_cancer.assoc
