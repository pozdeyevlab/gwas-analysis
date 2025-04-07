library(optparse)
library(stringr)
library(data.table)
library(pROC)
library(tidyr)
library(dplyr)


# Define command-line arguments
option_list <- list(
  make_option(c("--manual"), type="character", help="Manual Result"),
  make_option(c("--prscs"), type="character", help="PRSCS Results"),
  make_option(c("--output"), type="character", help="Output File")
)

# Parse command-line arguments
opt <- parse_args(OptionParser(option_list=option_list))

print(opt$manual)
print(opt$prscs)
print(opt$output)

# Load data
manual<- fread(opt$manual, sep='\t', header=TRUE)
prscs <- fread(opt$prscs, sep='\t', header=TRUE)

# Get a list of all overlapping column headers ending in prs
cols_df1 <- grep("_prs$", colnames(manual), value = TRUE)
cols_df2 <- grep("_prs$", colnames(prscs), value = TRUE)
prs_cols <- intersect(cols_df1, cols_df2)

# Get a list of all overlapping column headers NOT ending in prs
cols_df1 <- grep("_prs$", colnames(manual), value = TRUE, invert = TRUE)
cols_df2 <- grep("_prs$", colnames(prscs), value = TRUE, invert = TRUE)

# Exclude "IID" from the column names
cols_df1 <- setdiff(cols_df1, "IID")
cols_df2 <- setdiff(cols_df2, "IID")
pheno_cols <- intersect(cols_df1, cols_df2)

# For every overlapping PRS calculate a p-value
df_total = NULL

for (prs in prs_cols){
  for (pheno in pheno_cols){
    tmp <- manual %>% select(pheno, prs) %>% drop_na()
	  tmp_p <- prscs %>% select(pheno, prs) %>% drop_na()

    roc_manual <- roc(tmp[[pheno]], tmp[[prs]])
    roc_prscs <- roc(tmp_p[[pheno]], tmp_p[[prs]])

    p_value <- roc.test(roc_manual, roc_prscs, method = "delong")$p.value
    auc_m <- auc(roc_manual)
    auc_p <- auc(roc_prscs)
    predicting=pheno

    method_1 = str_glue('manual_{prs} ({round(auc_m, 4)})')
    method_2 = str_glue('prscs_{prs} ({round(auc_p, 4)})')

    delong=p_value

    df = data.frame(method_1,method_2,predicting,delong)
    df_total=rbind(df_total,df)
  }
}
print(df_total)

# If available compare thyroid cancer from GBMI I vs GBMI II in both data frames
if (('thyroid_cancer_GBMI_I_prs' %in% colnames(prscs))& ('thyroid_cancer_GBMI_II_prs' %in% colnames(prscs))) {
  for (pheno in pheno_cols){
    tmp_1 <- prscs %>% select(pheno, 'thyroid_cancer_GBMI_I_prs') %>% drop_na()
	  tmp_2 <- prscs %>% select(pheno, 'thyroid_cancer_GBMI_II_prs') %>% drop_na()

    roc_1 <- roc(tmp_1[[pheno]], tmp_1[['thyroid_cancer_GBMI_I_prs']])
    roc_2 <- roc(tmp_2[[pheno]], tmp_2[['thyroid_cancer_GBMI_II_prs']])

    p_value <- roc.test(roc_1, roc_2, method = "delong")$p.value
    auc_1 <- auc(roc_1)
    auc_2 <- auc(roc_2)
    predicting=pheno

    method_1 = str_glue('prscs_thyroid_cancer_GBMI_I_prs ({round(auc_1, 4)})')
    method_2 = str_glue('prscs_thyroid_cancer_GBMI_II_prs ({round(auc_2, 4)})')

    delong=p_value

    df = data.frame(method_1,method_2,predicting,delong)
    df_total=rbind(df_total,df)
  }
}

if (('thyroid_cancer_GBMI_I_prs' %in% colnames(manual))& ('thyroid_cancer_GBMI_II_prs' %in% colnames(manual))) {
  for (pheno in pheno_cols){
    tmp_1 <- manual %>% select(pheno, 'thyroid_cancer_GBMI_I_prs') %>% drop_na()
	  tmp_2 <- manual %>% select(pheno, 'thyroid_cancer_GBMI_II_prs') %>% drop_na()

    roc_1 <- roc(tmp_1[[pheno]], tmp_1[['thyroid_cancer_GBMI_I_prs']])
    roc_2 <- roc(tmp_2[[pheno]], tmp_2[['thyroid_cancer_GBMI_II_prs']])

    p_value <- roc.test(roc_1, roc_2, method = "delong")$p.value
    auc_1 <- auc(roc_1)
    auc_2 <- auc(roc_2)
    predicting=pheno

    method_1 = str_glue('manual_thyroid_cancer_GBMI_I_prs ({round(auc_1, 4)})')
    method_2 = str_glue('manual_thyroid_cancer_GBMI_II_prs ({round(auc_2, 4)})')

    delong=p_value

    df = data.frame(method_1,method_2,predicting,delong)
    df_total=rbind(df_total,df)
  }
}
print(df_total)

write.table(df_total, opt$output, row.names=FALSE, sep='\t')