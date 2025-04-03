library(optparse)
library(data.table)
library(dplyr)
library(tidyr)
library(glue)

# Define command-line arguments
option_list <- list(
  make_option(c("--input"), type="character", help="TSV with PRS, covariates, and risk categories"),
  make_option(c("--output"), type="character", help="Output file"),
  make_option(c("--prs"), type="character", help="Name of the PRS column"),
  make_option(c("--invnorm"), type="logical", default=FALSE, help="Binary of T/F ro inv normalize continuous risk strata")
)

# Parse command-line arguments
opt <- parse_args(OptionParser(option_list=option_list))

# Load data
df = fread(glue(opt$input), sep='\t')

# Scale PRS & Format Risk
df$PRS <- scale(df$PRS)
# T-Stage
df$T_stage[df$T_stage == "T1"] <- 1
df$T_stage[df$T_stage == "T1a"] <- 1
df$T_stage[df$T_stage == "T2"] <- 2
df$T_stage[df$T_stage == "T3"] <- 3
df$T_stage[df$T_stage == "T3a"] <- 3
df$T_stage[df$T_stage == "T4a"] <- 4
df$T_stage[df$T_stage == "T4b"] <- 4
df$T_stage <- as.numeric(df$T_stage)
# N-Stage
df$N_stage[df$N_stage == "N0"] <- 0
df$N_stage[df$N_stage == "N1a"] <- 1
df$N_stage[df$N_stage == "N1"] <- 1
df$N_stage[df$N_stage == "N1b"] <- 1

# M-Stage
df$M_stage[df$M_stage == "M0"] <- 0
df$M_stage[df$M_stage == "M1"] <- 1

# Tumor Focality
df$tumor_focality[df$tumor_focality == "Multifocal"] <- 1
df$tumor_focality[df$tumor_focality == "Unifocal"] <- 0

# Vascular Invasion
df$vascular_invasion[df$vascular_invasion == "Vascular invasion present"] <- 1
df$vascular_invasion[df$vascular_invasion == "Vascular invasion not present"] <- 0

# Extra-Thyroidal Extension
df$extra_thyroidal_extension[df$extra_thyroidal_extension == "Extra-thyroidal extension present"] <- 1
df$extra_thyroidal_extension[df$extra_thyroidal_extension == "Extra-thyroidal extension not present"] <- 0

# Lymphatic Invasion
df$lymphatic_invasion[df$lymphatic_invasion == "Lymphatic invasion present"] <- 1
df$lymphatic_invasion[df$lymphatic_invasion == "Lymphatic invasion not present"] <- 0

# Surgical Margins
df$surgical_margins[df$surgical_margins == "Tumor extends beyond surgical margins"] <- 1
df$surgical_margins[df$surgical_margins == "Tumor does not extends beyond surgical margins"] <- 0

# Extranodal Extension
df$extranodal_extension[df$extranodal_extension == "Extranodal extension present"] <- 1
df$extranodal_extension[df$extranodal_extension == "Extranodal extension not present"] <- 0

# Categorical Recurrence
df$risk_of_recurrence_categorical[df$risk_of_recurrence_categorical == 'Low Risk'] <- 0
df$risk_of_recurrence_categorical[df$risk_of_recurrence_categorical == 'Intermediate Risk'] <- 0
df$risk_of_recurrence_categorical[df$risk_of_recurrence_categorical == 'High Risk'] <- 1

# Inv-Norm if required
# Inv-Norm continuous outputs
invnormal <- function(x){qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))}
if (opt$invnorm) {
  df$age_at_diag_years <- invnormal(df$age_at_diag_years)
  df$risk_of_recurrence_continuous <- invnormal(df$risk_of_recurrence_continuous)
  df$tumor_max_dimension <- invnormal(df$tumor_max_dimension)
  df$lymph_node_metastasis_max_dimension <- invnormal(df$lymph_node_metastasis_max_dimension)
  df$number_of_metastatic_lymph_nodes <- invnormal(df$number_of_metastatic_lymph_nodes)
}

# Linear Regresison Analysis Function
run_lm_analysis <- function(df, prs, risk, out, covariates) {
  # Subset the data and remove NAs
  tmp <- df %>%
    dplyr::select(all_of(risk), all_of(covariates), PRS) %>%
    na.omit() %>%
    filter(risk != '') %>%
    filter(risk != 'NA')
    
  # Create the formula dynamically for the linear model
  formula <- as.formula(paste(risk, "~ ."))

  # Fit the linear model
  m <- lm(formula, data = tmp)

  # Get the summary of the model
  model_summary <- summary(m)

  # Extract coefficients and p-values
  coefficients <- data.frame(model_summary$coefficients)

  # Add Bonferroni correction for p-value
  coefficients$bonferroni <- 0.05 / (ncol(tmp) - 1)

  # Add PRS and risk information
  coefficients$prs <- prs
  coefficients$risk <- risk
    
  # Add the covariate names
  coefficients$covariates <- rownames(coefficients)

  # Extract p-values
  coefficients$pval <- coefficients$`Pr(>|t|)`

  # Bind the results to the output
  out <- rbind(out, coefficients)

  # Return the final output
  return(out)
}

# Logistic Regresison Analysis Function
run_glm_analysis <- function(df, prs, risk, out, covariates) {
  # Subset the data and remove NAs
  df[[risk]] <- as.numeric(df[[risk]])
  tmp <- df %>%
    dplyr::select(all_of(risk), all_of(covariates), PRS) %>%
    na.omit() %>%
    filter(risk != '') %>%
    filter(risk != 'NA')
  
  # Create the formula dynamically for the linear model
  formula <- as.formula(paste(risk, "~ ."))

  # Fit the logistic model
  m <- glm(formula, data = tmp, family='binomial')

  # Get the summary of the model
  model_summary <- summary(m)

  # Extract coefficients and p-values
  coefficients <- data.frame(model_summary$coefficients)

  # Add Bonferroni correction for p-value
  coefficients$bonferroni <- 0.05 / (ncol(tmp) - 1)

  # Add PRS and risk information
  coefficients$prs <- prs
  coefficients$risk <- risk
    
  # Add the covariate names
  coefficients$covariates <- rownames(coefficients)

  # Extract p-values
  coefficients$pval <- coefficients$`Pr(>|t|)`

  # Conform names 
  names(coefficients)[names(coefficients) == 'z.value'] <- 't.value'
  names(coefficients)[names(coefficients) == 'Pr...z..'] <- 'Pr...t..'

  # Bind the results to the output
  out <- rbind(out, coefficients)

  # Return the final output
  return(out)
}


# Call Logistic Regression
# Exclude 'age' as a covariate in 'age_At_diag_years'
out=NULL
covariates <- c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10','sex','Genotyping_Batch')
result <- run_lm_analysis(df, opt$prs, 'age_at_diag_years', out, covariates)
out=result
# Base-line covariates
covariates = c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10','sex','age','Genotyping_Batch')

# Loop through logistic risks
logistic_risks <- c('N_stage', 'M_stage', 'risk_of_recurrence_categorical', 'tumor_focality', 'extra_thyroidal_extension', 'vascular_invasion', 'lymphatic_invasion', 'surgical_margins', 'extranodal_extension')

for (risk in logistic_risks) {
  result <- run_glm_analysis(df, opt$prs, risk, out, covariates)
  out=result
}

linear_risks <- c('T_stage','risk_of_recurrence_continuous', 'lymph_node_metastasis_max_dimension', 'number_of_metastatic_lymph_nodes', 'tumor_max_dimension')

for (risk in linear_risks) {
  result <- run_lm_analysis(df, opt$prs, risk, out, covariates)
  out=result
}

# Save results
write.table(out, opt$output, row.names=FALSE, sep='\t', quote = FALSE, append = FALSE)
