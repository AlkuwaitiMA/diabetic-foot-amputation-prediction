# =============================================================================
# DIABETIC FOOT AMPUTATION STUDY — QASSIM REGION, SAUDI ARABIA (N=241)
# MASTER RUN SCRIPT
# =============================================================================
# Author  : Analysis pipeline
# Purpose : Execute all analysis steps in sequence
# Requires: R >= 4.3.3
# =============================================================================

cat("================================================================\n")
cat("  DIABETIC FOOT AMPUTATION STUDY — FULL ANALYSIS PIPELINE\n")
cat("================================================================\n\n")

# Set working directory to project root
setwd("/home/claude/diabetic_foot_analysis")

# Record start time
start_time <- Sys.time()

# ---- Run all steps in order ------------------------------------------------

cat(">> STEP 1: Data Cleaning & Recoding\n")
source("R/01_data_cleaning.R")
cat("   DONE\n\n")

cat(">> STEP 2: Descriptive Statistics\n")
source("R/02_descriptive_statistics.R")
cat("   DONE\n\n")

cat(">> STEP 3: Missing Data Imputation (MICE)\n")
source("R/03_mice_imputation.R")
cat("   DONE\n\n")

cat(">> STEP 4: Temporal Trend Analysis (Joinpoint-equivalent)\n")
source("R/04_temporal_trend.R")
cat("   DONE\n\n")

cat(">> STEP 5: LASSO Variable Selection\n")
source("R/05_lasso_selection.R")
cat("   DONE\n\n")

cat(">> STEP 6: Multivariable Regression (POM vs MLR vs PPOM)\n")
source("R/06_multivariable_regression.R")
cat("   DONE\n\n")

cat(">> STEP 7: Binary Logistic Models + Bootstrap Validation + Nomogram\n")
source("R/07_binary_logistic.R")
cat("   DONE\n\n")

cat(">> STEP 8: Machine Learning Models + SHAP\n")
source("R/08_machine_learning.R")
cat("   DONE\n\n")

cat(">> STEP 9: Decision Curve Analysis\n")
source("R/09_decision_curve.R")
cat("   DONE\n\n")

cat(">> STEP 10: Sensitivity Analyses\n")
source("R/10_sensitivity_analyses.R")
cat("   DONE\n\n")

# ---- Report total runtime --------------------------------------------------
end_time <- Sys.time()
cat("================================================================\n")
cat(sprintf("  PIPELINE COMPLETE — Total runtime: %.1f minutes\n",
            as.numeric(difftime(end_time, start_time, units="mins"))))
cat("  All outputs saved to: /home/claude/diabetic_foot_analysis/outputs/\n")
cat("================================================================\n")
