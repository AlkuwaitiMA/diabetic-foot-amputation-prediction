# =============================================================================
# STEP 3: MISSING DATA IMPUTATION — MICE
# =============================================================================
# Input : outputs/df_clean.rds
# Output: outputs/mice_imputed.rds     — mids object (20 imputed datasets)
#         outputs/mice_diagnostics.csv — convergence summary
# =============================================================================
# Variables with missing data:
#   hba1c       : 9  missing (3.7%) — predictive mean matching
#   dm_duration : 51 missing (21%) — polytomous logistic
#   cause       : 23 missing (9.5%)— polytomous logistic
#   chronic_disease: 5 missing (2%) — logistic
# =============================================================================

suppressPackageStartupMessages({
  library(mice)
  library(dplyr)
})

set.seed(123)

df <- readRDS("outputs/df_clean.rds")

cat(sprintf("  Loaded clean data: %d rows\n", nrow(df)))

# ---- 3.1  Select variables for imputation model -----------------------------
# Include all variables that will enter the regression models + the outcome
# Exclude: derived/redundant columns (amp_group_ord, amp_group_nom, major_binary
#          are all derived from amp_group — keep only amp_group_nom for imputation)
# Exclude: level-specific columns (minor_level, major_level — not used in models)
# Exclude: log transforms (derived from wbc/creatinine — impute originals, then derive)

imp_vars <- c(
  "amp_group_nom",          # outcome
  "gender", "age_group", "nationality",
  "dm_type", "dm_duration", "oha_use", "insulin_use",
  "chronic_disease", "cvd", "htn", "dm_neuro", "dm_nephro", "ckd", "malignancy",
  "dry_gangrene", "wet_gangrene", "cause",
  "site_admission", "hosp_stay", "icu_admission", "dm_admissions",
  "osteomyelitis", "prev_surgery",
  "hba1c", "wbc", "creatinine", "hb",
  "year"
)

df_imp <- df %>% select(all_of(imp_vars))

# Confirm missingness
miss_counts <- colSums(is.na(df_imp))
cat("\n  Missing values in imputation dataset:\n")
print(miss_counts[miss_counts > 0])

# ---- 3.2  Set up imputation methods per variable ----------------------------

# Default method matrix from mice
ini <- mice(df_imp, maxit = 0, print = FALSE)
meth <- ini$method

# Assign methods explicitly
# hba1c: continuous — predictive mean matching (pmm)
meth["hba1c"]        <- "pmm"
# dm_duration: ordered factor — polytomous logistic regression
meth["dm_duration"]  <- "polyreg"
# cause: nominal factor — polytomous logistic regression
meth["cause"]        <- "polyreg"
# chronic_disease: binary — logistic regression
meth["chronic_disease"] <- "logreg"

# Variables with no missing — set to "" (don't impute, use as predictors)
meth[miss_counts == 0] <- ""

cat("\n  Imputation methods assigned:\n")
print(meth[meth != ""])

# ---- 3.3  Set predictor matrix ----------------------------------------------
# All variables predict each other (default); 
# exclude outcome from predicting itself

pred <- ini$predictorMatrix

# Outcome predicts covariates but we keep it as predictor for covariates
# (fully conditional specification — default behaviour is fine)

# ---- 3.4  Run MICE ----------------------------------------------------------

cat("\n  Running MICE (m=20 imputations, maxit=20 iterations)...\n")
cat("  This may take 1-2 minutes.\n")

mice_out <- mice(
  data            = df_imp,
  m               = 20,        # 20 imputed datasets
  maxit           = 20,        # 20 iterations per imputation
  method          = meth,
  predictorMatrix = pred,
  seed            = 123,
  print           = FALSE
)

cat("  MICE completed.\n")

# ---- 3.5  Convergence diagnostics -------------------------------------------

cat("\n  Checking convergence (mean and SD of imputed values per iteration):\n")

# Extract imputed values for hba1c across iterations
hba1c_means <- apply(mice_out$chainMean["hba1c", , ], 2, mean)
hba1c_sds   <- apply(mice_out$chainVar["hba1c",  , ], 2, mean)

cat(sprintf("  HbA1C — imputed mean range across iterations: %.2f to %.2f\n",
            min(hba1c_means), max(hba1c_means)))
cat(sprintf("  HbA1C — imputed SD   range across iterations: %.2f to %.2f\n",
            min(sqrt(hba1c_sds)), max(sqrt(hba1c_sds))))

# Check convergence: variance of chain means should be small relative to
# within-chain variance (R-hat analog)
diag_df <- data.frame(
  variable  = c("hba1c"),
  iter_mean_min = min(hba1c_means),
  iter_mean_max = max(hba1c_means),
  iter_sd_min   = min(sqrt(hba1c_sds)),
  iter_sd_max   = max(sqrt(hba1c_sds)),
  converged     = ifelse(
    (max(hba1c_means) - min(hba1c_means)) / mean(hba1c_means) < 0.05,
    "Yes", "Check"
  )
)

cat("\n  Convergence summary:\n")
print(diag_df)

# ---- 3.6  Add log-transformed variables to each imputed dataset -------------

mice_out_extended <- mice_out   # copy

# We need to add log_wbc and log_creatinine to each completed dataset
# Do this by completing all 20 datasets and adding derived variables

imputed_list <- lapply(1:20, function(i) {
  d <- complete(mice_out, action = i)
  d$log_wbc        <- log10(d$wbc)
  d$log_creatinine <- log10(d$creatinine)
  # Re-add ordered/nominal outcome variants
  d$amp_group_ord  <- factor(d$amp_group_nom,
                              levels = c("Minor","Both","Major"), ordered=TRUE)
  d$major_binary   <- factor(
    ifelse(d$amp_group_nom %in% c("Major","Both"), "Yes","No"),
    levels = c("No","Yes")
  )
  d
})

cat(sprintf("\n  Imputed datasets created: %d\n", length(imputed_list)))
cat(sprintf("  Each dataset: %d rows x %d columns\n",
            nrow(imputed_list[[1]]), ncol(imputed_list[[1]])))

# ---- 3.7  Inspect imputed values vs observed --------------------------------

cat("\n  Imputed vs observed HbA1C comparison:\n")
obs_hba1c  <- df$hba1c[!is.na(df$hba1c)]
imp_hba1c  <- unlist(lapply(imputed_list, function(d) {
  d$hba1c[is.na(df$hba1c)]
}))

cat(sprintf("  Observed  — Mean: %.2f, SD: %.2f, Range: %.1f–%.1f\n",
            mean(obs_hba1c), sd(obs_hba1c),
            min(obs_hba1c), max(obs_hba1c)))
cat(sprintf("  Imputed   — Mean: %.2f, SD: %.2f, Range: %.1f–%.1f\n",
            mean(imp_hba1c), sd(imp_hba1c),
            min(imp_hba1c), max(imp_hba1c)))

# ---- 3.8  Save outputs ------------------------------------------------------

saveRDS(mice_out,      "outputs/mice_imputed.rds")
saveRDS(imputed_list,  "outputs/imputed_list.rds")   # list of 20 completed dfs
write.csv(diag_df,     "outputs/mice_diagnostics.csv", row.names = FALSE)

cat("\n  Saved: outputs/mice_imputed.rds\n")
cat("  Saved: outputs/imputed_list.rds (list of 20 completed datasets)\n")
cat("  Saved: outputs/mice_diagnostics.csv\n")
