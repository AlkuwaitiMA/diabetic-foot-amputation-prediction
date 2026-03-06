# =============================================================================
# STEP 5: LASSO VARIABLE SELECTION
# =============================================================================
# Input : outputs/imputed_list.rds
# Output: outputs/lasso_selected_vars.rds  — character vector of selected vars
#         outputs/lasso_lambda_results.csv  — lambda, deviance, vars retained
#         outputs/lasso_coefficients.csv    — non-zero coefficients at lambda.1se
# =============================================================================
# Method: 10-fold cross-validated multinomial LASSO via glmnet.
#         Applied to ALL 20 imputed datasets; a variable is retained if it
#         has a non-zero coefficient in >= 50% of imputed datasets (Rubin rule
#         analog for variable selection — conservative and stable).
#         lambda.1se used (more parsimonious than lambda.min).
# =============================================================================

suppressPackageStartupMessages({
  library(glmnet)
  library(dplyr)
})

set.seed(42)

imputed_list <- readRDS("outputs/imputed_list.rds")

cat(sprintf("  Loaded %d imputed datasets\n", length(imputed_list)))

# ---- 5.1  Define predictor set for LASSO ------------------------------------

# These are the candidate predictors entering LASSO
# (outcome: amp_group_nom — 3 levels: Minor, Both, Major)

predictor_vars <- c(
  # Demographics
  "gender", "age_group", "nationality",
  # DM
  "dm_type", "dm_duration", "oha_use", "insulin_use",
  # Comorbidities
  "cvd", "htn", "dm_neuro", "dm_nephro", "ckd", "malignancy",
  # Gangrene / clinical
  "dry_gangrene", "wet_gangrene", "cause",
  # Admission
  "site_admission",  "dm_admissions",
  # History
  "osteomyelitis", "prev_surgery",
  # Labs (log-transformed continuous)
  "hba1c", "log_wbc", "log_creatinine", "hb"
)

outcome_var <- "amp_group_nom"

# ---- 5.2  Helper: build model matrix from one imputed dataset ---------------

build_matrix <- function(data, pred_vars, outcome) {
  # Convert outcome to factor
  y <- factor(data[[outcome]], levels = c("Minor","Both","Major"))

  # Build model matrix (dummy-codes factors, scales continuous vars)
  form  <- as.formula(paste("~", paste(pred_vars, collapse = " + ")))
  X_raw <- model.matrix(form, data = data)[, -1]   # drop intercept

  # Standardise continuous columns
  cont_cols <- c("hba1c","log_wbc","log_creatinine","hb")
  for (col in colnames(X_raw)) {
    if (any(sapply(cont_cols, function(cc) grepl(cc, col)))) {
      X_raw[, col] <- scale(X_raw[, col])
    }
  }

  list(X = X_raw, y = y)
}

# ---- 5.3  Run LASSO on each imputed dataset ---------------------------------

cat("\n  Running 10-fold CV LASSO on 20 imputed datasets...\n")

lasso_results <- lapply(seq_along(imputed_list), function(i) {

  dat <- imputed_list[[i]]
  dat <- dat[complete.cases(dat[, c(predictor_vars, outcome_var)]), ]

  mat <- build_matrix(dat, predictor_vars, outcome_var)

  # 10-fold cross-validated LASSO (multinomial)
  cv_fit <- cv.glmnet(
    x             = mat$X,
    y             = mat$y,
    family        = "multinomial",
    alpha         = 1,            # LASSO
    nfolds        = 10,
    type.measure  = "class",      # misclassification error
    standardize   = FALSE         # already standardised above
  )

  # Extract non-zero coefficients at lambda.1se
  lambda_use <- cv_fit$lambda.1se

  coefs <- coef(cv_fit, s = "lambda.1se")
  # coefs is a list (one per outcome class); variable is selected if non-zero
  # in ANY class comparison
  nonzero_vars <- unique(unlist(lapply(coefs, function(co) {
    rn <- rownames(co)[co[,1] != 0]
    rn[rn != "(Intercept)"]
  })))

  list(
    imputation   = i,
    lambda_min   = cv_fit$lambda.min,
    lambda_1se   = lambda_use,
    cvm_min      = min(cv_fit$cvm),
    n_vars_min   = length(unique(unlist(lapply(
      coef(cv_fit, s="lambda.min"), function(co)
        rownames(co)[co[,1] != 0 & rownames(co) != "(Intercept)"])))),
    n_vars_1se   = length(nonzero_vars),
    selected     = nonzero_vars
  )
})

# ---- 5.4  Variable selection across imputations (Rubin analog) --------------

# Count how many imputed datasets each variable was selected in
all_selected <- unlist(lapply(lasso_results, function(r) r$selected))
var_counts   <- sort(table(all_selected), decreasing = TRUE)

cat("\n  Variable selection frequency across 20 imputations:\n")
print(var_counts)

# Retain variables selected in >= 50% of imputations (>= 10 out of 20)
n_imp        <- length(imputed_list)
threshold    <- ceiling(n_imp * 0.50)   # 10

selected_vars <- names(var_counts[var_counts >= threshold])

cat(sprintf("\n  Variables selected (>= %d/%d imputations = >= 50%%):\n",
            threshold, n_imp))
print(selected_vars)

# ---- 5.5  Extract coefficients from median lambda model ---------------------

# Re-fit on first imputed dataset at lambda.1se to get coefficient magnitudes
dat1 <- imputed_list[[1]]
dat1 <- dat1[complete.cases(dat1[, c(predictor_vars, outcome_var)]), ]
mat1 <- build_matrix(dat1, predictor_vars, outcome_var)

median_lambda <- median(sapply(lasso_results, function(r) r$lambda_1se))

fit1 <- glmnet(mat1$X, mat1$y, family = "multinomial", alpha = 1,
               lambda = median_lambda)

# Collect non-zero coefficients per class
coef_list <- coef(fit1)
coef_df <- lapply(names(coef_list), function(cls) {
  co <- as.matrix(coef_list[[cls]])
  df_co <- data.frame(
    variable  = rownames(co),
    class     = cls,
    coef      = round(co[,1], 4)
  )
  df_co[df_co$variable != "(Intercept)" & df_co$coef != 0, ]
}) %>% bind_rows()

cat("\n  Non-zero LASSO coefficients (median lambda.1se, imputation 1):\n")
print(coef_df)

# ---- 5.6  Lambda selection summary ------------------------------------------

lambda_summary <- data.frame(
  imputation  = sapply(lasso_results, function(r) r$imputation),
  lambda_min  = round(sapply(lasso_results, function(r) r$lambda_min), 5),
  lambda_1se  = round(sapply(lasso_results, function(r) r$lambda_1se), 5),
  n_vars_1se  = sapply(lasso_results, function(r) r$n_vars_1se),
  cvm_min     = round(sapply(lasso_results, function(r) r$cvm_min), 4)
)

cat("\n  Lambda selection summary:\n")
print(lambda_summary)
cat(sprintf("\n  Median lambda.1se: %.5f\n", median_lambda))
cat(sprintf("  Median vars at 1se: %d\n",
            median(lambda_summary$n_vars_1se)))

# ---- 5.7  Final selected variable list --------------------------------------

# Map dummy-coded names back to original variable names
# e.g., "genderMale" → "gender", "causeInfection" → "cause"
original_selected <- unique(sapply(selected_vars, function(v) {
  # Try to match to original predictor names
  matched <- predictor_vars[sapply(predictor_vars, function(p) grepl(paste0("^", p), v))]
  if (length(matched) > 0) matched[1] else v
}))

cat("\n  Final original variable names selected by LASSO:\n")
print(original_selected)

# ---- 5.8  Save outputs ------------------------------------------------------

write.csv(lambda_summary, "outputs/lasso_lambda_results.csv", row.names = FALSE)
write.csv(coef_df,        "outputs/lasso_coefficients.csv",   row.names = FALSE)
write.csv(
  data.frame(variable = names(var_counts), n_selected = as.integer(var_counts),
             pct_selected = round(100*as.integer(var_counts)/n_imp, 1)),
  "outputs/lasso_selection_frequency.csv", row.names = FALSE
)

saveRDS(original_selected,   "outputs/lasso_selected_vars.rds")
saveRDS(selected_vars,       "outputs/lasso_selected_dummies.rds")
saveRDS(lasso_results,       "outputs/lasso_all_results.rds")

cat("\n  Saved: outputs/lasso_selected_vars.rds (original variable names)\n")
cat("  Saved: outputs/lasso_selected_dummies.rds (dummy-coded names)\n")
cat("  Saved: outputs/lasso_lambda_results.csv\n")
cat("  Saved: outputs/lasso_coefficients.csv\n")
cat("  Saved: outputs/lasso_selection_frequency.csv\n")
