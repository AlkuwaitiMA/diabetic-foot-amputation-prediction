# =============================================================================
# STEP 7: BINARY LOGISTIC MODELS — BOOTSTRAP VALIDATION + NOMOGRAM
# =============================================================================
# Input : outputs/imputed_list.rds
#         outputs/lasso_selected_vars.rds
# Output: outputs/binary_model_major.rds       — rms lrm fit
#         outputs/binary_model_icu.rds          — rms lrm fit
#         outputs/bootstrap_validation.csv      — optimism-corrected stats
#         outputs/calibration_data_major.rds    — for calibration plot
#         outputs/nomogram_major.rds             — rms nomogram object
#         outputs/binary_coefficients_major.csv
#         outputs/binary_coefficients_icu.csv
# =============================================================================

suppressPackageStartupMessages({
  library(rms)
  library(pROC)
  library(dplyr)
  library(mice)
})

set.seed(2024)

# Required by rms for ordered factor variables
options(contrasts = c("contr.treatment", "contr.treatment"))

imputed_list  <- readRDS("outputs/imputed_list.rds")
selected_vars <- readRDS("outputs/lasso_selected_vars.rds")

cat(sprintf("  Loaded %d imputed datasets\n", length(imputed_list)))

# ---- 7.1  Prepare pooled dataset for rms (use Rubin's rules manually) -------
# rms::lrm does not directly support mice; we fit on each imputation and pool

# Use imputed dataset 1 for rms (representative) + bootstrap validation
# Then pool coefficients across all 20 imputations via Rubin's rules

dat1 <- imputed_list[[1]]

# ---- 7.2  MODEL A: Major Amputation (binary) --------------------------------
# Outcome: major_binary (Yes = Major or Both; No = Minor)
# Predictors: LASSO-selected variables

cat("\n  === MODEL A: Major Amputation (binary) ===\n")

# Check events per variable
n_events <- sum(dat1$major_binary == "Yes", na.rm=TRUE)
n_pred   <- length(selected_vars)
cat(sprintf("  Events: %d (EPV = %.1f)\n", n_events, n_events/n_pred))

# Set rms data distribution for imputation 1
dd <- datadist(dat1)
options(datadist = "dd")

form_major <- as.formula(paste(
  "major_binary == 'Yes' ~",
  paste(selected_vars, collapse = " + ")
))

# ---- 7.2a  Fit lrm on all 20 imputations; pool coefficients ----------------

cat("  Fitting lrm on 20 imputations for pooling...\n")

lrm_fits <- lapply(imputed_list, function(d) {
  options(contrasts = c("contr.treatment", "contr.treatment"))
  assign("dd_lrm", datadist(d), envir=.GlobalEnv)
  options(datadist = "dd_lrm")
  tryCatch(
    lrm(form_major, data=d, x=TRUE, y=TRUE),
    error=function(e) NULL
  )
})

success <- !sapply(lrm_fits, is.null)
cat(sprintf("  Converged: %d/20\n", sum(success)))

# Rubin's rules for binary model
pool_binary <- function(fits_list) {
  coef_mat <- do.call(rbind, lapply(fits_list, coef))
  vcov_list <- lapply(fits_list, vcov)
  m   <- length(fits_list)
  Q   <- colMeans(coef_mat)
  U   <- Reduce("+", lapply(vcov_list, diag)) / m
  B   <- apply(coef_mat, 2, var)
  T_v <- U + (1 + 1/m) * B
  list(coef=Q, se=sqrt(T_v), var=T_v)
}

pooled_major <- pool_binary(lrm_fits[success])

major_coef_df <- data.frame(
  term      = names(pooled_major$coef),
  estimate  = round(pooled_major$coef, 4),
  std.error = round(pooled_major$se, 4),
  OR        = round(exp(pooled_major$coef), 3),
  OR_lo95   = round(exp(pooled_major$coef - 1.96*pooled_major$se), 3),
  OR_hi95   = round(exp(pooled_major$coef + 1.96*pooled_major$se), 3),
  z_stat    = round(pooled_major$coef / pooled_major$se, 3),
  p_value   = round(2*pnorm(abs(pooled_major$coef/pooled_major$se),
                            lower.tail=FALSE), 4)
) %>%
  filter(term != "Intercept")

cat("\n  Pooled Major Amputation Model — ORs:\n")
print(major_coef_df)

# ---- 7.2b  Bootstrap internal validation (imputation 1) --------------------
cat("\n  Bootstrap internal validation (B=500, imputation 1)...\n")

lrm_fit1 <- lrm(form_major, data=dat1, x=TRUE, y=TRUE)
val_major <- validate(lrm_fit1, B=500, method="boot")

cat("  Bootstrap validation results:\n")
print(val_major)

# Extract optimism-corrected c-statistic
# Dxy = 2*(c - 0.5); c = Dxy/2 + 0.5
Dxy_corrected  <- val_major["Dxy", "index.corrected"]
c_corrected    <- round(Dxy_corrected/2 + 0.5, 3)
Brier_apparent <- val_major["B", "index.orig"]
Brier_corrected<- val_major["B", "index.corrected"]
R2_corrected   <- val_major["R2", "index.corrected"]

cat(sprintf("  Optimism-corrected c-statistic: %.3f\n", c_corrected))
cat(sprintf("  Optimism-corrected R2: %.3f\n", R2_corrected))
cat(sprintf("  Brier score (corrected): %.4f\n", Brier_corrected))

# ---- 7.2c  Calibration curve ------------------------------------------------
cat("\n  Computing calibration curve (B=500)...\n")
cal_major <- calibrate(lrm_fit1, B=500, method="boot")

# ---- 7.2d  Apparent AUROC with 95% CI (pROC) --------------------------------
pred_probs <- predict(lrm_fit1, type="fitted")
roc_major  <- roc(dat1$major_binary, pred_probs,
                   levels=c("No","Yes"), direction="<", quiet=TRUE)
auc_major  <- round(auc(roc_major), 3)
auc_ci     <- round(ci.auc(roc_major), 3)

cat(sprintf("  Apparent AUROC: %.3f (95%% CI: %.3f–%.3f)\n",
            auc_major, auc_ci[1], auc_ci[3]))

# ---- 7.2e  Nomogram ---------------------------------------------------------
cat("  Building nomogram...\n")

nom_major <- nomogram(lrm_fit1,
                       fun        = plogis,
                       fun.at     = c(0.05, 0.10, 0.20, 0.30,
                                      0.40, 0.50, 0.60, 0.70, 0.80),
                       funlabel   = "P(Major Amputation)")

# ---- 7.3  MODEL B: ICU Admission (binary) -----------------------------------

cat("\n  === MODEL B: ICU Admission (binary) ===\n")

n_icu <- sum(dat1$icu_admission == "Yes", na.rm=TRUE)
cat(sprintf("  ICU events: %d\n", n_icu))

# For ICU model — use clinically plausible predictors regardless of LASSO
icu_vars <- c("log_wbc", "log_creatinine", "hb", "major_binary",
              "cause", "cvd", "htn", "ckd", "dm_nephro",
              "osteomyelitis", "prev_surgery")

# Keep only those available
icu_vars <- intersect(icu_vars, names(dat1))

form_icu <- as.formula(paste(
  "icu_admission == 'Yes' ~",
  paste(icu_vars, collapse=" + ")
))

lrm_fits_icu <- lapply(imputed_list, function(d) {
  options(contrasts = c("contr.treatment", "contr.treatment"))
  assign("dd_lrm", datadist(d), envir=.GlobalEnv)
  options(datadist = "dd_lrm")
  tryCatch(
    lrm(form_icu, data=d, x=TRUE, y=TRUE),
    error=function(e) NULL
  )
})

success_icu <- !sapply(lrm_fits_icu, is.null)
cat(sprintf("  ICU model converged: %d/20\n", sum(success_icu)))

pooled_icu <- pool_binary(lrm_fits_icu[success_icu])

icu_coef_df <- data.frame(
  term      = names(pooled_icu$coef),
  estimate  = round(pooled_icu$coef, 4),
  std.error = round(pooled_icu$se, 4),
  OR        = round(exp(pooled_icu$coef), 3),
  OR_lo95   = round(exp(pooled_icu$coef - 1.96*pooled_icu$se), 3),
  OR_hi95   = round(exp(pooled_icu$coef + 1.96*pooled_icu$se), 3),
  z_stat    = round(pooled_icu$coef / pooled_icu$se, 3),
  p_value   = round(2*pnorm(abs(pooled_icu$coef/pooled_icu$se),
                            lower.tail=FALSE), 4)
) %>%
  filter(term != "Intercept")

cat("\n  Pooled ICU Admission Model — ORs:\n")
print(icu_coef_df)

# Bootstrap validation for ICU model
lrm_fit_icu1 <- lrm(form_icu, data=dat1, x=TRUE, y=TRUE)
val_icu      <- validate(lrm_fit_icu1, B=500, method="boot")

Dxy_icu    <- val_icu["Dxy","index.corrected"]
c_icu      <- round(Dxy_icu/2 + 0.5, 3)
cat(sprintf("  ICU optimism-corrected c-statistic: %.3f\n", c_icu))

# ---- 7.4  Hosmer-Lemeshow goodness-of-fit -----------------------------------

hl_test <- function(fit, data, outcome_col, g=10) {
  pred  <- predict(fit, type="fitted")
  obs   <- as.integer(data[[outcome_col]] == "Yes")
  groups <- cut(pred, breaks=quantile(pred, probs=seq(0,1,1/g)),
                include.lowest=TRUE)
  tab   <- data.frame(pred=pred, obs=obs, grp=groups) %>%
    group_by(grp) %>%
    summarise(n=n(), obs_events=sum(obs),
              exp_events=sum(pred), .groups="drop")
  hl_stat <- sum((tab$obs_events - tab$exp_events)^2 / tab$exp_events)
  hl_p    <- pchisq(hl_stat, df=g-2, lower.tail=FALSE)
  list(stat=round(hl_stat,3), df=g-2, p=round(hl_p,4))
}

hl_major <- hl_test(lrm_fit1, dat1, "major_binary")
hl_icu   <- hl_test(lrm_fit_icu1, dat1, "icu_admission")

cat(sprintf("\n  Hosmer-Lemeshow (Major): Chi2=%.3f, df=%d, p=%.4f\n",
            hl_major$stat, hl_major$df, hl_major$p))
cat(sprintf("  Hosmer-Lemeshow (ICU):   Chi2=%.3f, df=%d, p=%.4f\n",
            hl_icu$stat, hl_icu$df, hl_icu$p))

# ---- 7.5  Compile validation summary ----------------------------------------

validation_summary <- data.frame(
  Model             = c("Major Amputation", "ICU Admission"),
  N_events          = c(n_events, n_icu),
  Apparent_c        = c(round(auc_major, 3),
                        round(auc(roc(dat1$icu_admission,
                          predict(lrm_fit_icu1, type="fitted"),
                          levels=c("No","Yes"), direction="<", quiet=TRUE)), 3)),
  Corrected_c       = c(c_corrected, c_icu),
  Optimism          = c(round(auc_major - c_corrected, 3),
                        NA),
  Brier_corrected   = c(round(Brier_corrected, 4), NA),
  R2_corrected      = c(round(R2_corrected, 3), NA),
  HL_chi2           = c(hl_major$stat, hl_icu$stat),
  HL_p              = c(hl_major$p, hl_icu$p)
)

cat("\n  Validation Summary:\n")
print(validation_summary)

# ---- 7.6  Save outputs -------------------------------------------------------

saveRDS(lrm_fit1,     "outputs/binary_model_major.rds")
saveRDS(lrm_fit_icu1, "outputs/binary_model_icu.rds")
saveRDS(val_major,    "outputs/bootstrap_val_major.rds")
saveRDS(cal_major,    "outputs/calibration_data_major.rds")
saveRDS(nom_major,    "outputs/nomogram_major.rds")
saveRDS(roc_major,    "outputs/roc_major.rds")

write.csv(major_coef_df,      "outputs/binary_coefficients_major.csv", row.names=FALSE)
write.csv(icu_coef_df,        "outputs/binary_coefficients_icu.csv",   row.names=FALSE)
write.csv(validation_summary, "outputs/bootstrap_validation.csv",       row.names=FALSE)

cat("\n  Saved: all binary model outputs to outputs/\n")
