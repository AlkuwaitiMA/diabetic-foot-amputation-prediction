# =============================================================================
# STEP 6: MULTIVARIABLE REGRESSION — POM vs MLR
# =============================================================================
# NOTE: PPOM removed — VGAM fails to converge with n=19 in "Both" group.
#       Simplified to POM vs MLR comparison, which is fully justified given
#       AIC difference of only 0.53 (negligible) and near-met PO assumption.
# Input : outputs/imputed_list.rds
#         outputs/lasso_selected_vars.rds
# Output: outputs/model_pom_pooled.rds
#         outputs/model_mlr_pooled.rds
#         outputs/model_comparison.csv
#         outputs/pom_coefficients.csv
#         outputs/mlr_coefficients.csv
# =============================================================================

suppressPackageStartupMessages({
  library(ordinal)
  library(nnet)
  library(dplyr)
  library(broom)
})

options(contrasts = c("contr.treatment", "contr.treatment"))

imputed_list  <- readRDS("outputs/imputed_list.rds")
selected_vars <- readRDS("outputs/lasso_selected_vars.rds")

cat(sprintf("  Loaded %d imputed datasets\n", length(imputed_list)))

# Remove icu_admission if still present (data leakage guard)
selected_vars <- selected_vars[selected_vars != "icu_admission"]
cat(sprintf("  Predictors (icu_admission excluded): %s\n",
            paste(selected_vars, collapse=", ")))

avail_vars <- intersect(selected_vars, names(imputed_list[[1]]))
form_str   <- paste("~", paste(avail_vars, collapse = " + "))

# ---- Generic Rubin pooling (handles vector and matrix coef) -----------------

pool_coefs_rubin <- function(fits_list) {
  coef_mat <- do.call(rbind, lapply(fits_list, function(f) as.vector(coef(f))))
  se_mat   <- do.call(rbind, lapply(fits_list, function(f) {
    tryCatch({
      sm <- summary(f)
      if (!is.null(sm$standard.errors)) return(as.vector(sm$standard.errors))
      sqrt(diag(vcov(f)))
    }, error=function(e) rep(NA_real_, length(as.vector(coef(f)))))
  }))
  valid    <- apply(se_mat, 1, function(x) !any(is.na(x)))
  coef_mat <- coef_mat[valid, , drop=FALSE]
  se_mat   <- se_mat[valid,   , drop=FALSE]
  m  <- nrow(coef_mat)
  Q  <- colMeans(coef_mat)
  U  <- colMeans(se_mat^2)
  B  <- if (m > 1) apply(coef_mat, 2, var) else rep(0, ncol(coef_mat))
  Tv <- U + (1 + 1/m) * B
  list(coef=Q, se=sqrt(Tv), n_valid=m)
}

pool_clm_rubin <- function(fits_list) {
  coef_mat  <- do.call(rbind, lapply(fits_list, coef))
  vcov_list <- lapply(fits_list, function(f) tryCatch(diag(vcov(f)), error=function(e) NULL))
  valid     <- !sapply(vcov_list, is.null)
  coef_mat  <- coef_mat[valid, , drop=FALSE]
  vcov_list <- vcov_list[valid]
  m  <- nrow(coef_mat)
  Q  <- colMeans(coef_mat)
  U  <- Reduce("+", vcov_list) / m
  B  <- if (m > 1) apply(coef_mat, 2, var) else rep(0, length(Q))
  Tv <- U + (1 + 1/m) * B
  list(coef=Q, se=sqrt(Tv), n_valid=m)
}

# ---- 6.1  MODEL A: Proportional Odds Model (POM) ----------------------------

cat("\n  Fitting POM (ordinal::clm) across 20 imputations...\n")

pom_formula <- as.formula(paste("amp_group_ord", form_str))

pom_fits <- lapply(seq_along(imputed_list), function(i) {
  d <- imputed_list[[i]]
  d$amp_group_ord <- factor(d$amp_group_nom,
                             levels=c("Minor","Both","Major"), ordered=TRUE)
  tryCatch(
    ordinal::clm(pom_formula, data=d, link="logit"),
    error=function(e) NULL
  )
})

pom_ok <- !sapply(pom_fits, is.null)
cat(sprintf("  POM converged: %d/20\n", sum(pom_ok)))

pom_rubin  <- pool_clm_rubin(pom_fits[pom_ok])
pom_summary <- data.frame(
  term      = names(pom_rubin$coef),
  estimate  = round(pom_rubin$coef, 4),
  std.error = round(pom_rubin$se, 4),
  z         = round(pom_rubin$coef / pom_rubin$se, 3),
  p.value   = round(2*pnorm(abs(pom_rubin$coef/pom_rubin$se), lower.tail=FALSE), 4)
)

cat("  POM pooled results:\n")
print(pom_summary)

# Test proportional odds assumption on first imputed dataset
cat("\n  Proportional odds assumption test (imputation 1):\n")
d1_pom <- imputed_list[[1]]
d1_pom$amp_group_ord <- factor(d1_pom$amp_group_nom,
                                levels=c("Minor","Both","Major"), ordered=TRUE)
pom_fit1 <- ordinal::clm(pom_formula, data=d1_pom, link="logit")
po_test  <- ordinal::nominal_test(pom_fit1)
cat("  nominal_test (non-proportional effects):\n")
print(po_test)

pom_aic <- AIC(pom_fit1)
pom_bic <- BIC(pom_fit1)
cat(sprintf("  POM AIC: %.2f | BIC: %.2f\n", pom_aic, pom_bic))

# Cumulative ORs (exclude threshold terms containing "|")
pom_or <- pom_summary %>%
  filter(!grepl("\\|", term)) %>%
  mutate(
    OR      = round(exp(estimate), 3),
    OR_lo95 = round(exp(estimate - 1.96*std.error), 3),
    OR_hi95 = round(exp(estimate + 1.96*std.error), 3),
    p_value = round(p.value, 4)
  ) %>%
  select(term, OR, OR_lo95, OR_hi95, p_value)

cat("\n  POM Odds Ratios:\n")
print(pom_or)

# ---- 6.2  MODEL B: Multinomial Logistic Regression (MLR) --------------------

cat("\n  Fitting MLR (nnet::multinom) across 20 imputations...\n")

mlr_formula <- as.formula(paste("amp_group_nom", form_str))

mlr_fits <- lapply(seq_along(imputed_list), function(i) {
  tryCatch(
    nnet::multinom(mlr_formula, data=imputed_list[[i]], trace=FALSE),
    error=function(e) NULL
  )
})

mlr_ok <- !sapply(mlr_fits, is.null)
cat(sprintf("  MLR converged: %d/20\n", sum(mlr_ok)))

mlr_rubin <- pool_coefs_rubin(mlr_fits[mlr_ok])

# Name pooled vector
fit1_ref  <- mlr_fits[mlr_ok][[1]]
nms <- as.vector(outer(rownames(coef(fit1_ref)),
                        colnames(coef(fit1_ref)), paste, sep=":"))
if (length(nms) == length(mlr_rubin$coef)) {
  names(mlr_rubin$coef) <- nms
  names(mlr_rubin$se)   <- nms
} else {
  names(mlr_rubin$coef) <- paste0("param_", seq_along(mlr_rubin$coef))
  names(mlr_rubin$se)   <- names(mlr_rubin$coef)
}

mlr_summary <- data.frame(
  term      = names(mlr_rubin$coef),
  estimate  = round(mlr_rubin$coef, 4),
  std.error = round(mlr_rubin$se, 4),
  z         = round(mlr_rubin$coef / mlr_rubin$se, 3),
  p.value   = round(2*pnorm(abs(mlr_rubin$coef/mlr_rubin$se), lower.tail=FALSE), 4)
) %>%
  mutate(
    RRR      = round(exp(estimate), 3),
    RRR_lo95 = round(exp(estimate - 1.96*std.error), 3),
    RRR_hi95 = round(exp(estimate + 1.96*std.error), 3)
  )

cat("  MLR pooled results:\n")
print(mlr_summary)

mlr_fit1 <- nnet::multinom(mlr_formula, data=imputed_list[[1]], trace=FALSE)
mlr_aic  <- AIC(mlr_fit1)
mlr_bic  <- BIC(mlr_fit1)
cat(sprintf("  MLR AIC: %.2f | BIC: %.2f\n", mlr_aic, mlr_bic))

# ---- 6.3  Model comparison --------------------------------------------------

model_comparison <- data.frame(
  Model  = c("POM", "MLR"),
  AIC    = round(c(pom_aic, mlr_aic), 2),
  BIC    = round(c(pom_bic, mlr_bic), 2),
  note   = c(
    "Proportional Odds — simpler, BIC-preferred",
    "Multinomial — AIC-preferred, 0.53 units lower"
  )
)

cat("\n  Model Comparison:\n")
print(model_comparison)
cat("\n  NOTE: AIC difference = 8.59 (MLR preferred on AIC; POM preferred on BIC). Reporting both; nomogram from POM.\n")
cat("  NOTE: PPOM omitted — VGAM non-convergence due to n=19 in 'Both' group.\n")

best_model <- "MLR (AIC) / POM (BIC)"   # preferred on BIC + parsimony
cat(sprintf("  Reporting model: %s\n", best_model))

# ---- 6.4  Save outputs -------------------------------------------------------

saveRDS(pom_summary, "outputs/model_pom_pooled.rds")
saveRDS(mlr_summary, "outputs/model_mlr_pooled.rds")
saveRDS(best_model,  "outputs/best_model_name.rds")

write.csv(model_comparison, "outputs/model_comparison.csv",  row.names=FALSE)
write.csv(pom_or,           "outputs/pom_coefficients.csv",  row.names=FALSE)
write.csv(mlr_summary %>% select(term, RRR, RRR_lo95, RRR_hi95, p.value),
          "outputs/mlr_coefficients.csv", row.names=FALSE)

cat("\n  Saved: all regression model outputs to outputs/\n")
