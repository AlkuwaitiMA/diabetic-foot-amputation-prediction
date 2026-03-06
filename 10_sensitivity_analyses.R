# =============================================================================
# STEP 10: SENSITIVITY ANALYSES (FIXED)
# =============================================================================
# SA-1: Complete case vs MICE-imputed (lrm with contrasts fix)
# SA-2: Alternative outcome orderings (POM A vs B — fixed concordance)
# SA-3: Exclude 2024 partial year
# =============================================================================

suppressPackageStartupMessages({
  library(ordinal)
  library(nnet)
  library(rms)
  library(dplyr)
  library(pROC)
})

options(contrasts = c("contr.treatment", "contr.treatment"))
set.seed(99)

df_clean      <- readRDS("outputs/df_clean.rds")
df_complete   <- readRDS("outputs/df_complete.rds")
imputed_list  <- readRDS("outputs/imputed_list.rds")
selected_vars <- readRDS("outputs/lasso_selected_vars.rds")
selected_vars <- selected_vars[selected_vars != "icu_admission"]  # guard

main_coefs <- read.csv("outputs/binary_coefficients_major.csv")

cat("  Sensitivity analyses starting...\n")
cat(sprintf("  Predictors used: %s\n", paste(selected_vars, collapse=", ")))

avail_cc  <- intersect(selected_vars, names(df_complete))
avail_imp <- intersect(selected_vars, names(imputed_list[[1]]))

# ============================================================
# SA-1: Complete case vs MICE
# ============================================================
cat("\n  === SA-1: Complete Case vs MICE Imputation ===\n")

form_cc_ord <- as.formula(paste("amp_group_ord ~", paste(avail_cc, collapse=" + ")))
form_cc_bin <- as.formula(paste("major_binary == 'Yes' ~",
                                 paste(avail_cc, collapse=" + ")))

# Ordinal model on complete cases
cc_pom <- tryCatch(
  clm(form_cc_ord, data=df_complete, link="logit"),
  error=function(e) { cat("  CC POM error:", e$message, "\n"); NULL }
)

# Binary logistic on complete cases — with datadist in global env
assign("dd_cc", datadist(df_complete), envir=.GlobalEnv)
options(datadist="dd_cc")
cc_lrm <- tryCatch(
  lrm(form_cc_bin, data=df_complete, x=TRUE, y=TRUE),
  error=function(e) { cat("  CC LRM error:", e$message, "\n"); NULL }
)

if (!is.null(cc_lrm)) {
  cc_coef_df <- data.frame(
    term      = names(coef(cc_lrm)),
    CC_OR     = round(exp(coef(cc_lrm)), 3),
    CC_OR_lo  = round(exp(coef(cc_lrm) - 1.96*sqrt(diag(vcov(cc_lrm)))), 3),
    CC_OR_hi  = round(exp(coef(cc_lrm) + 1.96*sqrt(diag(vcov(cc_lrm)))), 3)
  ) %>% filter(term != "Intercept")

  comparison_cc <- merge(
    cc_coef_df,
    main_coefs[, c("term","OR","OR_lo95","OR_hi95")],
    by="term", all=TRUE
  ) %>% mutate(
    OR_ratio   = round(CC_OR / OR, 3),
    CI_overlap = (!is.na(CC_OR_lo) & !is.na(OR_hi95)) &
                 (CC_OR_lo <= OR_hi95) & (CC_OR_hi >= OR_lo95)
  )

  n_concordant_sa1 <- sum(comparison_cc$CI_overlap, na.rm=TRUE)
  n_compared_sa1   <- sum(!is.na(comparison_cc$CI_overlap))
  pct_sa1          <- round(100 * n_concordant_sa1 / n_compared_sa1, 1)

  cat(sprintf("  SA-1: %d/%d terms concordant (%.1f%%)\n",
              n_concordant_sa1, n_compared_sa1, pct_sa1))
  print(comparison_cc)
  write.csv(comparison_cc, "outputs/sensitivity_complete_case.csv", row.names=FALSE)
} else {
  cat("  CC LRM still failed — writing NA placeholder\n")
  comparison_cc <- data.frame(note="CC model failed to converge")
  n_concordant_sa1 <- NA; n_compared_sa1 <- NA; pct_sa1 <- NA
  write.csv(comparison_cc, "outputs/sensitivity_complete_case.csv", row.names=FALSE)
}

# ============================================================
# SA-2: Alternative outcome orderings
# ============================================================
cat("\n  === SA-2: Alternative Outcome Orderings ===\n")

dat1 <- imputed_list[[1]]
dat1$amp_ord_A <- factor(dat1$amp_group_nom,
                          levels=c("Minor","Both","Major"), ordered=TRUE)
dat1$amp_ord_B <- factor(dat1$amp_group_nom,
                          levels=c("Minor","Major","Both"), ordered=TRUE)

form_A <- as.formula(paste("amp_ord_A ~", paste(avail_imp, collapse=" + ")))
form_B <- as.formula(paste("amp_ord_B ~", paste(avail_imp, collapse=" + ")))

pom_A <- tryCatch(clm(form_A, data=dat1, link="logit"), error=function(e) NULL)
pom_B <- tryCatch(clm(form_B, data=dat1, link="logit"), error=function(e) NULL)

extract_clm_sig <- function(fit, label) {
  if (is.null(fit)) return(data.frame())
  co <- coef(fit)
  se <- sqrt(diag(vcov(fit)))
  keep <- !grepl("\\|", names(co))
  z <- co[keep] / se[keep]
  p <- 2 * pnorm(abs(z), lower.tail=FALSE)
  data.frame(
    term = names(co)[keep],
    OR   = round(exp(co[keep]), 3),
    p    = round(p, 4),
    sig  = p < 0.05,
    model = label
  )
}

res_A <- extract_clm_sig(pom_A, "A_Minor_Both_Major")
res_B <- extract_clm_sig(pom_B, "B_Minor_Major_Both")

if (nrow(res_A) > 0 & nrow(res_B) > 0) {
  sig_A <- res_A %>% select(term, sig_A = sig, OR_A = OR, p_A = p)
  sig_B <- res_B %>% select(term, sig_B = sig, OR_B = OR, p_B = p)

  concordance <- merge(sig_A, sig_B, by="term", all=TRUE) %>%
    mutate(
      concordant     = !is.na(sig_A) & !is.na(sig_B) & (sig_A == sig_B),
      direction_same = !is.na(OR_A) & !is.na(OR_B) &
                       (sign(log(OR_A)) == sign(log(OR_B)))
    )

  n_conc  <- sum(concordance$concordant, na.rm=TRUE)
  n_dir   <- sum(concordance$direction_same, na.rm=TRUE)
  n_terms <- nrow(concordance)

  cat(sprintf("  Significance concordance: %d/%d terms (%.1f%%)\n",
              n_conc, n_terms, 100*n_conc/n_terms))
  cat(sprintf("  Direction concordance:    %d/%d terms (%.1f%%)\n",
              n_dir, n_terms, 100*n_dir/n_terms))
  print(concordance)

  write.csv(concordance, "outputs/sensitivity_ordering_concordance.csv", row.names=FALSE)

  ordering_df <- bind_rows(res_A, res_B)
  write.csv(ordering_df, "outputs/sensitivity_outcome_ordering.csv", row.names=FALSE)

  n_concordant_sa2 <- n_conc
  n_compared_sa2   <- n_terms
  pct_sa2          <- round(100 * n_conc / n_terms, 1)
  dir_sa2          <- round(100 * n_dir  / n_terms, 1)
} else {
  cat("  One or both ordering models failed\n")
  n_concordant_sa2 <- NA; n_compared_sa2 <- NA
  pct_sa2 <- NA; dir_sa2 <- NA
}

# ============================================================
# SA-3: Exclude 2024 (partial year)
# ============================================================
cat("\n  === SA-3: Excluding 2024 (Partial Year) ===\n")

imputed_no2024 <- lapply(imputed_list, function(d) d[d$year != 2024, ])
dat1_no2024    <- imputed_no2024[[1]]

cat(sprintf("  N without 2024: %d (excluded %d)\n",
            nrow(dat1_no2024),
            nrow(imputed_list[[1]]) - nrow(dat1_no2024)))
cat(sprintf("  Outcome: %s\n", paste(table(dat1_no2024$major_binary), collapse="/")))

form_bin_ex <- as.formula(paste("major_binary == 'Yes' ~",
                                 paste(avail_imp, collapse=" + ")))

lrm_fits_ex <- lapply(imputed_no2024, function(d) {
  options(contrasts = c("contr.treatment","contr.treatment"))
  assign("dd_lrm", datadist(d), envir=.GlobalEnv)
  options(datadist="dd_lrm")
  tryCatch(lrm(form_bin_ex, data=d, x=TRUE, y=TRUE), error=function(e) NULL)
})

success_ex <- !sapply(lrm_fits_ex, is.null)
cat(sprintf("  Converged: %d/20\n", sum(success_ex)))

if (sum(success_ex) >= 1) {
  # Rubin's rules
  coef_mat  <- do.call(rbind, lapply(lrm_fits_ex[success_ex], coef))
  if (is.null(dim(coef_mat))) coef_mat <- matrix(coef_mat, nrow=1)
  vcov_list <- lapply(lrm_fits_ex[success_ex], function(f)
    tryCatch(diag(vcov(f)), error=function(e) NULL))
  valid     <- !sapply(vcov_list, is.null)
  coef_mat  <- coef_mat[valid, , drop=FALSE]
  vcov_list <- vcov_list[valid]
  m  <- nrow(coef_mat)
  Q  <- colMeans(coef_mat)
  U  <- Reduce("+", vcov_list) / m
  B  <- if (m > 1) apply(coef_mat, 2, var) else rep(0, length(Q))
  Tv <- U + (1 + 1/m)*B

  ex2024_coefs <- data.frame(
    term      = names(Q),
    OR_ex2024 = round(exp(Q), 3),
    OR_lo     = round(exp(Q - 1.96*sqrt(Tv)), 3),
    OR_hi     = round(exp(Q + 1.96*sqrt(Tv)), 3)
  ) %>% filter(term != "Intercept")

  comparison_ex <- merge(
    ex2024_coefs,
    main_coefs[, c("term","OR","OR_lo95","OR_hi95")],
    by="term", all=TRUE
  ) %>% mutate(
    OR_ratio   = round(OR_ex2024 / OR, 3),
    CI_overlap = (!is.na(OR_lo) & !is.na(OR_hi95)) &
                 (OR_lo <= OR_hi95) & (OR_hi >= OR_lo95)
  )

  n_concordant_sa3 <- sum(comparison_ex$CI_overlap, na.rm=TRUE)
  n_compared_sa3   <- sum(!is.na(comparison_ex$CI_overlap))
  pct_sa3          <- round(100*n_concordant_sa3/n_compared_sa3, 1)

  cat(sprintf("  SA-3: %d/%d terms concordant (%.1f%%)\n",
              n_concordant_sa3, n_compared_sa3, pct_sa3))
  print(comparison_ex)
  write.csv(comparison_ex, "outputs/sensitivity_excl_2024.csv", row.names=FALSE)
} else {
  cat("  No models converged for excl-2024\n")
  n_concordant_sa3 <- NA; n_compared_sa3 <- NA; pct_sa3 <- NA
}

# ============================================================
# Robustness summary
# ============================================================

robustness_summary <- data.frame(
  Sensitivity_Analysis = c(
    "SA-1: Complete case vs MICE",
    "SA-2: Outcome ordering (significance concordance)",
    "SA-2: Outcome ordering (direction concordance)",
    "SA-3: Exclude 2024 partial year"
  ),
  N_compared    = c(n_compared_sa1,   n_compared_sa2, n_compared_sa2, n_compared_sa3),
  N_concordant  = c(n_concordant_sa1, n_concordant_sa2, ifelse(!is.na(dir_sa2), round(dir_sa2/100*n_compared_sa2), NA), n_concordant_sa3),
  Pct_concordant= c(pct_sa1, pct_sa2, dir_sa2, pct_sa3),
  Conclusion    = c(
    ifelse(!is.na(pct_sa1) & pct_sa1 >= 80,
           "Robust — CIs overlap in ≥80% of terms",
           ifelse(is.na(pct_sa1), "CC model failed", "Divergence — check")),
    ifelse(!is.na(pct_sa2) & pct_sa2 >= 80,
           "Robust — significant predictors stable across orderings",
           ifelse(is.na(pct_sa2), "Could not assess", "Moderate — some ordering sensitivity")),
    ifelse(!is.na(dir_sa2) & dir_sa2 >= 80,
           "Robust — effect directions stable across orderings",
           ifelse(is.na(dir_sa2), "Could not assess", "Some direction reversals")),
    ifelse(!is.na(pct_sa3) & pct_sa3 >= 80,
           "Robust — 2024 exclusion does not change conclusions",
           ifelse(is.na(pct_sa3), "Could not assess", "2024 influences some estimates"))
  )
)

cat("\n  === ROBUSTNESS SUMMARY ===\n")
print(robustness_summary)
write.csv(robustness_summary, "outputs/sensitivity_summary.csv", row.names=FALSE)

cat("\n  Saved: all sensitivity outputs to outputs/\n")
