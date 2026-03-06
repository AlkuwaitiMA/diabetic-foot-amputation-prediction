# =============================================================================
# STEP 11: SENSITIVITY — BOTH+MAJOR COLLAPSE
# =============================================================================
# Purpose : Test whether the three-group ordinal structure is driving results,
#           or whether collapsing Both (n=19) into Major yields the same
#           clinical conclusions as the primary binary model.
#
# Approach: Define collapsed_binary = Major+Both (n=73) vs Minor (n=168)
#           — identical to major_binary in the primary analysis.
#           THEN also reframe as a strict binary:
#           strict_major = Major only (n=54) vs Minor+Both (n=187)
#           — tests whether Both patients are clinically more like Minor.
#
# Both comparisons use:
#   - Same 6 LASSO-selected predictors
#   - Same rms::lrm + bootstrap validation (B=500)
#   - Same Hosmer-Lemeshow calibration test
#   - Pooled across 20 MICE imputed datasets
#
# Output: outputs/sensitivity_collapse_both.csv      — coefficient comparison
#         outputs/sensitivity_collapse_validation.csv — performance metrics
#         outputs/sensitivity_collapse_summary.rds
# =============================================================================

suppressPackageStartupMessages({
  library(rms)
  library(pROC)
  library(dplyr)
})

set.seed(2024)
options(contrasts = c("contr.treatment", "contr.treatment"))

imputed_list  <- readRDS("outputs/imputed_list.rds")
selected_vars <- readRDS("outputs/lasso_selected_vars.rds")
selected_vars <- selected_vars[selected_vars != "icu_admission"]

main_coefs <- read.csv("outputs/binary_coefficients_major.csv")
main_val   <- read.csv("outputs/bootstrap_validation.csv")

cat("================================================================\n")
cat("  SA-4: Both+Major Collapse Sensitivity Analysis\n")
cat("================================================================\n\n")
cat(sprintf("  Predictors: %s\n", paste(selected_vars, collapse=", ")))
cat(sprintf("  N imputed datasets: %d\n\n", length(imputed_list)))

avail_vars <- intersect(selected_vars, names(imputed_list[[1]]))
form_str   <- paste(avail_vars, collapse=" + ")

# ---- Helper: Rubin pooling for lrm ------------------------------------------

pool_lrm_rubin <- function(fits_list) {
  coef_mat  <- do.call(rbind, lapply(fits_list, coef))
  if (is.null(dim(coef_mat))) coef_mat <- matrix(coef_mat, nrow=1)
  vcov_list <- lapply(fits_list, function(f)
    tryCatch(diag(vcov(f)), error=function(e) NULL))
  valid     <- !sapply(vcov_list, is.null)
  coef_mat  <- coef_mat[valid, , drop=FALSE]
  vcov_list <- vcov_list[valid]
  m  <- nrow(coef_mat)
  Q  <- colMeans(coef_mat)
  U  <- Reduce("+", vcov_list) / m
  B  <- if (m > 1) apply(coef_mat, 2, var) else rep(0, length(Q))
  Tv <- U + (1 + 1/m) * B
  data.frame(
    term      = names(Q),
    estimate  = round(Q, 4),
    std.error = round(sqrt(Tv), 4),
    OR        = round(exp(Q), 3),
    OR_lo95   = round(exp(Q - 1.96*sqrt(Tv)), 3),
    OR_hi95   = round(exp(Q + 1.96*sqrt(Tv)), 3),
    z         = round(Q / sqrt(Tv), 3),
    p_value   = round(2*pnorm(abs(Q/sqrt(Tv)), lower.tail=FALSE), 4)
  ) %>% filter(term != "Intercept")
}

# ---- Helper: HL test --------------------------------------------------------

hl_test <- function(fit, data, outcome_expr, g=10) {
  pred <- predict(fit, type="fitted")
  obs  <- as.integer(eval(parse(text=outcome_expr), data))
  grp  <- cut(pred, breaks=quantile(pred, probs=seq(0,1,1/g)),
              include.lowest=TRUE)
  tab  <- data.frame(pred=pred, obs=obs, grp=grp) %>%
    group_by(grp) %>%
    summarise(n=n(), obs_ev=sum(obs), exp_ev=sum(pred), .groups="drop")
  hl   <- sum((tab$obs_ev - tab$exp_ev)^2 / (tab$exp_ev+1e-8))
  list(stat=round(hl,3), df=g-2, p=round(pchisq(hl, df=g-2, lower.tail=FALSE), 4))
}

# =============================================================
# MODEL 1 (reference): Primary binary — Major+Both vs Minor
# Already run; reload for comparison
# =============================================================

cat("  [Primary model loaded for reference]\n")
cat(sprintf("  Primary: Apparent c=%.3f | Corrected c=%.3f | Optimism=%.3f | Brier=%.4f\n\n",
            main_val$Apparent_c[1], main_val$Corrected_c[1],
            main_val$Optimism[1],   main_val$Brier_corrected[1]))

# =============================================================
# MODEL 2: Strict Major only vs Minor+Both
# Clinical question: are "Both" patients more like Minor?
# =============================================================

cat("  --- MODEL 2: Strict Major (n=54) vs Minor+Both (n=187) ---\n")

# Add strict_major to each imputed dataset
imputed_strict <- lapply(imputed_list, function(d) {
  d$strict_major <- factor(
    ifelse(d$amp_group_nom == "Major", "Yes", "No"),
    levels=c("No","Yes")
  )
  d
})

n_strict <- sum(imputed_strict[[1]]$strict_major == "Yes")
cat(sprintf("  Strict Major events: %d / %d (%.1f%%)\n",
            n_strict, nrow(imputed_strict[[1]]),
            100*n_strict/nrow(imputed_strict[[1]])))
cat(sprintf("  EPV: %.1f\n\n", n_strict / length(avail_vars)))

form_strict <- as.formula(paste('strict_major == "Yes" ~', form_str))

lrm_strict <- lapply(imputed_strict, function(d) {
  options(contrasts=c("contr.treatment","contr.treatment"))
  assign("dd_lrm", datadist(d), envir=.GlobalEnv)
  options(datadist="dd_lrm")
  tryCatch(lrm(form_strict, data=d, x=TRUE, y=TRUE), error=function(e) NULL)
})

ok_strict <- !sapply(lrm_strict, is.null)
cat(sprintf("  Converged: %d/20\n", sum(ok_strict)))

strict_coefs <- pool_lrm_rubin(lrm_strict[ok_strict])
cat("  Pooled ORs (Strict Major vs Minor+Both):\n")
print(strict_coefs)

# Bootstrap validation on first imputed dataset
d1_strict  <- imputed_strict[[1]]
assign("dd_lrm", datadist(d1_strict), envir=.GlobalEnv)
options(datadist="dd_lrm")
lrm_strict1 <- lrm(form_strict, data=d1_strict, x=TRUE, y=TRUE)

cat("\n  Bootstrap internal validation (B=500)...\n")
val_strict   <- validate(lrm_strict1, B=500, method="boot")
Dxy_s        <- val_strict["Dxy","index.corrected"]
c_strict_app <- round(val_strict["Dxy","index.orig"]/2 + 0.5, 3)
c_strict_cor <- round(Dxy_s/2 + 0.5, 3)
opt_strict   <- round(c_strict_app - c_strict_cor, 3)
brier_strict <- round(val_strict["B","index.corrected"], 4)
r2_strict    <- round(val_strict["R2","index.corrected"], 3)

cat(sprintf("  Apparent c:   %.3f\n", c_strict_app))
cat(sprintf("  Corrected c:  %.3f\n", c_strict_cor))
cat(sprintf("  Optimism:     %.3f\n", opt_strict))
cat(sprintf("  Brier (corr): %.4f\n", brier_strict))

roc_strict  <- roc(d1_strict$strict_major,
                    predict(lrm_strict1, type="fitted"),
                    levels=c("No","Yes"), direction="<", quiet=TRUE)
auc_strict  <- round(auc(roc_strict), 3)

hl_strict   <- hl_test(lrm_strict1, d1_strict, 'strict_major == "Yes"')
cat(sprintf("  AUROC (apparent): %.3f\n", auc_strict))
cat(sprintf("  Hosmer-Lemeshow:  Chi2=%.3f, df=%d, p=%.4f\n\n",
            hl_strict$stat, hl_strict$df, hl_strict$p))

# =============================================================
# MODEL 3: Confirm primary — Major+Both vs Minor (refit cleanly)
# =============================================================

cat("  --- MODEL 3: Confirm Primary — Major+Both (n=73) vs Minor (n=168) ---\n")

form_primary <- as.formula(paste('major_binary == "Yes" ~', form_str))

lrm_primary <- lapply(imputed_list, function(d) {
  options(contrasts=c("contr.treatment","contr.treatment"))
  assign("dd_lrm", datadist(d), envir=.GlobalEnv)
  options(datadist="dd_lrm")
  tryCatch(lrm(form_primary, data=d, x=TRUE, y=TRUE), error=function(e) NULL)
})

ok_primary    <- !sapply(lrm_primary, is.null)
primary_coefs <- pool_lrm_rubin(lrm_primary[ok_primary])

d1_primary <- imputed_list[[1]]
assign("dd_lrm", datadist(d1_primary), envir=.GlobalEnv)
options(datadist="dd_lrm")
lrm_primary1 <- lrm(form_primary, data=d1_primary, x=TRUE, y=TRUE)
val_primary  <- validate(lrm_primary1, B=500, method="boot")

c_prim_app <- round(val_primary["Dxy","index.orig"]/2 + 0.5, 3)
c_prim_cor <- round(val_primary["Dxy","index.corrected"]/2 + 0.5, 3)
opt_prim   <- round(c_prim_app - c_prim_cor, 3)
brier_prim <- round(val_primary["B","index.corrected"], 4)

# =============================================================
# Coefficient comparison across outcome definitions
# =============================================================

cat("\n  --- Coefficient Comparison Across Outcome Definitions ---\n")

coef_compare <- merge(
  primary_coefs %>% select(term, OR_primary=OR, p_primary=p_value),
  strict_coefs  %>% select(term, OR_strict=OR,  p_strict=p_value),
  by="term", all=TRUE
) %>%
  mutate(
    direction_consistent = !is.na(OR_primary) & !is.na(OR_strict) &
                           sign(log(OR_primary)) == sign(log(OR_strict)),
    significance_consistent = !is.na(p_primary) & !is.na(p_strict) &
                               (p_primary < 0.05) == (p_strict < 0.05),
    OR_ratio = round(OR_strict / OR_primary, 2)
  )

cat("\n  Direction consistent: ",
    sum(coef_compare$direction_consistent, na.rm=TRUE), "/",
    sum(!is.na(coef_compare$direction_consistent)), "\n")
cat("  Significance consistent: ",
    sum(coef_compare$significance_consistent, na.rm=TRUE), "/",
    sum(!is.na(coef_compare$significance_consistent)), "\n")

print(coef_compare)

# =============================================================
# Where do the "Both" patients sit clinically?
# =============================================================

cat("\n  --- Clinical Profile of 'Both' Patients (n=19) ---\n")

df_clean <- readRDS("outputs/df_clean.rds")

both_profile <- df_clean %>%
  group_by(amp_group) %>%
  summarise(
    n              = n(),
    median_hb      = round(median(hb,         na.rm=TRUE), 2),
    median_wbc     = round(median(wbc,         na.rm=TRUE), 2),
    median_creat   = round(median(creatinine,  na.rm=TRUE), 1),
    pct_dry_gang   = round(100*mean(dry_gangrene=="Yes", na.rm=TRUE), 1),
    pct_icu        = round(100*mean(icu_admission=="Yes", na.rm=TRUE), 1),
    median_dm_adm  = median(as.integer(dm_admissions), na.rm=TRUE),
    .groups        = "drop"
  )

cat("\n  Summary by group:\n")
print(both_profile)

cat("\n  Interpretation: Are 'Both' patients clinically closer to Major or Minor?\n")
for (v in c("median_hb","median_wbc","median_creat","pct_dry_gang","pct_icu")) {
  vals <- setNames(both_profile[[v]], both_profile$amp_group)
  if (all(c("Minor","Both","Major") %in% names(vals))) {
    range_mn_mj <- abs(vals["Major"] - vals["Minor"])
    both_to_minor <- abs(vals["Both"] - vals["Minor"])
    both_to_major <- abs(vals["Both"] - vals["Major"])
    closer <- ifelse(both_to_minor <= both_to_major, "Minor", "Major")
    cat(sprintf("  %-20s — Both closer to: %s (dist_minor=%.2f, dist_major=%.2f)\n",
                v, closer, both_to_minor, both_to_major))
  }
}

# =============================================================
# Performance summary table
# =============================================================

perf_summary <- data.frame(
  Model = c(
    "Primary: Major+Both vs Minor",
    "Strict:  Major only vs Minor+Both"
  ),
  N_events    = c(73, 54),
  EPV         = round(c(73, 54) / length(avail_vars), 1),
  Apparent_c  = c(c_prim_app, c_strict_app),
  Corrected_c = c(c_prim_cor, c_strict_cor),
  Optimism    = c(opt_prim,   opt_strict),
  Brier_corr  = c(brier_prim, brier_strict),
  HL_p        = c(
    hl_test(lrm_primary1, d1_primary, 'major_binary == "Yes"')$p,
    hl_strict$p
  )
)

cat("\n  === PERFORMANCE SUMMARY ===\n")
print(perf_summary)

# =============================================================
# Save outputs
# =============================================================

write.csv(coef_compare,  "outputs/sensitivity_collapse_both.csv",       row.names=FALSE)
write.csv(perf_summary,  "outputs/sensitivity_collapse_validation.csv",  row.names=FALSE)
write.csv(both_profile,  "outputs/both_group_clinical_profile.csv",      row.names=FALSE)

saveRDS(list(
  primary_coefs = primary_coefs,
  strict_coefs  = strict_coefs,
  coef_compare  = coef_compare,
  perf_summary  = perf_summary,
  both_profile  = both_profile,
  val_primary   = val_primary,
  val_strict    = val_strict
), "outputs/sensitivity_collapse_summary.rds")

cat("\n  Saved: outputs/sensitivity_collapse_both.csv\n")
cat("  Saved: outputs/sensitivity_collapse_validation.csv\n")
cat("  Saved: outputs/both_group_clinical_profile.csv\n")
cat("  Saved: outputs/sensitivity_collapse_summary.rds\n")

cat("\n================================================================\n")
cat("  SA-4 COMPLETE\n")
cat("================================================================\n")
