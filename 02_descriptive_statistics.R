# =============================================================================
# STEP 2: DESCRIPTIVE STATISTICS
# =============================================================================
# Input : outputs/df_clean.rds
# Output: outputs/descriptive_table1.csv
#         outputs/normality_tests.csv
#         outputs/posthoc_pairwise.csv
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(coin)
  library(multcomp)
})

df <- readRDS("outputs/df_clean.rds")

cat(sprintf("  Loaded clean data: %d rows\n", nrow(df)))

# ---- 2.1  Normality testing for continuous variables ------------------------

cat("\n  Normality tests (Shapiro-Wilk, full sample and by group):\n")

cont_vars <- c("hba1c", "wbc", "creatinine", "hb", "log_wbc", "log_creatinine")

normality_results <- lapply(cont_vars, function(v) {
  x <- df[[v]]
  x <- x[!is.na(x)]
  sw <- shapiro.test(x)
  data.frame(
    variable  = v,
    n         = length(x),
    mean      = round(mean(x), 3),
    sd        = round(sd(x), 3),
    median    = round(median(x), 3),
    iqr_low   = round(quantile(x, 0.25), 3),
    iqr_high  = round(quantile(x, 0.75), 3),
    sw_W      = round(sw$statistic, 4),
    sw_p      = round(sw$p.value, 4),
    normal    = ifelse(sw$p.value > 0.05, "Yes", "No")
  )
}) %>% bind_rows()

print(normality_results)

# ---- 2.2  Summary statistics by amputation group ----------------------------

cat("\n  Descriptive statistics by amputation group:\n")

# Continuous: median [IQR] for non-normal; mean (SD) for normal
summarise_continuous <- function(data, var, group_var = "amp_group") {
  data %>%
    group_by(.data[[group_var]]) %>%
    summarise(
      n         = sum(!is.na(.data[[var]])),
      mean_val  = round(mean(.data[[var]], na.rm=TRUE), 2),
      sd_val    = round(sd(.data[[var]], na.rm=TRUE), 2),
      median_val = round(median(.data[[var]], na.rm=TRUE), 2),
      q25       = round(quantile(.data[[var]], 0.25, na.rm=TRUE), 2),
      q75       = round(quantile(.data[[var]], 0.75, na.rm=TRUE), 2),
      .groups   = "drop"
    ) %>%
    mutate(variable = var)
}

cont_summary <- lapply(c("hba1c","wbc","creatinine","hb"), function(v) {
  summarise_continuous(df, v)
}) %>% bind_rows()

print(cont_summary)

# Categorical: counts and percentages
summarise_categorical <- function(data, var, group_var = "amp_group") {
  data %>%
    count(.data[[group_var]], .data[[var]], name = "n") %>%
    group_by(.data[[group_var]]) %>%
    mutate(
      total = sum(n),
      pct   = round(100 * n / total, 1)
    ) %>%
    ungroup() %>%
    rename(group = all_of(group_var), level = all_of(var)) %>%
    mutate(variable = var, level = as.character(level)) %>%
    dplyr::select(variable, level, group, n, pct) %>%
    arrange(variable, group, level)
}

cat_vars <- c("gender","age_group","nationality","dm_type","dm_duration",
              "oha_use","insulin_use","chronic_disease","cvd","htn",
              "dm_neuro","dm_nephro","ckd","malignancy",
              "dry_gangrene","wet_gangrene","cause",
              "site_admission","hosp_stay","icu_admission","dm_admissions",
              "osteomyelitis","prev_surgery")

cat_summary <- lapply(cat_vars, function(v) {
  summarise_categorical(df, v)
}) %>% bind_rows()

# ---- 2.3  Statistical tests: Kruskal-Wallis for continuous ------------------

cat("\n  Kruskal-Wallis tests (continuous variables by amputation group):\n")

kw_results <- lapply(c("hba1c","wbc","creatinine","hb"), function(v) {
  sub <- df %>% filter(!is.na(.data[[v]]), !is.na(amp_group))
  # Use base kruskal.test (coin requires numeric groups)
  kt <- kruskal.test(as.formula(paste(v, "~ amp_group")), data=sub)
  data.frame(
    variable  = v,
    chi_sq    = round(kt$statistic, 3),
    df        = kt$parameter,
    p_value   = round(kt$p.value, 4)
  )
}) %>% bind_rows()

print(kw_results)

# ---- 2.4  Post-hoc pairwise Wilcoxon with Bonferroni correction -------------

cat("\n  Post-hoc pairwise Wilcoxon tests (Bonferroni corrected):\n")

sig_vars <- kw_results %>% filter(p_value < 0.05) %>% pull(variable)
groups   <- c("Minor", "Major", "Both")
pairs    <- combn(groups, 2, simplify = FALSE)

posthoc_results <- lapply(sig_vars, function(v) {
  lapply(pairs, function(pair) {
    sub <- df %>%
      filter(amp_group %in% pair, !is.na(.data[[v]]))
    wt  <- wilcox.test(as.formula(paste(v, "~ amp_group")), data = sub,
                       exact = FALSE)
    data.frame(
      variable   = v,
      group1     = pair[1],
      group2     = pair[2],
      W          = round(wt$statistic, 1),
      p_raw      = round(wt$p.value, 4)
    )
  }) %>% bind_rows()
}) %>% bind_rows() %>%
  group_by(variable) %>%
  mutate(p_bonferroni = round(p.adjust(p_raw, method = "bonferroni"), 4)) %>%
  ungroup()

print(posthoc_results)

# ---- 2.5  Chi-squared tests for categorical variables -----------------------

cat("\n  Chi-squared tests (categorical variables by amputation group):\n")

chisq_results <- lapply(cat_vars, function(v) {
  sub <- df %>% filter(!is.na(.data[[v]]), !is.na(amp_group))
  tbl <- table(sub[[v]], sub$amp_group)

  # Use Fisher exact if any expected cell < 5
  expected_ok <- all(chisq.test(tbl, correct=FALSE)$expected >= 5)

  if (expected_ok) {
    ct <- chisq.test(tbl, correct = FALSE)
    data.frame(variable = v, test = "Chi-sq",
               statistic = round(ct$statistic, 3),
               df = ct$parameter,
               p_value = round(ct$p.value, 4))
  } else {
    ft <- tryCatch(
      fisher.test(tbl, simulate.p.value = TRUE, B = 10000),
      error = function(e) list(p.value = NA)
    )
    data.frame(variable = v, test = "Fisher",
               statistic = NA, df = NA,
               p_value = round(ft$p.value, 4))
  }
}) %>% bind_rows()

print(chisq_results)

# ---- 2.6  Overall sample descriptives (Table 1 raw data) --------------------

overall_cont <- data.frame(
  variable = c("HbA1C (%)","WBCs (10^3/uL)","Creatinine (umol/L)","Hb (g/dL)"),
  n        = sapply(c("hba1c","wbc","creatinine","hb"),
                    function(v) sum(!is.na(df[[v]]))),
  mean_sd  = sapply(c("hba1c","wbc","creatinine","hb"),
                    function(v) sprintf("%.2f ± %.2f",
                                        mean(df[[v]], na.rm=TRUE),
                                        sd(df[[v]], na.rm=TRUE))),
  median_iqr = sapply(c("hba1c","wbc","creatinine","hb"),
                      function(v) sprintf("%.2f [%.2f–%.2f]",
                                          median(df[[v]], na.rm=TRUE),
                                          quantile(df[[v]], 0.25, na.rm=TRUE),
                                          quantile(df[[v]], 0.75, na.rm=TRUE)))
)

cat("\n  Overall continuous variable summary:\n")
print(overall_cont)

# ---- 2.7  Save outputs ------------------------------------------------------

write.csv(normality_results, "outputs/normality_tests.csv",     row.names = FALSE)
write.csv(cont_summary,      "outputs/cont_summary_by_group.csv", row.names = FALSE)
write.csv(cat_summary,       "outputs/cat_summary_by_group.csv",  row.names = FALSE)
write.csv(kw_results,        "outputs/kruskal_wallis_tests.csv", row.names = FALSE)
write.csv(posthoc_results,   "outputs/posthoc_pairwise.csv",     row.names = FALSE)
write.csv(chisq_results,     "outputs/chisquared_tests.csv",     row.names = FALSE)
write.csv(overall_cont,      "outputs/overall_cont_summary.csv", row.names = FALSE)

# Save full table1 components as RDS for figure generation later
saveRDS(list(
  cont_summary  = cont_summary,
  cat_summary   = cat_summary,
  kw_results    = kw_results,
  chisq_results = chisq_results,
  posthoc       = posthoc_results
), "outputs/table1_components.rds")

cat("\n  Saved: all descriptive outputs to outputs/\n")
