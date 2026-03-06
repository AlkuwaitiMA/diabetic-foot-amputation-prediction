# =============================================================================
# STEP 4: TEMPORAL TREND ANALYSIS — JOINPOINT EQUIVALENT
# =============================================================================
# Input : outputs/df_clean.rds
# Output: outputs/temporal_apc_results.csv   — APC/AAPC per stratum
#         outputs/temporal_counts.csv         — annual counts
#         outputs/joinpoint_models.rds         — fitted segmented models
# =============================================================================
# Method: segmented package implements piecewise linear regression on log-
#         transformed annual counts, equivalent to NCI Joinpoint software.
#         APC = (exp(slope) - 1) × 100 per year.
#         Max 1 joinpoint for 5 data points (per NCI guidelines).
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(segmented)
})

df <- readRDS("outputs/df_clean.rds")

# ---- 4.1  Compute annual counts by stratum ----------------------------------

# Stratum 1: All amputations
# Stratum 2: Minor only
# Stratum 3: Major only (Major + Both)
# Stratum 4: Male patients
# Stratum 5: Female patients

annual_counts <- df %>%
  group_by(year) %>%
  summarise(
    all_amp   = n(),
    minor_amp = sum(amp_group == "Minor", na.rm=TRUE),
    major_amp = sum(amp_group %in% c("Major","Both"), na.rm=TRUE),
    male_amp  = sum(gender == "Male", na.rm=TRUE),
    female_amp= sum(gender == "Female", na.rm=TRUE),
    .groups   = "drop"
  ) %>%
  arrange(year)

# NOTE: 2024 is a partial year (Jan–Sep = 9 months)
# Annualise 2024 counts to 12-month equivalent: multiply by (12/9)
annual_counts <- annual_counts %>%
  mutate(across(-year, ~ ifelse(year == 2024, round(.x * (12/9)), .x)))

cat("  Annual counts (2024 annualised to 12 months):\n")
print(annual_counts)

# ---- 4.2  Joinpoint function -------------------------------------------------
# Fits log-linear model, then applies segmented() for breakpoint detection
# Returns: APC per segment, AAPC, breakpoint location (if any)

fit_joinpoint <- function(years, counts, stratum_name, max_breaks = 1) {

  # Guard: need at least 5 points; skip if all zeros
  if (length(years) < 4 || all(counts == 0)) {
    return(NULL)
  }

  # Log-linear base model
  df_t  <- data.frame(y = log(counts), x = years)
  lm0   <- lm(y ~ x, data = df_t)

  # Compute overall AAPC from simple linear model
  slope_lin <- coef(lm0)[["x"]]
  aapc      <- round((exp(slope_lin) - 1) * 100, 2)
  aapc_se   <- round(sqrt(vcov(lm0)[2,2]) * exp(slope_lin) * 100, 2)
  aapc_lo   <- round((exp(slope_lin - 1.96*sqrt(vcov(lm0)[2,2])) - 1)*100, 2)
  aapc_hi   <- round((exp(slope_lin + 1.96*sqrt(vcov(lm0)[2,2])) - 1)*100, 2)
  aapc_p    <- round(summary(lm0)$coefficients["x","Pr(>|t|)"], 4)

  result_aapc <- data.frame(
    stratum      = stratum_name,
    type         = "AAPC",
    segment      = "Overall",
    start_year   = min(years),
    end_year     = max(years),
    apc          = aapc,
    apc_95lo     = aapc_lo,
    apc_95hi     = aapc_hi,
    p_value      = aapc_p,
    joinpoint_yr = NA_real_
  )

  # Try fitting 1 joinpoint
  seg_model <- tryCatch({
    # Initial psi guess = midpoint
    psi_init <- median(years)
    segmented(lm0, seg.Z = ~x, psi = psi_init,
              control = seg.control(it.max = 100, display = FALSE))
  }, error = function(e) NULL)

  if (is.null(seg_model)) {
    cat(sprintf("    [%s] No significant joinpoint found; using linear trend\n",
                stratum_name))
    return(result_aapc)
  }

  # Test if breakpoint is significant (Davies test)
  davies_p <- tryCatch({
    dt <- davies.test(lm0, seg.Z = ~x, k = 4)
    dt$p.value
  }, error = function(e) NA)

  cat(sprintf("    [%s] Davies test p = %.4f\n", stratum_name,
              ifelse(is.na(davies_p), -1, davies_p)))

  if (!is.na(davies_p) && davies_p >= 0.05) {
    cat(sprintf("    [%s] Joinpoint not significant (p>=0.05); using AAPC only\n",
                stratum_name))
    return(result_aapc)
  }

  # Extract breakpoint
  bp <- seg_model$psi[,"Est."]

  # Extract slopes per segment
  slp <- slope(seg_model)$x
  n_seg <- nrow(slp)

  # Build segment breakpoints
  seg_breaks <- c(min(years), round(bp), max(years))

  result_segs <- lapply(1:n_seg, function(i) {
    apc_i    <- round((exp(slp[i,"Est."]) - 1)*100, 2)
    apc_lo_i <- round((exp(slp[i,"Est."] - 1.96*slp[i,"Std. Error"]) - 1)*100, 2)
    apc_hi_i <- round((exp(slp[i,"Est."] + 1.96*slp[i,"Std. Error"]) - 1)*100, 2)

    data.frame(
      stratum      = stratum_name,
      type         = "APC",
      segment      = paste0("Seg", i),
      start_year   = seg_breaks[i],
      end_year     = seg_breaks[i+1],
      apc          = apc_i,
      apc_95lo     = apc_lo_i,
      apc_95hi     = apc_hi_i,
      p_value      = NA_real_,
      joinpoint_yr = ifelse(i == 1, round(bp, 1), NA_real_)
    )
  }) %>% bind_rows()

  bind_rows(result_aapc, result_segs)
}

# ---- 4.3  Run Joinpoint for each stratum ------------------------------------

years  <- annual_counts$year
strata <- list(
  "All amputations"  = annual_counts$all_amp,
  "Minor amputation" = annual_counts$minor_amp,
  "Major amputation" = annual_counts$major_amp,
  "Male patients"    = annual_counts$male_amp,
  "Female patients"  = annual_counts$female_amp
)

cat("\n  Fitting joinpoint models:\n")

apc_results <- lapply(names(strata), function(nm) {
  fit_joinpoint(years, strata[[nm]], nm)
}) %>% bind_rows()

cat("\n  APC/AAPC Results:\n")
print(apc_results)

# ---- 4.4  Summary interpretation --------------------------------------------

cat("\n  Interpretation summary:\n")
aapc_summary <- apc_results %>%
  filter(type == "AAPC") %>%
  mutate(
    direction   = ifelse(apc > 0, "Increasing", "Decreasing"),
    significant = ifelse(p_value < 0.05, "Yes (p<0.05)", "No"),
    summary     = sprintf("%s: AAPC = %+.1f%% [%.1f, %.1f], p = %.3f (%s)",
                          stratum, apc, apc_95lo, apc_95hi,
                          p_value, significant)
  )

for (s in aapc_summary$summary) cat(" ", s, "\n")

# ---- 4.5  Save outputs ------------------------------------------------------

write.csv(annual_counts, "outputs/temporal_counts.csv",     row.names = FALSE)
write.csv(apc_results,   "outputs/temporal_apc_results.csv", row.names = FALSE)
saveRDS(apc_results,     "outputs/joinpoint_models.rds")

cat("\n  Saved: outputs/temporal_counts.csv\n")
cat("  Saved: outputs/temporal_apc_results.csv\n")
cat("  Saved: outputs/joinpoint_models.rds\n")
