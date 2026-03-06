# =============================================================================
# STEP 1: DATA CLEANING & RECODING
# =============================================================================
# Input : data/DF.xlsx
# Output: outputs/df_clean.rds     — clean data frame (no imputation yet)
#         outputs/df_complete.rds  — complete cases only (for sensitivity)
#         outputs/missingness_report.csv
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
})

# ---- 1.1  Load raw data -----------------------------------------------------

df_raw <- read_excel(
  path  = "data/DF.xlsx",
  sheet = 1,
  col_types = "text"   # read everything as text first; we recode explicitly
)

cat(sprintf("  Raw data loaded: %d rows x %d columns\n", nrow(df_raw), ncol(df_raw)))

# ---- 1.2  Standardise missing-value strings ---------------------------------
# All of these mean "missing" in the raw file

MISSING_STRINGS <- c("NA", "N/A", "NA.", "#NULL!", "no", "No data", "")

df <- df_raw %>%
  mutate(across(everything(), ~ {
    x <- trimws(as.character(.x))
    ifelse(x %in% MISSING_STRINGS | is.na(x), NA_character_, x)
  }))

# ---- 1.3  Primary outcome: amputation group ---------------------------------
# groups22: 1 = Minor only, 2 = Major only, 3 = Both
# Ordered as ordinal severity: Minor(1) < Both(3) < Major(2)
# NOTE: clinical ordering used in POM = Minor < Both < Major

df <- df %>%
  mutate(
    amp_group = case_when(
      groups22 == "1" ~ "Minor",
      groups22 == "2" ~ "Major",
      groups22 == "3" ~ "Both",
      TRUE            ~ NA_character_
    ),
    # Ordinal factor for proportional odds modelling
    amp_group_ord = factor(amp_group,
                           levels = c("Minor", "Both", "Major"),
                           ordered = TRUE),
    # Nominal factor for multinomial modelling (reference = Minor)
    amp_group_nom = factor(amp_group,
                           levels = c("Minor", "Both", "Major")),
    # Binary: major amputation yes/no (Major or Both vs Minor)
    major_binary = factor(
      ifelse(amp_group %in% c("Major", "Both"), "Yes", "No"),
      levels = c("No", "Yes")
    )
  )

# ---- 1.4  Demographics ------------------------------------------------------

df <- df %>%
  mutate(
    gender = factor(Gender, levels = c("Female", "Male")),

    age_group = factor(
      case_when(
        Age == "< 30ys"      ~ "<30",
        Age == "30-40 years" ~ "30-40",
        Age == "41-50 years" ~ "41-50",
        Age == "51-60 years" ~ "51-60",
        Age == "> 60ys"      ~ ">60",
        TRUE                 ~ NA_character_
      ),
      levels = c("<30", "30-40", "41-50", "51-60", ">60"),
      ordered = TRUE
    ),

    nationality = factor(
      ifelse(Nationality == "Saudi", "Saudi", "Non-Saudi"),
      levels = c("Saudi", "Non-Saudi")
    )
  )

# ---- 1.5  Diabetes characteristics ------------------------------------------

df <- df %>%
  mutate(
    dm_type = factor(
      case_when(
        Typeofdiabetes == "Type 1" ~ "Type1",
        Typeofdiabetes == "Type 2" ~ "Type2",
        TRUE                       ~ NA_character_
      ),
      levels = c("Type2", "Type1")   # Type2 as reference
    ),

    dm_duration = factor(
      case_when(
        Durationofdiabetes == "< 10ys"  ~ "<10",
        Durationofdiabetes == "10-20ys" ~ "10-20",
        Durationofdiabetes == "21-30ys" ~ "21-30",
        Durationofdiabetes == "31-40ys" ~ "31-40",
        Durationofdiabetes == "> 40ys"  ~ ">40",
        TRUE                            ~ NA_character_  # 51 true missing
      ),
      levels = c("<10", "10-20", "21-30", "31-40", ">40"),
      ordered = TRUE
    ),

    oha_use     = factor(Useoralhypoglycemicagents, levels = c("No", "Yes")),
    insulin_use = factor(UseInsulintreatment,       levels = c("No", "Yes"))
  )

# ---- 1.6  Comorbidities (all binary 0/1 numeric → factor) -------------------

df <- df %>%
  mutate(
    chronic_disease = factor(
      case_when(
        CHRONIC == "1" ~ "Yes",
        CHRONIC == "2" ~ "No",        # coded as 2 in raw = No
        TRUE           ~ NA_character_
      ),
      levels = c("No", "Yes")
    ),
    cvd          = factor(ifelse(CVD           == "1", "Yes", "No"), levels = c("No","Yes")),
    htn          = factor(ifelse(HTN           == "1", "Yes", "No"), levels = c("No","Yes")),
    dm_neuro     = factor(ifelse(DM_NEUROPATHY == "1", "Yes", "No"), levels = c("No","Yes")),
    dm_nephro    = factor(ifelse(DM_NEPHRO     == "1", "Yes", "No"), levels = c("No","Yes")),
    ckd          = factor(ifelse(CKD           == "1", "Yes", "No"), levels = c("No","Yes")),
    malignancy   = factor(ifelse(MALIGNNCY     == "1", "Yes", "No"), levels = c("No","Yes"))
  )

# ---- 1.7  Gangrene, osteomyelitis, surgery history --------------------------

df <- df %>%
  mutate(
    dry_gangrene   = factor(HistoryofDrygangrene,    levels = c("No", "Yes")),
    wet_gangrene   = factor(HistoryofWetgangrene,    levels = c("No", "Yes")),
    osteomyelitis  = factor(Presenceofosteomyelitis, levels = c("No", "Yes")),
    prev_surgery   = factor(Historyofprevioussurgeries, levels = c("No", "Yes"))
  )

# ---- 1.8  Amputation details ------------------------------------------------

df <- df %>%
  mutate(
    # Level of minor amputation — collapse rare categories
    minor_level = case_when(
      Ifyeswhatsthelevelofminoramputation %in%
        c("AMPUTATION OF BIG ,2ND & 3RD TOES", "big and second toe",
          "Index and middle finger") ~ "Other toes",
      Ifyeswhatsthelevelofminoramputation == "Great toe or first ray" ~ "Great toe/first ray",
      Ifyeswhatsthelevelofminoramputation == "Transmetatarsal"        ~ "Transmetatarsal",
      Ifyeswhatsthelevelofminoramputation == "Other toes"             ~ "Other toes",
      TRUE                                                             ~ NA_character_
    ),
    minor_level = factor(minor_level,
                         levels = c("Other toes", "Great toe/first ray", "Transmetatarsal")),

    # Level of major amputation
    major_level = factor(
      Ifyeswhatsthelevelofmajoramputation,
      levels = c("Below knee", "Above knee")
    ),

    # Cause of amputation — consolidate mixed categories
    cause = case_when(
      Causeofamputation == "Infection"                      ~ "Infection",
      Causeofamputation == "Critical ischemia"              ~ "Critical ischemia",
      Causeofamputation == "Trauma"                         ~ "Trauma",
      Causeofamputation %in% c("Critical ischemia & Infection",
        "Major amputation: Trauma, years ago. Minor amputation: Infection")
                                                            ~ "Mixed",
      TRUE                                                  ~ NA_character_  # 23 NAs
    ),
    cause = factor(cause, levels = c("Infection", "Critical ischemia", "Trauma", "Mixed"))
  )

# ---- 1.9  Clinical / admission variables ------------------------------------

df <- df %>%
  mutate(
    site_admission = factor(Siteofadmission, levels = c("ER", "Outpatient", "Referral")),

    hosp_stay = factor(
      Durationofhospitalstay,
      levels = c("<5 days", "5-10 days", "11-15 days", ">15 days"),
      ordered = TRUE
    ),

    icu_admission = factor(ICUadmission, levels = c("No", "Yes")),

    dm_admissions = factor(
      Numberofadmissionsrelatedtodiabetesstatus,
      levels = c("1", "2", "3", "4", ">4"),
      ordered = TRUE
    )
  )

# ---- 1.10  Continuous laboratory variables ----------------------------------
# Apply log10 transformation to right-skewed variables for modelling

df <- df %>%
  mutate(
    hba1c      = as.numeric(HbA1C),           # 9 NAs — will be imputed
    wbc        = as.numeric(WBCs103uL),
    creatinine = as.numeric(CreatinineumolL),
    hb         = as.numeric(HbgdL),

    # Log10-transformed versions (for regression models)
    log_wbc        = log10(wbc),
    log_creatinine = log10(creatinine),

    # Year as numeric
    year = as.integer(Yearofamputation)
  )

# ---- 1.11  Select and rename final clean columns ----------------------------

df_clean <- df %>%
  select(
    # Primary outcome
    amp_group, amp_group_ord, amp_group_nom, major_binary,
    # Demographics
    gender, age_group, nationality,
    # DM
    dm_type, dm_duration, oha_use, insulin_use,
    # Comorbidities
    chronic_disease, cvd, htn, dm_neuro, dm_nephro, ckd, malignancy,
    # Gangrene / clinical
    dry_gangrene, wet_gangrene, cause,
    # Amputation details
    minor_level, major_level,
    # Admission
    site_admission, hosp_stay, icu_admission, dm_admissions,
    # History
    osteomyelitis, prev_surgery,
    # Labs (raw + log)
    hba1c, wbc, creatinine, hb, log_wbc, log_creatinine,
    # Time
    year
  )

cat(sprintf("  Clean data: %d rows x %d columns\n", nrow(df_clean), ncol(df_clean)))

# ---- 1.12  Missingness report -----------------------------------------------

miss_report <- data.frame(
  variable    = names(df_clean),
  n_missing   = sapply(df_clean, function(x) sum(is.na(x))),
  pct_missing = round(sapply(df_clean, function(x) mean(is.na(x)) * 100), 1)
) %>%
  filter(n_missing > 0) %>%
  arrange(desc(n_missing))

cat("\n  Missingness summary (variables with any missing):\n")
print(miss_report)

# ---- 1.13  Complete-case subset for sensitivity analysis --------------------

# Variables used in primary regression model (excludes level details)
model_vars <- c("amp_group_ord", "amp_group_nom", "major_binary",
                "gender", "age_group", "oha_use", "insulin_use",
                "cvd", "htn", "dm_neuro", "dm_nephro", "ckd",
                "dry_gangrene", "wet_gangrene", "cause",
                "site_admission", "icu_admission", "dm_admissions",
                "osteomyelitis", "prev_surgery",
                "hba1c", "log_wbc", "log_creatinine", "hb")

df_complete <- df_clean %>%
  select(all_of(model_vars)) %>%
  drop_na()

cat(sprintf("\n  Complete cases (model variables): %d / %d (%.1f%%)\n",
            nrow(df_complete), nrow(df_clean),
            100 * nrow(df_complete) / nrow(df_clean)))

# ---- 1.14  Outcome distribution check ---------------------------------------

cat("\n  Outcome distribution (amp_group):\n")
print(table(df_clean$amp_group))
cat("\n  Major binary:\n")
print(table(df_clean$major_binary))

# ---- 1.15  Save outputs -----------------------------------------------------

dir.create("outputs", showWarnings = FALSE)

saveRDS(df_clean,    "outputs/df_clean.rds")
saveRDS(df_complete, "outputs/df_complete.rds")
write.csv(miss_report, "outputs/missingness_report.csv", row.names = FALSE)

cat("\n  Saved: outputs/df_clean.rds\n")
cat("  Saved: outputs/df_complete.rds\n")
cat("  Saved: outputs/missingness_report.csv\n")
