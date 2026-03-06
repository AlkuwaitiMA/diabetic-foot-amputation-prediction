# Predictors of Amputation Severity in Diabetic Foot Ulcers

Retrospective cohort study, Sheikh Fahad Alawidah Diabetic Foot Center,
Qassim, Saudi Arabia, 2020–2024. n = 241 patients.

## Study Overview
- Outcome: Three-category amputation severity (Minor / Both / Major)
- Primary model: Binary logistic regression (Major + Both vs Minor)
- Internal validation: Bootstrap (B = 500), corrected AUROC = 0.870
- Methods: MICE imputation, LASSO selection, POM, ML comparison, DCA

## Analysis Pipeline
| Script | Description |
|--------|-------------|
| 00_master_run.R | Runs all steps in sequence |
| 01_data_cleaning.R | Data import, recoding, variable creation |
| 02_descriptive_statistics.R | Table 1, KW tests, post-hoc comparisons |
| 03_mice_imputation.R | Multiple imputation (m = 20) |
| 04_temporal_trend.R | Joinpoint regression, AAPC |
| 05_lasso_selection.R | Cross-validated multinomial LASSO |
| 06_multivariable_regression.R | POM and MLR |
| 07_binary_logistic.R | Prediction model, bootstrap validation |
| 08_machine_learning.R | RF, GBM, XGBoost, DT, L2 logistic |
| 09_decision_curve.R | Net benefit analysis |
| 10_sensitivity_analyses.R | SA-1, SA-2, SA-3 |
| 11_sensitivity_collapse_both.R | SA-4: outcome redefinition |

## Requirements
R version 4.3.3. Key packages: rms, mice, glmnet, ordinal,
randomForest, xgboost, gbm, pROC, dcurves, nnet, ggplot2.

## Data Availability
Raw patient data are not publicly available due to ethical restrictions.
Available from the corresponding author on reasonable request.

## Citation
[To be added upon publication]
