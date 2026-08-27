#################################################################
##                       Soundarya S                           ##
##                            August 2026                      ##
##                            Inter-correlation and regression ##
##                           CKD-CC                           ##
#################################################################

# CKD-CC
# # Stage 1: Inter-correlation of the 12 temperature predictors -> select ~3.
# Stage 2: Log-linear within-state fixed-effects panel regression (+ year FE),
#          mirroring the CVD pipeline (pglm, Gaussian, BFGS).
# Requires: 00_build_panel.R already run (reads _rawData/derived/ckd_temp_panel.rds)
# Outcomes (counts): Deaths, Incidence, Prevalence x {all-CKD (589), CKD-U (593)}

# Author: Soundarya S
# Version: date

# Packages
pacman::p_load(tidyverse, here, pglm, corrplot, knitr)

panel <- readRDS(here("2_derivedData", "ckd_temp_panel.rds"))
dim(panel)


pred_cols <- c("mean_2_5","mean_97_5","median_2_5","median_97_5",
               "mode_2_5","mode_97_5","max_2_5","max_97_5",
               "min_2_5","min_97_5","range_2_5","range_97_5")

pred_labels <- c(
  mean_2_5    = "Mean temp (2.5th pct)",      mean_97_5   = "Mean temp (97.5th pct)",
  median_2_5  = "Median temp (2.5th pct)",    median_97_5 = "Median temp (97.5th pct)",
  mode_2_5    = "Mode temp (2.5th pct)",      mode_97_5   = "Mode temp (97.5th pct)",
  max_2_5     = "Max temp (2.5th pct)",       max_97_5    = "Max temp (97.5th pct)",
  min_2_5     = "Min temp (2.5th pct)",       min_97_5    = "Min temp (97.5th pct)",
  range_2_5   = "Temp variability (2.5th pct)", range_97_5 = "Temp variability (97.5th pct)"
)

outcomes <- c("deaths_allCKD","incidence_allCKD","prevalence_allCKD",
              "deaths_CKDu","incidence_CKDu","prevalence_CKDu")

# Pre-build log outcome columns (mirrors CVD pipeline; avoids inline-log quirks in pglm)
panel <- panel |> mutate(across(all_of(outcomes), log, .names = "log_{.col}"))

outcome_labels <- c(
  deaths_allCKD = "Deaths (all-CKD)",  incidence_allCKD = "Incidence (all-CKD)",
  prevalence_allCKD = "Prevalence (all-CKD)",
  deaths_CKDu = "Deaths (CKD-U)",      incidence_CKDu = "Incidence (CKD-U)",
  prevalence_CKDu = "Prevalence (CKD-U)"
)

fig_dir <- here("4_docs", "figures")
out_dir <- here("4_docs", "outputs")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# STAGE 1 — Inter-correlation of the 12 temperature predictors
# =============================================================================
cormat <- panel |> select(all_of(pred_cols)) |> cor(use = "complete.obs")

cat("\n===== Pearson correlation matrix (12 temperature predictors) =====\n")
print(round(cormat, 2))
write.csv(round(cormat, 3), file.path(out_dir, "aim1_correlation_matrix.csv"))

png(file.path(fig_dir, "aim1_correlation_matrix.png"), width = 1500, height = 1500, res = 200)
corrplot(cormat, method = "color", type = "upper", addCoef.col = "black",
         tl.col = "black", tl.srt = 45, number.cex = 0.6,
         title = "Inter-correlation of temperature predictors (1994-2023)",
         mar = c(0,0,2,0))
dev.off()
cat("Saved correlation figure -> figures/aim1_correlation_matrix.png\n")

# --- Automated selection: drop any predictor correlated >= cutoff with another,
#     Reproduces the |r| < 0.7 rule used in the CVD predictor screen.
cor_cutoff <- 0.70   # EDIT to test 0.8 vs 0.7

select_uncorrelated <- function(cmat, cutoff) {
  cm <- abs(cmat); diag(cm) <- 0
  repeat {
    if (max(cm) < cutoff) break
    # highest remaining pair
    idx <- which(cm == max(cm), arr.ind = TRUE)[1, ]
    i <- idx[1]; j <- idx[2]
    # drop whichever of the pair is more correlated with everything else
    drop <- if (mean(cm[i, ]) >= mean(cm[j, ])) i else j
    cm <- cm[-drop, -drop, drop = FALSE]
  }
  rownames(cm)
}

# Manual override: set to a character vector to bypass auto-selection; else NULL.
manual_metrics <- NULL

auto_metrics <- select_uncorrelated(cormat, cor_cutoff)
dropped      <- setdiff(pred_cols, auto_metrics)
selected_metrics <- if (is.null(manual_metrics)) auto_metrics else manual_metrics

cat("\n===== Predictor selection (cutoff |r| >=", cor_cutoff, ") =====\n")
cat("DROPPED (redundant):", if (length(dropped)) paste(dropped, collapse = ", ") else "none", "\n")
cat("RETAINED for final models:", paste(selected_metrics, collapse = ", "), "\n")
cat("Inter-correlation among retained (all should be <", cor_cutoff, "):\n")
print(round(cormat[selected_metrics, selected_metrics], 2))

# =============================================================================
# STAGE 2 — Fixed-effects panel regressions
# =============================================================================
# Log-linear within-state FE + year FE:  log(Y_it) = b1*Temp_it + g_t + a_i + e
# NOTE: outcomes are GBD absolute counts (Number). State FE absorbs time-invariant
#   population scale; year FE absorbs the common national trend. If the statistician
#   prefers per-100,000 rates (as in the CVD paper), re-pull the GBD Rate metric and
#   this script works unchanged on rate columns.

fit_fe <- function(outcome, predictor, data) {
  # mirrors CVD pipeline exactly: pglm Gaussian within + year FE, BFGS
  # pglm re-evaluates the formula by name internally, so a formula stored in a
  # local variable fails ("object 'form' not found"). do.call() inlines the
  # actual formula object into the call, which fixes it. Data forced to df too.
  form <- as.formula(paste0("log_", outcome, " ~ ", predictor, " + factor(year)"))
  do.call(pglm::pglm, list(
    formula = form,
    data    = as.data.frame(data),
    model   = "within",
    family  = "gaussian",
    index   = c("state_name", "year"),
    method  = "BFGS"
  ))
}

extract_pct <- function(model, predictor, outcome) {
  est <- as.numeric(coef(model)[predictor])
  se  <- as.numeric(sqrt(vcov(model)[predictor, predictor]))
  z   <- qnorm(0.975)
  tibble(
    Outcome    = outcome_labels[[outcome]],
    Predictor  = pred_labels[[predictor]],
    pred_key   = predictor,
    out_key    = outcome,
    Pct_change = (exp(est) - 1) * 100,
    CI_low     = (exp(est - z * se) - 1) * 100,
    CI_high    = (exp(est + z * se) - 1) * 100,
    P_value    = 2 * pnorm(-abs(est / se))
  )
}

run_grid <- function(metrics) {
  tidyr::crossing(outcome = outcomes, predictor = metrics) |>
    pmap_dfr(function(outcome, predictor) {
      m <- tryCatch(fit_fe(outcome, predictor, panel),
                    error = function(e) { message("FAIL ", outcome, " ~ ", predictor,
                                                   " : ", conditionMessage(e)); NULL })
      if (is.null(m)) return(tibble(Outcome = outcome_labels[[outcome]],
                                    Predictor = pred_labels[[predictor]],
                                    pred_key = predictor, out_key = outcome,
                                    Pct_change = NA, CI_low = NA, CI_high = NA,
                                    P_value = NA))
      extract_pct(m, predictor, outcome)
    })
}

fmt_tab <- function(df) {
  df |>
    mutate(
      `% change per 1°C (95% CI)` = sprintf("%.2f (%.2f, %.2f)", Pct_change, CI_low, CI_high),
      `p`   = ifelse(P_value < 0.001, "<0.001", sprintf("%.3f", P_value)),
      Sig   = ifelse(P_value < 0.05, "*", "")
    ) |>
    select(Outcome, Predictor, `% change per 1°C (95% CI)`, p, Sig)
}

# ---- FINAL models: retained predictors x 6 outcomes (each predictor separate)
cat("\n===== FINAL models (retained predictors only) =====\n")
res_final <- run_grid(selected_metrics)
write.csv(res_final, file.path(out_dir, "aim1_results_final.csv"), row.names = FALSE)

final_tab <- fmt_tab(res_final)
print(final_tab |> as.data.frame(), row.names = FALSE)
cat("\nSignificant associations (p<0.05):",
    sum(res_final$P_value < 0.05, na.rm = TRUE), "of", nrow(res_final), "\n")

cat("\nSaved: outputs/aim1_results_final.csv\n")
cat("       outputs/aim1_correlation_matrix.csv, figures/aim1_correlation_matrix.png\n")
cat("=================================================================\n")

