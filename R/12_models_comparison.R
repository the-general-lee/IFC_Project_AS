# ================================================================
# MODELS COMPARISON — LASSO benchmark vs Linear Mixed Models
# ----------------------------------------------------------------
# This script does not fit anything. It loads the CV-summary CSVs
# produced by:
#   - R/10_final_benchmark.R  -> outputs/final_benchmark/repeated_cv_summary.csv
#   - R/11_mixed_models.R     -> outputs/mixed_models/mixed_models_cv_summary.csv
# and produces a single comparison table + bar chart with error bars.
#
# Goal: answer the question "does adding an explicit hierarchy
# (region / province / GRINS macroclass) as a random effect improve
# generalisation beyond the LASSO benchmark?"
# ================================================================

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2)
})

lasso_file <- "outputs/final_benchmark/repeated_cv_summary.csv"
lmm_file   <- "outputs/mixed_models/mixed_models_cv_summary.csv"
out_dir    <- "outputs/mixed_models"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(lasso_file)) stop("Missing LASSO benchmark summary: ", lasso_file)
if (!file.exists(lmm_file))   stop("Missing LMM summary: ", lmm_file)

lasso <- read.csv(lasso_file, stringsAsFactors = FALSE)
lmm   <- read.csv(lmm_file,   stringsAsFactors = FALSE)

# Unify column names so the two tables can be stacked
lasso_std <- data.frame(
  model     = paste0("LASSO ", lasso$model),
  R2_mean   = lasso$R2_mean,    R2_sd   = lasso$R2_sd,
  RMSE_mean = lasso$RMSE_mean,  RMSE_sd = lasso$RMSE_sd,
  MAE_mean  = lasso$MAE_mean,   MAE_sd  = lasso$MAE_sd
)
lmm_std <- data.frame(
  model     = lmm$model,
  R2_mean   = lmm$R2_mean,    R2_sd   = lmm$R2_sd,
  RMSE_mean = lmm$RMSE_mean,  RMSE_sd = lmm$RMSE_sd,
  MAE_mean  = lmm$MAE_mean,   MAE_sd  = lmm$MAE_sd
)

comparison <- bind_rows(lasso_std, lmm_std) %>%
  mutate(class = ifelse(grepl("^LASSO", model), "LASSO", "Mixed model")) %>%
  arrange(desc(R2_mean))

cat("--- LASSO vs Mixed Models comparison (panel-aware CV) ---\n")
print(comparison, row.names = FALSE, digits = 4)
write.csv(comparison, file.path(out_dir, "models_comparison.csv"), row.names = FALSE)

# ----------------------------------------------------------------
# Plot: R^2 mean +/- sd, ordered, colour-coded by model class
# ----------------------------------------------------------------
p <- ggplot(comparison,
            aes(x = reorder(model, R2_mean), y = R2_mean, fill = class)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = R2_mean - R2_sd, ymax = R2_mean + R2_sd),
                width = 0.2) +
  coord_flip() +
  labs(title = "Out-of-fold R² — LASSO benchmark vs Mixed Models",
       x = NULL, y = "R² (mean ± SD across folds)",
       fill = "Model class") +
  scale_fill_manual(values = c("LASSO" = "#377eb8", "Mixed model" = "#e41a1c")) +
  theme_minimal()

ggsave(file.path(out_dir, "models_comparison_R2.png"),
       p, width = 9, height = 5, dpi = 150)

cat("\nWritten:\n  ",
    file.path(out_dir, "models_comparison.csv"), "\n  ",
    file.path(out_dir, "models_comparison_R2.png"), "\n")

# Verdict text shown to the console
best_lasso <- max(lasso_std$R2_mean)
best_lmm   <- max(lmm_std$R2_mean)
delta      <- best_lmm - best_lasso
cat(sprintf("\nBest LASSO R^2: %.4f\nBest LMM   R^2: %.4f\nDelta (LMM - LASSO): %+.4f\n",
            best_lasso, best_lmm, delta))
if (abs(delta) < 0.005) {
  cat("=> LMM and LASSO are statistically indistinguishable: the GRINS\n",
      "   indicators already absorb the hierarchical structure.\n")
} else if (delta > 0) {
  cat("=> LMM adds explanatory power: residual heterogeneity between groups\n",
      "   beyond what the fixed effects capture.\n")
} else {
  cat("=> LASSO outperforms LMM, possibly because the random-effect groups\n",
      "   are too coarse or fixed effects already over-explain the structure.\n")
}
