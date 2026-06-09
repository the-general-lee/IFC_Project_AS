# ================================================================
# LINEAR MIXED MODELS — geographic and taxonomic hierarchies
# ----------------------------------------------------------------
# Question: does explicitly modelling region/province/GRINS hierarchy
# as a random effect add information beyond what the LASSO benchmark
# (10_final_benchmark.R) already captures with its fixed effects?
#
# Three model variants are fit on the same feature matrix that the
# LASSO benchmark selected, so the comparison is "apples to apples":
#
#   M_geo        : IFC ~ X + (1 | region)
#   M_geo_nested : IFC ~ X + (1 | region/province)
#   M_tax        : IFC ~ X + (1 | GRINS_macroclass)
#
# Reported per model:
#   - R^2 marginal  (fixed effects only)
#   - R^2 conditional (fixed + random)
#   - ICC (share of residual variance that lies between groups)
#   - 5-fold panel-aware test R^2 / RMSE / MAE (same split logic as
#     the LASSO benchmark, for direct comparison)
# ================================================================

suppressPackageStartupMessages({
  library(readxl);  library(dplyr);  library(tidyr)
  library(lme4);    library(lmerTest); library(MuMIn)
})

set.seed(2026)

ifc_file       <- "data/raw/ifc/final_analysis_sorted.xlsx"
grins_file     <- "data/processed/grins_v3/comunale_v3.rds"
tassonomia_file<- "data/raw/grins/tassonomia_grins.xlsx"
benchmark_coef <- "outputs/final_benchmark/final_coefficients_standardised.csv"
out_dir        <- "outputs/mixed_models"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# OneDrive placeholder workaround
if (file.exists(ifc_file)) {
  tmp <- file.path(tempdir(), "IFC_tmp.xlsx")
  file.copy(ifc_file, tmp, overwrite = TRUE); ifc_file <- tmp
}

# ================================================================
# 1) LOAD IFC + GRINS (same as benchmark)
# ================================================================
ifc_long <- read_excel(ifc_file) %>%
  select(PRO_COM, IFC_2019, IFC_2021) %>%
  pivot_longer(c(IFC_2019, IFC_2021), names_to = "year", values_to = "IFC") %>%
  mutate(year = as.integer(sub("IFC_", "", year)),
         PRO_COM = as.numeric(PRO_COM)) %>%
  filter(!is.na(IFC), !is.na(PRO_COM))

grins <- readRDS(grins_file) %>%
  filter(anno %in% c(2019, 2021)) %>%
  mutate(codice_comune = as.numeric(codice_comune))

# Extract grouping variables for random effects
grins_groups <- grins %>%
  select(codice_comune, anno, nome_regione, nome_provincia)

# Load GRINS taxonomy for macroclass random effect
tassonomia <- read_excel(tassonomia_file, sheet = 1) %>%
  rename(PRO_COM = 1, GRINS_macroclass = 2, GRINS_class = 3) %>%
  mutate(PRO_COM = as.numeric(PRO_COM)) %>%
  filter(!is.na(PRO_COM)) %>%
  select(PRO_COM, GRINS_macroclass)

# Numeric predictors only
num_cols <- setdiff(names(grins)[sapply(grins, is.numeric)],
                    c("codice_comune", "anno"))

data <- ifc_long %>%
  inner_join(grins %>% select(all_of(c("codice_comune", "anno", num_cols))),
             by = c("PRO_COM" = "codice_comune", "year" = "anno")) %>%
  left_join(grins_groups,
            by = c("PRO_COM" = "codice_comune", "year" = "anno")) %>%
  left_join(tassonomia, by = "PRO_COM") %>%
  filter(!is.na(nome_regione), !is.na(nome_provincia), !is.na(GRINS_macroclass))

cat("Rows after merge + group filter:", nrow(data), "\n")
cat("Regions:",    length(unique(data$nome_regione)),
    " | Provinces:", length(unique(data$nome_provincia)),
    " | GRINS macroclasses:", length(unique(data$GRINS_macroclass)), "\n")

pop <- data$Popolazione
pop[is.na(pop) | pop <= 0] <- median(pop[pop > 0], na.rm = TRUE)
y <- data$IFC

# ================================================================
# 2) FEATURE ENGINEERING (identical to benchmark)
# ================================================================
classify <- function(nm) {
  ln <- tolower(nm)
  if (grepl("_anno_x$|_anno_y$|^cod_|^den_|sigla|nome_|^stringa|backcast|recovery|tipo_na", ln)) return("skip")
  if (ln %in% c("popolazione", "superficie_totale_kmq_formattato", "superfici kmq")) return("size")
  rp <- c("indice","incidenza","mobilità","mobilita","percentuale","pro capite","procapite","pro-capite",
          "per_addetto","per_dipendente","valori_percentuali","_media_","_media$","^eta_","^anzianita_",
          "media_donne","media_uomini","media_media","coverage","contribuenti_su_pop","reddito_medio",
          "reddito_pc","lacc_mean","employee_services","employee_concentration","degree_stem",
          "degree_concentration","^no2_","^pm10_","^pm25_","^o3_","_media_valori_annuali",
          "distanza_","^elezion","regionali","verde urbano","suolo consumato","acqua potabile",
          "produzione pro","sau","ricettività","ricettivita","densità")
  if (any(sapply(rp, function(p) grepl(p, ln)))) return("rate")
  if (grepl("saldo", ln)) return("count_neg")
  "count"
}

build_X <- function(df_num, pop) {
  cls <- sapply(names(df_num), classify)
  keep <- names(df_num)[cls != "skip"]; X <- df_num[, keep]; cls <- cls[keep]
  na_prop <- sapply(X, function(z) mean(is.na(z)))
  X <- X[, names(na_prop[na_prop <= 0.20])]; cls <- cls[names(X)]
  for (j in seq_along(X)) if (any(is.na(X[[j]]))) X[[j]][is.na(X[[j]])] <- median(X[[j]], na.rm = TRUE)
  for (nm in names(X)) {
    cc <- cls[nm]
    if (cc == "count") {
      X[[nm]] <- X[[nm]] / pop; X[[nm]][!is.finite(X[[nm]])] <- 0
      X[[nm]] <- log1p(pmax(X[[nm]], 0))
    } else if (cc == "count_neg") {
      X[[nm]] <- X[[nm]] / pop; X[[nm]][!is.finite(X[[nm]])] <- 0
    } else if (cc == "size") {
      X[[nm]] <- log1p(pmax(X[[nm]], 0))
    }
  }
  wins <- function(z) { q <- quantile(z, c(0.01, 0.99), na.rm = TRUE); pmin(pmax(z, q[1]), q[2]) }
  X[] <- lapply(X, wins)
  X[, sapply(X, var) > 0]
}

X_fe <- build_X(data %>% select(all_of(num_cols)), pop)
X_fe$year2021 <- as.integer(data$year == 2021)
names(X_fe) <- make.names(names(X_fe), unique = TRUE)
cat("Feature matrix:", dim(X_fe), "\n")

# ================================================================
# 3) RESTRICT TO THE LASSO-SELECTED PREDICTORS
# ----------------------------------------------------------------
# Fitting lmer on the full 326-column matrix would be slow and would
# also conflate "is the hierarchy useful?" with "is variable selection
# useful?". By restricting the fixed-effects design to the predictors
# that LASSO already kept, we isolate the contribution of the random
# effects. This is the standard "post-selection" use of LMM.
# ================================================================
sel_df <- read.csv(benchmark_coef, stringsAsFactors = FALSE)
selected <- sel_df$variable
selected <- intersect(selected, names(X_fe))
cat("Predictors carried over from LASSO benchmark:", length(selected), "\n")

X_sel <- X_fe[, selected, drop = FALSE]
# Standardise predictors so that lmer coefficients are comparable
X_sel <- scale(as.matrix(X_sel))

model_data <- data.frame(
  IFC              = y,
  region           = factor(data$nome_regione),
  province         = factor(data$nome_provincia),
  GRINS_macroclass = factor(data$GRINS_macroclass),
  PRO_COM          = data$PRO_COM,
  year             = data$year
)
model_data <- cbind(model_data, X_sel)

# Build the fixed-effects part of the formula
X_terms <- paste(selected, collapse = " + ")

# ================================================================
# 4) FIT THE FOUR MIXED MODELS ON FULL DATA
# ----------------------------------------------------------------
# Helper that captures any convergence warnings emitted by lme4
# so we can flag silently misbehaving fits at the end of the script.
# ================================================================
fit_lmm_safe <- function(formula_str, data) {
  warns <- character()
  fit <- withCallingHandlers(
    lmer(as.formula(formula_str), data = data, REML = FALSE,
         control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-4))),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  attr(fit, "warnings") <- warns
  attr(fit, "singular") <- isSingular(fit, tol = 1e-4)
  fit
}

cat("\n[Fitting M_geo]\n")
M_geo <- fit_lmm_safe(paste("IFC ~", X_terms, "+ (1 | region)"), model_data)

cat("[Fitting M_geo_nested]\n")
M_geo_nested <- fit_lmm_safe(paste("IFC ~", X_terms, "+ (1 | region/province)"), model_data)

cat("[Fitting M_tax]\n")
M_tax <- fit_lmm_safe(paste("IFC ~", X_terms, "+ (1 | GRINS_macroclass)"), model_data)

# Improvement #2: combined geography + taxonomy
# Tests whether GRINS_macroclass still contributes anything once region
# is already in the model. If the taxonomy is redundant with the
# geography (as we expect from the M_tax-alone result), the macroclass
# variance component should collapse to ~0 here too.
cat("[Fitting M_geo_plus_tax]\n")
M_geo_plus_tax <- fit_lmm_safe(
  paste("IFC ~", X_terms, "+ (1 | region) + (1 | GRINS_macroclass)"),
  model_data
)

# Improvement #4: convergence sanity check
convergence_report <- data.frame(
  model = c("M_geo", "M_geo_nested", "M_tax", "M_geo_plus_tax"),
  singular = sapply(list(M_geo, M_geo_nested, M_tax, M_geo_plus_tax),
                    function(f) attr(f, "singular")),
  n_warnings = sapply(list(M_geo, M_geo_nested, M_tax, M_geo_plus_tax),
                      function(f) length(attr(f, "warnings")))
)
cat("\n--- Convergence check ---\n")
print(convergence_report, row.names = FALSE)
write.csv(convergence_report,
          file.path(out_dir, "convergence_report.csv"), row.names = FALSE)

# ================================================================
# 5) IN-SAMPLE METRICS: R^2 marginal / conditional + ICC + variances
# ================================================================
extract_metrics <- function(m, name) {
  r2 <- tryCatch(r.squaredGLMM(m), error = function(e) c(R2m = NA, R2c = NA))
  vc <- as.data.frame(VarCorr(m))[, c("grp", "vcov")]
  total_var <- sum(vc$vcov)
  random_var <- sum(vc$vcov[vc$grp != "Residual"])
  data.frame(
    model        = name,
    R2_marginal  = r2[1],
    R2_conditional = r2[2],
    ICC          = random_var / total_var,
    sigma2_resid = vc$vcov[vc$grp == "Residual"],
    sigma2_random = random_var,
    AIC          = AIC(m),
    BIC          = BIC(m)
  )
}

insample <- bind_rows(
  extract_metrics(M_geo,          "M_geo (1|region)"),
  extract_metrics(M_geo_nested,   "M_geo_nested (1|region/province)"),
  extract_metrics(M_tax,          "M_tax (1|GRINS_macroclass)"),
  extract_metrics(M_geo_plus_tax, "M_geo_plus_tax (1|region)+(1|macroclass)")
)
print(insample, row.names = FALSE, digits = 4)
write.csv(insample, file.path(out_dir, "in_sample_summary.csv"), row.names = FALSE)

# ================================================================
# 6) PANEL-AWARE 5-FOLD CV (single repeat, for runtime reasons)
# ----------------------------------------------------------------
# Each municipality stays in a single fold (same split logic as
# the LASSO benchmark in script 10).
# ================================================================
municipalities <- unique(model_data$PRO_COM)
shuffled <- sample(municipalities)
fold_of  <- (seq_along(shuffled) - 1) %% 5 + 1

cv_records <- data.frame()
for (k in 1:5) {
  test_munis <- shuffled[fold_of == k]
  test_idx   <- which(model_data$PRO_COM %in% test_munis)
  train_idx  <- setdiff(seq_len(nrow(model_data)), test_idx)

  for (cfg in list(
    list(name = "M_geo",          formula = paste("IFC ~", X_terms, "+ (1 | region)")),
    list(name = "M_geo_nested",   formula = paste("IFC ~", X_terms, "+ (1 | region/province)")),
    list(name = "M_tax",          formula = paste("IFC ~", X_terms, "+ (1 | GRINS_macroclass)")),
    list(name = "M_geo_plus_tax", formula = paste("IFC ~", X_terms, "+ (1 | region) + (1 | GRINS_macroclass)"))
  )) {
    fit <- tryCatch(
      lmer(as.formula(cfg$formula), data = model_data[train_idx, ], REML = FALSE,
           control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-4))),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    # allow.new.levels handles regions/provinces present only in test
    pred <- tryCatch(
      predict(fit, newdata = model_data[test_idx, ], allow.new.levels = TRUE),
      error = function(e) rep(NA_real_, length(test_idx))
    )
    yt <- model_data$IFC[test_idx]
    rmse <- sqrt(mean((yt - pred)^2, na.rm = TRUE))
    mae  <- mean(abs(yt - pred), na.rm = TRUE)
    r2   <- 1 - sum((yt - pred)^2, na.rm = TRUE) / sum((yt - mean(yt))^2)
    cv_records <- rbind(cv_records,
      data.frame(model = cfg$name, fold = k, RMSE = rmse, MAE = mae, R2 = r2))
    cat(sprintf("  fold %d / %-15s done (R2=%.3f)\n", k, cfg$name, r2))
  }
}

cv_summary <- cv_records %>%
  group_by(model) %>%
  summarise(
    R2_mean   = mean(R2),   R2_sd   = sd(R2),
    RMSE_mean = mean(RMSE), RMSE_sd = sd(RMSE),
    MAE_mean  = mean(MAE),  MAE_sd  = sd(MAE),
    n_folds   = n(), .groups = "drop"
  )

cat("\n--- CV summary (5-fold panel-aware, 1 repeat) ---\n")
print(as.data.frame(cv_summary), row.names = FALSE, digits = 4)
write.csv(cv_summary,
          file.path(out_dir, "mixed_models_cv_summary.csv"), row.names = FALSE)
write.csv(cv_records,
          file.path(out_dir, "mixed_models_cv_per_fold.csv"), row.names = FALSE)

# ================================================================
# 7) RANDOM-EFFECT ESTIMATES (BLUPs) for the best model
# ----------------------------------------------------------------
# These are interpretable: how much each region/province/macroclass
# differs from the overall mean once X has been controlled for.
# ================================================================
saveRDS(M_geo,          file.path(out_dir, "fit_M_geo.rds"))
saveRDS(M_geo_nested,   file.path(out_dir, "fit_M_geo_nested.rds"))
saveRDS(M_tax,          file.path(out_dir, "fit_M_tax.rds"))
saveRDS(M_geo_plus_tax, file.path(out_dir, "fit_M_geo_plus_tax.rds"))

# Regional offsets from M_geo
ranef_geo <- as.data.frame(ranef(M_geo)$region)
ranef_geo$region <- rownames(ranef_geo)
names(ranef_geo)[1] <- "intercept_offset"
ranef_geo <- ranef_geo[order(-abs(ranef_geo$intercept_offset)), ]
write.csv(ranef_geo, file.path(out_dir, "ranef_region.csv"), row.names = FALSE)
cat("\nTop 10 regional offsets (M_geo) — higher = more fragile:\n")
print(head(ranef_geo, 10), row.names = FALSE, digits = 3)

# Improvement #3: province offsets from M_geo_nested (the best model)
# lme4 stores nested random effects under the name "province:region".
# Each row is "<province>:<region>", with the offset relative to the
# regional baseline (i.e. additional fragility on top of the region).
ranef_prov_raw <- ranef(M_geo_nested)
prov_key <- grep("province", names(ranef_prov_raw), value = TRUE)
ranef_prov <- as.data.frame(ranef_prov_raw[[prov_key]])
ranef_prov$province_region <- rownames(ranef_prov)
names(ranef_prov)[1] <- "intercept_offset"

# Polish B: fix UTF-8 mojibake in region/province names (Vallée d'Aoste etc.)
fix_encoding <- function(x) {
  # The source RDS stores some names with double-encoded UTF-8
  # (UTF-8 bytes mis-interpreted as Latin1 then re-encoded). Reversing
  # the Latin1->UTF-8 step recovers the original characters.
  out <- iconv(x, from = "UTF-8", to = "latin1")
  out <- iconv(out, from = "UTF-8", to = "UTF-8", sub = "?")
  Encoding(out) <- "UTF-8"
  out
}
ranef_prov$province_region <- fix_encoding(ranef_prov$province_region)

# Polish A: present BOTH "largest absolute offsets" AND "signed-sorted"
#  - sort by |offset| -> most extreme provinces (in either direction)
#  - sort by raw offset, top of head() and tail() -> most fragile / least fragile
ranef_prov_by_abs <- ranef_prov[order(-abs(ranef_prov$intercept_offset)), ]
ranef_prov_by_sign <- ranef_prov[order(-ranef_prov$intercept_offset), ]
write.csv(ranef_prov_by_abs, file.path(out_dir, "ranef_province.csv"), row.names = FALSE)

cat("\nTop 10 province offsets — most extreme |offset| (M_geo_nested):\n")
print(head(ranef_prov_by_abs, 10), row.names = FALSE, digits = 3)

cat("\nTop 10 MOST FRAGILE provinces (signed, beyond regional baseline):\n")
print(head(ranef_prov_by_sign, 10), row.names = FALSE, digits = 3)

cat("\nTop 10 LEAST FRAGILE provinces (signed, beyond regional baseline):\n")
print(head(ranef_prov_by_sign[order(ranef_prov_by_sign$intercept_offset), ], 10),
      row.names = FALSE, digits = 3)

# Same encoding fix on the regional table
ranef_geo$region <- fix_encoding(ranef_geo$region)
write.csv(ranef_geo, file.path(out_dir, "ranef_region.csv"), row.names = FALSE)

cat("\nOutputs in:", out_dir, "\n")
