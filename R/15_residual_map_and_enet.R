# ================================================================
# RESIDUAL MAP + ELASTIC NET CHECK on the parsimonious model
# ----------------------------------------------------------------
# Two diagnostics on top of the 25-covariate parsimonious model
# produced by R/14_parsimonious_model.R:
#
# 1) Residual map:
#    Where in Italy does the parsimonious model under- or over-
#    predict IFC? The municipality-level residual is averaged across
#    2019/2021 and plotted on the ISTAT 2021 shapefile with a
#    diverging colour scale. This visualises the geographic pattern
#    of unexplained fragility -- and is the natural visual companion
#    to the LISA maps produced by another team member.
#
# 2) Elastic Net cross-check:
#    The professor asked whether using a "LASSO with a ridge penalty"
#    (i.e. Elastic Net, alpha < 1) would change the selection. With
#    the VIF iterative pruning we have already removed multicollinear
#    predictors, so we expect Elastic Net to converge to a similar
#    set. We test alpha values in {0.1, 0.25, 0.5, 0.75, 1.0} and
#    report selection size and CV error for each.
# ================================================================

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr)
  library(ggplot2); library(sf); library(glmnet)
})

set.seed(2026)

ifc_file       <- "data/raw/ifc/final_analysis_sorted.xlsx"
grins_file     <- "data/processed/grins_v3/comunale_v3.rds"
benchmark_coef <- "outputs/final_benchmark/final_coefficients_standardised.csv"
shape_file     <- "data/raw/istat/shapefile/Com2021.shp"
ranking_file   <- "outputs/parsimonious_model/vif_clean_ranking.csv"
out_dir        <- "outputs/parsimonious_model"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (file.exists(ifc_file)) {
  tmp <- file.path(tempdir(), "IFC_tmp.xlsx")
  file.copy(ifc_file, tmp, overwrite = TRUE); ifc_file <- tmp
}

# ================================================================
# 1) Reproduce the parsimonious data set
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

num_cols <- setdiff(names(grins)[sapply(grins, is.numeric)],
                    c("codice_comune", "anno"))

data <- ifc_long %>%
  inner_join(grins %>% select(all_of(c("codice_comune", "anno", num_cols))),
             by = c("PRO_COM" = "codice_comune", "year" = "anno"))

pop <- data$Popolazione
pop[is.na(pop) | pop <= 0] <- median(pop[pop > 0], na.rm = TRUE)
y <- data$IFC

classify <- function(nm) {
  ln <- tolower(nm)
  if (grepl("_anno_x$|_anno_y$|^cod_|^den_|sigla|nome_|^stringa|backcast|recovery|tipo_na", ln)) return("skip")
  if (ln %in% c("popolazione","superficie_totale_kmq_formattato","superfici kmq")) return("size")
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
    if (cc == "count")     { X[[nm]] <- X[[nm]]/pop; X[[nm]][!is.finite(X[[nm]])] <- 0; X[[nm]] <- log1p(pmax(X[[nm]],0)) }
    else if (cc == "count_neg") { X[[nm]] <- X[[nm]]/pop; X[[nm]][!is.finite(X[[nm]])] <- 0 }
    else if (cc == "size") { X[[nm]] <- log1p(pmax(X[[nm]],0)) }
  }
  wins <- function(z){q <- quantile(z,c(0.01,0.99),na.rm=TRUE); pmin(pmax(z,q[1]),q[2])}
  X[] <- lapply(X, wins)
  X[, sapply(X, var) > 0]
}

X_fe <- build_X(data %>% select(all_of(num_cols)), pop)
X_fe$year2021 <- as.integer(data$year == 2021)
names(X_fe) <- make.names(names(X_fe), unique = TRUE)

# Take the 25 top-ranked parsimonious predictors (from the VIF-clean ranking)
ranking <- read.csv(ranking_file, stringsAsFactors = FALSE)
parsimonious_vars <- ranking$variable[1:25]
parsimonious_vars <- intersect(parsimonious_vars, names(X_fe))
cat("Parsimonious predictors used:", length(parsimonious_vars), "\n")

X_parsimonious <- as.data.frame(scale(as.matrix(X_fe[, parsimonious_vars])))
model_df <- data.frame(IFC = y, PRO_COM = data$PRO_COM, X_parsimonious)

# ================================================================
# 2) Fit final parsimonious model and compute residuals
# ================================================================
fit_pars <- lm(as.formula(paste("IFC ~", paste(parsimonious_vars, collapse = " + "))),
               data = model_df)
res <- residuals(fit_pars)
cat("Residual SD:", round(sd(res), 3),
    " | mean:", round(mean(res), 4), "\n")

# Aggregate to municipality level (average across 2019/2021)
muni_res <- model_df %>%
  mutate(residual = res) %>%
  group_by(PRO_COM) %>%
  summarise(mean_residual = mean(residual), .groups = "drop")

# ================================================================
# 3) RESIDUAL MAP
# ================================================================
cat("\nBuilding residual map...\n")
shapes <- st_read(shape_file, quiet = TRUE)
shapes$PRO_COM_num <- as.numeric(as.character(shapes$PRO_COM))
shapes_res <- shapes %>%
  inner_join(muni_res, by = c("PRO_COM_num" = "PRO_COM"))
cat("Matched municipalities for map:", nrow(shapes_res), "/", nrow(shapes), "\n")

# Cap residuals at the 1st / 99th percentile for the colour scale
q <- quantile(shapes_res$mean_residual, c(0.01, 0.99), na.rm = TRUE)
shapes_res$residual_capped <- pmin(pmax(shapes_res$mean_residual, q[1]), q[2])
lim <- max(abs(q))

p_map <- ggplot(shapes_res) +
  geom_sf(aes(fill = residual_capped), colour = NA) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0, limits = c(-lim, lim),
                       name = "Residual\n(observed -\npredicted IFC)") +
  labs(title = "Residuals of the parsimonious 25-covariate model",
       subtitle = "Red = under-predicted (model too optimistic) | Blue = over-predicted") +
  theme_void() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 10))

ggsave(file.path(out_dir, "residual_map.png"), p_map,
       width = 8, height = 10, dpi = 200)
cat("Residual map saved.\n")

# Save the data underlying the map
write.csv(muni_res, file.path(out_dir, "municipal_residuals.csv"), row.names = FALSE)

# ================================================================
# 4) ELASTIC NET CHECK across alpha values
# ----------------------------------------------------------------
# We use the 157 VIF-pruned predictors as the input matrix. For
# each alpha in the grid we run cv.glmnet (10-fold) and report:
#   - chosen lambda.1se
#   - number of selected predictors at lambda.1se
#   - CV MSE at lambda.1se
# A smaller selection at the same MSE would mean Elastic Net is
# more efficient; a similar selection means LASSO is fine.
# ================================================================
cat("\n--- Elastic Net comparison across alpha ---\n")
ranking_full <- read.csv(ranking_file, stringsAsFactors = FALSE)
vif_clean_vars <- intersect(ranking_full$variable, names(X_fe))
X_vif_clean <- as.matrix(scale(X_fe[, vif_clean_vars]))

lambda_grid <- 10^seq(2, -3, length.out = 60)
enet_summary <- data.frame()
for (a in c(0.10, 0.25, 0.50, 0.75, 1.00)) {
  cv_fit <- cv.glmnet(X_vif_clean, y, alpha = a, nfolds = 10, lambda = lambda_grid)
  lam <- cv_fit$lambda.1se
  fit  <- glmnet(X_vif_clean, y, alpha = a, lambda = lambda_grid)
  cf   <- predict(fit, s = lam, type = "coefficients")
  n_sel <- sum(as.vector(cf) != 0) - 1   # subtract intercept
  cv_err <- cv_fit$cvm[cv_fit$lambda == lam]
  enet_summary <- rbind(enet_summary,
    data.frame(alpha = a, n_selected = n_sel,
               lambda_1se = lam, CV_MSE = cv_err))
  cat(sprintf("alpha = %.2f | n_selected = %3d | lambda.1se = %.4f | CV MSE = %.4f\n",
              a, n_sel, lam, cv_err))
}
write.csv(enet_summary,
          file.path(out_dir, "elastic_net_summary.csv"), row.names = FALSE)

# ================================================================
# 5) Save a one-page summary plot: trade-off curve + alpha comparison
# ================================================================
p_alpha <- ggplot(enet_summary, aes(x = factor(alpha))) +
  geom_col(aes(y = n_selected), fill = "steelblue") +
  geom_text(aes(y = n_selected, label = n_selected), vjust = -0.4) +
  labs(title = "Elastic Net selection size at lambda.1se",
       subtitle = "alpha = 1 is pure LASSO, alpha = 0 is pure Ridge (kept here in [0.1, 1.0])",
       x = "alpha (mixing parameter)",
       y = "# selected covariates") +
  theme_minimal()
ggsave(file.path(out_dir, "elastic_net_alpha.png"),
       p_alpha, width = 8, height = 5, dpi = 150)

cat("\nOutputs in:", out_dir, "\n")
