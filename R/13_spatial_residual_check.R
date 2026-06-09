# ================================================================
# SPATIAL RESIDUAL DIAGNOSTIC for the best LMM (M_geo_nested)
# ----------------------------------------------------------------
# Question: after explaining IFC with (a) GRINS predictors and
# (b) a region/province nested random intercept, is there still
# spatial structure left in the residuals? If yes, neighbouring
# municipalities have correlated unexplained fragility — which
# would call for an explicit spatial model (CAR, SAR, or a Gaussian
# process) as the next step.
#
# Method: compute Moran's I on the BLUP-augmented residuals,
#         aggregated to the municipality level (averaging 2019/2021).
#
# Reference shapefile (Limiti ISTAT 2021) is expected at
#   data/raw/istat/shapefile/Com2021.shp
# and is NOT versioned in the repository (see README).
# ================================================================

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr)
  library(lme4); library(sf); library(spdep)
})

set.seed(2026)

ifc_file       <- "data/raw/ifc/final_analysis_sorted.xlsx"
grins_file     <- "data/processed/grins_v3/comunale_v3.rds"
tassonomia_file<- "data/raw/grins/tassonomia_grins.xlsx"
fit_file       <- "outputs/mixed_models/fit_M_geo_nested.rds"
shape_file     <- "data/raw/istat/shapefile/Com2021.shp"
out_dir        <- "outputs/mixed_models"

if (!file.exists(fit_file)) {
  stop("Missing LMM fit. Run R/11_mixed_models.R first.")
}
if (!file.exists(shape_file)) {
  stop("Missing ISTAT shapefile. See README for download instructions.")
}

# ================================================================
# 1) Recompute the data exactly as in script 11 (same FE pipeline)
#    so we can locate the rows the LMM saw.
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

grins_groups <- grins %>%
  select(codice_comune, anno, nome_regione, nome_provincia)

tassonomia <- read_excel(tassonomia_file, sheet = 1) %>%
  rename(PRO_COM = 1, GRINS_macroclass = 2, GRINS_class = 3) %>%
  mutate(PRO_COM = as.numeric(PRO_COM)) %>%
  filter(!is.na(PRO_COM)) %>%
  select(PRO_COM, GRINS_macroclass)

num_cols <- setdiff(names(grins)[sapply(grins, is.numeric)],
                    c("codice_comune", "anno"))

data <- ifc_long %>%
  inner_join(grins %>% select(all_of(c("codice_comune", "anno", num_cols))),
             by = c("PRO_COM" = "codice_comune", "year" = "anno")) %>%
  left_join(grins_groups,
            by = c("PRO_COM" = "codice_comune", "year" = "anno")) %>%
  left_join(tassonomia, by = "PRO_COM") %>%
  filter(!is.na(nome_regione), !is.na(nome_provincia), !is.na(GRINS_macroclass))

cat("Rows considered (matches script 11 N):", nrow(data), "\n")

# ================================================================
# 2) Load saved fit and compute residuals
# ================================================================
fit <- readRDS(fit_file)
res <- residuals(fit)
stopifnot(length(res) == nrow(data))

# Aggregate residuals to municipality level (average across 2019/2021)
muni_res <- data %>%
  mutate(residual = res) %>%
  group_by(PRO_COM) %>%
  summarise(mean_residual = mean(residual, na.rm = TRUE), .groups = "drop")

cat("Unique municipalities with residual:", nrow(muni_res), "\n")
cat("Mean residual:", round(mean(muni_res$mean_residual), 4),
    " | SD:", round(sd(muni_res$mean_residual), 4), "\n")

# ================================================================
# 3) Attach residuals to the shapefile geometry
# ================================================================
cat("\nReading shapefile...\n")
shapes <- st_read(shape_file, quiet = TRUE)

# ISTAT 2021 shapefile has a PRO_COM_T or similar code column; locate it
candidate <- intersect(c("PRO_COM", "PRO_COM_T", "COD_COM", "COD_ISTAT"),
                       names(shapes))
if (length(candidate) == 0) {
  candidate <- names(shapes)[
    sapply(shapes, function(v) is.numeric(v) || is.character(v))
  ]
  candidate <- candidate[grepl("PRO|COM|CODE", toupper(candidate))][1]
}
pcol <- candidate[1]
cat("Using shapefile code column:", pcol, "\n")

shapes$PRO_COM_num <- as.numeric(as.character(shapes[[pcol]]))
shapes_with_res <- shapes %>%
  inner_join(muni_res, by = c("PRO_COM_num" = "PRO_COM"))
cat("Shapefile rows matched:", nrow(shapes_with_res),
    "/", nrow(shapes), "\n")

# ================================================================
# 4) Build spatial weights matrix (k-nearest-neighbours, k=5)
# ----------------------------------------------------------------
# We use kNN rather than queen/rook contiguity because the latter
# would create islands (sea, lakes, enclaves) that have no neighbours
# and would be dropped from the test, biasing it.
# ================================================================
cat("\nBuilding kNN spatial weights (k=5)...\n")
centroids <- st_centroid(st_geometry(shapes_with_res), of_largest_polygon = TRUE)
coords    <- st_coordinates(centroids)

k_nn <- knearneigh(coords, k = 5)
nb   <- knn2nb(k_nn)
lw   <- nb2listw(nb, style = "W")

# ================================================================
# 5) Moran's I on the LMM residuals
# ================================================================
cat("\nComputing Moran's I (1000 permutations)...\n")
mc <- moran.mc(shapes_with_res$mean_residual, lw, nsim = 999,
               alternative = "greater")
cat("\n--- Moran's I on residuals of M_geo_nested ---\n")
print(mc)

# Standard analytical test as well (gives interpretable I + p-value)
mt <- moran.test(shapes_with_res$mean_residual, lw, alternative = "greater")
cat("\nAnalytic Moran test:\n")
print(mt)

result <- data.frame(
  test          = c("Monte Carlo (999 perm)", "Analytic"),
  statistic_I   = c(mc$statistic, unname(mt$estimate["Moran I statistic"])),
  expectation_I = c(NA, unname(mt$estimate["Expectation"])),
  variance_I    = c(NA, unname(mt$estimate["Variance"])),
  p_value       = c(mc$p.value, mt$p.value)
)
write.csv(result, file.path(out_dir, "moran_residuals.csv"), row.names = FALSE)
saveRDS(shapes_with_res, file.path(out_dir, "shapes_with_residuals.rds"))

# ================================================================
# 6) Quick verdict
# ================================================================
I_obs <- as.numeric(mc$statistic)
cat(sprintf("\nObserved Moran's I: %.4f\n", I_obs))
if (mc$p.value < 0.05 && I_obs > 0) {
  cat("=> Significant positive spatial autocorrelation IS still present.\n",
      "   The LMM with regional/provincial random effects has not absorbed\n",
      "   all of the spatial structure: neighbouring municipalities have\n",
      "   correlated unexplained fragility. A next step could be an\n",
      "   explicit spatial model (CAR/SAR or a Gaussian process).\n")
} else if (mc$p.value >= 0.05) {
  cat("=> No significant spatial autocorrelation in residuals.\n",
      "   The (1|region/province) hierarchy has effectively absorbed\n",
      "   the spatial structure; the LMM benchmark is adequate.\n")
} else {
  cat("=> Significant NEGATIVE spatial autocorrelation (unusual).\n",
      "   May indicate over-smoothing by the random effects; investigate.\n")
}

cat("\nOutputs in:", out_dir, "\n")
