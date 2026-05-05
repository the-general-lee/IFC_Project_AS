# ============================================================
# REGRESSION TASK (R VERSION)
# IFC Decile ~ GRINS Macroclass * Year
# Years: 2018, 2019, 2021, 2022
# ============================================================

# -----------------------------
# 1. Libraries
# -----------------------------
library(readxl)
library(dplyr)
library(ggplot2)
library(writexl)

# -----------------------------
# 2. File paths
# -----------------------------
ifc_file <- "/Users/shuvroahmed/Desktop/Imaging Program Academics/POLIMI/sem 2/applied stat/projects/Composite_fragility_index_Tutti_Anni.xlsx"

grins_file <- "/Users/shuvroahmed/Desktop/Imaging Program Academics/POLIMI/sem 2/applied stat/projects/Tassonomia GRINS_com2021_v2023-11-16_SOLO CLASSI GRINS(Classi GRINS_com2021).csv"

years <- c(2018, 2019, 2021, 2022)

# -----------------------------
# 3. Read IFC data for all years
# -----------------------------
all_data <- list()

for (year in years) {
  
  df <- read_excel(ifc_file, sheet = as.character(year))
  
  # Remove second header row
  df <- df[-1, ]
  
  # Rename columns
  colnames(df)[1] <- "PRO_COM"
  colnames(df)[2] <- "Municipality"
  colnames(df)[3] <- "IFC_decile"
  
  df <- df %>%
    select(PRO_COM, Municipality, IFC_decile)
  
  df <- df %>%
    mutate(
      PRO_COM = as.numeric(PRO_COM),
      IFC_decile = as.numeric(IFC_decile),
      Year = year
    ) %>%
    filter(!is.na(PRO_COM), !is.na(IFC_decile))
  
  all_data[[as.character(year)]] <- df
}

ifc_all <- bind_rows(all_data)

# -----------------------------
# 4. Read GRINS taxonomy
# -----------------------------
grins <- read.csv(
  grins_file,
  sep = ";",
  stringsAsFactors = FALSE
)

print(colnames(grins))
print(head(grins))

# Rename columns safely
grins <- grins %>%
  rename(
    PRO_COM = 1,
    GRINS_macroclass = 2,
    GRINS_class = 3
  ) %>%
  mutate(PRO_COM = as.numeric(PRO_COM)) %>%
  filter(!is.na(PRO_COM))

# -----------------------------
# 5. Merge datasets
# -----------------------------
data <- ifc_all %>%
  inner_join(grins, by = "PRO_COM")

# Pre/Post COVID variable
data <- data %>%
  mutate(
    Period = ifelse(Year %in% c(2018, 2019), "Pre-COVID", "Post-COVID"),
    Year = as.factor(Year),
    GRINS_macroclass = as.factor(GRINS_macroclass),
    Period = as.factor(Period)
  )

cat("Merged rows:", nrow(data), "\n")

# -----------------------------
# 6. Summary table
# -----------------------------
summary_table <- data %>%
  group_by(Year, Period, GRINS_macroclass) %>%
  summarise(
    count = n(),
    mean = mean(IFC_decile, na.rm = TRUE),
    sd = sd(IFC_decile, na.rm = TRUE)
  )

print(summary_table)

# -----------------------------
# 7. Regression 1: Year model
# -----------------------------
model_year <- lm(
  IFC_decile ~ GRINS_macroclass * Year,
  data = data
)

cat("\n===== YEAR REGRESSION =====\n")
summary(model_year)

# -----------------------------
# 8. Regression 2: Pre/Post COVID
# -----------------------------
model_period <- lm(
  IFC_decile ~ GRINS_macroclass * Period,
  data = data
)

cat("\n===== PRE/POST COVID REGRESSION =====\n")
summary(model_period)

# -----------------------------
# 9. Predictions by Year
# -----------------------------
pred_year <- data %>%
  select(GRINS_macroclass, Year) %>%
  distinct()

pred_year$predicted_IFC_decile <- predict(model_year, newdata = pred_year)

print(pred_year)

# -----------------------------
# 10. Predictions by Period
# -----------------------------
pred_period <- data %>%
  select(GRINS_macroclass, Period) %>%
  distinct()

pred_period$predicted_IFC_decile <- predict(model_period, newdata = pred_period)

print(pred_period)

# -----------------------------
# 11. Save outputs
# -----------------------------
write.csv(data, "/Users/shuvroahmed/Desktop/Imaging Program Academics/POLIMI/sem 2/applied stat/projects/regression_output/merged_IFC_GRINS_all_years.csv", row.names = FALSE)
write.csv(summary_table, "/Users/shuvroahmed/Desktop/Imaging Program Academics/POLIMI/sem 2/applied stat/projects/regression_output/summary_IFC_by_year_macroclass.csv", row.names = FALSE)
write.csv(pred_year, "/Users/shuvroahmed/Desktop/Imaging Program Academics/POLIMI/sem 2/applied stat/projects/regression_output/predicted_IFC_by_year_macroclass.csv", row.names = FALSE)
write.csv(pred_period, "/Users/shuvroahmed/Desktop/Imaging Program Academics/POLIMI/sem 2/applied stat/projects/regression_output/predicted_IFC_by_period_macroclass.csv", row.names = FALSE)

# Save coefficients
write.csv(summary(model_year)$coefficients, "/Users/shuvroahmed/Desktop/Imaging Program Academics/POLIMI/sem 2/applied stat/projects/regression_output/regression_coefficients_year.csv")
write.csv(summary(model_period)$coefficients, "/Users/shuvroahmed/Desktop/Imaging Program Academics/POLIMI/sem 2/applied stat/projects/regression_output/regression_coefficients_period.csv")

# -----------------------------
# 12. Plot
# -----------------------------
ggplot(plot_data, aes(
  x = Year, 
  y = mean_IFC, 
  group = GRINS_macroclass, 
  color = GRINS_macroclass   # <-- THIS FIXES IT
)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Mean IFC Decile by GRINS Macroclass Over Time",
    x = "Year",
    y = "Mean IFC Decile",
    color = "GRINS Macroclass"   # legend name
  ) +
  theme_minimal()

ggsave("/Users/shuvroahmed/Desktop/Imaging Program Academics/POLIMI/sem 2/applied stat/projects/regression_output/mean_IFC_decile_plot.png")