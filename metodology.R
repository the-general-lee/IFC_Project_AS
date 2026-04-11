# 1. Load required libraries
options(repos = c(CRAN = "https://cloud.r-project.org"))

load_or_install <- function(package_name) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    install.packages(package_name)
  }
  library(package_name, character.only = TRUE)
}

load_or_install("readxl")
load_or_install("writexl")
load_or_install("tidyverse")

# --- CONFIGURATION ---
file_path <- "Composite_fragility_index_Tutti_Anni.xlsx"
sheets_to_read <- c("2018", "2019", "2021")

# 2. Load and deep-clean function
load_clean_data <- function(s) {
  message(paste("Processing sheet:", s))
  
  # Read everything as text to avoid Excel format errors
  df <- read_excel(file_path, sheet = s, skip = 2, col_names = FALSE, col_types = "text")
  
  # Select only the first 15 columns (Code, Territory, Decile + 12 indicators)
  df <- df[, 1:15]
  colnames(df) <- c("PRO_COM", "Territory", "MFI_Dec", paste0("Ind", 1:12))
  
  # Data cleaning:
  df_clean <- df %>%
    # Remove fully empty rows or rows without municipality code
    filter(!is.na(PRO_COM), !is.na(Territory)) %>%
    mutate(MFI_Dec = as.integer(readr::parse_number(MFI_Dec))) %>%
    # Convert text to numbers (handling Italian commas and spaces)
    mutate(across(4:15, ~ as.numeric(gsub(",", ".", gsub("\\s+", "", .x))))) %>%
    mutate(Year = as.character(s)) %>%
    # Remove possible duplicates (same municipality in the same year)
    distinct(PRO_COM, Territory, Year, .keep_all = TRUE)
  
  return(df_clean)
}

# 3. Load data for all years
data_all <- map_df(sheets_to_read, load_clean_data)

# 4. Polarity inversion (according to Istat methodology)
# Invert indicators where "High = Good" to make them "High = Fragile"
# Ind3: Protected Areas, Ind9: Employment, Ind10: Migration, Ind11: Business Density
cols_invert <- c("Ind3", "Ind9", "Ind10", "Ind11")
all_indicators <- paste0("Ind", 1:12)

data_transformed <- data_all %>%
  mutate(across(all_of(cols_invert), ~ - .x))

# 5. Goalposts calculation (Italy 2018 reference = 100)
check_2018 <- data_transformed %>% filter(Year == "2018")
if (nrow(check_2018) == 0) stop("Error: Unable to find 2018 data for reference calculation!")

# NOTE: The official methodology uses Italy-wide 2018 indicator values as Ref_xj.
# In this reconstruction, Ref_xj is approximated with the 2018 municipal mean because
# the national reference table is not available in the workbook.
ref_2018_proxy <- check_2018 %>%
  summarise(across(all_of(all_indicators), ~ mean(.x, na.rm = TRUE)))

limits <- data_transformed %>%
  summarise(across(all_of(all_indicators), list(
    min = ~ min(.x, na.rm = TRUE), 
    max = ~ max(.x, na.rm = TRUE)
  )))

# 6. AMPI normalization (transform to base-100 indices)
data_norm <- data_transformed
for(i in all_indicators) {
  ref <- as.numeric(ref_2018_proxy[[i]])
  mi <- as.numeric(limits[[paste0(i, "_min")]])
  ma <- as.numeric(limits[[paste0(i, "_max")]])
  
  if(!is.na(ref)) {
    delta <- (ma - mi) / 2
    # AMPI+ formula: centers the distribution on 100 with width 60
    data_norm[[i]] <- ((data_norm[[i]] - (ref - delta)) / (2 * delta)) * 60 + 70
  }
}

# 7. AMPI+ aggregation function (Mean + Imbalance penalty)
calc_ampi_plus <- function(row_vals) {
  # Remove NA values from the row
  row_vals <- row_vals[!is.na(row_vals)]
  if(length(row_vals) < 2) return(NA) # Need at least mean and standard deviation
  
  m <- mean(row_vals)
  s <- sd(row_vals)
  if(is.na(s) || m == 0) return(m)
  
  # Penalty = Standard Deviation * Coefficient of Variation
  return(m + (s * (s/m)))
}

build_decile_thresholds <- function(ifc_values) {
  clean_values <- ifc_values[!is.na(ifc_values)]
  if(length(clean_values) == 0) stop("Error: Unable to derive 2018 decile thresholds without valid IFC values!")

  as.numeric(quantile(
    clean_values,
    probs = seq(0.1, 1, 0.1),
    na.rm = TRUE,
    names = FALSE,
    type = 7
  ))
}

assign_deciles_from_thresholds <- function(ifc_values, upper_bounds) {
  if(length(upper_bounds) != 10) stop("Error: Expected 10 upper bounds for decile assignment!")

  assigned_deciles <- rep(NA_integer_, length(ifc_values))
  non_missing <- !is.na(ifc_values)
  assigned_deciles[non_missing] <- findInterval(ifc_values[non_missing], upper_bounds[1:9], left.open = TRUE) + 1L

  pmin(assigned_deciles, 10L)
}

# 8. Final index calculation and table reshaping
ifc_results <- data_norm %>%
  rowwise() %>%
  mutate(IFC = calc_ampi_plus(c_across(all_of(all_indicators)))) %>%
  ungroup() %>%
  select(PRO_COM, Territory, Year, MFI_Dec, IFC)

decile_thresholds <- ifc_results %>%
  filter(Year == "2018") %>%
  pull(IFC) %>%
  build_decile_thresholds()

final_result <- ifc_results %>%
  mutate(
    Calc_Decile = assign_deciles_from_thresholds(IFC, decile_thresholds),
    Decile_Match = if_else(!is.na(Calc_Decile) & !is.na(MFI_Dec), Calc_Decile == MFI_Dec, NA)
  ) %>%
  filter(Year != "2018") %>%
  mutate(IFC = round(IFC, 2)) %>%
  # Transform from vertical (Year) to horizontal with year-specific IFC and decile checks
  pivot_wider(
    names_from = Year,
    values_from = c(IFC, Calc_Decile, MFI_Dec, Decile_Match),
    names_glue = "{.value}_{Year}",
    values_fn = list(
      IFC = ~ mean(.x, na.rm = TRUE),
      Calc_Decile = dplyr::first,
      MFI_Dec = dplyr::first,
      Decile_Match = dplyr::first
    )
  ) %>%
  select(
    PRO_COM,
    Territory,
    IFC_2019,
    Calc_Decile_2019,
    MFI_Dec_2019,
    Decile_Match_2019,
    IFC_2021,
    Calc_Decile_2021,
    MFI_Dec_2021,
    Decile_Match_2021
  ) %>%
  # ORDERING: From the smallest municipality code to the largest
  arrange(as.numeric(PRO_COM))

# 9. Final export
write_xlsx(final_result, "IFC_Final_Analysis_Sorted.xlsx")

message("---------------------------------------------------------")
message("OPERATION COMPLETED!")
message("Generated file: IFC_Final_Analysis_Sorted.xlsx")
message("Data are sorted by PRO_COM code (from smallest to largest).")
message("---------------------------------------------------------")