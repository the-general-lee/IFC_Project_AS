# 1. Load required libraries
if (!require("readxl")) install.packages("readxl")
if (!require("writexl")) install.packages("writexl")
if (!require("tidyverse")) install.packages("tidyverse")

library(readxl)
library(writexl)
library(tidyverse)

# --- CONFIGURATION ---
file_path <- "Composite_fragility_index_Tutti_Anni.xlsx"
sheets_to_read <- c("2018", "2019", "2021", "2022")

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

ref_2018 <- check_2018 %>%
  summarise(across(all_of(all_indicators), ~ mean(.x, na.rm = TRUE)))

limits <- data_transformed %>%
  summarise(across(all_of(all_indicators), list(
    min = ~ min(.x, na.rm = TRUE), 
    max = ~ max(.x, na.rm = TRUE)
  )))

# 6. AMPI normalization (transform to base-100 indices)
data_norm <- data_transformed
for(i in all_indicators) {
  ref <- as.numeric(ref_2018[[i]])
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

# 8. Final index calculation and table reshaping
final_result <- data_norm %>%
  rowwise() %>%
  mutate(IFC = calc_ampi_plus(c_across(all_of(all_indicators)))) %>%
  ungroup() %>%
  select(PRO_COM, Territory, Year, IFC) %>%
  mutate(IFC = round(IFC, 2)) %>%
  # Transform from vertical (Year) to horizontal (columns 2018, 2019...)
  pivot_wider(
    names_from = Year, 
    values_from = IFC,
    values_fn = list(IFC = ~ mean(.x, na.rm = TRUE)) # Anti-duplicate safeguard
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