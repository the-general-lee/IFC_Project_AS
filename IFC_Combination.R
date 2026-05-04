# install the packages needed.
if (!require("readxl")) install.packages("readxl")
if (!require("dplyr")) install.packages("dplyr")
if (!require("writexl")) install.packages("writexl")

library(readxl)
library(dplyr)
library(writexl)

# 1. set the path
file_2019_path <- "Composite fragility index - all municipalities 2019.xlsx"
file_2021_path <- "Composite fragility index - all municipalities 2021.xlsx"
file_final_path <- "IFC_Final_Analysis_Sorted.xlsx"

# 2. read the data
df_2019_raw <- read_excel(file_2019_path, skip = 4)
df_2021_raw <- read_excel(file_2021_path, skip = 4)
df_final <- read_excel(file_final_path)

colnames(df_2019_raw)[1] <- "PRO_COM"
colnames(df_2021_raw)[1] <- "PRO_COM"

ifc_mapping_2019 <- df_final %>% select(PRO_COM, IFC_2019)
ifc_mapping_2021 <- df_final %>% select(PRO_COM, IFC_2021)

df_2019_updated <- df_2019_raw %>%
  left_join(ifc_mapping_2019, by = "PRO_COM") %>%
  relocate(IFC_2019, .after = PRO_COM)

df_2021_updated <- df_2021_raw %>%
  left_join(ifc_mapping_2021, by = "PRO_COM") %>%
  relocate(IFC_2021, .after = PRO_COM)

write_xlsx(df_2019_updated, "Updated_Fragility_Index_2019.xlsx")
write_xlsx(df_2021_updated, "Updated_Fragility_Index_2021.xlsx")