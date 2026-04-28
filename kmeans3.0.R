###---------------------------------------------------###
### K-Means - Colori Alto Contrasto e K Dinamico      ###
###---------------------------------------------------###

# --- 1. Librerie ------------------------------------------------------------
library(readxl)     
library(cluster)    
library(factoextra) 
library(ggplot2)    
library(dplyr)      
library(openxlsx)   
library(gridExtra)  
library(sf)         
library(grid) 

# --- 2. Parametri Globali ---------------------------------------------------
FILE_PATH  <- "Composite_fragility_index_Tutti_Anni.xlsx"
PERCORSO_SHAPEFILE <- "ISTAT shape file/Limiti2021/Com2021/Com2021.shp"
ANNI       <- c("2018", "2019", "2021", "2022")
K_MAX      <- 10    
NSTART     <- 25    
SET_SEED   <- 2026  
CARTELLE_OUTPUT <- "Output_Grafici"

if (!dir.exists(CARTELLE_OUTPUT)) dir.create(CARTELLE_OUTPUT)

# --- 3. DEFINIZIONE COLORI AD ALTO CONTRASTO --------------------------------
# Colori visivamente molto distanti tra loro per distinguere bene i confini dei cluster
colori_fissi <- c(
  "Livello 1" = "#00429d", # Blu scuro (Meno fragile in assoluto)
  "Livello 2" = "#2e8ad8", # Azzurro acceso
  "Livello 3" = "#00fa9a", # Verde smeraldo
  "Livello 4" = "#b8df29", # Verde lime
  "Livello 5" = "#ffeb3b", # Giallo intenso
  "Livello 6" = "#ff9800", # Arancione forte
  "Livello 7" = "#f44336", # Rosso acceso
  "Livello 8" = "#c51b7d", # Fucsia/Magenta
  "Livello 9" = "#4a148c", # Viola scuro
  "Livello 10"= "#000000"  # Nero (Massima fragilità)
)

# I 3 gruppi di variabili
gruppi_variabili <- list(
  "Economia" = c("Employment_rate_20_64", "Density_local_units_ind_serv", "Persons_low_productivity_units"),
  "Sociale"  = c("Accessibility_essential_services", "Population_dependency_index", "Pop_low_education_25_64", "Incidence_net_migration"),
  "Ambiente" = c("Motorisation_rate_high_emissions", "Undiff_waste_generated", "Protected_natural_areas", "Areas_risk_landslides", "Land_consumption")
)

plot_list <- list()
profili_cluster <- list()

for (g in names(gruppi_variabili)) {
  plot_list[[g]] <- list(Elbow = list(), Sil = list(), Mappa = list())
  profili_cluster[[g]] <- list()
}

# --- 4. Caricamento Mappa Base ----------------------------------------------
cat("\nCaricamento base cartografica...\n")
mappa_italia <- st_read(PERCORSO_SHAPEFILE, quiet = TRUE) %>%
  mutate(PRO_COM = as.character(PRO_COM))

# --- 5. Funzioni di Supporto ------------------------------------------------
carica_foglio <- function(file_path, anno) {
  df <- read_excel(file_path, sheet = anno, skip = 2, col_names = FALSE)
  colnames(df) <- c("PRO_COM", "Territory", "MFI_Decile", "Motorisation_rate_high_emissions", 
                    "Undiff_waste_generated", "Protected_natural_areas", "Areas_risk_landslides", 
                    "Land_consumption", "Accessibility_essential_services", "Population_dependency_index",
                    "Pop_low_education_25_64", "Employment_rate_20_64", "Incidence_net_migration",
                    "Density_local_units_ind_serv", "Persons_low_productivity_units")
  return(df)
}

prepara_dati <- function(df) {
  indici <- df %>% select(MFI_Decile:Persons_low_productivity_units) %>% mutate(across(everything(), as.numeric))
  indici_puliti <- na.omit(indici)
  meta <- df[complete.cases(indici), c("PRO_COM", "Territory")] %>% mutate(PRO_COM = as.character(PRO_COM))
  list(indici = indici_puliti, meta = meta)
}

# --- 6. ESECUZIONE ANALISI --------------------------------------------------
fogli_risultati <- list()

for (anno in ANNI) {
  cat(paste0("\nElaborazione Anno: ", anno, "...\n"))
  df_raw <- carica_foglio(FILE_PATH, anno)
  dati <- prepara_dati(df_raw)
  risultati_anno <- dati$meta
  
  for (nome_gruppo in names(gruppi_variabili)) {
    vars <- gruppi_variabili[[nome_gruppo]]
    dati_scalati <- scale(dati$indici[, vars])
    
    # Ricerca K dinamico (Minimo 3)
    w <- numeric(K_MAX); sil <- numeric(K_MAX)
    set.seed(SET_SEED)
    for (k in 1:K_MAX) {
      res <- kmeans(dati_scalati, centers = k, nstart = NSTART)
      w[k] <- sum(res$withinss) / res$totss
      if (k >= 2) sil[k] <- mean(silhouette(res$cluster, dist(dati_scalati))[, "sil_width"])
    }
    
    k_ottimale <- which.max(sil[3:K_MAX]) + 2 
    cat(paste0("  -> ", nome_gruppo, ": K scelto = ", k_ottimale, "\n"))
    
    set.seed(SET_SEED)
    mod_kmeans <- kmeans(dati_scalati, centers = k_ottimale, nstart = NSTART)
    
    # Ordinamento dei cluster in base alla Fragilità Media (MFI)
    mfi_temp <- aggregate(dati$indici$MFI_Decile, by = list(Cluster_Num = mod_kmeans$cluster), FUN = mean)
    mfi_temp$Livello <- rank(mfi_temp$x, ties.method = "first") 
    
    mapping_cluster <- data.frame(
      Cluster_Num = mfi_temp$Cluster_Num,
      Nome_Cluster = paste("Livello", mfi_temp$Livello)
    )
    
    risultati_gruppo <- data.frame(
      PRO_COM = dati$meta$PRO_COM,
      Cluster_Num = mod_kmeans$cluster
    ) %>%
      left_join(mapping_cluster, by = "Cluster_Num")
    
    risultati_anno[[paste0("Cluster_", nome_gruppo)]] <- risultati_gruppo$Nome_Cluster
    
    # Profili per analisi
    dati_reali_gruppo <- dati$indici[, vars, drop = FALSE]
    medie_cluster <- aggregate(dati_reali_gruppo, by = list(Tipologia_Cluster = risultati_gruppo$Nome_Cluster), FUN = function(x) round(mean(x), 2))
    medie_cluster$MFI_Medio <- round(aggregate(dati$indici$MFI_Decile, by = list(risultati_gruppo$Nome_Cluster), FUN = mean)$x, 2)
    profili_cluster[[nome_gruppo]][[anno]] <- medie_cluster
    
    # --- CREAZIONE PLOT ---
    plot_list[[nome_gruppo]]$Elbow[[anno]] <- ggplot(data.frame(k=1:K_MAX, wss=w), aes(x=k, y=wss)) +
      geom_line(color="steelblue") + geom_point() + labs(title=anno) + theme_minimal()
    
    plot_list[[nome_gruppo]]$Sil[[anno]] <- ggplot(data.frame(k=2:K_MAX, s=sil[2:K_MAX]), aes(x=k, y=s)) +
      geom_line(color="darkorange") + geom_point() + labs(title=anno) + theme_minimal() +
      geom_vline(xintercept = k_ottimale, linetype="dashed", color="red")
    
    mappa_cluster <- mappa_italia %>%
      left_join(risultati_gruppo, by = "PRO_COM") %>%
      filter(!is.na(Nome_Cluster))
    
    plot_list[[nome_gruppo]]$Mappa[[anno]] <- ggplot(mappa_cluster) +
      geom_sf(aes(fill = Nome_Cluster), color = NA) +
      # Applichiamo la NUOVA scala colori
      scale_fill_manual(values = colori_fissi, name = "Grado di Fragilità") + 
      labs(title = paste(anno, "| K =", k_ottimale)) +
      theme_void() + 
      theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, size = 11, face = "bold"))
  }
  fogli_risultati[[anno]] <- risultati_anno
}

# --- 7. SALVATAGGIO IMMAGINI E STAMPA DATI ----------------------------------
cat("\nSalvataggio grafici in '", CARTELLE_OUTPUT, "'...\n")

for (nome_gruppo in names(gruppi_variabili)) {
  g_elbow <- arrangeGrob(grobs = plot_list[[nome_gruppo]]$Elbow, ncol = 2, top = textGrob(paste("Elbow -", nome_gruppo), gp=gpar(fontsize=14, font=2)))
  ggsave(filename = file.path(CARTELLE_OUTPUT, paste0("Elbow_", nome_gruppo, ".png")), plot = g_elbow, width = 10, height = 8)
  
  g_sil <- arrangeGrob(grobs = plot_list[[nome_gruppo]]$Sil, ncol = 2, top = textGrob(paste("Silhouette -", nome_gruppo), gp=gpar(fontsize=14, font=2)))
  ggsave(filename = file.path(CARTELLE_OUTPUT, paste0("Silhouette_", nome_gruppo, ".png")), plot = g_sil, width = 10, height = 8)
  
  g_mappa <- arrangeGrob(grobs = plot_list[[nome_gruppo]]$Mappa, ncol = 2, top = textGrob(paste("Mappe -", nome_gruppo), gp=gpar(fontsize=14, font=2)))
  ggsave(filename = file.path(CARTELLE_OUTPUT, paste0("Mappe_", nome_gruppo, ".png")), plot = g_mappa, width = 10, height = 12)
  
  cat(paste0("\n--- PROFILI ORDINATI: ", toupper(nome_gruppo), " ---\n"))
  for (anno in ANNI) {
    cat(paste0("\n> Anno: ", anno, "\n"))
    print(profili_cluster[[nome_gruppo]][[anno]])
  }
}

write.xlsx(fogli_risultati, file = "Risultati_KMeans_ColoriFissi.xlsx")
cat("\nAnalisi conclusa con successo! Controlla i grafici.\n")