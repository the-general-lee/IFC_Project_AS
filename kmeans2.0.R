###---------------------------------------------------###
### K-Means Clustering - Composite Fragility Index    ###
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
library(viridis)    
library(grid) # Caricata esplicitamente per evitare errori con textGrob

# --- 2. Parametri Globali ---------------------------------------------------
FILE_PATH  <- "Composite_fragility_index_Tutti_Anni.xlsx"
PERCORSO_SHAPEFILE <- "ISTAT shape file/Limiti2021/Com2021/Com2021.shp"
ANNI       <- c("2018", "2019", "2021", "2022")
K_MAX      <- 10    
NSTART     <- 25    
SET_SEED   <- 2026  

# Liste per salvare i grafici
lista_plot_elbow   <- list()
lista_plot_sil     <- list()
lista_plot_cluster <- list()
lista_plot_mappe   <- list()

# --- 3. Caricamento Mappa Base ----------------------------------------------
cat("\nCaricamento base cartografica...\n")
mappa_italia <- st_read(PERCORSO_SHAPEFILE, quiet = TRUE) %>%
  mutate(PRO_COM = as.character(PRO_COM))

# --- 4. Funzioni di Supporto ------------------------------------------------
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

# --- 5. ESECUZIONE ANALISI --------------------------------------------------
fogli_risultati <- list()

for (anno in ANNI) {
  cat(paste0("Elaborando ", anno, "...\n"))
  
  df_raw <- carica_foglio(FILE_PATH, anno)
  dati <- prepara_dati(df_raw)
  dati_scalati <- scale(dati$indici)
  
  # Scelta k (Migliore k > 2)
  w <- numeric(K_MAX); sil <- numeric(K_MAX)
  set.seed(SET_SEED)
  for (k in 1:K_MAX) {
    res <- kmeans(dati_scalati, centers = k, nstart = NSTART)
    w[k] <- sum(res$withinss) / res$totss
    if (k >= 2) sil[k] <- mean(silhouette(res$cluster, dist(dati_scalati))[, "sil_width"])
  }
  k_ottimale <- which.max(sil[3:K_MAX]) + 2 
  
  # K-Means finale
  set.seed(SET_SEED)
  mod_kmeans <- kmeans(dati_scalati, centers = k_ottimale, nstart = NSTART)
  
  # Profiling per nomi cluster
  centers <- as.data.frame(mod_kmeans$centers)
  centers$ClusterID <- 1:nrow(centers)
  id_urbani <- centers$ClusterID[which.min(centers$MFI_Decile)]
  id_socio  <- centers$ClusterID[which.max(centers$MFI_Decile)]
  id_territoriale <- setdiff(1:3, c(id_urbani, id_socio))
  
  # Assegnazione NOMI ABBREVIATI per la legenda
  dati$meta <- dati$meta %>%
    mutate(Cluster_Num = mod_kmeans$cluster) %>%
    mutate(Nome_Cluster = case_when(
      Cluster_Num == id_urbani       ~ "Poli Urbani",
      Cluster_Num == id_socio        ~ "Frag. Socio-Econ.",
      Cluster_Num == id_territoriale ~ "Frag. Territoriale",
      TRUE                           ~ "Altro"
    ))
  
  # --- Grafici Diagnostici ---
  lista_plot_elbow[[anno]] <- ggplot(data.frame(k=1:K_MAX, wss=w), aes(x=k, y=wss)) +
    geom_line(color="steelblue") + geom_point() + labs(title=anno) + theme_minimal()
  
  lista_plot_sil[[anno]] <- ggplot(data.frame(k=2:K_MAX, s=sil[2:K_MAX]), aes(x=k, y=s)) +
    geom_line(color="darkorange") + geom_point() + labs(title=anno) + theme_minimal() +
    geom_vline(xintercept = k_ottimale, linetype="dashed", color="red")
  
  # PCA Plot
  lista_plot_cluster[[anno]] <- fviz_cluster(mod_kmeans, data = dati_scalati, geom = "point", 
                                             palette = "jco", ggtheme = theme_minimal(), main = anno) +
    theme(legend.position = "none") # Nascondo la legenda qui per pulizia
  
  # --- MAPPA CON LEGENDA OTTIMIZZATA ---
  mappa_cluster <- mappa_italia %>%
    left_join(dati$meta, by = "PRO_COM") %>%
    filter(!is.na(Nome_Cluster))
  
  lista_plot_mappe[[anno]] <- ggplot(mappa_cluster) +
    geom_sf(aes(fill = Nome_Cluster), color = NA) +
    scale_fill_viridis_d(option = "plasma", name = "Legenda:") +
    labs(title = paste("Mappa", anno)) +
    theme_void() + 
    theme(
      legend.position = "bottom", 
      legend.text = element_text(size = 7),        # Testo più piccolo
      legend.title = element_text(size = 8, face = "bold"),
      legend.key.size = unit(0.4, "cm"),           # Quadratini colori più piccoli
      plot.title = element_text(hjust = 0.5, size = 10, face = "bold")
    ) +
    guides(fill = guide_legend(nrow = 1))          # Tutti i nomi su una riga
  
  fogli_risultati[[anno]] <- dati$meta
}

# --- STAMPA FINALE ---
cat("\nVisualizzazione griglie finali...\n")

# Griglia Elbow
grid.arrange(grobs = lista_plot_elbow, ncol = 2, top = textGrob("Elbow Method", gp=gpar(fontsize=14, font=2)))

# Griglia Silhouette
grid.arrange(grobs = lista_plot_sil, ncol = 2, top = textGrob("Silhouette Score", gp=gpar(fontsize=14, font=2)))

# Griglia Mappe (Legenda sotto ogni mappa)
grid.arrange(grobs = lista_plot_mappe, ncol = 2, top = textGrob("Mappe della Fragilità", gp=gpar(fontsize=14, font=2)))

# Salvataggio
write.xlsx(fogli_risultati, file = "Risultati_Legenda_OK.xlsx")