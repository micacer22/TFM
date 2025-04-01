# ----------------------------------------------------------------------------
# 1) PREPARACIÓN: librerías y lista de comparaciones
# ----------------------------------------------------------------------------
# Instalar paquetes si no están instalados (descomentar si te hace falta)
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("DESeq2")
# install.packages("biomartr")
# BiocManager::install("clusterProfiler")
# BiocManager::install("org.At.tair.db")

library(DESeq2)
library(tidyverse)
library(EnhancedVolcano)
library(pheatmap)
library(biomartr)
library(clusterProfiler)
library(org.At.tair.db)
library(gridExtra)  # Para arrange de plots

# Ajustar directorio de trabajo según corresponda
setwd("/Users/miriamcaballerocervero/Documents/Miriam/24-25/TFM/Gráficas_Expr_Dif/Comparaciones")

# Ya generaste los archivos CSV de DESeq2 en esta carpeta.
# Ahora, solo necesitas leerlos y hacer el enriquecimiento GO.

# Lista de comparaciones con nombres de archivos (cada sublista tiene un solo elemento)
comparaciones <- list(
  c("DESeq2_results_WT_responsive.csv"),
  c("DESeq2_results_sdg8_responsive.csv"),
  c("DESeq2_results_clf_responsive.csv"),
  c("DESeq2_results_atx1_responsive.csv"),
  c("DESeq2_results_sdg8_mock_primed.csv"),
  c("DESeq2_results_clf_mock_primed.csv"),
  c("DESeq2_results_atx1_mock_primed.csv"),
  c("DESeq2_results_sdg8_infected_primed.csv"),
  c("DESeq2_results_clf_infected_primed.csv"),
  c("DESeq2_results_atx1_infected_primed.csv")
)

# ----------------------------------------------------------------------------
# 2) DEFINIR FUNCIÓN DE ENRIQUECIMIENTO GO
# ----------------------------------------------------------------------------
# Esta función:
#  - Lee un archivo CSV con resultados DESeq2 (con columna "X" = geneID),
#  - Filtra genes up/down con padj <= 0.05,
#  - Anota con biomartr,
#  - Corre enrichGO para BP, CC y MF,
#  - Retorna una lista con los resultados de cada ontología.

enriquecer_comparacion <- function(csv_file, out_prefix = "GO") {
  
  # Leer resultados DESeq2
  res <- read.csv(csv_file, header = TRUE)
  
  # 1) Define el cutoff de significancia
  padj_cutoff <- 0.05
  
  # 2) Filtra genes Up (log2FoldChange > 0) con padj <= 0.05
  genes_up <- unique(as.character(
    res$X[res$padj <= padj_cutoff & res$log2FoldChange > 0],
    na.rm = TRUE
  ))
  
  # 3) Filtra genes Down (log2FoldChange < 0) con padj <= 0.05
  genes_down <- unique(as.character(
    res$X[res$padj <= padj_cutoff & res$log2FoldChange < 0],
    na.rm = TRUE
  ))
  
  cat("\nArchivo:", csv_file, "\n")
  cat("Número de genes Up   =", length(genes_up), "\n")
  cat("Número de genes Down =", length(genes_down), "\n")
  
  # 4) Definir el universo de genes (todos los presentes en el CSV)
  #    Si deseas, puedes usar rownames(dds) u otra lista más completa.
  all_genes <- as.character(res$X)
  
  # 5) Anotar con biomartr (IDs ENTREZ)
  GO.attributes <- c("tair_symbol", "uniprotswissprot", "entrezgene_id")
  options(timeout = 100000000)
  
  ann_universe <- biomart(
    genes     = all_genes,
    mart      = "plants_mart",
    dataset   = "athaliana_eg_gene",
    attributes = GO.attributes,
    filters   = "ensembl_gene_id"
  )
  ann_universe$entrezgene_id <- as.character(ann_universe$entrezgene_id)
  universe_entrez <- unique(ann_universe$entrezgene_id)
  
  ann_up <- biomart(
    genes     = genes_up,
    mart      = "plants_mart",
    dataset   = "athaliana_eg_gene",
    attributes = GO.attributes,
    filters   = "ensembl_gene_id"
  )
  ann_down <- biomart(
    genes     = genes_down,
    mart      = "plants_mart",
    dataset   = "athaliana_eg_gene",
    attributes = GO.attributes,
    filters   = "ensembl_gene_id"
  )
  
  up_entrez   <- unique(ann_up$entrezgene_id)
  down_entrez <- unique(ann_down$entrezgene_id)
  
  # 6) Enriquecer con clusterProfiler
  enrich_genes <- function(genes, ont) {
    enrichGO(
      gene          = genes,
      universe      = universe_entrez,
      OrgDb         = org.At.tair.db,
      keyType       = "ENTREZID",
      ont           = ont,
      pAdjustMethod = "BH",
      qvalueCutoff  = 0.05,
      readable      = TRUE
    )
  }
  
  res_bp_up   <- simplify(enrich_genes(up_entrez,   "BP"))
  res_bp_down <- simplify(enrich_genes(down_entrez, "BP"))
  res_cc_up   <- simplify(enrich_genes(up_entrez,   "CC"))
  res_cc_down <- simplify(enrich_genes(down_entrez, "CC"))
  res_mf_up   <- simplify(enrich_genes(up_entrez,   "MF"))
  res_mf_down <- simplify(enrich_genes(down_entrez, "MF"))
  
  # 7) Convertir a data.frame
  df_bp_up   <- as.data.frame(res_bp_up)
  df_bp_down <- as.data.frame(res_bp_down)
  df_cc_up   <- as.data.frame(res_cc_up)
  df_cc_down <- as.data.frame(res_cc_down)
  df_mf_up   <- as.data.frame(res_mf_up)
  df_mf_down <- as.data.frame(res_mf_down)
  
  # 8) Retornar en una lista
  list(
    BP_up   = df_bp_up,
    BP_down = df_bp_down,
    CC_up   = df_cc_up,
    CC_down = df_cc_down,
    MF_up   = df_mf_up,
    MF_down = df_mf_down
  )
}

# ----------------------------------------------------------------------------
# 3) APLICAR FUNCIÓN A CADA COMPARACIÓN
# ----------------------------------------------------------------------------
resultados_enriquecimiento <- list()

for (comp in comparaciones) {
  # Aquí, comp es algo como c("DESeq2_results_WT_responsive.csv")
  # comp[[1]] extrae la cadena "DESeq2_results_WT_responsive.csv"
  archivo_csv <- comp[[1]]
  
  # Crear un nombre identificativo
  nombre_corto <- gsub("^DESeq2_results_", "", archivo_csv)
  nombre_corto <- gsub("\\.csv$", "", nombre_corto)
  
  cat("\nProcesando archivo:", archivo_csv, "\n")
  res_go <- enriquecer_comparacion(archivo_csv)
  
  # Guardar en la lista
  resultados_enriquecimiento[[nombre_corto]] <- res_go
  
  # (Opcional) Exportar cada tabla a CSV
  # write.csv(res_go$BP_up,   paste0("GO_BP_up_",   nombre_corto, ".csv"), row.names = FALSE)
  # write.csv(res_go$BP_down, paste0("GO_BP_down_", nombre_corto, ".csv"), row.names = FALSE)
  # ...
}

# ----------------------------------------------------------------------------
# 4) EJEMPLO DE VISUALIZACIÓN
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# Suponiendo que tu objeto 'resultados_enriquecimiento' tiene la forma:
# resultados_enriquecimiento[["WT_responsive"]]$BP_up
# resultados_enriquecimiento[["WT_responsive"]]$BP_down
# resultados_enriquecimiento[["WT_responsive"]]$CC_up
# ...
# Y así para cada comparación (WT_responsive, sdg8_responsive, etc.).
# ----------------------------------------------------------------------------

# 1) Definimos una función auxiliar para crear y guardar el dotplot
#plot_enrichment_dotplot <- function(df_enrich, comp_name, sub_name) {
  # Filtramos términos con p.adjust < 0.05
  #df_sig <- subset(df_enrich, p.adjust < 0.05)
  
  # Si no hay términos significativos, avisamos y salimos
  #if (nrow(df_sig) == 0) {
   # cat("Sin términos significativos (p.adjust < 0.05) para", sub_name, "en", comp_name, "\n")
   # return(NULL)
  #}
  
  # Creamos el gráfico
 # p <- ggplot(df_sig, aes(x = Count, y = reorder(Description, Count))) +
    #geom_point(aes(size = Count, color = p.adjust)) +
   # scale_color_gradient(low = "red", high = "yellow") +
   # theme_minimal() +
   # labs(
     # title = paste(sub_name, "-", comp_name),
   #   x     = "Gene Count",
    #  y     = "GO Term",
   #   color = "p.adjust",
   #   size  = "Gene Count"
   # )
  
  # Nombre de archivo, por ejemplo: "dotplot_BP_up_WT_responsive.pdf"
  #file_name <- paste0("dotplot_", sub_name, "_", comp_name, ".pdf")
 # ggsave(file_name, p)
  
  # Retornamos el objeto ggplot por si quieres manipularlo en R
  #return(p)
#}

# 2) Generamos los dotplots para cada comparación y cada subcategoría
#    (BP_up, BP_down, CC_up, CC_down, MF_up, MF_down) de tu lista.

# Vector con los nombres de sub-listas que tienes en cada comparación
#ontologias <- c("BP_up", "BP_down", "CC_up", "CC_down", "MF_up", "MF_down")

# Iterar sobre cada comparación en 'resultados_enriquecimiento'
#for (comp_name in names(resultados_enriquecimiento)) {
  
 # cat("\nProcesando comparación:", comp_name, "\n")
  
  # Extraer la lista de enriquecimiento para esta comparación
 # sub_res <- resultados_enriquecimiento[[comp_name]]
  
  # Iterar sobre cada ontología/dirección (BP_up, BP_down, etc.)
#  for (ont_name in ontologias) {
  #  if (!ont_name %in% names(sub_res)) {
 #     cat("No existe", ont_name, "en la comparación", comp_name, "\n")
 #     next
  #  }
    
    # Extraer el data frame de enriquecimiento
  #  df_enrich <- sub_res[[ont_name]]
    
    # Generar y guardar el dotplot
   # plot_enrichment_dotplot(df_enrich, comp_name, ont_name)
 # }
#}

# ----------------------------------------------------------------------------
# 4) VISUALIZACIÓN: GRÁFICO COMBINADO (Up y Down en dos columnas)
# ----------------------------------------------------------------------------
plot_combined_dotplot <- function(df_up, df_down, comp_name, ontologia, top_n = NULL) {
  # Procesamiento de df_up y df_down
  if (!is.null(df_up) && nrow(df_up) > 0 && "p.adjust" %in% colnames(df_up)) {
    df_up$p.adjust <- as.numeric(as.character(df_up$p.adjust))
    df_up$Direction <- "Up"
  } else {
    df_up <- data.frame()
  }
  if (!is.null(df_down) && nrow(df_down) > 0 && "p.adjust" %in% colnames(df_down)) {
    df_down$p.adjust <- as.numeric(as.character(df_down$p.adjust))
    df_down$Direction <- "Down"
  } else {
    df_down <- data.frame()
  }
  
  if(nrow(df_up) > 0) {
    df_up <- subset(df_up, p.adjust < 0.05)
  }
  if(nrow(df_down) > 0) {
    df_down <- subset(df_down, p.adjust < 0.05)
  }
  
  if(!is.null(top_n)) {
    if(nrow(df_up) > 0) {
      df_up <- df_up %>% arrange(desc(Count)) %>% head(top_n)
    }
    if(nrow(df_down) > 0) {
      df_down <- df_down %>% arrange(desc(Count)) %>% head(top_n)
    }
  }
  
  df_combined <- dplyr::bind_rows(df_up, df_down)
  df_combined$dummy <- 1
  df_combined$Direction <- factor(df_combined$Direction, levels = c("Up", "Down"))
  # Cambiar etiquetas a "upregulated" y "downregulated"
  levels(df_combined$Direction) <- c("upregulated", "downregulated")
  
  df_combined$Description <- str_wrap(df_combined$Description, width = 40)
  
  # Ordenar los GO terms por Count descendente (los de mayor Count aparecerán primero)
  df_order <- df_combined %>%
    dplyr::group_by(Description) %>%
    dplyr::summarise(max_Count = max(Count, na.rm = TRUE)) %>%
    dplyr::arrange(desc(max_Count))
  df_combined$Description <- factor(df_combined$Description, levels = rev(df_order$Description))
  
  p <- ggplot(df_combined, aes(x = dummy, y = Description, fill = p.adjust, size = Count)) +
    geom_point(shape = 21, color = "black", stroke = 0.5) +
    facet_grid(. ~ Direction, scales = "free_y", space = "free_y") +
    scale_fill_viridis_c(option = "magma") +
    # Aumentar el rango de tamaño para que las bolitas sean más grandes
    scale_size_continuous(range = c(4, 12)) +
    theme_classic() +
    labs(
      title = paste0(ontologia, " - ", comp_name,
                     ifelse(is.null(top_n), " (Todos los términos)", paste0(" (Top ", top_n, " terms)"))),
      x = "",
      y = "GO Term",
      fill = "p.adjust",
      size = "Gene Count"
    ) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      # Reducir el espacio entre paneles y márgenes para disminuir el espacio en blanco
      panel.spacing = unit(0.1, "lines"),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"),
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.text.y = element_text(size = 10),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 10)
    )
  
  # Aumentar la altura para que la gráfica sea más larga y se use mejor el espacio
  file_name <- paste0("dotplot_combined_", ontologia, "_", comp_name, 
                      ifelse(is.null(top_n), "_all", paste0("_top", top_n)), ".pdf")
  ggsave(file_name, p, width = 7, height = 14)
  
  return(p)
}

# Definir vector con las ontologías a procesar
ontologias <- c("BP", "CC", "MF")

# Iterar sobre cada comparación en 'resultados_enriquecimiento'
for (comp_name in names(resultados_enriquecimiento)) {
  cat("\nGenerando gráficos combinados para la comparación:", comp_name, "\n")
  sub_res <- resultados_enriquecimiento[[comp_name]]
  
  # Para cada ontología, generar gráfico combinando Up y Down en columnas separadas
  for (ont in ontologias) {
    cat("Procesando ontología:", ont, "\n")
    df_up   <- sub_res[[paste0(ont, "_up")]]
    df_down <- sub_res[[paste0(ont, "_down")]]
    
    # Generar y guardar el gráfico (top 15 términos)
    plot_combined_dotplot(
      df_up   = df_up,
      df_down = df_down,
      comp_name = comp_name,
      ontologia = ont,
      top_n = 15
    )
  }
}


# ----------------------------------------------------------------------------
# FIN DEL SCRIPT
# ----------------------------------------------------------------------------
