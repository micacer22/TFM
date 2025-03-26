# Instalar paquetes si no están instalados
#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2")

#install.packages("tidyverse")
#install.packages("pheatmap")
#install.packages("EnhancedVolcano")
setwd("/Users/miriamcaballerocervero/Documents/Miriam/24-25/TFM/Gráficas_Expr_Dif/Comparaciones")

# Cargar librerías necesarias
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(EnhancedVolcano)

# Cargar los datos de expresión
counts_data <- read.delim("/Users/miriamcaballerocervero/Documents/Miriam/24-25/TFM/Gráficas_Expr_Dif/counts_DESeq2_ready.txt", 
                          header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
counts_data <- as.matrix(counts_data)

# Leer la metadata
metadata <- read.csv("/Users/miriamcaballerocervero/Documents/Miriam/24-25/TFM/Gráficas_Expr_Dif/metadata_conditions.csv", row.names = 1)

# Asegurar alineación de nombres
colnames(counts_data) <- trimws(colnames(counts_data))
rownames(metadata) <- trimws(rownames(metadata))
metadata <- metadata[colnames(counts_data), , drop = FALSE]

if (!identical(colnames(counts_data), rownames(metadata))) {
  stop("Los nombres de las muestras en counts_data y metadata no están alineados.")
}

# Filtrar solo las condiciones de interés
condiciones_interes <- c("WT_mock", "WT_infected", "sdg8_mock", "sdg8_infected", 
                         "atx1_mock", "atx1_infected", "clf_mock", "clf_infected")

metadata_filtered <- metadata %>% filter(condition %in% condiciones_interes)

# Filtrar counts_data para mantener solo estas muestras
counts_data_filtered <- counts_data[, colnames(counts_data) %in% rownames(metadata_filtered)]

# Asegurar alineación nuevamente
metadata_filtered <- metadata_filtered[colnames(counts_data_filtered), , drop = FALSE]

if (!identical(colnames(counts_data_filtered), rownames(metadata_filtered))) {
  stop("Los nombres de las muestras en counts_data y metadata no están alineados.")
}

# Crear objeto DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts_data_filtered, 
                              colData = metadata_filtered, 
                              design = ~ condition)

# Filtrar genes con baja expresión
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

# Función para realizar comparaciones, guardar resultados (log2FC = mytreatment / mycontrol)
comparar_condiciones <- function(mycontrol, mytreatment, nombre_archivo) {
  # Definir contraste: log2FoldChange =  mytreatment / mycontrol
  res <- results(dds, contrast = c("condition", mytreatment, mycontrol))
  res <- res[order(res$pvalue), ]
  
  write.csv(as.data.frame(res), nombre_archivo)
}

# La lista de comparaciones se interpreta así:
# c("WT_mock", "WT_infected", "DESeq2_results_WT_responsive.csv")
# significa log2FoldChange = WT_infected / WT_mock
comparaciones <- list(
  c("WT_mock",     "WT_infected",  "DESeq2_results_WT_responsive.csv"),
  c("sdg8_mock",   "sdg8_infected","DESeq2_results_sdg8_responsive.csv"),
  c("clf_mock",    "clf_infected", "DESeq2_results_clf_responsive.csv"),
  c("atx1_mock",   "atx1_infected","DESeq2_results_atx1_responsive.csv"),
  c("WT_mock",     "sdg8_mock",    "DESeq2_results_sdg8_mock_primed.csv"),
  c("WT_mock",     "clf_mock",     "DESeq2_results_clf_mock_primed.csv"),
  c("WT_mock",     "atx1_mock",    "DESeq2_results_atx1_mock_primed.csv"),
  c("WT_infected", "sdg8_infected","DESeq2_results_sdg8_infected_primed.csv"),
  c("WT_infected", "clf_infected", "DESeq2_results_clf_infected_primed.csv"),
  c("WT_infected", "atx1_infected","DESeq2_results_atx1_infected_primed.csv")
)

# Ejecutar el bucle para generar los CSV
for (comp in comparaciones) {
  comparar_condiciones(comp[1], comp[2], comp[3])
}

# 🔹 Guardar datos normalizados
rld <- rlog(dds, blind = TRUE)
write.csv(assay(rld), "normalized_counts.csv")

# --- Después de haber corrido el análisis DESeq2 y definido la lista de comparaciones ---

# Crear data frame vacío para almacenar los conteos
df_counts <- data.frame(Comparacion = character(), 
                        Regulation = character(), 
                        Count = integer(), 
                        stringsAsFactors = FALSE)

# Iterar por cada comparación y contar los genes
for (comp in comparaciones) {
  file_name <- comp[3]
  # Extraer el nombre de la comparación usando anclas para el inicio y final
  condition_name <- gsub("^DESeq2_results_", "", file_name)
  condition_name <- gsub("\\.csv$", "", condition_name)
  condition_name <- gsub("_primed", " primed", condition_name)
  
  print(paste("Comparando:", condition_name))
  
  # Obtener resultados para la comparación 
  res <- results(dds, contrast = c("condition", comp[1], comp[2]))
  res <- res[order(res$pvalue), ]
  
  # Usar padj <= 0.05 para definir significancia
  padj_cutoff <- 0.05
  up_count <- sum(res$padj <= padj_cutoff & res$log2FoldChange > 0, na.rm = TRUE)
  down_count <- sum(res$padj <= padj_cutoff & res$log2FoldChange < 0, na.rm = TRUE)
  
  # Agregar los conteos al data frame
  df_counts <- rbind(df_counts,
                     data.frame(Comparacion = condition_name, Regulation = "Up", Count = up_count, stringsAsFactors = FALSE),
                     data.frame(Comparacion = condition_name, Regulation = "Down", Count = down_count, stringsAsFactors = FALSE))
}

# Verificar el data frame generado
print(df_counts)

# --- Graficar los conteos con ggplot2 en barras apiladas ---
# Reordena para que Down sea el primer nivel (abajo) y Up el segundo (arriba)
df_counts$Regulation <- factor(df_counts$Regulation, levels = c("Down", "Up"))

p <- ggplot(df_counts, aes(x = Comparacion, y = Count, fill = Regulation)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = c("Up" = "#E69F00", "Down" = "#999999")) +
  theme_minimal() +
  labs(title = "Number of Differentially Expressed Genes Per Condition",
       x = "Condition",
       y = "Number of Differentially Expressed Genes",
       fill = "Regulation") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p)

# Guardar la gráfica en diferentes formatos
ggsave("Grafica_DEGs.png", plot = p, width = 8, height = 6, dpi = 300)  # Guardar como PNG
ggsave("Grafica_DEGs.pdf", plot = p, width = 8, height = 6)  # Guardar como PDF

# 🔹 PCA Plot
pca_plot <- plotPCA(rld, intgroup = "condition") + theme_minimal()
ggsave("PCA_plot.pdf", pca_plot)

# 🔹 Heatmap de los 100 genes más variables
norm_counts <- assay(rld)
variances <- apply(norm_counts, 1, var)
top_genes <- order(variances, decreasing = TRUE)[1:100]
norm_counts_top_genes <- norm_counts[top_genes, ]

pdf("Heatmap.pdf", height = 12)
pheatmap(norm_counts_top_genes, 
         scale = "row", 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete", 
         annotation_col = metadata_filtered[, "condition", drop = FALSE], 
         main = "Heatmap of the 100 most variable genes")
dev.off()

# 🔹 Volcano Plot para cada comparación
crear_volcano_plot <- function(archivo_csv, nombre_pdf, titulo) {
  res <- read.csv(archivo_csv, row.names = 1)
  enhanced_volcano_plot <- EnhancedVolcano(res,
                                           lab = rownames(res),
                                           x = 'log2FoldChange',
                                           y = 'pvalue',
                                           title = titulo,
                                           pCutoff = 0.05,
                                           FCcutoff = 1.5)
  ggsave(nombre_pdf, enhanced_volcano_plot)
}

crear_volcano_plot("DESeq2_results_WT_responsive.csv", "Volcano_WT_mock_vs_WT_infected.pdf", "WT Mock vs WT Infected")
crear_volcano_plot("DESeq2_results_sdg8_responsive.csv", "Volcano_sdg8_mock_vs_sdg8_infected.pdf", "sdg8 Mock vs sdg8 Infected")
crear_volcano_plot("DESeq2_results_atx1_responsive.csv", "Volcano_atx1_mock_vs_atx1_infected.pdf", "atx1 Mock vs atx1 Infected")
crear_volcano_plot("DESeq2_results_clf_responsive.csv", "Volcano_clf_mock_vs_clf_infected.pdf", "clf Mock vs clf Infected")

crear_volcano_plot("DESeq2_results_sdg8_mock_primed.csv", "Volcano_WT_mock_vs_sdg8_mock.pdf", "WT Mock vs sdg8 Mock")
crear_volcano_plot("DESeq2_results_clf_mock_primed.csv", "Volcano_WT_mock_vs_clf_mock.pdf", "WT Mock vs clf Mock")
crear_volcano_plot("DESeq2_results_atx1_mock_primed.csv", "Volcano_WT_mock_vs_atx1_mock.pdf", "WT Mock vs atx1 Mock")

crear_volcano_plot("DESeq2_results_sdg8_infected_primed.csv", "Volcano_WT_infected_vs_sdg8_infected.pdf", "WT Infected vs sdg8 Infected")
crear_volcano_plot("DESeq2_results_clf_infected_primed.csv", "Volcano_WT_infected_vs_clf_infected.pdf", "WT Infected vs clf Infected")
crear_volcano_plot("DESeq2_results_atx1_infected_primed.csv", "Volcano_WT_infected_vs_atx1_infected.pdf", "WT Infected vs atx1 Infected")

# 🔹 Dendrograma de clustering
pdf("Dendrogram.pdf")
dist_matrix <- dist(t(norm_counts))
hc <- hclust(dist_matrix, method = "complete")
plot(hc, main = "Clustering Dendrogram", xlab = "", sub = "")
dev.off()

print("Análisis DESeq2 completado con éxito.")
