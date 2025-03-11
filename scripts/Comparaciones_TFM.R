# Instalar paquetes si no están instalados
#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2")

#install.packages("tidyverse")
#install.packages("pheatmap")
#install.packages("EnhancedVolcano")

# Cargar librerías necesarias
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(EnhancedVolcano)

# Cargar los datos de expresión
counts_data <- read.delim("/Users/miriamcaballerocervero/Desktop/counts_DESeq2_ready.txt", 
                          header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
counts_data <- as.matrix(counts_data)

# Leer la metadata
metadata <- read.csv("/Users/miriamcaballerocervero/Desktop/metadata_conditions.csv", row.names = 1)

# Asegurar alineación de nombres
colnames(counts_data) <- trimws(colnames(counts_data))
rownames(metadata) <- trimws(rownames(metadata))
metadata <- metadata[colnames(counts_data), , drop = FALSE]

if (!identical(colnames(counts_data), rownames(metadata))) {
  stop("Los nombres de las muestras en counts_data y metadata no están alineados.")
}

# 🔹 Filtrar solo las condiciones de interés
metadata_filtered <- metadata %>% 
  filter(condition %in% c("WT_mock", "WT_infected", "sdg8_mock", "sdg8_infected", "atx1_mock", "atx1_infected", "clf_mock", "clf_infected"))

# Filtrar counts_data para mantener solo estas muestras
counts_data_filtered <- counts_data[, colnames(counts_data) %in% rownames(metadata_filtered)]

# Asegurar alineación
metadata_filtered <- metadata_filtered[colnames(counts_data_filtered), , drop = FALSE]

if (!identical(colnames(counts_data_filtered), rownames(metadata_filtered))) {
  stop("Los nombres de las muestras en counts_data y metadata no están alineados.")
}

# 🔹 Crear objeto DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts_data_filtered, 
                              colData = metadata_filtered, 
                              design = ~ condition)

# Filtrar genes con baja expresión
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

# Función para realizar comparaciones y guardar resultados
comparar_condiciones <- function(grupo1, grupo2, nombre_archivo) {
  res <- results(dds, contrast = c("condition", grupo1, grupo2))
  res <- res[order(res$pvalue), ]
  write.csv(as.data.frame(res), nombre_archivo)
}

# 🔹 Realizar comparaciones y guardar resultados
comparar_condiciones("WT_mock", "WT_infected", "DESeq2_results_WT_mock_vs_WT_infected.csv")

comparar_condiciones("sdg8_mock", "sdg8_infected", "DESeq2_results_sdg8_mock_vs_sdg8_infected.csv")
comparar_condiciones("clf_mock", "clf_infected", "DESeq2_results_clf_mock_vs_clf_infected.csv")
comparar_condiciones("atx1_mock", "atx1_infected", "DESeq2_results_atx1_mock_vs_atx1_infected.csv")

comparar_condiciones("WT_mock", "sdg8_mock", "DESeq2_results_WT_mock_vs_sdg8_mock.csv")
comparar_condiciones("WT_mock", "clf_mock", "DESeq2_results_WT_mock_vs_clf_mock.csv")
comparar_condiciones("WT_mock", "atx1_mock", "DESeq2_results_WT_mock_vs_atx1_mock.csv")

comparar_condiciones("WT_infected", "sdg8_infected", "DESeq2_results_WT_infected_vs_sdg8_infected.csv")
comparar_condiciones("WT_infected", "clf_infected", "DESeq2_results_WT_infected_vs_clf_infected.csv")
comparar_condiciones("WT_infected", "atx1_infected", "DESeq2_results_WT_infected_vs_atx1_infected.csv")

# 🔹 Guardar datos normalizados
rld <- rlog(dds, blind = TRUE)
write.csv(assay(rld), "normalized_counts.csv")

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

crear_volcano_plot("DESeq2_results_WT_mock_vs_WT_infected.csv", "Volcano_WT_mock_vs_WT_infected.pdf", "WT Mock vs WT Infected")

crear_volcano_plot("DESeq2_results_sdg8_mock_vs_sdg8_infected.csv", "Volcano_sdg8_mock_vs_sdg8_infected.pdf", "sdg8 Mock vs sdg8 Infected")
crear_volcano_plot("DESeq2_results_atx1_mock_vs_atx1_infected.csv", "Volcano_atx1_mock_vs_atx1_infected.pdf", "atx1 Mock vs atx1 Infected")
crear_volcano_plot("DESeq2_results_clf_mock_vs_clf_infected.csv", "Volcano_clf_mock_vs_clf_infected.pdf", "clf Mock vs clf Infected")

crear_volcano_plot("DESeq2_results_WT_mock_vs_sdg8_mock.csv", "Volcano_WT_mock_vs_sdg8_mock.pdf", "WT Mock vs sdg8 Mock")
crear_volcano_plot("DESeq2_results_WT_mock_vs_atx1_mock.csv", "Volcano_WT_mock_vs_atx1_mock.pdf", "WT Mock vs atx1 Mock")
crear_volcano_plot("DESeq2_results_WT_mock_vs_clf_mock.csv", "Volcano_WT_mock_vs_clf_mock.pdf", "WT Mock vs clf Mock")

crear_volcano_plot("DESeq2_results_WT_infected_vs_sdg8_infected.csv", "Volcano_WT_infected_vs_sdg8_infected.pdf", "WT Infected vs sdg8 Infected")
crear_volcano_plot("DESeq2_results_WT_infected_vs_atx1_infected.csv", "Volcano_WT_infected_vs_atx1_infected.pdf", "WT Infected vs atx1 Infected")
crear_volcano_plot("DESeq2_results_WT_infected_vs_clf_infected.csv", "Volcano_WT_infected_vs_clf_infected.pdf", "WT Infected vs clf Infected")

# 🔹 Dendrograma de clustering
pdf("Dendrogram.pdf")
dist_matrix <- dist(t(norm_counts))
hc <- hclust(dist_matrix, method = "complete")
plot(hc, main = "Clustering Dendrogram", xlab = "", sub = "")
dev.off()

print("Análisis DESeq2 completado con éxito.")
