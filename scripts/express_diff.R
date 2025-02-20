# Instalar pheatmap
#install.packages("pheatmap") # Descomentar si no está instalado

# Cargar librerías necesarias
library(DESeq2)
library(tidyverse)
library(conflicted)
library(gplots)
library(RColorBrewer)
library(pheatmap)

# Cargar los datos de expresión
counts_data <- read.delim("/Users/miriamcaballerocervero/Desktop/counts_DESeq2_ready.txt",
                          header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

# Convertir a matriz
counts_data <- as.matrix(counts_data)

# Verificar estructura
print(dim(counts_data))
print(head(counts_data))

# Leer la metadata
metadata <- read.csv("/Users/miriamcaballerocervero/Desktop/metadata_conditions.csv", row.names = 1)

# Verificar que los nombres coincidan
print(colnames(counts_data))
print(rownames(metadata))

# Corregir nombres si hay espacios en blanco
colnames(counts_data) <- trimws(colnames(counts_data))
rownames(metadata) <- trimws(rownames(metadata))

# Verificar si hay diferencias
base::setdiff(colnames(counts_data), rownames(metadata))
base::setdiff(rownames(metadata), colnames(counts_data))

# Ordenar metadata según el orden de las columnas de counts
metadata <- metadata[colnames(counts_data), , drop = FALSE]

# Asegurar que están alineadas
if (!identical(colnames(counts_data), rownames(metadata))) {
  stop("Los nombres de las muestras en counts_data y metadata no están alineados.")
}

# Crear objeto DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts_data, 
                              colData = metadata, 
                              design = ~ condition)

print(dim(dds))  # Confirmar que el objeto se ha creado correctamente

# Filtrar genes con baja expresión
dds <- dds[rowSums(counts(dds)) >= 10, ]

# Realizar análisis de DESeq2
dds <- DESeq(dds)

# Obtener resultados
res <- results(dds)
res <- res[order(res$pvalue), ]  # Ordenar por p-valor

# Guardar resultados
write.csv(as.data.frame(res), "DESeq2_results.csv")

# Graficar MA plot
pdf("MA_plot.pdf")
plotMA(res, main="MA Plot", ylim=c(-2,2))
dev.off()

# PCA Plot
rld <- rlog(dds, blind = TRUE)
pca_plot <- plotPCA(rld, intgroup = "condition") + theme_minimal()
ggsave("PCA_plot.pdf", pca_plot)

# Guardar datos normalizados
write.csv(assay(rld), "normalized_counts.csv")

# Extraer datos normalizados
norm_counts <- assay(rld)

# Seleccionar los genes más variables
# Usamos la varianza de cada gen para seleccionar los que más varían entre las muestras
variances <- apply(norm_counts, 1, var)
top_genes <- order(variances, decreasing = TRUE)[1:100]  # Seleccionamos los 100 genes con mayor varianza

# Subconjunto de los datos normalizados para los genes más variables
norm_counts_top_genes <- norm_counts[top_genes, ]

# Crear el heatmap
pdf("Heatmap.pdf", height = 12)  # Para guardar el gráfico en formato PDF
pheatmap(norm_counts_top_genes, 
         scale = "row",  # Escala por fila para resaltar la variabilidad relativa
         clustering_distance_rows = "euclidean",  # Método de distancia para las filas
         clustering_distance_cols = "euclidean",  # Método de distancia para las columnas
         clustering_method = "complete",  # Método de agrupamiento
         annotation_col = metadata[, "condition", drop = FALSE],  # Colorear las columnas por condición
         main = "Heatmap of the 100 most variable genes")
dev.off()  # Cerrar el archivo PDF


