# Instalar pheatmap si no está instalado
# install.packages("pheatmap")
#install.packages("BiocManager")
#BiocManager::install("EnhancedVolcano")

# Cargar librerías necesarias
library(DESeq2)
library(tidyverse)
library(conflicted)
library(gplots)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
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

# Crear objeto DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts_data, 
                              colData = metadata, 
                              design = ~ condition)

dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

# Obtener resultados
res <- results(dds)
res <- res[order(res$pvalue), ]
write.csv(as.data.frame(res), "DESeq2_results.csv")

# MA plot
pdf("MA_plot.pdf")
plotMA(res, main="MA Plot", ylim=c(-2,2))
dev.off()

# PCA Plot
rld <- rlog(dds, blind = TRUE)
pca_plot <- plotPCA(rld, intgroup = "condition") + theme_minimal()
ggsave("PCA_plot.pdf", pca_plot)

# Guardar datos normalizados
write.csv(assay(rld), "normalized_counts.csv")

# Heatmap de los 100 genes más variables
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
         annotation_col = metadata[, "condition", drop = FALSE], 
         main = "Heatmap of the 100 most variable genes")
dev.off()

# Volcano Plot
enhanced_volcano_plot <- EnhancedVolcano(res,
                                         lab = rownames(res),
                                         x = 'log2FoldChange',
                                         y = 'pvalue',
                                         title = 'Volcano Plot',
                                         pCutoff = 0.05,
                                         FCcutoff = 1.5)
ggsave("Volcano_plot.pdf", enhanced_volcano_plot)

# Boxplot de expresión normalizada
# Convertir rownames a columna para la unión
metadata <- metadata %>%
  rownames_to_column(var = "Sample")  

# Transformar datos en formato largo y unir con metadata
norm_counts_long <- as.data.frame(norm_counts) %>%
  rownames_to_column(var = "Gene") %>%
  pivot_longer(-Gene, names_to = "Sample", values_to = "Expression") %>%
  left_join(metadata, by = "Sample")

# Crear el boxplot
boxplot_expr <- ggplot(norm_counts_long, aes(x = condition, y = Expression, fill = condition)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Normalized expression boxplot")
ggsave("Boxplot_expression.pdf", boxplot_expr, width = 12, height = 8)  


# Dendrograma de clustering
pdf("Dendrogram.pdf")
dist_matrix <- dist(t(norm_counts))
hc <- hclust(dist_matrix, method = "complete")
plot(hc, main = "Clustering Dendrogram", xlab = "", sub = "")
dev.off()

# Seleccionar los 10 genes más significativos por p-valor
top_genes_counts <- rownames(res[order(res$pvalue), ])[1:10]

# Extraer datos de conteo normalizados solo para estos genes
counts_subset <- as.data.frame(t(norm_counts[top_genes_counts, ])) %>%
  rownames_to_column(var = "Sample") %>%
  pivot_longer(-Sample, names_to = "Gene", values_to = "Normalized_Counts") %>%
  left_join(metadata, by = "Sample")  

# Crear el Counts Plot
counts_plot <- ggplot(counts_subset, aes(x = condition, y = Normalized_Counts, fill = condition)) +
  geom_boxplot(alpha = 1.7) +
  geom_jitter(position = position_jitter(0.2), size = 1.5, alpha = 1.5) +
  facet_wrap(~ Gene, scales = "free_x", ncol=10) +  
  theme_minimal() +
  labs(title = "Normalized Counts for Top 10 Differentially Expressed Genes",
       x = "Condition",
       y = "Normalized Counts") +
  theme(legend.position = "none",
        panel.spacing = unit(2, "lines"), 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        strip.text = element_text(size = 10))  

ggsave("Counts_Plot.pdf", counts_plot, width = 40, height = 20)  



