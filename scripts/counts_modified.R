# Cargar la tabla de conteos
counts <- read.table("/Users/miriamcaballerocervero/Desktop/counts_genes.txt", sep="\t", header=TRUE, comment.char="#")

# Eliminar columnas innecesarias
counts <- counts[, c(1, 7:ncol(counts))]  # Mantiene "Geneid" y los conteos

# Renombrar columnas eliminando sufijos
colnames(counts) <- gsub("_merged_Ath_R60.sorted.bam", "", colnames(counts))

# Guardar el archivo limpio
write.table(counts, "counts_DESeq2_ready.txt", sep="\t", quote=FALSE, row.names=FALSE)