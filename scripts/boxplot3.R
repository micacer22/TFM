# Cargar librerías necesarias
library(ggplot2)
library(data.table)

# Función para extraer el porcentaje de la fila 12
extract_properly_paired_percentage <- function(file_path) {
  flagstat_data <- fread(file_path, header = FALSE, fill = TRUE)
  paired_row <- flagstat_data[12, V6]
  percentage_str <- sub(".*\\((\\d+\\.\\d+)%.*", "\\1", paired_row)
  if (grepl("^[0-9.]+$", percentage_str)) {
    as.numeric(percentage_str)
  } else {
    NA
  }
}

# Lista de archivos .flagstat.txt
files <- list.files(path = "/Users/miriamcaballerocervero/Desktop/Calidad_mapeo/Stats/Virus_stats", 
                    pattern = "_merged_TuMV.flagstat.txt", full.names = TRUE)

# Crear tabla de datos
mapping_data <- data.table(
  Library = gsub("(.*)_merged_TuMV.flagstat.txt", "\\1", basename(files)),
  Mapped_Percentage = sapply(files, extract_properly_paired_percentage)
)

# Asignar condición basado en el nombre
mapping_data$Condition <- ifelse(grepl("I", mapping_data$Library), "Infected", 
                                 ifelse(grepl("M", mapping_data$Library), "Mock", NA))

# Agrupar réplicas eliminando números finales
mapping_data$Genotype_Prefix <- gsub("[0-9]+$", "", mapping_data$Library)

# Definir la paleta de colores cbPalette
cbPalette <- c("Mock" = "#56B4E9", "Infected" = "#D55E00")  # Azul y rojo de la paleta cbPalette

# Crear el gráfico
myp1 <- ggplot(mapping_data, aes(x = Genotype_Prefix, y = Mapped_Percentage, fill = Condition)) +
  geom_boxplot(outlier.shape = 16, width = 0.5, color = "black") +  # Dibujar la caja combinada
  geom_jitter(position = position_jitter(width = 0.1, height = 0), size = 3, aes(color = Condition), alpha = 0.6) +  # Puntos individuales
  scale_fill_manual(values = cbPalette) +  # Usar la paleta cbPalette para las cajas
  scale_color_manual(values = cbPalette) +  # Usar la paleta cbPalette para los puntos
  theme_minimal() +
  ggtitle("TuMV Load by Genotype and Condition") +
  xlab("Genotype") +
  ylab("Virus Percentage Mapped to TuMV") +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  guides(fill = guide_legend(title = "Condition"), color = guide_legend(title = "Condition"))

# Guardar el gráfico
ggsave("TuMV_load_by_genotype_condition.pdf", plot = myp1, width = 10, height = 7)
