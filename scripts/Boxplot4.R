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

# Extraer el prefijo genotípico (letra más condición)
mapping_data$Genotype_Prefix <- gsub("[0-9]+$", "", mapping_data$Library)

# Asignar condición basado en el prefijo
mapping_data$Condition <- ifelse(grepl("I$", mapping_data$Genotype_Prefix), "Infected",
                                 ifelse(grepl("M$", mapping_data$Genotype_Prefix), "Mock", NA))

# Convertir 'Condition' a factor con niveles explícitos
mapping_data$Condition <- factor(mapping_data$Condition, levels = c("Mock", "Infected"))

# Definir la paleta de colores
cbPalette <- c("Mock" = "#56B4E9", "Infected" = "#D55E00")  # Azul y rojo de la paleta

# Crear el gráfico principal
myp1 <- ggplot(mapping_data, aes(x = Genotype_Prefix, y = Mapped_Percentage, fill = Condition)) +
  geom_boxplot(outlier.shape = 16, width = 0.5, color = "black") +
  geom_jitter(position = position_jitter(width = 0.1, height = 0), size = 3, aes(color = Condition), alpha = 0.6) +
  scale_fill_manual(values = cbPalette) +
  scale_color_manual(values = cbPalette) +
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

# Guardar el gráfico principal
ggsave("TuMV_load_by_genotype_condition.pdf", plot = myp1, width = 10, height = 7)

# Filtrar datos para las librerías que empiezan con W, A y S
data_WAS <- mapping_data[grepl("^(W|A|S)", Genotype_Prefix)]

# Crear el gráfico para W, A y S
plot_WAS <- ggplot(data_WAS, aes(x = Genotype_Prefix, y = Mapped_Percentage, fill = Condition)) +
  geom_boxplot(outlier.shape = 16, width = 0.5, color = "black") +
  geom_jitter(position = position_jitter(width = 0.1, height = 0), size = 3, aes(color = Condition), alpha = 0.6) +
  scale_fill_manual(values = cbPalette) +
  scale_color_manual(values = cbPalette) +
  theme_minimal() +
  ggtitle("TuMV Load by Genotype (WT/atx1/sdg8)") +
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

# Guardar el gráfico para W, A y S
ggsave("TuMV_load_WAS.pdf", plot = plot_WAS, width = 10, height = 7)

# Filtrar datos para las librerías que empiezan con W y C
data_WC <- mapping_data[grepl("^(W|C)", Genotype_Prefix)]

# Crear el gráfico para W y C
plot_WC <- ggplot(data_WC, aes(x = Genotype_Prefix, y = Mapped_Percentage, fill = Condition)) +
  geom_boxplot(outlier.shape = 16, width = 0.5, color = "black") +
  geom_jitter(position = position_jitter(width = 0.1, height = 0), size = 3, aes(color = Condition), alpha = 0.6) +
  scale_fill_manual(values = cbPalette) +
  scale_color_manual(values = cbPalette) +
  theme_minimal() +
  ggtitle("TuMV Load by Genotype (WT/clf)") +
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

# Guardar el gráfico para W y C
ggsave("TuMV_load_WC.pdf", plot = plot_WC, width = 10, height = 7)
