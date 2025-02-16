#!/bin/bash
# Lista de archivos BAM
bam_files=(
    AI1_merged_Ath_R60.sorted.bam AI2_merged_Ath_R60.sorted.bam AI3_merged_Ath_R60.sorted.bam
    AM1_merged_Ath_R60.sorted.bam AM2_merged_Ath_R60.sorted.bam AM3_merged_Ath_R60.sorted.bam
    CI1_merged_Ath_R60.sorted.bam CI2_merged_Ath_R60.sorted.bam CI3_merged_Ath_R60.sorted.bam
    CM1_merged_Ath_R60.sorted.bam CM2_merged_Ath_R60.sorted.bam CM3_merged_Ath_R60.sorted.bam
    SI1_merged_Ath_R60.sorted.bam SI2_merged_Ath_R60.sorted.bam SI3_merged_Ath_R60.sorted.bam
    SM1_merged_Ath_R60.sorted.bam SM2_merged_Ath_R60.sorted.bam SM3_merged_Ath_R60.sorted.bam
    WI1_merged_Ath_R60.sorted.bam WI2_merged_Ath_R60.sorted.bam WI3_merged_Ath_R60.sorted.bam
    WM1_merged_Ath_R60.sorted.bam WM2_merged_Ath_R60.sorted.bam WM3_merged_Ath_R60.sorted.bam
)

# Archivo de anotaciones
annotation="/home/cacermi/miriam/2024/InSilico/proyecto_epigenetica/Data/SRR/combined_output/Annotations/atRTD3_TS_21Feb22_transfix.gtf"

# Lista de genes de interés
genes=("AT2G31650" "AT1G77300" "AT2G23380") # ATX1, SDG8, CLF

# Archivo de salida
output_file="counts_genes_filtered.txt"

featureCounts -a "$annotation" \
              -o counts_genes.txt \
              -T 8 \
              -s 0 \
              -p \
              --countReadPairs \
              -t exon \
              -g gene_id \
              "${bam_files[@]}"


# Filtrar los genes de interés y guardar en un nuevo archivo
grep -E "$(IFS="|"; echo "${genes[*]}")" counts_genes.txt > "$output_file"

echo "Proceso completado. Resultados filtrados en $output_file"
