#!/bin/bash
SBATCH --job-name=hisat2_alignment
SBATCH --ntasks=1
SBATCH --cpus-per-task=8
SBATCH --mem=32G
SBATCH --time=48:00:00
SBATCH --partition=global  
SBATCH --qos=medium                   

# Ruta de trabajo principal
wdir="/home/cacermi/miriam/2024/TFM/Data/X204SC24104429-Z01-F001/"
raw_dir="$wdir/01.RawData"
trimmed_dir="$wdir/1_trimmed"
qc_dir="$wdir/2_qc"
multiqc_dir="$wdir/3_multiqc"
hisat2_dir="$wdir/4_hisat2"
stats_dir="$wdir/5_stats"
virus_stats_dir="$wdir/6_virus_stats"

# Ruta al genoma de TuMV
tuvm_genome="/home/cacermi/miriam/2024/TFM/Data/X204SC24104429-Z01-F001/TuMV_Anc.fa"

# Directorios de salida
mkdir -p $raw_dir $trimmed_dir $qc_dir $multiqc_dir $hisat2_dir $stats_dir $virus_stats_dir

# Lista de bibliotecas
libs=("SI1_FKRN240341962-1A_22JCNWLT4_L1"
      "SI1_FKRN240341962-1A_22JCNWLT4_L1"
      "SI1_FKRN240341962-1A_22JKMCLT4_L8"
      "SI1_FKRN240341962-1A_22JKMCLT4_L8"
      "SI2_FKRN240341963-1A_22JCNWLT4_L2"
      "SI2_FKRN240341963-1A_22JCNWLT4_L2"
      "SI3_FKRN240341964-1A_22JCNWLT4_L2"
      "SI3_FKRN240341964-1A_22JCNWLT4_L2"
      "SI3_FKRN240341964-1A_22JKMCLT4_L8"
      "SI3_FKRN240341964-1A_22JKMCLT4_L8"
      "WM1_FKRN240341953-1A_22JCNWLT4_L7"
      "WM1_FKRN240341953-1A_22JKMCLT4_L8"
      "WM1_FKRN240341953-1A_22JCNWLT4_L7"
      "WM1_FKRN240341953-1A_22JKMCLT4_L8"
      "WM2_FKRN240341954-1A_22JCNWLT4_L4"
      "WM2_FKRN240341954-1A_22JCNWLT4_L4"
      "WM2_FKRN240341954-1A_22JKMCLT4_L8"
      "WM2_FKRN240341954-1A_22JKMCLT4_L8"
      "WM3_FKRN240341955-1A_22JCNWLT4_L8"
      "WM3_FKRN240341955-1A_22JCNWLT4_L8"
      "WM3_FKRN240341955-1A_22JKMCLT4_L8"
      "WM3_FKRN240341955-1A_22JKMCLT4_L8"
      "AI1_FKRN240341968-1A_22JCNWLT4_L8"
      "AI1_FKRN240341968-1A_22JCNWLT4_L8"
      "AI1_FKRN240341968-1A_22JCNWLT4_L8"
      "AI1_FKRN240341968-1A_22JCNWLT4_L8"
      "AI2_FKRN240341969-1A_22JCNWLT4_L3"
      "AI2_FKRN240341969-1A_22JCNWLT4_L3"
      "AI2_FKRN240341969-1A_22JKMCLT4_L8"
      "AI2_FKRN240341969-1A_22JKMCLT4_L8"
      "AI3_FKRN240341970-1A_22JCNWLT4_L2"
      "AI3_FKRN240341970-1A_22JCNWLT4_L2"
      "AI3_FKRN240341970-1A_22JKMCLT4_L8"
      "AI3_FKRN240341970-1A_22JKMCLT4_L8"
      "CI1_FKRN240341974-1A_22JCNWLT4_L6"
      "CI1_FKRN240341974-1A_22JCNWLT4_L6"
      "CI2_FKRN240341975-1A_22JCNWLT4_L1"
      "CI2_FKRN240341975-1A_22JCNWLT4_L1"
      "CI2_FKRN240341975-1A_22JKMCLT4_L8"
      "CI2_FKRN240341975-1A_22JKMCLT4_L8"
      "CI3_FKRN240341976-1A_22JCNWLT4_L1"
      "CI3_FKRN240341976-1A_22JCNWLT4_L1"
      "CI3_FKRN240341976-1A_22JKMCLT4_L8"
      "CI3_FKRN240341976-1A_22JKMCLT4_L8"
      "AM1_FKRN240341965-1A_22JCNWLT4_L3"
      "AM1_FKRN240341965-1A_22JCNWLT4_L3"
      "AM1_FKRN240341965-1A_22JKMCLT4_L8"
      "AM1_FKRN240341965-1A_22JKMCLT4_L8"
      "AM2_FKRN240341966-1A_22JCNWLT4_L5"
      "AM2_FKRN240341966-1A_22JCNWLT4_L5"
      "AM2_FKRN240341966-1A_22JKMCLT4_L8"
      "AM2_FKRN240341966-1A_22JKMCLT4_L8"
      "AM3_FKRN240341967-1A_22JCNWLT4_L5"
      "AM3_FKRN240341967-1A_22JCNWLT4_L5"
      "CM1_FKRN240341971-1A_22JCNWLT4_L3"
      "CM1_FKRN240341971-1A_22JCNWLT4_L3"
      "CM1_FKRN240341971-1A_22JKMCLT4_L8"
      "CM1_FKRN240341971-1A_22JKMCLT4_L8"
      "CM2_FKRN240341972-1A_22JCNWLT4_L6"
      "CM2_FKRN240341972-1A_22JCNWLT4_L6"
      "CM3_FKRN240341973-1A_22JCNWLT4_L4"
      "CM3_FKRN240341973-1A_22JCNWLT4_L4"
      "WI1_FKRN240341956-1A_22JCNWLT4_L4"
      "WI1_FKRN240341956-1A_22JCNWLT4_L4"
      "WI1_FKRN240341956-1A_22JKMCLT4_L8"
      "WI1_FKRN240341956-1A_22JKMCLT4_L8"
      "WI2_FKRN240341957-1A_22JCNWLT4_L5"
      "WI2_FKRN240341957-1A_22JCNWLT4_L5"
      "WI3_FKRN240341958-1A_22JCNWLT4_L1"
      "WI3_FKRN240341958-1A_22JCNWLT4_L1"
      "WI3_FKRN240341958-1A_22JKMCLT4_L8"
      "WI3_FKRN240341958-1A_22JKMCLT4_L8"
      "SM1_FKRN240341959-1A_22JCNWLT4_L4"
      "SM1_FKRN240341959-1A_22JCNWLT4_L4"
      "SM2_FKRN240341960-1A_22JCNWLT4_L5"
      "SM2_FKRN240341960-1A_22JCNWLT4_L5"
      "SM3_FKRN240341961-1A_22JCNWLT4_L4"
      "SM3_FKRN240341961-1A_22JCNWLT4_L4")

# 1: Verificación de calidad con FastQC
echo "Ejecutando FastQC..."
for lib in "${libs[@]}"; do
    # Determinar la carpeta 
    folder=$(echo "$lib" | cut -d'_' -f1)

    echo "Verificando calidad para $lib en carpeta $folder"
    fastqc -o $qc_dir "$raw_dir/$folder/${lib}_1.fq.gz" "$raw_dir/$folder/${lib}_2.fq.gz"
done


 2: Cortar secuencias con Trim Galore
echo "Ejecutando Trim Galore..."
for lib in "${libs[@]}"; do
    folder=$(echo "$lib" | cut -d'_' -f1)

    # Crear la carpeta salida de trimmed
    mkdir -p "$trimmed_dir/$folder"

    # Recortar las bases con una calidad menor a 20, eliminar primers y recortar los primeros 10 nt
    echo "Recortando secuencias para $lib en carpeta $folder con calidad menor a 20"
    trim_galore --paired \
        --quality 20 \
        --fastqc \
        --stringency 1 \
        --clip_R1 10 --clip_R2 10 \
        --output_dir "$trimmed_dir/$folder" \
        "$raw_dir/$folder/${lib}_1.fq.gz" "$raw_dir/$folder/${lib}_2.fq.gz"
done

# 3: Generación de informe MultiQC
echo "Ejecutando MultiQC..."
multiqc $trimmed_dir -o $multiqc_dir
echo "Resultados de MultiQC disponibles en: $multiqc_dir"

 4: Preparación de las referencias(Arabidopsis thaliana y TuMV)

# Descargar el genoma de Arabidopsis thaliana
echo "Descargando el genoma de Arabidopsis thaliana..."
wget -q -P $wdir https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-60/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.chromosome.1.fa.gz
wget -q -P $wdir https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-60/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.chromosome.2.fa.gz
wget -q -P $wdir https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-60/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.chromosome.3.fa.gz
wget -q -P $wdir https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-60/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.chromosome.4.fa.gz
wget -q -P $wdir https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-60/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.chromosome.5.fa.gz
wget -q -P $wdir https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-60/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.chromosome.Mt.fa.gz
wget -q -P $wdir https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-60/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.chromosome.Pt.fa.gz

# Descomprimir los archivos del genoma
echo "Descomprimiendo los archivos del genoma de Arabidopsis..."
gunzip $wdir*.fa.gz

# Concatenar los cromosomas de Arabidopsis thaliana
echo "Concatenando los cromosomas de Arabidopsis thaliana..."
if ls $wdir/Arabidopsis_thaliana.TAIR10.dna.chromosome.*.fa 1> /dev/null 2>&1; then
    cat $wdir/Arabidopsis_thaliana.TAIR10.dna.chromosome.*.fa > $wdir/Ath_R60.fa
else
    echo "Error: No se encontraron los cromosomas descargados."
    exit 1
fi


# 5: Indexación de los genomas con HISAT2
echo "Indexando el genoma de Arabidopsis thaliana..."
hisat2-build $wdir/Ath_R60.fa $hisat2_dir/Ath_R60

echo "Indexando el genoma de TuMV..."
hisat2-build $tuvm_genome $hisat2_dir/TuMV_genome

# 6: Mapeo de las lecturas con HISAT2

# Mapeo para Arabidopsis thaliana
echo "Mapeando las lecturas a Arabidopsis thaliana..."
for lib in "${libs[@]}"; do
    folder=$(echo "$lib" | cut -d'_' -f1)

    if [[ -f "$trimmed_dir/$folder/${lib}_1_val_1.fq.gz" && -f "$trimmed_dir/$folder/${lib}_2_val_2.fq.gz" ]]; then
        echo "Mapeando $lib desde carpeta $folder"
        
        # Ejecutar hisat2 con el número de hilos especificado
        hisat2 -x $hisat2_dir/Ath_R60 \
               -1 "$trimmed_dir/$folder/${lib}_1_val_1.fq.gz" \
               -2 "$trimmed_dir/$folder/${lib}_2_val_2.fq.gz" \
               -S "$hisat2_dir/${lib}_Ath_R60.sam" \
               -p 8

        # Convertir SAM a BAM con samtools
        echo "Convirtiendo archivo SAM a BAM..."
        samtools view -bS "$hisat2_dir/${lib}_Ath_R60.sam" > "$hisat2_dir/${lib}_Ath_R60.unsorted.bam"

        # Ordenar archivo BAM
        echo "Ordenando archivo BAM..."
        samtools sort "$hisat2_dir/${lib}_Ath_R60.unsorted.bam" -o "$hisat2_dir/${lib}_Ath_R60.sorted.bam"

        # Crear índice .bai para el archivo BAM
        echo "Creando índice BAI para el archivo BAM..."
        samtools index "$hisat2_dir/${lib}_Ath_R60.sorted.bam"
    fi
done

# Fusionar BAMs por muestra (agrupar réplicas técnicas) y calcular media de reads mapeadas
echo "Fusionando réplicas técnicas de Arabidopsis y calculando media de reads mapeadas..."
declare -A sample_bams
for bam in $hisat2_dir/*_Ath_R60.sorted.bam; do
    sample_name=$(basename "$bam" | cut -d'_' -f1)  # Extraer el identificador base
    sample_bams[$sample_name]="${sample_bams[$sample_name]} $bam"
done

for sample in "${!sample_bams[@]}"; do
    merged_bam="$hisat2_dir/${sample}_merged_Ath_R60.bam"
    echo "Fusionando BAMs para $sample..."
    samtools merge "$merged_bam" ${sample_bams[$sample]}
    samtools sort -o "$sorted_bam" "$merged_bam"
    samtools index "$sorted_bam"
    samtools flagstat "$sorted_bam" > "$hisat2_dir/${sample}_merged_Ath_R60.flagstat.txt"
    
    # Extraer porcentaje de reads mapeadas y calcular media
    mapped_reads=$(grep "properly paired" "$hisat2_dir/${sample}_merged_Ath_R60.flagstat.txt" | awk '{print $1}')
    total_reads=$(grep "in total" "$hisat2_dir/${sample}_merged_Ath_R60.flagstat.txt" | awk '{print $1}')
    percentage=$(echo "scale=2; ($mapped_reads/$total_reads)*100" | bc)
    echo "$sample $percentage" >> "$hisat2_dir/Ath_mapped_percentages.txt"
done

# Fusionar todos los BAMs de Arabidopsis en un solo archivo
echo "Fusionando todos los BAMs de Arabidopsis en un solo archivo..."
samtools merge "$hisat2_dir/All_Ath_R60.bam" $hisat2_dir/*_merged_Ath_R60.bam
samtools index "$hisat2_dir/All_Ath_R60.bam"
     

# Mapeo para el genoma del virus TuMV
echo "Mapeando las lecturas al genoma de TuMV..."
for lib in "${libs[@]}"; do
    folder=$(echo "$lib" | cut -d'_' -f1)

    if [[ -f "$trimmed_dir/$folder/${lib}_1_val_1.fq.gz" && -f "$trimmed_dir/$folder/${lib}_2_val_2.fq.gz" ]]; then
        echo "Mapeando $lib desde carpeta $folder"
        hisat2 -x $hisat2_dir/TuMV_genome \
               -1 "$trimmed_dir/$folder/${lib}_1_val_1.fq.gz" \
               -2 "$trimmed_dir/$folder/${lib}_2_val_2.fq.gz" \
               -S "$virus_stats_dir/${lib}_TuMV.sam" \
               -p 8 

        # Convertir SAM a BAM
        echo "Convirtiendo archivo SAM a BAM..."
        samtools view -bS "$virus_stats_dir/${lib}_TuMV.sam" > "$virus_stats_dir/${lib}_TuMV.unsorted.bam"

        # Ordenar archivo BAM
        echo "Ordenando archivo BAM..."
        samtools sort "$virus_stats_dir/${lib}_TuMV.unsorted.bam" -o "$virus_stats_dir/${lib}_TuMV.sorted.bam"

        # Crear índice .bai para el archivo BAM
        echo "Creando índice BAI para el archivo BAM..."
        samtools index "$virus_stats_dir/${lib}_TuMV.sorted.bam"
    fi
done

# Fusionar BAMs de TuMV por muestra y calcular media de reads mapeadas
echo "Fusionando réplicas técnicas de TuMV y calculando media de reads mapeadas..."
declare -A sample_bams_virus
for bam in $virus_stats_dir/*_TuMV.sorted.bam; do
    sample_name=$(basename "$bam" | cut -d'_' -f1)  # Extraer el identificador base
    sample_bams_virus[$sample_name]="${sample_bams_virus[$sample_name]} $bam"
done

for sample in "${!sample_bams_virus[@]}"; do
    merged_bam="$virus_stats_dir/${sample}_merged_TuMV.bam"
    echo "Fusionando BAMs para $sample..."
    samtools merge "$merged_bam" ${sample_bams_virus[$sample]}
    samtools index "$merged_bam"
    samtools flagstat "$merged_bam" > "$virus_stats_dir/${sample}_merged_TuMV.flagstat.txt"
    
    # Extraer porcentaje de reads mapeadas y calcular media
    mapped_reads=$(grep "properly paired" "$virus_stats_dir/${sample}_merged_TuMV.flagstat.txt" | awk '{print $1}')
    total_reads=$(grep "in total" "$virus_stats_dir/${sample}_merged_TuMV.flagstat.txt" | awk '{print $1}')
    percentage=$(echo "scale=2; ($mapped_reads/$total_reads)*100" | bc)
    echo "$sample $percentage" >> "$virus_stats_dir/TuMV_mapped_percentages.txt"#done

# Fusionar todos los BAMs de TuMV en un solo archivo
echo "Fusionando todos los BAMs de TuMV en un solo archivo..."
samtools merge "$virus_stats_dir/All_TuMV.bam" $virus_stats_dir/*_merged_TuMV.bam
samtools index "$virus_stats_dir/All_TuMV.bam"

# Generar gráfico con gnuplot
echo "Generando gráfico de porcentaje de virus en muestras..."
cat << EOF > plot_virus_load.gnu
set terminal pdfcairo enhanced font "Arial,12" size 7,5
set output "TuMV_virus_load_plot.pdf"
set title "TuMV Load"
set xlabel "Samples"
set ylabel "Virus percentage in RNAseq"
set xtics rotate by -45
set style data boxplot
set style fill solid
set boxwidth 0.5
plot "$virus_stats_dir/TuMV_mapped_percentages.txt" using 2:xtic(1) with boxes notitle
EOF
gnuplot plot_virus_load.gnu

echo "Gráfico guardado como TuMV_virus_load_plot.pdf"
echo "Análisis completado."

