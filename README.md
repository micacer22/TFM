# **Integrating Transcriptomics and Genome Occupancy Data to Identify Epigenetically-Regulated Genes in Plant-Virus Interactions**  

## **Description**  
This Master's Thesis aims to identify epigenetically regulated genes in plant-virus interactions using Arabidopsis thaliana infected with Turnip mosaic virus (TuMV). This research explores how epigenetic modifications, such as histone marks, may influence plant responses to viral infections. Introduction: Plant-virus interactions are governed by complex mechanisms, with epigenetic regulation playing a crucial role. This study investigates how modifications to histone marks H3K4 and H3K36 affect gene expression in Arabidopsis during viral infection, using mutants with defects in histone-modifying proteins (sdg8 and atx1). These mutants have shown differences in symptom severity and viral load, suggesting a potential role for epigenetic regulation in antiviral defense. Objectives: The main objectives of this project are: 1) To analyze ChIP-seq data from infected and control epigenetic mutants, 2) To perform transcriptomic profiling to identify differentially expressed genes, 3) To select genes potentially regulated by epigenetic mechanisms, and 4) To verify the involvement of transposable elements in gene regulation. Methodology: The study follows a work plan that includes collecting samples from mutant and wild-type plants, both infected and uninfected. ChIP-seq analysis will be conducted to study histone modifications (H3K4 and H3K36), and RNA-seq will be used to evaluate gene expression. Bioinformatics tools will be applied to integrate these data sets to identify genes whose epigenetic control is associated with viral infection. Additionally, the role of transposable elements in gene regulation will be investigated.
## **Project Structure**  
📂 **scripts/** – Bash and R scripts for data analysis.  
📄 **README.md** – This document.  

## **Requirements**  
- Operating system: Linux or macOS (HPC with SLURM recommended).  
- Required bioinformatics tools:  
  - FastQC  
  - Trim Galore  
  - MultiQC  
  - Hisat2  
  - Samtools  
  - IGV  
  - featureCounts  

## **Usage Instructions**  

### **1. Download Data**  
Run the script to download sRNA sequencing data from GEO/SRA:  
```bash
bash scripts/Download_script_sRNAs.sh
```

### **2. Quality Control and Preprocessing**  
```bash
bash scripts/data_preprocessing.sh
```

### **3. Sequence Alignment**  
```bash
bash scripts/alignment_sRNAs.sh
```

### **4. Differential Expression Analysis**    
```bash
bash scripts/featureCounts_sRNAs.sh
```

### **5. Visualization in IGV**  
Load the BAM files into IGV to explore expression patterns.  

## **Expected Results**  
- Identification of differentially expressed sRNAs in epigenetic mutants.  
- Analysis of the relationship between TEs and virus-responsive genes.  
- Integration of methylation and expression data.  

## **Authors and Credits**  
This project was developed by **Miriam Caballero Cerveró**, in collaboration with the *Epigenetic Complexes in Plant Immunity* group at I2SysBio.  

## **Contact**  
📧 miriam.caballero@csic.es  
