# **INTEGRATING TRANSCRIPTOMICS AND GENOME OCCUPANCY DATA TO IDENTIFY EPIGENETICALLY-REGULATED GENES IN PLANT-VIRUS INTERACTIONS**  

## **Description**  
This project investigates the role of small RNAs (sRNAs) in the epigenetic regulation of virus-responsive genes in *Arabidopsis thaliana*. It analyzes the differential expression of sRNAs in epigenetic mutants and their relationship with nearby transposable elements (TEs).  

## **Project Structure**  
📂 **scripts/** – Bash scripts for data analysis.  
📄 **README.md** – This document.  

## **Requirements**  
- Operating system: Linux or macOS (HPC with SLURM recommended).  
- Required bioinformatics tools:  
  - FastQC  
  - Trim Galore  
  - MultiQC  
  - ShortStack  
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
