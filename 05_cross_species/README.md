Here's a clear and well-structured GitHub README description for your script:  

---

# Cross-Species Cell Type Similarity Analysis  

## Overview  

This repository provides scripts for assessing the transcriptional similarity of cell clusters across species using **MetaNeighbor**. The analysis calculates **AUROC scores** to identify highly similar clusters between species. The key steps include data preprocessing, pseudo-cell generation, cross-species correlation analysis, and visualization.  

## Computational Principle  

We assessed the pairwise transcriptional similarity of cell clusters across species using **MetaNeighbor (v1.10.0)**. The AUROC score was calculated based on the expression patterns of **orthologous marker genes** between species.  

Key steps in the computation:  
- **Marker gene selection:** Orthologous genes identified as markers for each cluster (log2FC > 1.25, FDR < 0.05) were used for the analysis.  
- **Pseudo-cell generation:** A pseudo-cell was created by randomly selecting **ten cells** from a cluster (without replacement) and summing their UMI counts by gene.  
- **Normalization:** Marker gene expression was normalized as **CP10K** (UMI counts per 10,000) and standardized using **z-score transformation**.  
- **Similarity assessment:** Pearson’s correlation coefficients were used to compare transcriptional similarity between clusters from different species.  

## Pipeline Overview  

The analysis consists of three main steps:  

### **Step 1: Data Preprocessing (01tidySCE.R)**  
This script formats the dataset into **SingleCellExperiment (SCE)** format.  

**Input requirements:**  
1. A **gene (row) × cell (column)** raw UMI expression matrix.  
2. Metadata describing cell attributes.  

**Key processing steps:**  
- Generate pseudo-cell UMI matrices by summing the UMI counts of **ten randomly selected cells** per cluster.  
- Construct metadata for pseudo-cells.  
- Identify **orthologous genes** between species to create a common pseudo-cell UMI matrix.  
- Compute **marker genes** for each cell cluster (using z-score normalized pseudo-cell matrices).  

### **Step 2: AUROC Calculation (02auroc_cal.R)**  
This script computes **AUROC similarity scores** between cell clusters from different species.  

**Processing steps:**  
- Merge pseudo-cell UMI matrices from both species.  
- Perform **z-score normalization**.  
- Compute **AUROC scores** to quantify transcriptional similarity.  

**Dependencies:**  
This step requires two additional scripts for AUROC calculations:  
```r
source("2017-08-28-runMN-US.R")  
source("2017-08-28-runMN-US.pearson.R")  
```  

### **Step 3: Visualization (03heatmap_pheat.R)**  
This script organizes the AUROC scores and generates a heatmap for visualizing the transcriptional similarity between cell clusters across species.  

## Customization  

The scripts can be modified to fit different datasets. For example:  
- We included **gene symbol to gene ID conversion** for standardizing gene names across species. If not needed, these steps can be removed.  

## Usage  

Users can modify the scripts to accommodate their specific data and adjust parameters accordingly.  

---  