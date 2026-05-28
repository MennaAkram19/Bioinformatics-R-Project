# Bioinformatics R Project

This project demonstrates a complete bioinformatics workflow implemented in R, covering transcriptomic, viral genomics, and proteomics data analysis.

## 📊 Project Overview

The project is divided into three main parts:

### 1. Transcriptomic Pipeline & Quality Control
- Generation of simulated RNA-seq count data (5,000 genes × 6 samples)
- Noise filtering based on low-expression genes
- Missing value imputation using group-wise mean
- Normalization using CPM (Counts Per Million)
- Visualization of raw vs normalized data

### 2. Viral Variant & Sequence Integrity
- DNA sequence alignment (Wild Type vs Variant) using global alignment
- Transcription (DNA → RNA conversion)
- Translation (RNA → protein sequence)
- GC content calculation
- Mutation and stop codon detection

### 3. Proteomics Sparsity & Fold-Change Analysis
- Simulated protein intensity matrix (800 proteins × 6 samples)
- Quality filtering based on replicate detection
- Missing value imputation using minimum value strategy
- Z-score normalization
- Fold-change and Log2 fold-change analysis
- Visualization of protein expression changes

## 🧪 Tools & Libraries Used
- Biostrings
- pwalign
- tidyverse

## 📈 Key Features
- End-to-end omics data simulation and analysis
- Data cleaning and normalization techniques
- Sequence alignment and mutation analysis
- Protein expression analysis with visualization

## 🚀 Goal
This project was created to practice bioinformatics data analysis workflows using R and simulate real-world omics pipelines.

---

⭐ If you like this project, feel free to star the repository!
