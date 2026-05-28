# Bioinformatics R Project

## 📌 Project Description

This project is a complete simulated bioinformatics pipeline implemented using R.  
It demonstrates key concepts in **transcriptomics, viral genomics, and proteomics data analysis**.

The project was built for educational purposes to practice real-world bioinformatics workflows including:

- Data generation and preprocessing
- Quality control and filtering
- Missing value handling
- Normalization techniques
- Sequence alignment and mutation analysis
- Protein expression analysis and visualization

---

# 🧬 Project Structure

The project is divided into three main biological modules:

---

# 🧪 Problem 1: Transcriptomic Pipeline & Quality Control

## Objective
To simulate RNA-seq count data and perform preprocessing and normalization.

## Steps performed:

### 1. Data Simulation
- A synthetic gene expression dataset was generated using Poisson distribution.
- Matrix size: **5000 genes × 6 samples**
- Samples:
  - 3 Control samples
  - 3 Treated samples

### 2. Noise Filtering
- Genes with total counts less than 15 across all samples were removed.
- This step ensures low-expression genes do not affect downstream analysis.

### 3. Missing Value Handling
- Random NA values were introduced to simulate real biological noise.
- Missing values were replaced using **group-wise mean imputation**:
  - Control group (S1–S3)
  - Treated group (S4–S6)

### 4. Normalization (CPM)
- Counts Per Million (CPM) normalization was applied:
  
  CPM = (gene count / total sample counts) × 1e6

- This allows comparison between samples with different sequencing depths.

### 5. Visualization
- Boxplots were generated to compare:
  - Raw log2(counts + 1)
  - Normalized CPM values

---

# 🧬 Problem 2: Viral Variant & Sequence Integrity

## Objective
To analyze DNA mutations and their biological impact.

## Steps performed:

### 1. Sequence Alignment
- Global alignment was performed between:
  - Wild Type DNA sequence
  - Variant DNA sequence
- Tool used: `pwalign`

### 2. Transcription
- DNA sequence was converted into RNA by replacing:
  
  T → U

### 3. Translation
- RNA sequence was translated into amino acid sequence using genetic code.

### 4. GC Content Analysis
- GC content was calculated using:

  GC% = (G + C) / total nucleotides × 100

- GC content is important for understanding DNA stability.

### 5. Mutation Analysis
- Protein sequence was checked for:
  - Premature stop codons (*)
- Output:
  - “Truncated Protein Detected” if mutation affects protein length
  - “Full Length Protein” otherwise

---

# 🧫 Problem 3: Proteomics Sparsity & Fold-Change Analysis

## Objective
To simulate protein expression data and analyze differential expression.

## Steps performed:

### 1. Data Simulation
- Protein intensity matrix generated:
  - 800 proteins × 6 samples
- Values include random intensities and missing values (NA)

### 2. Filtering Strategy
- Proteins retained only if:
  - At least 2 valid values in either:
    - Group A (A1–A3)
    - Group B (B1–B3)

### 3. Missing Value Imputation
- Missing values replaced with:

  minimum detected value / 5

- This avoids bias from high-value imputation.

### 4. Z-score Normalization
- Each protein was standardized:

  Z = (x - mean) / standard deviation

- Ensures comparability across proteins.

### 5. Fold Change Analysis
- Fold change calculated:

  Fold Change = Mean(Group B) / Mean(Group A)

- Log2 transformation applied for interpretation.

### 6. Visualization
- Scatter plot of Log2 Fold Change:
  - Upregulated proteins (Log2FC > 1)
  - Downregulated proteins (Log2FC < -1)
  - No major change proteins

---

# 🧰 Libraries Used

- Biostrings → DNA/RNA sequence handling
- pwalign → sequence alignment
- tidyverse → data manipulation and visualization

---

# 🎯 Key Skills Demonstrated

✔ Bioinformatics data simulation  
✔ RNA-seq preprocessing pipeline  
✔ CPM normalization  
✔ DNA alignment and mutation analysis  
✔ Protein expression analysis  
✔ Statistical scaling (Z-score normalization)  
✔ Fold-change and log2 analysis  
✔ Data visualization using ggplot2  

---

# 🚀 Project Goal

This project was built to simulate real bioinformatics workflows and strengthen understanding of:

- Omics data analysis
- Molecular sequence processing
- Differential expression analysis
- R programming for biological data science

---

# 📌 Author

Created as part of a Bioinformatics learning journey using R.

---

⭐ If you find this project useful, feel free to star the repository!
