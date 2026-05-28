#=========================
# Bioinformatics R Project 
#=========================
install.packages("BiocManager")
BiocManager::install("pwalign")
# importing Library
library(Biostrings)   # Biological sequence analysis
library(tidyverse)    # Data manipulation and visualization
library(pwalign)      # Sequence alignment

#================================================================
# Problem 1:Transcriptomic Pipeline & Quality Control
#================================================================
## 1. Create the count matrix ( (Random Date with function rpois()) with 5,000 genes and 6 samples (3 Control + 3 Treated)).

set.seed(123)
counts <- matrix(rpois(30000,lambda=7) , nrow = 5000 ,ncol = 6)
rownames(counts) <- paste0("Gene",1:5000)    # named rows
colnames(counts) <- c("S1_Control","S2_Control","S3_Control", "S4_Treated" , "S5_Treated", "S6_Treated") # named columns


## 2. Noise filtration — Remove genes where the total count across all samples is less than 15.
keep_geens <- rowSums(counts) > 15 
filtered_data <- counts[keep_geens,]


## 3. Missing value recovery

#A- Add 5 random cells to become NA.)
filtered_data[sample(length(filtered_data) , 5) ] <- NA

#B- replace each using the mean of the other replicates in the same group.
# control group
control_sampl <- filtered_data[,1:3]

# treatment group
treatment_sampl <- filtered_data[,4:6]


for(i in 1:nrow(filtered_data)){
  
  # mean of control of each row
  control_mean <- mean(control_sampl[i,], na.rm =TRUE)
  
  #mean of treatment of eact row
  treatment_mean <- mean( treatment_sampl[i,], na.rm =TRUE)
  
  #fill Missing values in control data                     
  for(j in 1:ncol(control_sampl)){
    
    if(is.na(control_sampl[i,j])){
      control_sampl[i ,j] <- control_mean
    }           
  }                
  
  #fill Missing values in treatment data                     
  for(j in 1:ncol(treatment_sampl)){
    
    if(is.na(treatment_sampl[i,j])){
      treatment_sampl[i, j]<- treatment_mean
    }
  }
}

# edit in original data
filtered_data <- cbind(control_sampl, treatment_sampl)

## 4.CMP

#Create empty matrix for CPM
set.seed(456)
CPM <- matrix(0,nrow = nrow(filtered_data), ncol = ncol(filtered_data))

#Add column names
colnames(CPM) <- colnames(filtered_data)
rownames(CPM) <- rownames(filtered_data)

#Calculate CPM manually
for(i in 1:ncol(filtered_data)){
  
  # Calculate total reads for current sample
  sample_sum <- sum(filtered_data[, i])
  
  # Apply CPM formula
  CPM[,i] <- (filtered_data[,i] / sample_sum) * 1e6
}
print(CPM)

## 5.Verification plot

# Create two plots side by side
par(mfrow = c(1 ,2))

#Plot A: boxplot of raw counts after log2(counts + 1)
boxplot(log2(filtered_data +1), 
        main="count",
        col = "pink",
        las =2
)

#Plot B: boxplot of normalized CPM counts
boxplot(CPM,
        main="CPM",
        col="skyblue",
        las=2
)


# ==========================================
# Problem 2: Viral Variant & Sequence Integrity
# ==========================================

###### Viral Variant & Sequence Integrity ########

# Wild Type DNA sequence
wt_seq <- DNAString(
  "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
) # Store WT viral DNA sequence

# Variant DNA sequence
variant_seq <- DNAString(
  "ATGGCCATTGTAATGTAGCCGCTGAAAGGGTGCCCGATAG"
) # Store mutated viral DNA sequence

#######################Step 1: Global Alignment###############################

# Perform global alignment between WT and Variant
alignment <- pairwiseAlignment(
  wt_seq,
  variant_seq,
  type = "global"
)

# Print alignment result
alignment

# Calculate alignment score
score(alignment)

################################ Step 2: Transcription#################################

# Convert DNA into RNA by replacing T with U
rna_seq <- RNAString(
  chartr("T", "U", as.character(wt_seq))
)

# Display RNA sequence
rna_seq

############################# Step 3: Translation###############################

# Translate Variant DNA into amino acid sequence
protein_seq <- translate(variant_seq)

# Display translated protein
protein_seq

############## Step 4: GC Content#############

# Count G and C nucleotides
gc_count <- letterFrequency(wt_seq, c("G", "C"))

# Calculate total GC count
total_gc <- sum(gc_count)

# Get total sequence length
seq_length <- nchar(as.character(wt_seq))

# Calculate GC percentage
gc_percent <- (total_gc / seq_length) * 100

# Print GC content percentage
gc_percent

#########Step 5: Mutation Report##################

# Convert protein sequence to character
protein_char <- as.character(protein_seq)

# Check for internal stop codon
if(grepl("\\*", substr(protein_char, 1, nchar(protein_char)-1))) {
  
  # Print warning if premature stop codon exists
  print("Truncated Protein Detected")
  
} else {
  
  # Print if protein is complete
  print("Full Length Protein")
  
}


# ==========================================
# Problem 3: Proteomics Sparsity & Fold-Change Analysis
# ==========================================

# ------------------------------------------
# STEP 1: Create Proteomics Intensity Matrix
# ------------------------------------------

# Set random seed for reproducibility
set.seed(1)

# Create matrix with random protein intensity values
proteins <- matrix(
  
  sample(c(NA, 100:5000),
         800 * 6,
         replace = TRUE),
  
  nrow = 800,
  ncol = 6
  
)

# Assign sample names
colnames(proteins) <- c("A1", "A2", "A3",
                        "B1", "B2", "B3")

# Assign protein names
rownames(proteins) <- paste0("Protein_", 1:800)

# View first rows
head(proteins)

# ------------------------------------------
# STEP 2: Valid Value Filtering
# Keep proteins detected in at least
# 2 out of 3 replicates in one group
# ------------------------------------------

valid_rows <- apply(proteins, 1, function(x){
  
  # Count valid values in Group A
  groupA_valid <- sum(!is.na(x[1:3])) >= 2
  
  # Count valid values in Group B
  groupB_valid <- sum(!is.na(x[4:6])) >= 2
  
  # Keep protein if condition is true
  groupA_valid | groupB_valid
  
})

# Filter proteins
filtered_proteins <- proteins[valid_rows, ]

# Check dimensions after filtering
dim(filtered_proteins)

# ------------------------------------------
# STEP 3: Low-Abundance Imputation
# Replace NA with minimum detected value / 5
# ------------------------------------------

# Find minimum detected value
min_value <- min(filtered_proteins, na.rm = TRUE)

# Calculate replacement value
replacement_value <- min_value / 5

# Replace missing values
filtered_proteins[is.na(filtered_proteins)] <- replacement_value

# Check remaining NA values
sum(is.na(filtered_proteins))

# ------------------------------------------
# STEP 4: Z-score Scaling
# Scale each protein to mean = 0 and sd = 1
# ------------------------------------------

scaled_proteins <- t(
  
  apply(filtered_proteins, 1, function(x){
    
    (x - mean(x)) / sd(x)
    
  })
  
)

# View first scaled proteins
head(scaled_proteins)

# ------------------------------------------
# STEP 5: Fold Change Calculation
# Fold Change = Group B Mean / Group A Mean
# ------------------------------------------

# Calculate Group A mean
groupA_mean <- rowMeans(filtered_proteins[,1:3])

# Calculate Group B mean
groupB_mean <- rowMeans(filtered_proteins[,4:6])

# Calculate fold change
fold_change <- groupB_mean / groupA_mean

# View first fold changes
head(fold_change)

# ------------------------------------------
# STEP 6: Better Visualization
# Log2 Fold Change Plot
# ------------------------------------------

# Calculate Log2 Fold Change
log2FC <- log2(fold_change)

# Create dataframe
plot_data <- data.frame(
  
  Protein = rownames(filtered_proteins),
  
  Log2FC = log2FC
  
)

# Add regulation category
plot_data$Regulation <- ifelse(
  
  plot_data$Log2FC > 1,
  "Upregulated",
  
  ifelse(
    plot_data$Log2FC < -1,
    "Downregulated",
    "No Major Change"
  )
)

# Plot
ggplot(plot_data,
       
       aes(
         x = 1:nrow(plot_data),
         y = Log2FC,
         color = Regulation
       )) +
  
  geom_point(size = 2) +
  
  geom_hline(yintercept = c(-1,1),
             linetype = "dashed") +
  
  labs(
    
    title = "Protein Expression Changes",
    
    x = "Proteins",
    
    y = "Log2 Fold Change"
    
  ) +
  
  theme_minimal()