# Glioblastoma Multi-Cohort Analysis Pipeline

## Overview
This pipeline performs differential gene expression analysis across three glioblastoma datasets (TCGA, CGGA, GLASS) to identify survival-associated genes. The analysis uses limma for differential expression and includes batch correction, quality control, and cross-dataset validation.

---

## Dataset Processing

### TCGA Dataset (`load_TCGA`)

**Input Files:**
- Expression data: `GBMLGG_EB_RmDiffFullGenesRanRmDup.csv`
- Clinical data: `Compiled data_TCGA_CGGA_GLASS_jh080424.xlsx` (sheet: 'TCGA')

**Processing Steps:**
1. **Platform filtering**: Retains only RNAseq and Agilent microarray platforms
2. **Sample filtering**: Keeps only samples present in both expression and clinical data
3. **Gene filtering**: Removes genes with all zeros or >50% zero values (`remove_zero_genes`)
4. **Normalization**: StandardScaler applied **separately** to each platform (RNAseq and Agilent) to reduce batch effects
5. **Quality control**: PCA plots generated before and after normalization to visualize platform effects

**Key Variables:**
- `Overall.Survival.Months`: Survival time
- `Diagnosis.Age`: Age at diagnosis
- `Sex`: Patient sex
- `case.vital.status`: Vital status (alive/dead)

**Platform Handling:**
- Critical: RNAseq and Agilent data are normalized separately before merging to account for platform-specific technical variation
- This reduces the separation between platforms visible in PCA

---

### CGGA Dataset (`load_CGGA`)

**Input Files:**
- Clinical 325: `CGGA.mRNAseq_325_clinical.20200506.txt`
- Clinical 693: `CGGA.mRNAseq_693_clinical.20200506.txt`
- Expression 325: `CGGA.mRNAseq_325.RSEM-genes.20200506.txt`
- Expression 693: `CGGA.mRNAseq_693.RSEM-genes.20200506.txt`

**Processing Steps:**
1. **Clinical filtering**: 
   - PRS_type == 'Primary'
   - Histology == 'GBM'
   - IDH_mutation_status == 'Wildtype'
2. **Survival conversion**: Days to months (divided by 30.44)
3. **Gene intersection**: Only genes present in both batches (325 and 693) are retained
4. **Gene filtering**: Removes genes with all zeros or >50% zero values
5. **Transformation**: Log1p transformation applied to expression values
6. **Outlier removal**: PCA-based outlier detection (samples with |PC1| > 100 or |PC2| > 100 removed)
7. **Batch normalization**: StandardScaler applied **separately** to batch 325 and batch 693
8. **Batch variable**: 'batch' column added ('batch325' or 'batch693') for downstream correction

**Key Variables:**
- `Overall.Survival.Months`: Converted from days
- `Diagnosis.Age`: Age at diagnosis
- `Sex` (Gender): Patient sex
- `case.vital.status`: 0=alive, 1=dead
- `batch`: Batch identifier (325 or 693)

**Important Notes:**
- Two separate sequencing batches are processed independently then combined
- Outlier removal based on extreme PCA values is performed before final normalization

---

### GLASS Dataset (`load_GLASS`)

**Input Files:**
- Expression: `gene_tpm_matrix_all_samples.tsv`
- Clinical: `Compiled data_TCGA_CGGA_GLASS_jh080424.xlsx` (sheet: 'GLASS')

**Processing Steps:**
1. **Sample name transformation**: Custom function `transform_sample_name_GLASS` standardizes sample IDs
2. **Gene filtering**: Removes genes with all zeros or >50% zero values
3. **Sample filtering**: Keeps only samples present in both expression and clinical data
4. **Transformation**: Log1p transformation with small constant (1e-5)
5. **Normalization**: StandardScaler applied to all expression data
6. **Deduplication**: Duplicate samples removed based on Sample.ID

**Key Variables:**
- `Overall.Survival.Months`: Survival time
- `Diagnosis.Age`: Age at diagnosis
- `Sex`: Patient sex
- `case.vital.status`: Vital status

**Important Notes:**
- Single batch dataset, simpler processing than TCGA/CGGA
- Expression values are TPM (transcripts per million)

---

## Differential Expression Analysis

### Limma Analysis (`limma` function)

**Method**: Linear modeling with empirical Bayes moderation using the R limma package

**Survival Group Definition:**

**Fixed Threshold Mode** (e.g., 24, 36, 48 months):
- **Low group**: Patients who died with OS ≤ threshold
- **High group**: Patients with OS > threshold (any vital status)
- Exclusion: Patients alive with OS ≤ threshold (censored before threshold)

**Range Threshold Mode** (e.g., 24-48, 24-60 months):
- **Low group**: Patients who died with OS ≤ low_threshold
- **High group**: Patients with OS > high_threshold (any vital status)
- Exclusion: Patients alive with OS ≤ low_threshold AND patients with low_threshold < OS ≤ high_threshold

**Design Matrix:**

For **CGGA** (has batch variable):
```R
~ group + Diagnosis.Age + Sex + batch
```

For **TCGA and GLASS** (no batch variable):
```R
~ group + Diagnosis.Age + Sex
```

**Covariates:**
- `group`: Survival group (low vs high) - **primary variable of interest**
- `Diagnosis.Age`: Continuous covariate
- `Sex`: Categorical covariate
- `batch`: Batch identifier (CGGA only)

**Statistical Testing:**
- Coefficient tested: `grouphigh` (high survival vs low survival)
- Multiple testing correction: FDR (Benjamini-Hochberg)
- Significance threshold: FDR < 0.1

**Output Files:**
- `all_genes_[consortium]_[threshold].csv`: All genes with statistics
- `significant_genes_[consortium]_[threshold].csv`: FDR < 0.1 genes only

---

## Key Analysis Functions

### `define_survival_groups`
Creates binary survival groups (Low/High) based on a single threshold. Used for limma analysis with fixed thresholds.

### `define_survival_groups_all`
Creates four survival categories for visualization:
- 0-24 months
- >24 and ≤48 months
- >48 and ≤60 months
- >60 months

### `pca_plot_glass`
Generates PCA plots colored by survival groups. Uses 2 components to visualize gene expression patterns.

### `boxplots_glass`
Creates 2×3 grids of boxplots comparing gene expression across survival timepoints. Includes stripplots for individual data points.

### `remove_zero_genes`
Quality control function that removes:
- Genes with sum = 0 across all samples
- Genes with >50% zero values

### `compute_log_fold_change`
Manual calculation of log fold change (mean_high - mean_low) for specified genes. Provides simple effect size estimation independent of limma.

### `coexpression_analysis`
Creates hierarchically clustered correlation heatmaps to identify co-expressed gene modules. Uses average linkage clustering.

### `plot_scatter_logfc_across_datasets`
Scatter plots comparing log fold changes between two datasets for the same genes. Helps identify consistent vs dataset-specific effects.

### `plot_bar_logfc_across_datasets`
Grouped bar plots showing log fold change across all three datasets simultaneously. Useful for cross-validation of findings.

### `plot_expression_subplots_before_after` & `plot_expression_grid_three_datasets`
Quality control plots showing expression distributions before and after normalization. Critical for verifying that normalization removed technical artifacts without over-correcting biological signal.

---

## Critical Points for Result Interpretation

### 1. Batch Effects and Normalization
- **TCGA**: Two platforms (RNAseq, Agilent) normalized separately
- **CGGA**: Two batches (325, 693) normalized separately, batch included as covariate in limma
- **GLASS**: Single batch, standard normalization
- **Implication**: Cross-dataset comparisons are valid only after within-dataset batch correction

### 2. Survival Group Definitions
- Patients alive with OS ≤ threshold are **excluded** (right-censored)
- High group includes both alive and dead patients (assumes alive patients will survive beyond threshold)
- Range thresholds exclude the middle interval to increase contrast
- **Implication**: Results apply to extreme survival groups, not intermediate patients

### 3. Sample Size Considerations
- Different datasets have different sample sizes after filtering
- IDH-wildtype GBM only (CGGA explicitly filtered, TCGA/GLASS assumed)
- Primary tumors only (CGGA explicitly filtered)
- **Implication**: Statistical power varies across datasets

### 4. Covariate Adjustment
- Age and sex are adjusted for in all analyses
- Batch is adjusted for in CGGA only
- **Implication**: Identified genes are associated with survival independent of these clinical variables

### 5. Gene Filtering
- Genes with >50% zeros removed
- Different gene sets retained across datasets due to different sequencing depths/platforms
- Union of significant genes used for cross-dataset comparisons
- **Implication**: Some genes may be dataset-specific due to technical rather than biological reasons

### 6. Expression Data Types
- **TCGA**: Mixed (RNAseq + microarray)
- **CGGA**: RSEM-normalized RNA-seq counts
- **GLASS**: TPM values
- **Implication**: Absolute expression levels not directly comparable across datasets; focus on fold changes and direction of effects

### 7. Multiple Testing
- FDR correction applied within each dataset separately
- No correction across multiple thresholds or datasets
- **Implication**: Genes significant at multiple thresholds or across multiple datasets are more robust findings

### 8. PCA Quality Control
- Before/after normalization plots verify technical correction
- Platform/batch separation should decrease after normalization
- Remaining separation may indicate biological heterogeneity
- **Implication**: Review PCA plots to ensure normalization was effective

---

## Output Structure

```
processed_data/
├── results/
│   ├── all_genes_[CONSORTIUM]_[THRESHOLD].csv
│   ├── significant_genes_[CONSORTIUM]_[THRESHOLD].csv
│   ├── [CONSORTIUM]_pca_before_norm.png
│   ├── [CONSORTIUM]_pca_after_norm.png
│   ├── coexpr_clustered_[CONSORTIUM]_[THRESHOLD].pdf
│   ├── scatter_logFC_[COMPARISON]_[THRESHOLD].png
│   └── logFC_bar_[THRESHOLD].png
```

---

## Dependencies

**Python packages:**
- pandas, numpy
- seaborn, matplotlib
- scikit-learn (StandardScaler, PCA)
- scipy (hierarchical clustering)
- rpy2 (R interface)

**R packages:**
- limma (differential expression)

---

## Reproducibility Notes

1. Random seed not set - PCA and clustering may show minor variations
2. File paths are hardcoded - update paths in `load_CGGA` and `load_GLASS` functions
3. StandardScaler uses default parameters (mean=0, std=1)
4. Log transformation uses log1p with small constant (1e-5) to avoid log(0)

---

## Recommended Usage

```python
# Run full pipeline
python script.py

# Typical thresholds tested
# Fixed: 24, 36, 48 months
# Range: (24,48), (24,60) months

# Output interpretation priority:
# 1. Check PCA plots for normalization quality
# 2. Review significant gene counts per dataset/threshold
# 3. Identify genes significant across multiple datasets (robust findings)
# 4. Examine log fold change consistency across datasets
# 5. Explore coexpression patterns for functional modules
```
