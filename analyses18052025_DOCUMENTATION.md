# Multi-Cohort Normalization Comparison Pipeline

## Overview
This pipeline systematically compares **z-score normalization** vs **quantile normalization + ComBat** for TCGA and CGGA datasets, with and without treatment filtering, and tests cross-cohort integration. The goal is to determine which normalization strategy yields the most robust survival-associated genes.

---

## 🔴 CRITICAL DIFFERENCES FROM PREVIOUS ANALYSES

### **Comparison to Analysis #1 (Separate Platform Analysis)**
| Aspect | Analysis #1 | Current Analysis |
|--------|-------------|------------------|
| **CGGA treatment filter** | Not applied | **Both & and OR filtering tested** |
| **TCGA sample filter** | Basic overlap | **Whitelist from treatment file** |
| **Normalization** | StandardScaler only | **Z-score vs Quantile+ComBat** |
| **Batch correction** | None (CGGA) or in model (CGGA) | **ComBat applied before limma** |
| **Cross-cohort** | Not performed | **TCGA+CGGA integrated analysis** |

### **Comparison to Analysis #2 (Integrated TCGA)**
| Aspect | Analysis #2 | Current Analysis |
|--------|-------------|------------------|
| **Focus** | TCGA RNA-seq + Agilent | **TCGA + CGGA cross-cohort** |
| **RNA-seq norm** | CPM → voom | **Log2 → z-score or quantile** |
| **Integration method** | ComBat (platform batch) | **ComBat (cohort batch)** |
| **Comparisons** | AR vs RA stacking | **Z-score vs Quantile methods** |
| **Treatment filter** | Not tested | **AND vs OR tested** |

---

## 🎯 EXPERIMENTAL DESIGN

This analysis tests **8 different preprocessing strategies** systematically:

### **TCGA Variants (4 total)**
1. **TCGA-z**: Z-score normalization, **treatment-filtered** samples
2. **TCGA-z_nc**: Z-score normalization, **no treatment filter** (broader inclusion)
3. **TCGA-qn**: Quantile normalization + ComBat, treatment-filtered
4. **(Implied) TCGA-qn_nc**: Quantile + ComBat, no filter (not explicitly run)

### **CGGA Variants (4 total)**
1. **CGGA-z**: Z-score normalization, **AND filter** (radio=1 AND chemo=1)
2. **CGGA-z_or**: Z-score normalization, **OR filter** (radio=1 OR chemo=1)
3. **CGGA-qn**: Quantile + ComBat, AND filter
4. **CGGA-qn_or**: Quantile + ComBat, OR filter

### **Integrated Cross-Cohort (2 total)**
1. **Integrated_z**: TCGA_z_nc + CGGA_z_or merged, z-score normalized
2. **Integrated_qn**: TCGA_qn + CGGA_qn_or merged, global quantile normalization

**Total comparisons**: 9 limma analyses (TCGA×3 + CGGA×4 + Integrated×2)

---

## KEY METHODOLOGICAL DECISIONS

### 1. Treatment Filtering Logic (NEW)

**TCGA Treatment Filter (`care=True`):**
```python
# Uses external whitelist from treatment file
id_data = pd.read_excel("TCGA_filter_treatment1.xlsx")
keep_ids = set(id_data['Case ID'])
```
**Purpose:** Restrict to patients with documented treatment (likely radiotherapy + chemotherapy)

**CGGA Treatment Filter:**

**AND filter (`care=True`):**
```python
c325 = c325.loc[
    (c325['Radio_status (treated=1;un-treated=0)'] == 1) &
    (c325['Chemo_status (TMZ treated=1;un-treated=0)'] == 1)
]
```
**Result:** Only patients who received **both** radiotherapy **and** TMZ chemotherapy

**OR filter (`care=False`):**
```python
c325 = c325.loc[
    (c325['Radio_status (treated=1;un-treated=0)'] == 1) |
    (c325['Chemo_status (TMZ treated=1;un-treated=0)'] == 1)
]
```
**Result:** Patients who received **either** radiotherapy **or** TMZ (more inclusive)

**CRITICAL IMPACT:**
- AND filter: Smaller N, more homogeneous treatment, genes specific to combined therapy
- OR filter: Larger N, more heterogeneous, genes generalizable to broader treatment

---

### 2. Two Normalization Pipelines (CORE COMPARISON)

#### **Pipeline A: Z-Score Normalization**

**TCGA (`load_TCGA`):**
```
Raw expression → 
Per-platform StandardScaler (separately for RNAseq & Agilent) →
Concatenate platforms
```

**CGGA (`CGGA_zscores`):**
```
RSEM counts →
log1p transform →
Per-batch StandardScaler (separately for 325 & 693) →
Concatenate batches
```

**Properties:**
- Gene-wise centering (mean=0) and scaling (std=1)
- **Does NOT equalize distributions** across samples
- Preserves biological variance structure
- No cross-sample assumptions

---

#### **Pipeline B: Quantile Normalization + ComBat**

**TCGA (`TCGA_QN`):**
```
Raw expression →
Quantile normalization (all samples together) →
ComBat batch correction (platform = batch) →
PCA before/after ComBat
```

**CGGA (`CGGA_QN_combat`):**
```
RSEM counts →
log1p transform →
Quantile normalization (all samples together) →
ComBat batch correction (batch325 vs batch693) →
PCA before/after ComBat
```

**Properties:**
- **Forces identical distributions** across all samples
- Removes technical variation more aggressively
- ComBat adjusts for location/scale shifts per batch
- Assumes batch effects are additive and parametric

---

### 3. ComBat Implementation (Scanpy)

**Function: `combat` (from helpers.py)**
```python
def combat(combined, batch):
    adata = ad.AnnData(combined.T)
    adata.X = adata.X.astype("float64", copy=False)
    adata.obs['batch'] = batch.astype('category')
    sc.pp.combat(adata, key='batch')
    corrected = pd.DataFrame(adata.X.T, ...)
    return corrected
```

**Applied to:**
- TCGA: Platform batch (RNAseq vs Agilent)
- CGGA: Sequencing batch (325 vs 693)
- Integrated: Cohort batch (TCGA vs CGGA)

**ComBat Assumptions:**
1. Batch effects are **additive**
2. Follow **parametric distributions** (location-scale family)
3. Biological groups are **balanced** across batches
4. Non-biological variance is primarily due to batch

**Risk:** If survival groups are unbalanced across batches, ComBat may remove real biology

---

### 4. Quantile Normalization Implementation

**Function: `quantile_normalize_combined` (from helpers.py)**
```python
def quantile_normalize_combined(combined):
    rmat = pandas2ri.py2rpy(combined)
    norm_r = limma.normalizeBetweenArrays(rmat, method="quantile")
    return pd.DataFrame(np.array(norm_r), ...)
```

**Applied:**
1. **Within-cohort**: Before ComBat (TCGA-qn, CGGA-qn)
2. **Cross-cohort**: After merging TCGA+CGGA (Integrated_qn)

**Effect:** Every sample has **identical empirical distribution**
- Same median, quartiles, extremes
- Removes all distributional differences (technical AND biological?)
- Most aggressive normalization

---

## Dataset Loading Functions

### `load_TCGA(care: bool)`

**Differences based on `care` flag:**

| `care=True` (Treatment-filtered) | `care=False` (Broader) |
|----------------------------------|------------------------|
| Load whitelist from `TCGA_filter_treatment1.xlsx` | Use all samples in outcome file |
| Strip `-01` suffix from sample IDs | Keep original IDs |
| Intersection with whitelist | Intersection with outcome file |
| Result: ~N samples (treatment-documented) | Result: Larger N (all available) |

**Common steps:**
1. Filter to RNAseq + Agilent platforms
2. Gene filtering (`remove_zero_genes`)
3. PCA before normalization (visualize platform separation)
4. **Separate normalization** per platform (StandardScaler)
5. PCA after normalization (verify correction)
6. Merge with clinical data
7. Survival group assignment (24-month threshold)
8. Sex encoding: Female=0, Male=1

**Output:** 
- `outcome_data`: Sample metadata with survival groups
- `df_norm`: Normalized expression (samples × genes)

---

### `load_CGGA(care: bool)`

**Treatment filtering:**

| `care=True` (AND logic) | `care=False` (OR logic) |
|-------------------------|-------------------------|
| Radio=1 **AND** Chemo=1 | Radio=1 **OR** Chemo=1 |
| Stricter: both treatments | Permissive: either treatment |
| Smaller N | Larger N |

**Common filters (both modes):**
- PRS_type == 'Primary'
- Histology == 'GBM'
- IDH_mutation_status == 'Wildtype'

**Processing:**
1. Load 325 and 693 batches separately
2. Transpose expression (genes as columns)
3. Convert OS from days to months (÷30.44)
4. Intersect samples between expression and clinical
5. Gene filtering (`remove_zero_genes`)
6. Find common genes across batches

**Output:**
- `df_325_t`, `df_693_t`: Expression DataFrames
- `c325`, `c693`: Clinical DataFrames

**DOES NOT merge batches yet** - returns separate for flexible processing

---

### `TCGA_QN(care: bool)` (NEW - Quantile pipeline)

**Major difference from `load_TCGA`:**
- Applies **quantile normalization** to ALL samples together (not per-platform)
- Then applies **ComBat** to remove platform batch effects
- Generates PCA plots before/after ComBat

**Pipeline:**
```
1. Load data (with or without treatment filter)
2. Gene filtering
3. Quantile normalize: genes × all_samples
4. ComBat correct: batch = platform
5. PCA visualization (2 plots)
6. Convert back to samples × genes
7. Merge with clinical, survival filtering
```

**Critical:** Platform effects corrected by ComBat, **not included in limma design**

---

## CGGA Processing Functions

### `PCA_cgga(batch1, batch2)` (Outlier removal)

**Purpose:** Visual QC + outlier detection before normalization

**Method:**
1. Concatenate batch325 + batch693
2. PCA on raw log-transformed data
3. Scatter plot colored by batch
4. **Remove outliers**: |PC1| > 100 OR |PC2| > 100
5. Filter original batches to exclude outliers

**Output:** Cleaned batch1, batch2 DataFrames

**Note:** Threshold of 100 is **hardcoded** - may need adjustment for different scales

---

### `CGGA_zscores(batch1, batch2)` (Z-score pipeline)

**Pipeline:**
```
1. Find common genes between batches
2. StandardScaler applied SEPARATELY to batch325 and batch693
3. Concatenate batches
4. PCA after normalization
5. Visualize batch separation (should be reduced)
```

**Output:** Normalized batch1, batch2 (still separate DataFrames)

**Design choice:** Per-batch normalization → each batch has mean=0, std=1 independently

---

### `CGGA_QN_combat(batch1, batch2)` (Quantile + ComBat pipeline)

**Pipeline:**
```
1. Extract gene columns from each batch
2. Transpose to genes × samples (per batch)
3. Concatenate: [batch325_samples, batch693_samples]
4. SINGLE quantile normalization on combined matrix
5. Create batch vector (325×N1 + 693×N2)
6. ComBat correction (batch as covariate)
7. PCA before ComBat (still has batch effects)
8. PCA after ComBat (batch effects removed)
9. Transpose back to samples × genes
10. Split back into batch1, batch2 DataFrames
```

**Output:** ComBat-corrected batch1, batch2

**Critical assertion:**
```python
assert expr_df.shape[0] == concat.shape[0], "Row mismatch"
```
Ensures sample counts match after transformations

---

## Differential Expression

### `run_limma` (Simplified from Analysis #2)

**Design matrix:**
```R
~ group + Diagnosis.Age + Sex + batch  # if batch column exists
~ group + Diagnosis.Age + Sex          # otherwise
```

**For quantile-normalized data:**
- **No additional normalization** in limma
- Goes straight to `lmFit` (assumes already normalized)

```R
fit <- lmFit(expr, design)
fit <- eBayes(fit)
tt  <- topTable(fit, coef=2, number=Inf, adjust.method='fdr')
```

**Coefficient tested:** `coef=2` → `grouphigh` (high survival vs low)

**Output:** 
- Returns DataFrame with all genes
- Prints count of significant genes (FDR < 0.1)
- **Does NOT write Excel directly** (collected later)

---

## Quality Control & Visualization

### `volcano_plot` (NEW)

**Simpler than Analysis #2:**
- Single dataset per plot
- Red = significant (FDR < 0.1)
- Grey = not significant
- Threshold lines: FDR=0.1, logFC=±1

**Purpose:** Quick visual summary of DE results

---

### `plot_pca_integration` (Cross-cohort QC)

**Purpose:** Assess batch effects when merging TCGA + CGGA

**Method:**
1. Find common genes between cohorts
2. Subset to common genes
3. Concatenate samples: TCGA + CGGA
4. **StandardScaler** on concatenated matrix (centers for PCA)
5. PCA (2 components)
6. Color by cohort (TCGA=green, CGGA=purple)

**Interpretation:**
- **Good**: Cohorts overlap → batch effects minimal
- **Bad**: Cohorts separate → batch effects dominate
- Compare z-score vs quantile versions

**Files generated:**
- `PCA_TCGA-CGGA_z.png`
- `PCA_TCGA-CGGA_qn.png`

---

## Main Pipeline Workflow

### **Phase 1: CGGA Analysis (4 versions)**

```
A. Treatment AND filter (care=True)
   1. Load + log transform
   2. Merge batches 325 + 693
   3. PCA + outlier removal
   4. Split into two normalization paths:
      → Z-score: CGGA-z
      → Quantile+ComBat: CGGA-qn

B. Treatment OR filter (care=False)
   1. Load + log transform
   2. Merge batches 325 + 693
   3. PCA + outlier removal
   4. Split into two normalization paths:
      → Z-score: CGGA-z_or
      → Quantile+ComBat: CGGA-qn_or

For each version:
   → Survival filtering (24 months)
   → Sex encoding (Female=0, Male=1)
   → Run limma
   → Volcano plot
```

---

### **Phase 2: TCGA Analysis (3 versions)**

```
A. Treatment-filtered (care=True)
   1. Load from whitelist
   2. Split into two normalization paths:
      → Z-score: TCGA-z
      → Quantile+ComBat: TCGA-qn

B. No treatment filter (care=False)
   1. Load all available samples
   2. Z-score normalization: TCGA-z_nc

For each version:
   → Survival filtering (24 months)
   → Sex encoding (Female=0, Male=1)
   → Run limma
   → Volcano plot
```

---

### **Phase 3: Cross-Cohort Integration (2 versions)**

**Z-Score Integration:**
```
1. TCGA-z_nc (broader TCGA) + CGGA-z_or (broader CGGA)
2. Find common genes
3. Concatenate samples (samples × common_genes)
4. Assign batch labels (CGGA vs TCGA)
5. Run limma with batch covariate
6. Volcano plot
7. PCA visualization (assess batch effects)
```

**Quantile Integration:**
```
1. TCGA-qn + CGGA-qn_or (both already ComBat-corrected internally)
2. Find common genes
3. Concatenate samples
4. GLOBAL quantile normalization on merged matrix
5. Assign batch labels
6. Run limma with batch covariate
7. Volcano plot
8. PCA visualization
```

**CRITICAL DIFFERENCE:**
- Z-score: No additional normalization after merge
- Quantile: **Second round** of quantile normalization on merged data

---

## Output Structure

### Excel File: `limma_significant_genes2.xlsx`

**Sheets (9 total):**
1. **TCGA_z**: Treatment-filtered, z-score normalized
2. **TCGA_z_nc**: No treatment filter, z-score normalized
3. **TCGA_qn**: Treatment-filtered, quantile+ComBat
4. **CGGA_z**: AND treatment filter, z-score
5. **CGGA_qn**: AND treatment filter, quantile+ComBat
6. **CGGA_z_or**: OR treatment filter, z-score
7. **CGGA_qn_or**: OR treatment filter, quantile+ComBat
8. **Integrated_z**: Cross-cohort, z-score
9. **Integrated_qn**: Cross-cohort, quantile+ComBat

**Each sheet contains:**
- Only significant genes (FDR < 0.1)
- Columns: Gene, logFC, AveExpr, t, P.Value, adj.P.Val, B

---

## Interpretation Strategy

### 1. **Within-Cohort Comparisons**

**TCGA: Z-score vs Quantile**
```python
compare_gene_lists(xlsx, "TCGA_z", "TCGA_qn")
```
**Question:** Does aggressive normalization (quantile+ComBat) change gene rankings?

**High overlap** → Robust to normalization choice  
**Low overlap** → Method-sensitive genes (caution!)

---

**CGGA: Z-score vs Quantile**
```python
compare_gene_lists(xlsx, "CGGA_z", "CGGA_qn")
```
**Question:** Same as TCGA, but with two batches

---

### 2. **Treatment Filter Sensitivity**

**CGGA: AND vs OR**
```python
compare_gene_lists(xlsx, "CGGA_z", "CGGA_z_or")
compare_gene_lists(xlsx, "CGGA_qn", "CGGA_qn_or")
```
**Question:** Do genes differ when we:
- Require **both** treatments (AND, N=smaller)
- Allow **either** treatment (OR, N=larger)

**High overlap** → Genes general to treated patients  
**Low overlap** → Genes specific to combined therapy

---

**TCGA: Treatment-filtered vs All**
```python
compare_gene_lists(xlsx, "TCGA_z", "TCGA_z_nc")
```
**Question:** Does restricting to treatment-documented samples matter?

---

### 3. **Cross-Cohort Validation**

**TCGA vs CGGA (same normalization)**
```python
compare_gene_lists(xlsx, "TCGA_z", "CGGA_z_or")
compare_gene_lists(xlsx, "TCGA_qn", "CGGA_qn_or")
```
**Question:** Which genes replicate across independent cohorts?

**These are the MOST ROBUST findings** → prioritize for validation

---

**Single-cohort vs Integrated**
```python
compare_gene_lists(xlsx, "TCGA_z", "Integrated_z")
compare_gene_lists(xlsx, "CGGA_z_or", "Integrated_z")
```
**Question:** Does integration (more power) find genes missed in single cohorts?

---

### 4. **Normalization Robustness**

**Within same samples, across methods:**
```python
compare_gene_lists(xlsx, "TCGA_z", "TCGA_qn")
compare_gene_lists(xlsx, "Integrated_z", "Integrated_qn")
```

**Gold standard:** Genes significant in **both** z-score **and** quantile methods

---

## Critical Interpretation Points

### 1. **Quantile Normalization May Over-Correct**
- Forces **all samples** to have identical distribution
- Can remove **biological variance** along with technical
- Especially risky if:
  - Survival groups have different expression profiles
  - Sample quality varies substantially

**Check:** If TCGA_qn finds far fewer genes than TCGA_z → possible over-correction

---

### 2. **ComBat Batch Assumptions**
Applied to:
- TCGA: Platform batch (RNAseq vs Agilent)
- CGGA: Sequencing batch (325 vs 693)
- Integrated: Cohort batch (TCGA vs CGGA)

**Assumes:**
- Batch effects are additive
- Survival groups balanced across batches

**Risk:** If low-survival patients concentrated in one batch, ComBat may remove signal

**Check PCA plots:**
- Before ComBat: batches should separate
- After ComBat: batches should overlap **but biological groups still separate**

---

### 3. **Treatment Filtering Changes the Question**

**AND filter (both treatments):**
- Smaller N, less power
- Genes specific to **combined** radio+chemo
- Clinically: "What predicts survival in fully-treated patients?"

**OR filter (either treatment):**
- Larger N, more power
- Genes general to **any treatment**
- Clinically: "What predicts survival in treated vs untreated?"

**No filter:**
- Largest N, maximum power
- Includes untreated patients (may confound survival)
- Clinically: "What predicts survival regardless of treatment?"

**Recommendation:** Report AND and OR separately; genes in both are most robust

---

### 4. **Sample Size Hierarchy**

**Expected ranking (largest to smallest N):**
1. TCGA_z_nc (no treatment filter)
2. TCGA_z (treatment filter)
3. CGGA_z_or (OR filter)
4. CGGA_z (AND filter)
5. Integrated_z (depends on which TCGA/CGGA used)

**Statistical power increases with N:**
- Larger cohorts detect smaller effect sizes
- But heterogeneity also increases (noise)

**Trade-off:** TCGA_z_nc has most power but most heterogeneity

---

### 5. **Limma Batch Covariate**

**Integrated analyses include batch:**
```R
~ group + Diagnosis.Age + Sex + batch
```
Where `batch = TCGA or CGGA`

**Purpose:** Account for cohort-level differences not removed by normalization

**Assumption:** Batch effects are **linear and additive** in the model

**Note:** This is **in addition to** ComBat (for quantile version)
- ComBat: removes global batch effects
- Batch covariate: adjusts for residual batch effects

**Potential issue:** Double-correction may be overly conservative

---

### 6. **PCA Interpretation for Integration**

**Z-score integration PCA:**
- Expect some separation (no aggressive correction)
- If TCGA and CGGA completely separate → batch effects dominate
- If they overlap → successful integration

**Quantile integration PCA:**
- Should show **more overlap** than z-score (due to ComBat)
- If still separated → ComBat failed or biology differs
- If perfectly overlapping → good correction (but check volcano for over-correction)

**Red flag:** If quantile PCA shows **perfect** overlap but finds far fewer DE genes than z-score → likely over-corrected

---

### 7. **Gene Count Expectations**

**Typical patterns:**

**Within normalization:**
- TCGA_z vs TCGA_qn: 60-80% overlap (some method-specific)
- CGGA_z vs CGGA_qn: 60-80% overlap

**Across cohorts:**
- TCGA_z vs CGGA_z: 30-50% overlap (cohort differences expected)
- Lower overlap normal due to different platforms, populations

**Integration boost:**
- Integrated should find **more** genes than single cohorts (more power)
- If Integrated finds **fewer** → over-correction or heterogeneity issues

---

## Recommended Analysis Sequence

### **Step 1: Within-Cohort Method Comparison**
```python
# Do z-score and quantile agree within same cohort?
compare_gene_lists(xlsx, "TCGA_z", "TCGA_qn")
compare_gene_lists(xlsx, "CGGA_z", "CGGA_qn")
```
**Goal:** Identify method-robust genes

---

### **Step 2: Treatment Filter Sensitivity**
```python
# Does treatment filtering change results?
compare_gene_lists(xlsx, "TCGA_z", "TCGA_z_nc")
compare_gene_lists(xlsx, "CGGA_z", "CGGA_z_or")
```
**Goal:** Identify treatment-independent genes

---

### **Step 3: Cross-Cohort Replication**
```python
# Do genes replicate across TCGA and CGGA?
compare_gene_lists(xlsx, "TCGA_z_nc", "CGGA_z_or")  # broadest
compare_gene_lists(xlsx, "TCGA_qn", "CGGA_qn_or")
```
**Goal:** Identify cohort-independent genes

---

### **Step 4: Integration Value**
```python
# Does integration find additional genes?
compare_gene_lists(xlsx, "TCGA_z_nc", "Integrated_z")
compare_gene_lists(xlsx, "CGGA_z_or", "Integrated_z")
```
**Goal:** Identify genes requiring combined power

---

### **Step 5: Final Robust Set**

**Criteria for MOST ROBUST genes:**
1. Significant in **both** z-score **and** quantile (method-independent)
2. Significant in **both** TCGA **and** CGGA (cohort-independent)
3. **Consistent direction** (logFC same sign) across all analyses
4. FDR < 0.05 (stricter threshold) in at least one analysis

**Conservative estimate:** Intersection of TCGA_z, TCGA_qn, CGGA_z_or, CGGA_qn_or

---

## Utility Function

### `compare_gene_lists`

**Purpose:** Quick summary of overlap between two analyses

**Usage:**
```python
stats = compare_gene_lists(
    xlsx_path="limma_significant_genes2.xlsx",
    sheet_a="TCGA_z",
    sheet_b="CGGA_z_or",
    gene_column="Gene"
)
print(stats)
```

**Output:**
```
                    Sheet  n_genes         comment
                   TCGA_z      150                
               CGGA_z_or      120                
  TCGA_z ∩ CGGA_z_or       45  genes in common
```

**Interpretation:**
- 45/150 (30%) of TCGA genes replicate in CGGA
- 45/120 (38%) of CGGA genes replicate in TCGA
- These 45 genes are prioritized for validation

---

## File Dependencies

**External:**
- `pipeline_CGGA.py`: Provides `remove_zero_genes`, `log_transform`, `make_surv_groups_int`
- `helpers.py`: Provides `quantile_normalize_combined`, `combat`, `pca_plot`, `pca_TCGA`

**Input files:**
- TCGA expression: `GBMLGG_EB_RmDiffFullGenesRanRmDup.csv`
- TCGA clinical: `Compiled data_TCGA_CGGA_GLASS_jh080424.xlsx`
- TCGA treatment whitelist: `TCGA_filter_treatment1.xlsx`
- CGGA 325 expression: `CGGA.mRNAseq_325.RSEM-genes.20200506.txt`
- CGGA 693 expression: `CGGA.mRNAseq_693.RSEM-genes.20200506.txt`
- CGGA 325 clinical: `CGGA.mRNAseq_325_clinical.20200506.txt`
- CGGA 693 clinical: `CGGA.mRNAseq_693_clinical.20200506.txt`

---

## Summary: Three-Way Analysis Comparison

| Question | Analysis #1 | Analysis #2 | Current Analysis |
|----------|-------------|-------------|------------------|
| **Goal** | Platform-specific | TCGA integration | Method comparison |
| **Datasets** | TCGA, CGGA, GLASS separate | TCGA (RNA+Agilent) | TCGA + CGGA |
| **Normalization** | StandardScaler | CPM/voom/quantile | Z-score vs Quantile |
| **Batch correction** | None or model-based | ComBat (platform) | ComBat (cohort) |
| **Treatment filter** | Not tested | Not tested | **AND vs OR tested** |
| **Integration** | None | Within-TCGA | **Cross-cohort** |
| **Best for** | Conservative | Max TCGA power | Robustness testing |

---

## When to Use This Analysis

**Use this analysis when:**
1. Want to test **robustness** of findings across normalization methods
2. Need to validate genes across **independent cohorts**
3. Concerned about **over-correction** from aggressive normalization
4. Want to assess **treatment-specific** vs **general** survival genes
5. Planning experimental validation (need highest-confidence genes)

**Use Analysis #1 when:**
- Want platform-specific insights
- Concerned about over-correction
- Need GLASS cohort included

**Use Analysis #2 when:**
- Focused on TCGA only
- Want maximum power from platform integration
- RNA-seq specific analysis needed

---

## Reproducibility Notes

1. **PCA outlier threshold**: Hardcoded at |PC| > 100 for CGGA
2. **Random seeds**: Not set (minor PCA variation)
3. **Quantile normalization**: Uses R's limma implementation (via rpy2)
4. **ComBat**: Uses scanpy implementation (Python)
5. **Gene filtering**: Applied independently to each batch before merging

---

## Key Takeaways

🎯 **Primary goal:** Identify genes robust to:
- Normalization method (z-score vs quantile)
- Treatment filtering (AND vs OR)
- Cohort (TCGA vs CGGA)

🔴 **Red flags to watch:**
- Quantile version finds <50% of z-score genes → over-correction
- TCGA and CGGA have <20% overlap → cohort-specific biology
- Integrated finds fewer genes than single cohorts → heterogeneity issue

✅ **High-confidence genes:**
- Significant in ≥4 analyses
- Same direction across all
- Replicate in TCGA and CGGA
- Robust to normalization method
