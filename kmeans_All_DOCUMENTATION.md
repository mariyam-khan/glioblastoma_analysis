# GLASS K-means Clustering & Unsupervised Subgroup Discovery

## Overview
This pipeline performs **unsupervised clustering** on the GLASS glioblastoma dataset to discover molecular subgroups without using survival information. It then tests whether these data-driven clusters associate with clinical outcomes, providing an unbiased approach to patient stratification.

---

## 🔴 CRITICAL DIFFERENCES FROM PREVIOUS ANALYSES

### **Paradigm Shift: Supervised vs Unsupervised**

| Aspect | Previous Analyses (#1-3) | Current Analysis |
|--------|-------------------------|------------------|
| **Approach** | **Supervised**: Pre-define groups by survival threshold | **Unsupervised**: Discover groups from expression data |
| **Group definition** | Low vs High survival (24 months) | K-means clusters (data-driven) |
| **Question** | "What genes differ between known survival groups?" | "Do expression-based clusters predict survival?" |
| **Bias** | Survival threshold is arbitrary (24, 36, 48 months) | No prior knowledge of outcomes |
| **Discovery** | Confirms known biology | Potentially reveals novel subgroups |
| **Clinical utility** | Requires outcome data | Can classify new patients immediately |

---

### **Key Methodological Innovations**

1. **Dual Clustering Strategy**:
   - **Gene-space clustering**: Find co-expressed gene modules
   - **Sample-space clustering**: Find patient subgroups

2. **Multi-group Comparison**:
   - Previous: Two groups (low vs high)
   - Current: K groups (2, 3, 4, 5, 6 clusters tested)
   - Statistical test: **F-test** across all groups (not t-test)

3. **Cluster Validation**:
   - Kruskal-Wallis test: Do clusters differ in survival?
   - Clinical annotation: Age, sex, survival characteristics per cluster
   - Heat-maps: Expression patterns that define clusters

4. **PCA Integration**:
   - Option to cluster in **reduced space** (first 2-3 PCs)
   - Reduces noise, focuses on major axes of variation

---

## EXPERIMENTAL DESIGN

### **Clustering Variants Tested (4 strategies)**

#### **1. Gene-Space Clustering (`kmeans_cluster_genes`)**
```
Genes × Samples matrix →
Transpose to Genes as observations →
K-means on full gene expression →
Clusters = co-expressed gene modules
```

**Purpose:** Identify genes with similar expression patterns across samples

**Use case:** Gene network analysis, pathway enrichment

**Not used for limma** (clustering genes, not samples)

---

#### **2. Gene-Space PCA Clustering (`kmeans_cluster_genes_pca`)**
```
Genes × Samples →
Transpose →
PCA (2-3 components) →
K-means on gene PCs →
Clusters = major gene modules
```

**Advantage:** Reduces dimensionality, focuses on dominant patterns

**Visualization:** Scatter plot of genes in PC1-PC2 space, colored by cluster

---

#### **3. Sample-Space Clustering (`kmeans_cluster_samples`)** ⭐ **MAIN ANALYSIS**
```
Samples × Genes matrix →
K-means on full expression →
Clusters = patient subgroups →
Run limma: gene expression ~ cluster + age + sex
```

**Purpose:** Discover molecular subtypes of patients

**This is the primary analysis for differential expression**

---

#### **4. Sample-Space PCA Clustering (`kmeans_cluster_samples_pca`)** ⭐ **DENOISED ANALYSIS**
```
Samples × Genes →
PCA (2-3 components) →
K-means on sample PCs →
Clusters = patient subgroups (noise-reduced) →
Run limma
```

**Advantage:** 
- Reduces curse of dimensionality
- Removes technical noise
- More stable clusters

**Trade-off:** May miss biology in lower PCs

---

## CORE FUNCTIONS

### `load_glass_data(filepath)`

**Purpose:** Load preprocessed GLASS data and separate metadata from genes

**Input:** CSV file from previous analyses (e.g., `GLASS_processed.csv`)

**Metadata columns (excluded from clustering):**
- `Sample.ID`
- `Overall.Survival.Months`
- `Diagnosis.Age`
- `Sex`
- `case.vital.status`

**Output:**
- `df`: Full DataFrame
- `gene_cols`: List of gene column names

**Note:** Assumes data is already normalized (from Analysis #1)

---

### `kmeans_cluster_samples(filepath, n_clusters=4)` ⭐ **PRIMARY**

**Pipeline:**
```python
1. Load data (samples × genes)
2. Extract expression matrix (samples × genes, numeric only)
3. Run K-means with k clusters
4. Assign cluster labels to samples
5. Add 'Cluster' column to DataFrame
6. Return: (df_with_clusters, kmeans_model)
```

**K-means Parameters:**
- `n_clusters`: Number of groups (typically 2-6)
- `random_state=42`: Reproducible clustering
- Distance metric: **Euclidean** (standard for K-means)

**Output:**
- `df['Cluster']`: Integer cluster assignment (0, 1, 2, ...)
- `kmeans`: Fitted KMeans object (for centroids, inertia)

---

### `kmeans_cluster_samples_pca(filepath, n_clusters=4, n_components=2)`

**Differences from basic K-means:**

1. **Dimensionality reduction first:**
```python
pca = PCA(n_components=2)
sample_pcs = pca.fit_transform(sample_data)
```

2. **Cluster in PC space:**
```python
kmeans.fit_predict(sample_pcs)  # Not on original genes
```

3. **Visualization included:**
- Scatter plot: PC1 vs PC2
- Colors: Cluster assignment
- Title indicates method

**When to use:**
- Many genes (>10,000): Reduces noise
- Batch effects: Major axes capture biology, minor axes capture technical
- Interpretability: PC1/PC2 often biologically meaningful

**When NOT to use:**
- Signal in lower PCs: Disease biology beyond major axes
- Small gene sets: Dimensionality not an issue

---

### `run_limma_multi(expr, clusters, covars, label)` ⭐ **CRITICAL**

**Purpose:** Multi-group differential expression (F-test across K clusters)

**Major difference from previous limma functions:**

**Previous analyses:**
```R
# Two-group t-test
~ group + Diagnosis.Age + Sex
coef = "grouphigh"  # high vs low
```

**Current analysis:**
```R
# Multi-group F-test
~ 0 + grp + Diagnosis.Age + Sex
coef = NULL  # F-test across all cluster levels
```

**Design matrix construction:**
```R
design <- model.matrix(~ 0 + grp)  # 0 = no intercept, one column per cluster
# If covariates provided:
design <- cbind(design, Diagnosis.Age, Sex)
```

**Statistical test:**
```R
topTable(fit, coef=NULL, number=Inf, sort.by='F')
```
- `coef=NULL` → Tests ALL cluster coefficients simultaneously
- Returns **F-statistic** and **F-test p-value**
- Null hypothesis: "Gene has same mean across all K clusters"
- Alternative: "At least one cluster differs"

**Output columns:**
- `Gene`: Gene symbol
- `AveExpr`: Average expression across all samples
- `F`: F-statistic (larger = more differential)
- `P.Value`: Raw p-value from F-test
- `adj.P.Val`: FDR-corrected p-value
- `(No logFC)`: F-test doesn't have a single fold-change

**Interpretation:**
- Significant gene → Expression differs across clusters (but doesn't say which pair differs)
- Top genes by F-statistic → Best cluster-defining genes

---

## VISUALIZATION FUNCTIONS

### `plot_volcano_and_heatmap(tt, expr, clusters, label, top_n=40)`

**Purpose:** Dual visualization of multi-cluster DE results

#### **Volcano Plot (adapted for F-test)**

**Challenge:** F-test has no single logFC (compares K groups, not 2)

**Solution:** Compute **maximum pairwise difference**
```python
cluster_means = expr.groupby(clusters).mean()  # mean per cluster
logfc_max = cluster_means.max() - cluster_means.min()  # range across clusters
```

**Plot:**
- X-axis: `logfc_max` (maximum difference between any two clusters)
- Y-axis: `-log10(P.Value)` from F-test
- Color: Red = significant (FDR < 0.1), Grey = not significant
- Threshold line: FDR = 0.1

**Interpretation:**
- Top-right: Large differences AND significant
- Top-left: Significant but small effect
- Bottom-right: Large differences but not significant (low power)

---

#### **Heat-map (cluster-defining genes)**

**Purpose:** Visualize expression patterns of top genes across clusters

**Method:**
1. Select top N genes (by F-statistic, default=40)
2. Extract expression matrix (samples × top genes)
3. **Z-score normalize** per gene (mean=0, std=1)
4. Hierarchical clustering of **samples** (not genes)
5. Color samples by cluster assignment
6. Display heat-map with cluster color bar

**Settings:**
- `row_cluster=True`: Reorder samples by similarity (within cluster)
- `col_cluster=False`: Keep genes in F-statistic order
- `cmap='vlag'`: Diverging colormap (blue–white–red)

**Interpretation:**
- Blocks of color → Cluster-specific expression patterns
- Genes in columns → Ordered by F-statistic (most differential at left)
- Row colors → Cluster assignment (should see grouping)

**Files saved:**
- `Volcano_{label}.png`
- `Heatmap_{label}.png`

---

### `plot_volcano_and_heatmap2` (Enhanced version)

**New features:**

1. **Dual row annotation bars:**
   - **Cluster**: Color by K-means cluster
   - **Sex**: Blue=Female (0), Pink=Male (1)

2. **Custom color scale:**
   - Symmetric around zero (white at zero)
   - `vmin = -vmax`, `vmax = ceiling(abs_max)`
   - Ensures equal color intensity for up/down

3. **Clinical summary printout:**
```
Cluster 0: n=35 | Age (median [IQR]) = 58.2 [51.3–64.7] | 
           Survival (median [IQR]) = 18.3 [9.2–26.1] | 
           Sex F/M = 12/23
Cluster 1: n=28 | ...
```

**Purpose:** Check if clusters differ in age, sex, or survival **before looking at DE genes**

**Critical check:**
- If clusters perfectly separate by sex → Sex confounds cluster assignment
- If clusters strongly differ in survival → Unsupervised method recovered prognostic groups!

---

### `overlay_clinical_data(df, case, n_cluster, cluster_col='Cluster')`

**Purpose:** Visual assessment of cluster-outcome associations

#### **1. Survival Boxplot**
```python
sns.boxplot(x='Cluster', y='Overall.Survival.Months', data=df)
sns.swarmplot(...)  # Individual points
```

**Interpretation:**
- Different medians → Clusters have prognostic value
- Overlapping boxes → Clusters don't predict survival well
- Outliers → Check if these are misclassified

---

#### **2. Vital Status Bar Plot**
```python
status_counts = pd.crosstab(df['Cluster'], df['case.vital.status'])
status_counts.plot(kind='bar', stacked=True)
```

**Shows:** Proportion alive vs dead per cluster

**Interpretation:**
- Cluster with high proportion dead → Poor prognosis subgroup
- Even distribution → Clusters not associated with outcome

---

### `kruskal_test_survival(df, cluster_col='Cluster', survival_col='Overall.Survival.Months')`

**Purpose:** Statistical test for survival differences across clusters

**Why Kruskal-Wallis (not ANOVA)?**
- **Non-parametric**: No normality assumption
- **Robust**: Handles skewed survival distributions
- **Multiple groups**: Compares K clusters simultaneously

**Null hypothesis:** All K clusters have same survival distribution

**Output:**
```
Kruskal-Wallis test result:
Statistic: 12.45, p-value: 0.006
```

**Interpretation:**
- p < 0.05 → Clusters differ significantly in survival ✅ 
- p > 0.05 → No evidence of survival differences ❌

**Limitation:** Doesn't say **which** clusters differ (just that at least one does)

**Follow-up:** Post-hoc pairwise tests (e.g., Dunn's test) if significant

---

## ELBOW PLOT ANALYSIS

### `elbow_plot_kmeans(X, case, max_clusters=10)`

**Purpose:** Determine optimal number of clusters (K)

**Method:**
1. Fit K-means for k = 1 to 10
2. Record **WCSS** (within-cluster sum of squares) = `kmeans.inertia_`
3. Compute **reduction in WCSS** between consecutive k:
```python
reductions[k] = WCSS[k-1] - WCSS[k]
```
4. Plot reduction vs k

**Elbow criterion:**
- Large reduction → Adding cluster explains a lot
- Small reduction → Diminishing returns
- **Elbow point**: Where reduction flattens

**Example interpretation:**
```
k=2: Reduction = 5000  (big gain)
k=3: Reduction = 3000  (still good)
k=4: Reduction = 800   (flattening) ← Elbow
k=5: Reduction = 500   (marginal)
```
**Recommendation:** Choose k=4

**Note:** Currently only runs for `n_clusters==2` in the code (can be enabled for all)

---

## MAIN PIPELINE WORKFLOW

### **Phase 1: Data Loading**
```python
filepath = "GLASS_processed.csv"
df, gene_cols = load_glass_data(filepath)
```

---

### **Phase 2: Clustering (Multiple K values)**
```python
for k in [2, 3, 4, 5, 6]:
    # Sample-space clustering
    df_clustered, kmeans = kmeans_cluster_samples(filepath, n_clusters=k)
    
    # (Optional) PCA-based clustering
    df_pca, kmeans_pca, pca = kmeans_cluster_samples_pca(
        filepath, n_clusters=k, n_components=2
    )
```

---

### **Phase 3: Clinical Association Testing**
```python
    # Visual inspection
    overlay_clinical_data(df_clustered, "sample_kmeans", k)
    
    # Statistical test
    stat, p_val = kruskal_test_survival(df_clustered)
    print(f"K={k}: Kruskal-Wallis p={p_val:.3f}")
```

**Decision point:** 
- If p < 0.05 → Clusters are prognostically relevant, proceed to DE
- If p > 0.05 → Clusters don't predict survival, reconsider K or method

---

### **Phase 4: Differential Expression (K=3 and K=4 only)**

**Why only K=3 and K=4?**
- K=2: Too simple, likely driven by one major axis
- K=3-4: Captures major subtypes without overfitting
- K>4: Risk of splitting clusters arbitrarily (overfitting)

```python
if k in (3, 4):
    # 1. Sex encoding
    df_clustered['Sex'] = df_clustered['Sex'].map({'female': 0, 'male': 1})
    
    # 2. Build expression matrix (genes × samples)
    gene_cols = [c for c in df.columns if c not in metadata_cols]
    expr = df_clustered.set_index('Sample.ID')[gene_cols].T
    
    # 3. Cluster labels
    clusters = df_clustered.set_index('Sample.ID')['Cluster']
    
    # 4. Covariates (adjust for age and sex)
    covars = df_clustered.set_index('Sample.ID')[['Diagnosis.Age', 'Sex']]
    
    # 5. Run multi-group limma
    tt = run_limma_multi(expr, clusters, covars, label=f"K{k}")
    
    # 6. Visualization
    plot_volcano_and_heatmap2(tt, expr.T, clusters, clinical, label=f"K{k}")
    
    # 7. Save significant genes
    tt_sig = tt[tt['adj.P.Val'] < 0.1]
    save_sheets_to_excel({f"K{k}": tt_sig}, EXCEL_OUT)
```

---

### **Phase 5: Cross-Comparison**
```python
# Compare gene lists between K=3 and K=4
compare_gene_lists(
    xlsx_path="GLASS_Kmeans_limma2.xlsx",
    sheet_a="K3",
    sheet_b="K4"
)
```

**Question:** Are the same genes differential regardless of K?

**High overlap** → Robust cluster-defining genes  
**Low overlap** → K-sensitive (clusters are different)

---

## OUTPUT FILES

### Excel File: `GLASS_Kmeans_limma2.xlsx`

**Sheets:**
- `K3`: Significant genes (FDR < 0.1) from 3-cluster analysis
- `K4`: Significant genes (FDR < 0.1) from 4-cluster analysis

**Columns:**
- `Gene`: Gene symbol
- `AveExpr`: Mean expression
- `F`: F-statistic
- `P.Value`: Raw p-value
- `adj.P.Val`: FDR-adjusted p-value

**Note:** No `logFC` column (F-test compares >2 groups)

---

### Visualization Files (per K value):

- `Volcano_K{k}.png`: Volcano plot (max pairwise difference vs F-test p)
- `Heatmap_K{k}.png`: Heat-map of top 40 genes with cluster colors
- `2Volcano_K{k}.png`: Enhanced volcano (if using `plot_volcano_and_heatmap2`)
- `2Heatmap_K{k}.png`: Enhanced heat-map with sex annotation

---

### Quality Control Plots:

- `PCA_plotsample_kmeans_{k}.png`: PCA of samples colored by cluster
- `PCA_plotsamplePCs_kmeans_{k}.png`: PCA clustering version
- `kmeans_survival_sample_kmeans_{k}.png`: Survival boxplot by cluster
- `kmeans_vital_sample_kmeans_{k}.png`: Vital status bar plot
- `elbow_plotsample_kmeans.png`: Elbow plot (if K=2)

---

## INTERPRETATION STRATEGY

### **Step 1: Validate Clustering Biologically**

**Check 1: Kruskal-Wallis test**
```
K=3: p=0.002 ✅  → Clusters predict survival
K=4: p=0.045 ✅  → Clusters predict survival  
K=5: p=0.31  ❌  → Too many clusters (overfitting)
```

**Check 2: Clinical summary**
```
Cluster 0: Median survival = 8 months   (Poor prognosis)
Cluster 1: Median survival = 24 months  (Intermediate)
Cluster 2: Median survival = 42 months  (Good prognosis)
```
→ Clusters recapitulate survival groups **without using survival in clustering**

**Check 3: Sex distribution**
```
Cluster 0: F/M = 15/20  (balanced)
Cluster 1: F/M = 8/22   (balanced)
```
→ Clusters not driven by sex (good)

---

### **Step 2: Inspect Cluster-Defining Genes**

**From heat-map:**
- Do blocks of genes show cluster-specific patterns?
- Are patterns consistent within clusters?
- Do high-survival clusters have different signature than low?

**From volcano:**
- How many significant genes? (K3 vs K4)
- What is the maximum difference (logfc_max)?

**Typical expectations:**
- K=3: 50-200 significant genes (broad differences)
- K=4: 100-400 genes (more refined)

---

### **Step 3: Compare K=3 vs K=4**

```python
compare_gene_lists(xlsx, "K3", "K4")
```

**High overlap (>70%):**
- Core subtype-defining genes robust to K
- Biology is stable
- Recommendation: Use K=3 (simpler)

**Moderate overlap (40-70%):**
- Some K-specific genes
- K=4 may reveal finer subtypes
- Recommendation: Report both, validate overlap

**Low overlap (<40%):**
- Clusters are unstable
- K-means may not be appropriate
- Consider: Hierarchical clustering, consensus clustering

---

### **Step 4: Biological Interpretation**

**Pathway enrichment on cluster-defining genes:**
```
Cluster 0 (poor prognosis): Proliferation, cell cycle, invasion genes UP
Cluster 1 (intermediate):  Immune response, inflammation genes UP  
Cluster 2 (good prognosis): Neuronal differentiation genes UP
```

**Compare to supervised analysis (Analyses #1-3):**
- Do cluster-defining genes overlap with survival-associated genes?
- Are there **novel genes** found only in unsupervised analysis?

**Validation in independent cohorts:**
- Can TCGA or CGGA patients be classified into these clusters?
- Do classifications predict survival in validation cohorts?

---

## CRITICAL INTERPRETATION POINTS

### 1. **Unsupervised ≠ Unbiased (Common Misconception)**

**Sources of bias:**
- Normalization method (z-score, quantile)
- Gene selection (all genes vs DE genes vs variable genes)
- Number of PCs (if using PCA clustering)
- Choice of K
- Random initialization (despite `random_state=42`)

**Recommendation:** 
- Run clustering on multiple gene sets (all genes, top variable genes)
- Test multiple K values
- Use consensus clustering for stability

---

### 2. **K-means Assumptions**

**Assumes:**
- **Spherical clusters**: Equal variance in all directions
- **Similar cluster sizes**: Roughly same N per cluster
- **Euclidean distance**: Appropriate for gene expression

**Violations:**
- If clusters are elongated → Consider Gaussian mixture models
- If cluster sizes very unequal → May split large cluster, miss small
- If batch effects present → May cluster by batch, not biology

**Check:** PCA plot should show roughly spherical, separated clusters

---

### 3. **Covariates in Multi-Group Limma**

**Design:**
```R
~ 0 + cluster + Diagnosis.Age + Sex
```

**Adjusting for Age and Sex:**
- Controls for confounding (like previous analyses)
- Ensures gene differences are due to cluster, not demographics

**No batch covariate:**
- GLASS is single-cohort (no batch)
- If using multi-cohort: Add batch

---

### 4. **F-test Interpretation Challenges**

**F-test says:** "Gene differs across K groups"

**F-test does NOT say:**
- Which clusters differ from which
- Direction of difference (up vs down)
- Magnitude of differences

**Post-hoc analyses needed:**
- Pairwise t-tests (Cluster 0 vs 1, 0 vs 2, 1 vs 2)
- Compute mean expression per cluster
- Identify cluster with highest/lowest expression

**Example:**
```python
# For a significant gene
gene_means = expr.T['GENE1'].groupby(clusters).mean()
print(gene_means)
# Cluster
# 0    5.2  ← Lowest
# 1    6.8
# 2    7.9  ← Highest
```

---

### 5. **Cluster Stability**

**Question:** Are clusters reproducible?

**Test 1: Different random seeds**
```python
for seed in [42, 123, 999]:
    kmeans = KMeans(n_clusters=4, random_state=seed)
    clusters_seed = kmeans.fit_predict(data)
    # Compare cluster assignments (adjusted for label permutation)
```

**Test 2: Bootstrap / cross-validation**
- Resample data, re-cluster
- Check: Are same samples consistently grouped?

**Test 3: Consensus clustering**
- Cluster 100 times with random subsets
- Build consensus matrix (how often are pairs co-clustered?)
- Stable clusters have high within-cluster consensus

**Recommendation:** If using for clinical decision-making, must validate stability

---

### 6. **Clinical Utility Assessment**

**For unsupervised clusters to be clinically useful:**

**Must satisfy:**
1. ✅ Clusters predict survival (Kruskal-Wallis p < 0.05)
2. ✅ Clusters differ in molecular profile (significant DE genes)
3. ✅ Clusters are stable (reproducible across seeds/methods)
4. ✅ Classification method is practical (can assign new patients)

**Bonus points:**
- Clusters respond differently to treatment
- Clusters have distinct pathway signatures
- Can be measured with clinical assay (e.g., qPCR panel)

---

### 7. **Comparison to Supervised Analysis**

**Supervised (Analyses #1-3):**
- Pre-define groups by survival cutoff
- Find genes that differ between groups
- Validates: "These genes associate with survival"

**Unsupervised (Current):**
- Discover groups from expression alone
- Test if groups differ in survival
- Validates: "Expression patterns predict survival"

**Complementary findings:**
- Genes in **both** approaches → Most robust biomarkers
- Genes only in unsupervised → Novel biology not captured by arbitrary cutoffs
- Genes only in supervised → May be confounded by survival-related effects (e.g., tumor burden)

---

## RECOMMENDED ANALYSIS SEQUENCE

### **Preliminary:**
1. Load GLASS data (already normalized from Analysis #1)
2. Run elbow plot (k=1 to 10) to guide K selection

### **Main analysis:**
3. Cluster with K=2,3,4,5,6
4. For each K:
   - Kruskal-Wallis test (survival association)
   - Clinical overlay plots (survival boxplot, vital status)
   - Print clinical summary

### **Focused DE analysis:**
5. Choose 1-2 K values with:
   - Significant Kruskal-Wallis (p < 0.05)
   - Biologically interpretable (3-4 clusters typically optimal)
   - Distinct survival curves

6. For chosen K:
   - Run multi-group limma (F-test)
   - Generate volcano + heat-map
   - Save significant genes

### **Validation:**
7. Compare gene lists across K values (overlap)
8. Pathway enrichment on cluster-defining genes
9. Cross-reference with supervised DE genes (Analyses #1-3)

### **Optional advanced:**
10. Consensus clustering for stability
11. Apply cluster model to TCGA/CGGA for validation
12. Post-hoc pairwise tests (which clusters differ for each gene)

---

## ADVANCED: Applying Clusters to New Data

**Goal:** Classify TCGA or CGGA samples into GLASS-derived clusters

**Method:**
```python
# Train on GLASS
kmeans = KMeans(n_clusters=4, random_state=42)
kmeans.fit(glass_expr)  # samples × genes

# Predict on TCGA (must have same genes)
tcga_clusters = kmeans.predict(tcga_expr)  # same gene columns

# Test: Do predicted clusters differ in survival in TCGA?
kruskal_test_survival(tcga_data_with_clusters)
```

**Interpretation:**
- If p < 0.05 in TCGA → Clusters are generalizable across cohorts ✅
- If p > 0.05 in TCGA → Clusters are GLASS-specific ❌

---

## COMPARISON TO PREVIOUS ANALYSES

| Question | Supervised (Analyses #1-3) | Unsupervised (Current) |
|----------|---------------------------|------------------------|
| **Starting point** | Survival groups (24mo cutoff) | Expression data only |
| **Genes found** | Survival-associated | Cluster-defining |
| **Statistical test** | t-test (2 groups) | F-test (K groups) |
| **Bias** | Arbitrary cutoff | K selection bias |
| **Discovery** | Limited to known groups | Can find novel subtypes |
| **Clinical use** | Requires follow-up | Immediate classification |
| **Validation** | Cross-cohort replication | Cluster stability |

---

## WHEN TO USE THIS ANALYSIS

**Use unsupervised clustering when:**
1. Exploring heterogeneity **without prior knowledge** of subgroups
2. Survival cutoffs are **arbitrary** or **unclear**
3. Want to discover **novel subtypes** beyond survival
4. Planning **prospective classification** of patients
5. Suspect **more than 2 subtypes** exist

**Use supervised analysis when:**
1. Have clear clinical grouping (e.g., responders vs non-responders)
2. Want to validate known biology
3. Need interpretable fold-changes (up/down)
4. Have limited sample size (supervised more powerful)

**Use BOTH when:**
- Maximum rigor required for biomarker discovery
- Genes robust to both approaches are highest confidence

---

## KEY TAKEAWAYS

🎯 **Main question:** Can molecular subtypes be discovered without using survival?

✅ **If Kruskal-Wallis p < 0.05:** Unsupervised clusters recovered prognostic groups

🔴 **Red flags:**
- Clusters differ by sex/age only → Confounded, not biological
- K=2 finds same split as survival cutoff → Just rediscovered supervised groups
- High K (>5) needed for significance → Overfitting

🌟 **Gold standard genes:**
- Significant in K=3 **AND** K=4 (robust to K)
- Also significant in supervised analyses (robust to approach)
- Consistent direction (up/down) across all tests
- Replicate in independent cohorts

---

## Reproducibility Notes

1. **K-means random seed:** Set to 42 (reproducible)
2. **PCA random state:** Should be set but not visible in code (minor variation possible)
3. **Gene order:** Affects K-means (uses Euclidean distance)
4. **Sex encoding:** female=0, male=1 (check data consistency)
5. **Elbow plot:** Only runs for K=2 (enable for all K to guide selection)
