# Cox Proportional Hazards Survival Analysis Pipeline

## Overview
This pipeline performs **Cox proportional hazards regression** for genome-wide survival analysis, identifying genes whose expression levels predict patient survival while adjusting for clinical covariates. This complements the supervised limma analyses by providing hazard ratios and handling censored survival data appropriately.

---

## 🔴 CRITICAL DIFFERENCES FROM LIMMA ANALYSES (Analyses #1-3)

### **Cox vs Limma: Fundamental Paradigm Differences**

| Aspect | Limma (Analyses #1-3) | Cox PH (Current) |
|--------|----------------------|------------------|
| **Outcome** | Binary groups (low/high survival) | **Continuous time-to-event** |
| **Censoring** | Excluded (alive ≤ threshold dropped) | **Properly handled** (uses all data) |
| **Statistic** | t-test / F-test (mean difference) | **Hazard ratio** (relative risk) |
| **Cutoff** | Arbitrary (24, 36, 48 months) | **None** (uses full time range) |
| **Covariates** | Adjusted in linear model | **Adjusted in Cox model** |
| **Interpretation** | "Gene differs between groups" | **"Gene predicts survival risk"** |
| **Clinical relevance** | Group comparison | **Individualized risk prediction** |

---

## 🎯 WHY USE COX REGRESSION?

### **1. Proper Handling of Censored Data**

**Censoring:** Patient alive at last follow-up (we don't know final survival time)

**Limma approach:**
```python
# Patients alive with OS ≤ 24 months are EXCLUDED
mask = (df['OS'] <= 24) & (df['vital_status'] == 'alive')
df = df[~mask]  # Loses information!
```

**Cox approach:**
```python
# ALL patients included, censoring status encoded
event = 1 if dead, 0 if alive (censored)
# Cox model uses partial likelihood accounting for censoring
```

**Impact:** Cox uses **100% of data**, limma loses 10-30% (censored patients)

---

### **2. No Arbitrary Survival Cutoff**

**Limma:**
- Requires choosing threshold (24, 36, 48 months)
- Different cutoffs → different genes
- Results depend on arbitrary decision

**Cox:**
- Uses **full survival time** (continuous)
- No cutoff needed
- More objective, less researcher degrees of freedom

---

### **3. Hazard Ratio Interpretation**

**Limma output:**
```
Gene X: logFC = 0.8, p = 0.001
Interpretation: Gene X is 0.8 log2-units higher in high-survival group
```

**Cox output:**
```
Gene X: HR = 1.5, p = 0.001
Interpretation: One unit increase in Gene X → 1.5× death risk
```

**Hazard ratio (HR):**
- HR = 1: No effect on survival
- HR > 1: **Higher expression = higher death risk** (poor prognosis marker)
- HR < 1: **Higher expression = lower death risk** (protective marker)
- HR = 2: Doubling of death risk

**More clinically interpretable than logFC**

---

## CORE FUNCTIONS

### `create_event_column(df, status_col='case.vital.status')`

**Purpose:** Convert vital status to binary event indicator for Cox model

**Input formats handled:**
1. **Numeric:** 0=alive, 1=dead (already coded)
2. **String:** 'alive'/'dead', 'LIVING'/'DECEASED', etc.

**Transformation:**
```python
# String → Numeric
'alive' / 'LIVING' / 0     → 0 (censored)
'dead' / 'DECEASED' / 1    → 1 (event occurred)
```

**Output:**
- Modified DataFrame with new column `'event'` (0/1)
- Column name `'event'` (for consistency)

**Critical:** Event=1 must be the outcome of interest (death), not censoring

---

### `run_cox_for_all_genes(df, time_col, event_col, covariates)` ⭐ **CORE FUNCTION**

**Purpose:** Genome-wide Cox regression adjusting for clinical covariates

**Model specification:**
```python
# For each gene:
hazard(t) = h₀(t) × exp(β₁×Gene + β₂×Age + β₃×Sex)
```

**Where:**
- `h₀(t)`: Baseline hazard (unspecified)
- `β₁`: Log hazard ratio for gene
- `exp(β₁)`: Hazard ratio (HR)
- Covariates: Age, Sex (adjust for confounding)

**Pipeline:**
```python
for each gene:
    1. Build DataFrame: [time, event, gene, Age, Sex]
    2. Drop missing values
    3. Fit Cox model: cph.fit(df, duration_col, event_col)
    4. Extract gene coefficient (β) and statistics
    5. Compute HR = exp(β)
    6. Store results
```

**Output DataFrame columns:**
- `gene`: Gene symbol
- `coef`: β coefficient (log hazard ratio)
- `p`: Wald test p-value
- `exp(coef)`: **Hazard ratio (HR)** = exp(β)
- `coef lower 95%`: 95% CI lower bound for β
- `coef upper 95%`: 95% CI upper bound for β
- `z`: z-statistic (β / SE)
- `gene_n`: Number of patients included (after dropna)

---

#### **Statistical Details:**

**Wald test:**
```
H₀: β = 0 (gene has no effect on survival)
H₁: β ≠ 0 (gene affects survival)

z = β / SE(β)
p = 2 × P(|Z| > |z|)  # Two-tailed
```

**Hazard ratio confidence interval:**
```
HR_lower = exp(β - 1.96 × SE(β))
HR_upper = exp(β + 1.96 × SE(β))
```

**Interpretation:**
- If CI excludes 1 → Significant effect
- CI includes 1 → No significant effect

---

#### **Covariate Adjustment:**

**Why adjust for Age and Sex?**

**Without adjustment:**
```
Gene X associated with survival → Could be confounded by age
(Gene X high in older patients → appears prognostic but it's just age)
```

**With adjustment:**
```
Gene X associated with survival AFTER controlling for age/sex
→ Independent prognostic value
```

**Model with covariates:**
```python
# Gene effect (β₁) is adjusted for Age (β₂) and Sex (β₃)
hazard = h₀(t) × exp(β₁×Gene + β₂×Age + β₃×Sex)
```

**Critical:** `coef` and `p` for the gene are **conditional on covariates**

---

### `fdr_correction(results_df, p_col='p', alpha=0.10)`

**Purpose:** Multiple testing correction via Benjamini-Hochberg FDR

**Why needed:**
- Testing 10,000-20,000 genes simultaneously
- Expected false positives = 0.05 × 20,000 = 1,000 at α=0.05
- FDR controls false discovery rate

**Method:** Benjamini-Hochberg procedure
```python
from statsmodels.stats.multitest import multipletests
corrected = multipletests(p_values, alpha=0.10, method='fdr_bh')
```

**Output:**
- `p_adj`: FDR-adjusted p-value
- `significant`: Boolean (True if p_adj < alpha)

**Typical thresholds:**
- α = 0.05: Stringent (5% FDR)
- α = 0.10: Standard (10% FDR)
- α = 0.20: Exploratory (20% FDR)

**Interpretation:**
- FDR = 0.10 → Among significant genes, expect ~10% false positives
- More lenient than Bonferroni (controls family-wise error rate)

---

### `plot_km_for_gene(df, gene, time_col, event_col, outpath)`

**Purpose:** Kaplan-Meier survival curve stratified by gene expression

**Method:**
1. **Dichotomize gene expression:**
   ```python
   median_expr = df[gene].median()
   low_group = df[df[gene] <= median_expr]
   high_group = df[df[gene] > median_expr]
   ```

2. **Fit KM curves separately:**
   ```python
   kmf_low.fit(durations=low_group['time'], 
               event_observed=low_group['event'])
   kmf_high.fit(durations=high_group['time'], 
                event_observed=high_group['event'])
   ```

3. **Plot both curves:**
   - X-axis: Time (months)
   - Y-axis: Survival probability S(t)
   - Two lines: Low expression vs High expression

**Interpretation:**

**Poor prognosis marker (HR > 1):**
```
High expression curve BELOW low expression curve
→ High expression = worse survival
```

**Protective marker (HR < 1):**
```
High expression curve ABOVE low expression curve
→ High expression = better survival
```

**No effect (HR ≈ 1):**
```
Curves overlap or cross
→ No clear survival difference
```

**Note:** KM plot uses median split (Cox model uses continuous expression)

---

### `analyze_dataset(dataset_csv, dataset_name, outdir, alpha_fdr)`

**Purpose:** Complete Cox analysis pipeline for one dataset

**Workflow:**
```
1. Load CSV
   ↓
2. Create event column (vital status → 0/1)
   ↓
3. Encode Sex (Male=1, Female=0)
   ↓
4. Run Cox for all genes (with Age, Sex covariates)
   ↓
5. Apply FDR correction (p_adj < alpha)
   ↓
6. Save results:
   - All genes: {dataset}_cox_results.csv
   - Significant: {dataset}_cox_significant_genes.csv
   ↓
7. Generate KM plots for significant genes
   - One plot per gene → {dataset}_KM_plots/{gene}_KM.png
```

**Output files:**
- `{dataset}_cox_results.csv`: All genes, sorted by p_adj
- `{dataset}_cox_significant_genes.csv`: FDR < alpha only
- `{dataset}_KM_plots/`: Directory with KM plots

---

### `summarize_threshold_overlaps` (Cross-validation function)

**Purpose:** Compare Cox-significant genes to limma-significant genes

**Question:** Do genes found by Cox also appear in limma analyses?

**Method:**
1. Load Cox significant genes (FDR < threshold)
2. For each limma threshold (24, 36, 48, 60, 24-48, 24-60):
   - Load limma significant genes
   - Compute intersection with Cox genes
3. Summarize: n_limma, n_cox, n_common

**Output DataFrame:**
```
consortium | threshold_label | n_limma | n_cox | n_common
TCGA       | 24-months      | 150     | 80    | 45
TCGA       | 36-months      | 120     | 80    | 52
TCGA       | 48-months      | 95      | 80    | 38
...
```

**Interpretation:**

**High overlap (>60%):**
- Cox and limma agree → Robust survival-associated genes
- Threshold-independent effects

**Moderate overlap (30-60%):**
- Some agreement
- Threshold-dependent effects

**Low overlap (<30%):**
- Cox finds different genes than limma
- Censoring handling or cutoff choice matters

**Genes in BOTH:**
- Highest confidence
- Validated by two statistical methods

---

## MAIN PIPELINE WORKFLOW

### **Phase 1: Cox Analysis (per dataset)**

```python
# TCGA
analyze_dataset(
    dataset_csv="TCGA_processed.csv",
    dataset_name="TCGA",
    outdir="cox_results",
    alpha_fdr=0.10  # 10% FDR
)

# CGGA (relaxed threshold - fewer samples)
analyze_dataset(
    dataset_csv="CGGA_processed.csv",
    dataset_name="CGGA",
    outdir="cox_results",
    alpha_fdr=0.50  # 50% FDR (exploratory)
)

# GLASS (relaxed threshold - fewer samples)
analyze_dataset(
    dataset_csv="GLASS_processed.csv",
    dataset_name="GLASS",
    outdir="cox_results",
    alpha_fdr=0.40  # 40% FDR
)
```

**Note:** Different alpha_fdr for each dataset
- TCGA: Large N → stringent (0.10)
- CGGA: Medium N → moderate (0.50)
- GLASS: Small N → relaxed (0.40)

**Rationale:** Smaller datasets have less power → need relaxed threshold to find anything

---

### **Phase 2: Cross-Method Comparison**

```python
# Compare Cox vs limma for multiple thresholds
for consortium in ["TCGA", "CGGA", "GLASS"]:
    summary = summarize_threshold_overlaps(
        consortium=consortium,
        cox_sig_genes_file=f"{consortium}_cox_significant_genes.csv",
        outdir="results",
        thresholds_fixed=[24, 36, 48, 60],
        thresholds_range=[(24,48), (24,60)]
    )
```

**Output:** `cox_limma_summary.csv`

---

## OUTPUT FILES

### **Per-Dataset Files:**

**1. `{dataset}_cox_results.csv`** (All genes)
```
gene     | coef  | p      | exp(coef) | coef lower 95% | coef upper 95% | z     | gene_n | p_adj | significant
EGFR     | 0.42  | 0.001  | 1.52      | 0.15           | 0.69           | 3.21  | 98     | 0.08  | True
TP53     | -0.18 | 0.12   | 0.84      | -0.45          | 0.09           | -1.55 | 98     | 0.45  | False
...
```

**2. `{dataset}_cox_significant_genes.csv`** (FDR < threshold only)
- Subset of above, filtered for significance
- Sorted by p_adj (ascending)

**3. `{dataset}_KM_plots/` directory**
- One PNG per significant gene
- Filename: `{gene}_KM.png`
- Kaplan-Meier curves (high vs low expression)

---

### **Summary File:**

**`cox_limma_summary.csv`** (Cross-method comparison)
```
consortium | threshold_label | n_limma | n_cox | n_common
TCGA       | 24-months       | 150     | 80    | 45
TCGA       | 36-months       | 120     | 80    | 52
CGGA       | 24-months       | 95      | 120   | 38
GLASS      | 24-48-months    | 62      | 75    | 28
...
```

**Columns:**
- `n_limma`: Number of significant genes in limma (that threshold)
- `n_cox`: Number of significant genes in Cox
- `n_common`: Intersection (genes in both)

---

## INTERPRETATION STRATEGY

### **Step 1: Examine Cox Results per Dataset**

**Load results:**
```python
tcga_cox = pd.read_csv("TCGA_cox_results.csv")
tcga_sig = tcga_cox[tcga_cox['significant']]
```

**Key questions:**
1. How many significant genes? (n_sig)
2. What's the median HR? (risk vs protective)
3. Are top genes biologically plausible?

**Typical expectations:**
- TCGA: 50-200 significant genes (large N, stringent FDR)
- CGGA: 100-300 genes (medium N, relaxed FDR)
- GLASS: 50-150 genes (small N, relaxed FDR)

---

### **Step 2: Interpret Hazard Ratios**

**For each significant gene:**

**Example 1: HR = 1.8, p_adj = 0.003**
```
Interpretation: 
- 1-unit increase in expression → 1.8× death risk
- Poor prognosis marker
- Strong statistical evidence (p_adj < 0.01)

Clinical implication:
- High expression = high risk
- Could be therapeutic target (inhibit gene)
```

**Example 2: HR = 0.65, p_adj = 0.04**
```
Interpretation:
- 1-unit increase in expression → 0.65× death risk
- Protective marker (35% risk reduction)
- Moderate evidence (p_adj = 0.04)

Clinical implication:
- High expression = low risk
- Could be positive prognostic marker
- Pathway might be therapeutic target (activate gene)
```

**Magnitude guidelines:**
- HR = 1.0-1.2 or 0.8-1.0: Weak effect
- HR = 1.2-1.5 or 0.65-0.8: Moderate effect
- HR > 1.5 or < 0.65: Strong effect
- HR > 2.0 or < 0.5: Very strong effect

---

### **Step 3: Visual Validation (KM Plots)**

**For top genes, inspect KM plots:**

**Good separation:**
```
High expression curve clearly below low expression curve
→ Validates HR > 1 (poor prognosis)
```

**Poor separation:**
```
Curves overlap or cross
→ HR may be statistically significant but clinically weak
```

**Recommendation:** Prioritize genes with:
- Significant p_adj (< 0.10)
- Moderate-to-strong HR (> 1.3 or < 0.75)
- Clear KM separation

---

### **Step 4: Cross-Method Validation**

**Load summary:**
```python
summary = pd.read_csv("cox_limma_summary.csv")
```

**For each dataset, examine overlap:**

**High overlap example:**
```
TCGA | 24-months | n_limma=150 | n_cox=80 | n_common=60
→ 60/80 = 75% of Cox genes also in limma ✅
→ 60/150 = 40% of limma genes also in Cox
```

**Interpretation:** 
- Cox genes are **robust** (75% replicate in limma)
- Limma finds more genes (possibly threshold-dependent)
- Focus on **n_common genes** (validated by both methods)

**Low overlap example:**
```
GLASS | 24-months | n_limma=62 | n_cox=75 | n_common=15
→ 15/75 = 20% of Cox genes in limma ❌
→ 15/62 = 24% of limma genes in Cox
```

**Interpretation:**
- Methods find different genes (concerning)
- Possible reasons:
  - Censoring handling (Cox includes censored, limma excludes)
  - Statistical power (small N in GLASS)
  - Biological heterogeneity

---

### **Step 5: Cross-Cohort Replication**

**Most robust genes:** Significant in Cox across multiple cohorts

**Example:**
```python
tcga_genes = set(tcga_sig['gene'])
cgga_genes = set(cgga_sig['gene'])
glass_genes = set(glass_sig['gene'])

# Triple overlap
robust_genes = tcga_genes & cgga_genes & glass_genes
print(f"Genes significant in all 3 cohorts: {len(robust_genes)}")
```

**Gold standard:** Genes in ≥2 cohorts with:
- Consistent HR direction (all > 1 or all < 1)
- p_adj < 0.10 in at least one cohort
- Also appear in limma analyses

---

## CRITICAL INTERPRETATION POINTS

### 1. **Proportional Hazards Assumption**

**Cox model assumes:** Hazard ratio is **constant over time**

**Assumption:**
```
HR = 1.5 at 6 months
HR = 1.5 at 24 months
HR = 1.5 at 48 months
```

**Violation:**
```
HR = 2.0 at 6 months  (strong early effect)
HR = 1.0 at 48 months (no late effect)
```

**Test:** Schoenfeld residuals (not implemented in current code)

**Practical impact:**
- If violated: HR estimate is **averaged over time**
- May miss time-varying effects
- Alternative: Time-dependent Cox model

**Recommendation:** For top genes, consider time-stratified analysis

---

### 2. **Censoring Assumptions**

**Cox model assumes:** Censoring is **non-informative**

**Non-informative censoring:**
- Patients censored for administrative reasons (end of study)
- Random loss to follow-up
- Censoring unrelated to prognosis

**Informative censoring (violates assumption):**
- Patients with poor prognosis leave study early (too sick)
- Patients with good prognosis censored more (healthier, mobile)

**Impact:** Biased HR estimates

**Mitigation:** Check censoring patterns across expression levels

---

### 3. **Covariate Adjustment**

**Current model:** Adjusts for Age and Sex only

**Missing confounders:**
- Treatment (radio, chemo)
- Performance status (KPS)
- Tumor size
- Extent of resection

**Impact:** Omitted variable bias
- Gene may appear prognostic but it's really a surrogate for treatment
- HR estimates are confounded

**Recommendation:** 
- For TCGA_care / CGGA_care: Add treatment variables if available
- Sensitivity analysis: Re-run Cox with/without treatment covariates

---

### 4. **Multiple Testing Trade-offs**

**Stringent FDR (α = 0.05):**
- Fewer false positives
- More false negatives (misses true signals)
- Conservative

**Relaxed FDR (α = 0.40):**
- More true positives
- More false positives (need validation)
- Exploratory

**Current choices:**
- TCGA: α = 0.10 (balanced)
- CGGA: α = 0.50 (exploratory)
- GLASS: α = 0.40 (exploratory)

**Recommendation:** 
- Use stringent (0.05-0.10) for biomarker discovery
- Use relaxed (0.20-0.50) for hypothesis generation
- **Always validate in independent cohort**

---

### 5. **Cox vs Limma Discordance**

**Why might they disagree?**

**1. Censoring:**
- Cox uses all patients
- Limma excludes censored (alive ≤ threshold)
- Cox may find genes associated with early/late survival
- Limma only captures differences at specific cutoff

**2. Statistical model:**
- Cox: Continuous hazard
- Limma: Mean difference between groups
- Different assumptions, different power

**3. Cutoff choice (limma):**
- Limma at 24 months ≠ 36 months ≠ 48 months
- Cox is cutoff-free

**Recommendation:** Prioritize genes in **both** Cox and limma

---

### 6. **Sample Size Dependency**

**TCGA (N ≈ 300-400):**
- High power
- Can detect small effects (HR = 1.2)
- Stringent FDR feasible

**CGGA (N ≈ 150-200):**
- Moderate power
- Detects moderate effects (HR = 1.3-1.5)
- Relaxed FDR needed

**GLASS (N ≈ 100-120):**
- Low power
- Only detects strong effects (HR > 1.5)
- Very relaxed FDR needed

**Implication:** Genes found in GLASS (small N) are likely **strong effects**

---

### 7. **Hazard Ratio Scale Dependency**

**Gene expression scale matters:**

**Z-scored expression (mean=0, std=1):**
```
HR = 1.5 per 1 SD increase
Interpretation: 1 standard deviation ↑ → 1.5× risk
```

**Log2 expression (typical range 4-14):**
```
HR = 1.1 per 1 unit increase
Interpretation: Doubling expression → 1.1^1 = 1.1× risk
```

**Current analysis:** Uses normalized expression from Analyses #1-3
- TCGA: Z-score (per platform) → HR per SD
- CGGA: Z-score (per batch) → HR per SD
- GLASS: Z-score → HR per SD

**Recommendation:** For interpretation, convert HR to **per SD** units

---

## COMPARISON TO LIMMA ANALYSES

### **Complementary Strengths:**

| Question | Limma | Cox |
|----------|-------|-----|
| **Uses censored data** | ❌ Excluded | ✅ Included |
| **Cutoff-free** | ❌ Requires threshold | ✅ Continuous time |
| **Hazard ratio** | ❌ Only logFC | ✅ Direct HR |
| **Clinical interpretation** | Moderate | ✅ Strong (risk) |
| **Statistical power** | High (two groups) | Moderate (continuous) |
| **Batch correction** | ✅ In linear model | ⚠️ As covariate only |
| **Multi-group** | ✅ F-test (K groups) | ❌ Pairwise only |

---

### **When to Use Which:**

**Use Limma when:**
- Clear clinical grouping (e.g., responders vs non-responders)
- Want interpretable fold-changes (up/down)
- Multi-group comparisons (K > 2)
- Already have censored patients excluded

**Use Cox when:**
- Want to use ALL data (including censored)
- Avoid arbitrary cutoffs
- Need hazard ratios for risk prediction
- Clinical endpoint is time-to-event

**Use BOTH when:**
- Maximum rigor for biomarker discovery
- Genes in **both** are highest confidence
- Cross-validate findings

---

## RECOMMENDED ANALYSIS SEQUENCE

### **Step 1: Run Cox for all datasets**
```python
for dataset in ["TCGA", "CGGA", "GLASS"]:
    analyze_dataset(f"{dataset}_processed.csv", dataset, "cox_results", alpha_fdr)
```

### **Step 2: Examine significant genes**
```python
for dataset in ["TCGA", "CGGA", "GLASS"]:
    sig_genes = pd.read_csv(f"{dataset}_cox_significant_genes.csv")
    print(f"{dataset}: {len(sig_genes)} significant genes")
    print(sig_genes.head(10))  # Top 10
```

### **Step 3: Cross-cohort overlap**
```python
tcga_genes = set(tcga_sig['gene'])
cgga_genes = set(cgga_sig['gene'])
glass_genes = set(glass_sig['gene'])

# Pairwise overlaps
print("TCGA ∩ CGGA:", len(tcga_genes & cgga_genes))
print("TCGA ∩ GLASS:", len(tcga_genes & glass_genes))
print("CGGA ∩ GLASS:", len(cgga_genes & glass_genes))

# Triple overlap
print("All 3 cohorts:", len(tcga_genes & cgga_genes & glass_genes))
```

### **Step 4: Cox-Limma comparison**
```python
summary = summarize_threshold_overlaps(...)
# Examine overlaps
```

### **Step 5: Prioritize genes**
**Criteria for highest-confidence genes:**
1. ✅ Significant in Cox (p_adj < 0.10) in ≥2 cohorts
2. ✅ Consistent HR direction across cohorts
3. ✅ Also significant in limma (≥1 threshold)
4. ✅ Moderate-to-strong HR (>1.3 or <0.75)
5. ✅ Clear KM separation

### **Step 6: Biological interpretation**
- Pathway enrichment on Cox-significant genes
- Literature review of top genes
- Compare to known GBM markers (EGFR, TP53, IDH1, etc.)

---

## OUTPUT INTERPRETATION EXAMPLES

### **Example 1: Successful Cross-Validation**

**TCGA Cox:**
```
Gene: EGFR
HR = 1.85, p_adj = 0.001
Interpretation: High EGFR → 1.85× death risk
```

**TCGA Limma (24 months):**
```
Gene: EGFR
logFC = 0.92, p_adj = 0.003
Interpretation: EGFR higher in low-survival group
```

**CGGA Cox:**
```
Gene: EGFR
HR = 1.62, p_adj = 0.008
```

**Conclusion:** EGFR is robust poor prognosis marker ✅
- Validated by Cox (2 cohorts)
- Validated by limma
- Consistent HR direction (>1)

---

### **Example 2: Discordant Results**

**GLASS Cox:**
```
Gene: GENEXY
HR = 0.55, p_adj = 0.02
Interpretation: Protective marker
```

**GLASS Limma (24 months):**
```
Gene: GENEXY
Not significant (p_adj = 0.45)
```

**Possible explanations:**
1. Cox detects continuous effect (limma misses with dichotomization)
2. Effect is early/late (Cox captures, limma misses with 24mo cutoff)
3. Censoring matters (Cox uses censored patients, limma excludes)
4. False positive in Cox (p_adj = 0.02 is borderline with α=0.40)

**Next steps:**
- Check if significant in TCGA/CGGA Cox
- Try different limma thresholds (36, 48 months)
- Examine KM plot (is separation clear?)

---

## KEY TAKEAWAYS

🎯 **Main advantage:** Uses **all survival data** (no censoring exclusions, no arbitrary cutoffs)

✅ **Use Cox when:**
- Want to maximize data utilization
- Avoid threshold selection bias
- Need hazard ratios for clinical interpretation
- Planning risk prediction models

🔴 **Limitations:**
- Assumes proportional hazards (may not hold)
- Requires non-informative censoring
- Less power for very small effects than limma

🌟 **Gold standard genes:**
- Significant in Cox (p_adj < 0.10) in ≥2 cohorts
- Consistent HR direction (all >1 or all <1)
- Also significant in limma (≥1 threshold)
- Moderate-to-strong effect size (HR >1.3 or <0.75)

💡 **Complementary to limma:**
- Cox and limma address slightly different questions
- Genes in **both** are most robust
- Use Cox for risk prediction, limma for mechanistic interpretation

---

## Future Enhancements

**High priority:**
1. Test proportional hazards assumption (Schoenfeld residuals)
2. Add treatment covariates (radio, chemo)
3. Time-dependent Cox (allow HR to vary over time)
4. Stratified Cox (separate baseline hazards per stratum)

**Medium priority:**
5. ROC curves for risk prediction (continuous expression → risk score)
6. Nomograms (combine gene + clinical variables)
7. C-index (concordance) for model performance
8. Bootstrap confidence intervals

**Nice to have:**
9. Multivariate Cox (multiple genes simultaneously)
10. Penalized Cox (LASSO, elastic net) for high-dimensional data
11. Competing risks (death vs other causes)
12. Landmark analysis (conditional survival)
