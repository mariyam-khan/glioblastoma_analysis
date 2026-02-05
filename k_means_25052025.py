
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations
from pathlib import Path
from scipy.stats import kruskal
from sklearn.cluster import KMeans

# rpy2 imports
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
pandas2ri.activate()

import warnings
warnings.filterwarnings(
    "ignore",
    message="The default value of `n_init` will change"
)
BASE_DIR = Path("/home/mkh062/Desktop/scratch/TCGA_project/processed_data/glass_clustering")
BASE_DIR.mkdir(parents=True, exist_ok=True)
TCGA_DIR    = BASE_DIR   # adjust if different
# Dataset names
tcga_sets = [
    "TCGA-zscore_non-care",
    "TCGA-raw_non-care",
    "TCGA-qn_non-care",
    "TCGA-zscore_care",
    "TCGA-raw_care",
    "TCGA-qn_care",
]

#OUTPUT_XLSX = BASE_DIR / "GLASS_cluster_DE_results.xlsx"
#-------------- your existing limma wrapper --------------
PADJ_CUTOFF = 0.1


def run_limma_multi(expr: pd.DataFrame,
                    clusters: pd.Series,
                    covars:   pd.DataFrame,
                    label:    str) -> pd.DataFrame:
    """
    Two‐group limma (via rpy2) for expr (genes×samples) contrasting
    grp == '1' vs grp == '0'. clusters must be 0/1.
    Returns topTable with logFC, P.Value, adj.P.Val, etc.
    """
    # prepare R
    r_expr = pandas2ri.py2rpy(expr)
    ro.globalenv['expr'] = r_expr

    # build grp factor in R: levels '0','1'
    grp = clusters.astype(int).astype(str)
    r_grp = ro.FactorVector(grp.values)
    ro.globalenv['grp'] = r_grp

    # design matrix: ~ 0 + grp (+ covars)
    design_cmds = ["design <- model.matrix(~ 0 + grp)"]
    if covars is not None:
        for cov in covars.columns:
            ro.globalenv[cov] = pandas2ri.py2rpy(covars[cov])
            design_cmds.append(f"design <- cbind(design, {cov})")

    ro.r("""
        suppressPackageStartupMessages(library(limma))
        {designs}
        fit <- lmFit(expr, design)
        fit <- eBayes(fit)
        # define contrast grp1 - grp0
        cont <- makeContrasts(grp1 - grp0, levels=design)
        fit2 <- contrasts.fit(fit, cont)
        fit2 <- eBayes(fit2)
        tt <- topTable(fit2, number=Inf, sort.by='P')
    """.format(designs="\n".join(design_cmds)))
    tt = ro.r("tt")
    tt.index.name = "Gene"
    tt.reset_index(inplace=True)
    return tt

#-------------- utilities --------------

def volcano_plot(tt: pd.DataFrame, out_png: Path, title: str):
    """Basic volcano: logFC vs -log10(adj.P.Val)"""
    plt.figure(figsize=(6,5))
    x = tt['logFC']
    y = -np.log10(tt['adj.P.Val'])
    plt.scatter(x, y, alpha=0.5, s=10)
    sig = tt['adj.P.Val'] < PADJ_CUTOFF
    plt.scatter(x[sig], y[sig], color='red', alpha=0.7, s=10)
    plt.axhline(-np.log10(PADJ_CUTOFF), color='grey', lw=1, ls='--')
    plt.xlabel('log₂ fold-change')
    plt.ylabel('-log₁₀ adj.P.Val')
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_png)
    plt.close()


def load_glass_data(filepath):
    """
    Load GLASS data, separate metadata and gene expression columns.
    """
    df = pd.read_csv(filepath)
    metadata_cols = ['Sample.ID', 'Overall.Survival.Months', 'Diagnosis.Age', 'Sex', 'case.vital.status']
    gene_cols = [col for col in df.columns if col not in metadata_cols]
    return df, gene_cols


def kmeans_cluster_samples(df, gene_cols, n_clusters=4):
    """
    Perform KMeans clustering on gene expression columns, assign cluster labels.
    """
    from sklearn.cluster import KMeans
    sample_data = df[gene_cols].values
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    df['Cluster'] = kmeans.fit_predict(sample_data)
    return df, kmeans

# Main analysis
df_glass = pd.read_csv(INPUT_FILE)
df_glass['Sex'] = df_glass['Sex'].map({'female':0,'male':1})
meta_glass = ['Sample.ID','Overall.Survival.Months','Diagnosis.Age','Sex','case.vital.status']
genes_glass = [c for c in df_glass if c not in meta_glass]

# --- Loop over each TCGA dataset ---
for ds in tcga_sets:
    tcga_file = TCGA_DIR / f"{ds}.csv"
    print("##############  tcga dataset ", ds)
    if not tcga_file.exists():
        print(f"⚠️  File not found: {tcga_file}")
        continue

    df_tcga = pd.read_csv(tcga_file)
    # drop platform if present
    if 'platform' in df_tcga.columns:
        df_tcga = df_tcga.drop(columns=['platform'])
    # map sex to numeric
    df_tcga['Sex'] = df_tcga['Sex'].map({'Female': 0, 'Male': 1})
    # define metadata vs. gene columns
    meta_tcga = ['Sample.ID','Overall.Survival.Months','Diagnosis.Age','Sex','case.vital.status']
    genes_tcga = [c for c in df_tcga.columns if c not in meta_tcga]

    print(f"\n{ds}: {len(df_tcga)} samples, {len(genes_tcga)} genes")

    #writer = pd.ExcelWriter(OUTPUT_XLSX, engine='openpyxl')

    #########################################################
    km = KMeans(n_clusters=3, random_state=42, n_init=10)
    df_glass['Cluster'] = km.fit_predict(df_glass[genes_glass].values)

    # — 1) compute centroids & log2(fold-change) —
    centroids_glass = df_glass.groupby('Cluster')[genes_glass].mean()
    fc_list_glass = []
    for c in centroids_glass.index:
        μ_c     = centroids_glass.loc[c]
        μ_other = centroids_glass.drop(c).mean(axis=0)
        log2fc  = np.log2((μ_c + 1e-9) / (μ_other + 1e-9))
        tmp = pd.DataFrame({
            'cluster': c,
            'log2FC':  log2fc.values
        })
        fc_list_glass.append(tmp)
    fc_df_glass = pd.concat(fc_list_glass, ignore_index=True)

    # — Compute mean/median stats —
    fc_stats_glass = fc_df_glass.groupby('cluster')['log2FC'].agg(['mean','median'])
    sv_stats = df_glass.groupby('Cluster')['Overall.Survival.Months'].agg(['mean','median'])
    print("     fc_stats_glass  ", fc_stats_glass)
    print("  sv_stats   ", sv_stats)
    # prepare colors
    cmap = {0:'C0', 1:'C1', 2:'C2'}

    # — 2) plotting with 2-line legend entries —
    fig, axes = plt.subplots(2, 2, figsize=(12,10))
    ax = axes[0,0]
    # 2a) log2FC scatter, legend includes μ & md on 2 lines
    for c, grp in fc_df_glass.groupby('cluster'):
        μ, md = fc_stats_glass.loc[c]
        label = f"Cluster {c}\nμ={μ:.2f}, md={md:.2f}"
        jitter = np.random.normal(scale=0.05, size=len(grp))
        x = grp['cluster'] + jitter
        ax.scatter(x, grp['log2FC'], s=5, c=cmap[c], alpha=0.6, label=label)

    ax.axhline(0, color='gray', linestyle='--', linewidth=1)
    ax.set_xlabel('Cluster')
    ax.set_ylabel('log$_2$ fold‐change')
    ax.set_title('Gene log$_2$FC vs. rest, by cluster GLASS')
    ax.set_xticks(fc_stats_glass.index)
    ax.legend(bbox_to_anchor=(1.05,1), loc='upper left')
    ax = axes[0,1]
    # 2b) overall survival scatter, legend includes μ & md on 2 lines
    for c, grp in df_glass.groupby('Cluster'):
        μ, md = sv_stats.loc[c]
        label = f"Cluster {c}\nμ={μ:.1f}, md={md:.1f}"
        jitter2 = np.random.normal(scale=0.05, size=len(grp))
        x2 = grp['Cluster'] + jitter2
        ax.scatter(x2, grp['Overall.Survival.Months'], s=20, c=cmap[c], alpha=0.7, label=label)

    ax.set_xlabel('Cluster')
    ax.set_ylabel('Overall Survival (months)')
    ax.set_title('Survival by cluster GLASS')
    ax.set_xticks(sv_stats.index)
    ax.legend(bbox_to_anchor=(1.05,1), loc='upper left')
    #???????????????????????????????????????????????????????????????????
    df_tcga['Cluster'] = km.fit_predict(df_tcga[genes_tcga].values)
    # — 1) compute centroids & log2(fold-change) —
    centroids_tcga = df_tcga.groupby('Cluster')[genes_tcga].mean()
    fc_list_tcga = []
    for c in centroids_tcga.index:
        μ_c     = centroids_tcga.loc[c]
        μ_other = centroids_tcga.drop(c).mean(axis=0)
        log2fc  = np.log2((μ_c + 1e-9) / (μ_other + 1e-9))
        tmp = pd.DataFrame({
            'cluster': c,
            'log2FC':  log2fc.values
        })
        fc_list_tcga.append(tmp)
    fc_df_tcga = pd.concat(fc_list_tcga, ignore_index=True)


    swap_map = {0: 2,
                1: 0,
                2: 1}

    # 1) remap in the main dataframe
    df_tcga['Cluster'] = df_tcga['Cluster'].map(swap_map)

    # 2) remap in the fold‐change dataframe
    fc_df_tcga['cluster'] = fc_df_tcga['cluster'].map(swap_map)

    # — Compute mean/median stats —
    fc_stats_tcga = fc_df_tcga.groupby('cluster')['log2FC'].agg(['mean','median'])
    sv_stats = df_tcga.groupby('Cluster')['Overall.Survival.Months'].agg(['mean','median'])
    print("     fc_stats_tcga  ", fc_stats_tcga)
    print("  sv_stats   ", sv_stats)
    # prepare colors
    cmap = {0:'C0', 1:'C1', 2:'C2'}

    ax = axes[1,0]
    # 2a) log2FC scatter, legend includes μ & md on 2 lines
    for c, grp in fc_df_tcga.groupby('cluster'):
        μ, md = fc_stats_tcga.loc[c]
        label = f"Cluster {c}\nμ={μ:.2f}, md={md:.2f}"
        jitter = np.random.normal(scale=0.05, size=len(grp))
        x = grp['cluster'] + jitter
        ax.scatter(x, grp['log2FC'], s=5, c=cmap[c], alpha=0.6, label=label)

    ax.axhline(0, color='gray', linestyle='--', linewidth=1)
    ax.set_xlabel('Cluster')
    ax.set_ylabel('log$_2$ fold‐change')
    ax.set_title(f"Gene log$_2$FC vs. rest, by cluster {ds}")
    ax.set_xticks(fc_stats_glass.index)
    ax.legend(bbox_to_anchor=(1.05,1), loc='upper left')
    ax = axes[1,1]
    # 2b) overall survival scatter, legend includes μ & md on 2 lines
    for c, grp in df_tcga.groupby('Cluster'):
        μ, md = sv_stats.loc[c]
        label = f"Cluster {c}\nμ={μ:.1f}, md={md:.1f}"
        jitter2 = np.random.normal(scale=0.05, size=len(grp))
        x2 = grp['Cluster'] + jitter2
        ax.scatter(x2, grp['Overall.Survival.Months'], s=20, c=cmap[c], alpha=0.7, label=label)

    ax.set_xlabel('Cluster')
    ax.set_ylabel('Overall Survival (months)')
    ax.set_title(f"Survival by cluster {ds}")
    ax.set_xticks(sv_stats.index)
    ax.legend(bbox_to_anchor=(1.05,1), loc='upper left')


    plt.tight_layout()
    plt.show()

    #########################################################
    km = KMeans(n_clusters=4, random_state=42, n_init=10)
    df_glass['Cluster'] = km.fit_predict(df_glass[genes_glass].values)

    # — 1) compute centroids & log2(fold-change) —
    centroids_glass = df_glass.groupby('Cluster')[genes_glass].mean()
    fc_list_glass = []
    for c in centroids_glass.index:
        μ_c     = centroids_glass.loc[c]
        μ_other = centroids_glass.drop(c).mean(axis=0)
        log2fc  = np.log2((μ_c + 1e-9) / (μ_other + 1e-9))
        tmp = pd.DataFrame({
            'cluster': c,
            'log2FC':  log2fc.values
        })
        fc_list_glass.append(tmp)
    fc_df_glass = pd.concat(fc_list_glass, ignore_index=True)

    # — Compute mean/median stats —
    fc_stats_glass = fc_df_glass.groupby('cluster')['log2FC'].agg(['mean','median'])
    sv_stats = df_glass.groupby('Cluster')['Overall.Survival.Months'].agg(['mean','median'])
    print("     fc_stats_glass  ", fc_stats_glass)
    print("  sv_stats   ", sv_stats)
    # prepare colors
    cmap = {0:'C0', 1:'C1', 2:'C2', 3:'C3'}

    # — 2) plotting with 2-line legend entries —
    fig, axes = plt.subplots(2, 2, figsize=(12,10))
    ax = axes[0,0]
    # 2a) log2FC scatter, legend includes μ & md on 2 lines
    for c, grp in fc_df_glass.groupby('cluster'):
        μ, md = fc_stats_glass.loc[c]
        label = f"Cluster {c}\nμ={μ:.2f}, md={md:.2f}"
        jitter = np.random.normal(scale=0.05, size=len(grp))
        x = grp['cluster'] + jitter
        ax.scatter(x, grp['log2FC'], s=5, c=cmap[c], alpha=0.6, label=label)

    ax.axhline(0, color='gray', linestyle='--', linewidth=1)
    ax.set_xlabel('Cluster')
    ax.set_ylabel('log$_2$ fold‐change')
    ax.set_title('Gene log$_2$FC vs. rest, by cluster')
    ax.set_xticks(fc_stats_glass.index)
    ax.legend(bbox_to_anchor=(1.05,1), loc='upper left')
    ax = axes[0,1]
    # 2b) overall survival scatter, legend includes μ & md on 2 lines
    for c, grp in df_glass.groupby('Cluster'):
        μ, md = sv_stats.loc[c]
        label = f"Cluster {c}\nμ={μ:.1f}, md={md:.1f}"
        jitter2 = np.random.normal(scale=0.05, size=len(grp))
        x2 = grp['Cluster'] + jitter2
        ax.scatter(x2, grp['Overall.Survival.Months'], s=20, c=cmap[c], alpha=0.7, label=label)

    ax.set_xlabel('Cluster')
    ax.set_ylabel('Overall Survival (months)')
    ax.set_title('Survival by cluster')
    ax.set_xticks(sv_stats.index)
    ax.legend(bbox_to_anchor=(1.05,1), loc='upper left')
    #???????????????????????????????????????????????????????????????????
    df_tcga['Cluster'] = km.fit_predict(df_tcga[genes_tcga].values)
    # — 1) compute centroids & log2(fold-change) —
    centroids_tcga = df_tcga.groupby('Cluster')[genes_tcga].mean()
    fc_list_tcga = []
    for c in centroids_tcga.index:
        μ_c     = centroids_tcga.loc[c]
        μ_other = centroids_tcga.drop(c).mean(axis=0)
        log2fc  = np.log2((μ_c + 1e-9) / (μ_other + 1e-9))
        tmp = pd.DataFrame({
            'cluster': c,
            'log2FC':  log2fc.values
        })
        fc_list_tcga.append(tmp)
    fc_df_tcga = pd.concat(fc_list_tcga, ignore_index=True)


    swap_map = {0: 1,
                1: 3,
                2: 0,
                3: 2}

    # 1) remap in the main dataframe
    df_tcga['Cluster'] = df_tcga['Cluster'].map(swap_map)

    # 2) remap in the fold‐change dataframe
    fc_df_tcga['cluster'] = fc_df_tcga['cluster'].map(swap_map)

    # — Compute mean/median stats —
    fc_stats_tcga = fc_df_tcga.groupby('cluster')['log2FC'].agg(['mean','median'])
    sv_stats = df_tcga.groupby('Cluster')['Overall.Survival.Months'].agg(['mean','median'])
    print("     fc_stats_tcga  ", fc_stats_tcga)
    print("  sv_stats   ", sv_stats)
    # prepare colors
    cmap = {0:'C0', 1:'C1', 2:'C2', 3:'C3'}

    ax = axes[1,0]
    # 2a) log2FC scatter, legend includes μ & md on 2 lines
    for c, grp in fc_df_tcga.groupby('cluster'):
        μ, md = fc_stats_tcga.loc[c]
        label = f"Cluster {c}\nμ={μ:.2f}, md={md:.2f}"
        jitter = np.random.normal(scale=0.05, size=len(grp))
        x = grp['cluster'] + jitter
        ax.scatter(x, grp['log2FC'], s=5, c=cmap[c], alpha=0.6, label=label)

    ax.axhline(0, color='gray', linestyle='--', linewidth=1)
    ax.set_xlabel('Cluster')
    ax.set_ylabel('log$_2$ fold‐change')
    ax.set_title('Gene log$_2$FC vs. rest, by cluster')
    ax.set_xticks(fc_stats_glass.index)
    ax.legend(bbox_to_anchor=(1.05,1), loc='upper left')
    ax = axes[1,1]
    # 2b) overall survival scatter, legend includes μ & md on 2 lines
    for c, grp in df_tcga.groupby('Cluster'):
        μ, md = sv_stats.loc[c]
        label = f"Cluster {c}\nμ={μ:.1f}, md={md:.1f}"
        jitter2 = np.random.normal(scale=0.05, size=len(grp))
        x2 = grp['Cluster'] + jitter2
        ax.scatter(x2, grp['Overall.Survival.Months'], s=20, c=cmap[c], alpha=0.7, label=label)

    ax.set_xlabel('Cluster')
    ax.set_ylabel('Overall Survival (months)')
    ax.set_title('Survival by cluster')
    ax.set_xticks(sv_stats.index)
    ax.legend(bbox_to_anchor=(1.05,1), loc='upper left')


    plt.tight_layout()
    plt.show()
