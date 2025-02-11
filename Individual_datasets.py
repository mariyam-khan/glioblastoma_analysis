import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")  # force non-GUI backend
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

# R-related imports
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
pandas2ri.activate()

###############################################################################
# 1. Helper functions
###############################################################################

def remove_zero_genes(df):
    """
    Removes genes (columns) that are all zeros or that have zeros in more than 50% of the rows.
    Assumes 'Sample ID' or 'Sample.ID' is the row-identifier column.
    """
    sample_col = 'Sample ID' if 'Sample ID' in df.columns else 'Sample.ID'
    # Ensure the sample column exists, otherwise raise an error
    if sample_col not in df.columns:
        raise ValueError("Neither 'Sample ID' nor 'Sample.ID' found in dataframe.")
    gene_cols = [c for c in df.columns if c != sample_col]
    plaform_col = 'platform'
    if plaform_col not in df.columns:
        gene_cols = [c for c in df.columns if c != sample_col and c != plaform_col]
    # Compute the sum of each gene column
    sum_mask = df[gene_cols].sum(axis=0) != 0
    # Compute the percentage of zeros in each column
    zero_percentage_mask = (df[gene_cols] == 0).mean(axis=0) <= 0.50
    # Find columns to keep by combining both conditions
    valid_columns = sum_mask.index.intersection(zero_percentage_mask.index)
    mask = sum_mask.loc[valid_columns] & zero_percentage_mask.loc[valid_columns]
    before = df.shape[1]
    df = df.loc[:, [sample_col] + mask.index[mask].tolist()]
    after = df.shape[1]
    print(f"Removed {before - after} columns with zero sum or more than 50% zeros.")
    return df

def transform_sample_name_GLASS(name):
    """
    Example transformation for GLASS sample IDs.
    Adjust as needed for your naming conventions.
    """
    if 'TP' in name:
        name = name.split('TP')[0] + 'TP'
    name = name.replace('.', '-')
    return name

###############################################################################
# 2. Plotting expression distributions
###############################################################################

def plot_expression_subplots_before_after(df_before, df_after, genes, dataset_label, outdir):
    """
    Creates 2x5 subplots (2 rows, 5 columns) for the given dataset:
      - row 1: boxplots (or violin) of each gene BEFORE normalization
      - row 2: the same genes AFTER normalization
    Saves the figure to disk.

    df_before, df_after must have the same samples and columns (besides scaling).
    The 'genes' list should be in the columns.
    """

    os.makedirs(outdir, exist_ok=True)
    fig, axes = plt.subplots(nrows=2, ncols=len(genes), figsize=(4*len(genes), 8))
    # Top row = BEFORE, bottom row = AFTER
    for i, gene in enumerate(genes):
        if gene not in df_before.columns or gene not in df_after.columns:
            continue
        # Row 0, col i => before
        sns.boxplot(y=df_before[gene], ax=axes[0,i], color='skyblue')
        axes[0,i].set_title(f"{gene} (Before)")
        axes[0,i].set_ylabel("Expr (Before)")

        # Row 1, col i => after
        sns.boxplot(y=df_after[gene], ax=axes[1,i], color='lightgreen')
        axes[1,i].set_title(f"{gene} (After)")
        axes[1,i].set_ylabel("Expr (After)")

    plt.suptitle(f"{dataset_label}: Expression Before vs After (Top=Before, Bottom=After)")
    plt.tight_layout(rect=[0,0,1,0.95])
    fname = os.path.join(outdir, f"{dataset_label}_2x5_before_after.png")
    plt.savefig(fname, dpi=120)
    plt.close()


def plot_expression_grid_three_datasets(df_before_dict, df_after_dict, genes, outdir):

    os.makedirs(outdir, exist_ok=True)
    dataset_list = list(df_before_dict.keys())  # e.g. ['TCGA','CGGA','GLASS']

    # 1) BEFORE
    fig, axes = plt.subplots(nrows=len(dataset_list), ncols=len(genes), figsize=(4*len(genes),4*len(dataset_list)))
    for r, ds_name in enumerate(dataset_list):
        df_bef = df_before_dict[ds_name]
        for c, gene in enumerate(genes):
            ax = axes[r,c]
            if gene in df_bef.columns:
                sns.boxplot(y=df_bef[gene], ax=ax, color='skyblue')
            ax.set_title(f"{ds_name} - {gene}")
            ax.set_ylabel("Expr (Before)")
    plt.suptitle("3x5 Subplots: Expression Before Normalization")
    plt.tight_layout(rect=[0,0,1,0.95])
    plt.savefig(os.path.join(outdir, "three_datasets_3x5_before_norm.png"), dpi=120)
    plt.close()

    # 2) AFTER
    fig, axes = plt.subplots(nrows=len(dataset_list), ncols=len(genes), figsize=(4*len(genes),4*len(dataset_list)))
    for r, ds_name in enumerate(dataset_list):
        df_aft = df_after_dict[ds_name]
        for c, gene in enumerate(genes):
            ax = axes[r,c]
            if gene in df_aft.columns:
                sns.boxplot(y=df_aft[gene], ax=ax, color='lightgreen')
            ax.set_title(f"{ds_name} - {gene}")
            ax.set_ylabel("Expr (After)")
    plt.suptitle("3x5 Subplots: Expression After Normalization")
    plt.tight_layout(rect=[0,0,1,0.95])
    plt.savefig(os.path.join(outdir, "three_datasets_3x5_after_norm.png"), dpi=120)
    plt.close()

###############################################################################
# 3. Limma in R
###############################################################################

def limma(data, threshold, consortium, threshold_type="fixed", outdir="./results"):
    """
    Runs a differential expression analysis in R using limma.
    """
    os.makedirs(outdir, exist_ok=True)
    data = data.copy()
    if pd.api.types.is_numeric_dtype(data['case.vital.status'] ):
        data['case.vital.status']  = data['case.vital.status'] .map({0: 'alive', 1: 'dead'}).astype(str)

    if threshold_type == "fixed":
        data = data[~((data['Overall.Survival.Months'] <= threshold) & (data['case.vital.status'] == 'alive'))]
        file_suffix = f"{consortium}_{threshold}-months"
        median_survival_r = threshold
    else:
        # range-based threshold
        low, high = threshold
        data = data[~((data['Overall.Survival.Months'] <= low) & (data['case.vital.status'] == 'alive'))]
        data = data[(data['Overall.Survival.Months'] <= low) | (data['Overall.Survival.Months'] > high)]
        file_suffix = f"{consortium}_{low}-{high}-months"
        median_survival_r = low

    data['group'] = np.where(data['Overall.Survival.Months'] > median_survival_r, 'high', 'low')
    n_low = sum(data['group']=='low')
    n_high= sum(data['group']=='high')
    print(f"{consortium} {file_suffix} => #Low group = {n_low}, #High group = {n_high}")

    data.drop(columns=['case.vital.status'], inplace=True)

    if 'batch' in data.columns:
        design_formula = "~ group + Diagnosis.Age + Sex + batch"
        exclude_cols_r = "c('Sample.ID','Overall.Survival.Months','Diagnosis.Age','Sex','batch','group')"
    else:
        design_formula = "~ group + Diagnosis.Age + Sex"
        exclude_cols_r = "c('Sample.ID','Overall.Survival.Months','Diagnosis.Age','Sex','group')"

    r_data = pandas2ri.py2rpy(data)
    r_script = f"""
    library(limma)
    data <- as.data.frame(data)
    gene_expression_data <- data[, !names(data) %in% {exclude_cols_r}]
    rownames(gene_expression_data) <- data$Sample.ID
    data$group <- factor(data$group, levels = c('low','high'))

    design <- model.matrix({design_formula}, data = data)
    rownames(design) <- data$Sample.ID

    fit <- lmFit(t(gene_expression_data), design)
    fit <- eBayes(fit)
    results <- topTable(fit, coef='grouphigh', number=Inf, adjust.method='fdr')
    out_all <- file.path("{outdir}", paste0("all_genes_{file_suffix}.csv"))
    write.csv(results, file=out_all, row.names=TRUE)
    sig <- results[results$adj.P.Val<0.1, ]
    out_sig <- file.path("{outdir}", paste0("significant_genes_{file_suffix}.csv"))
    write.csv(sig, file=out_sig, row.names=TRUE)
    """

    ro.globalenv['data'] = r_data
    ro.r(r_script)

    results_df = pd.read_csv(os.path.join(outdir, f"significant_genes_{file_suffix}.csv"), index_col=0)
    return results_df

###############################################################################
# 4. Manual log-fold-change + coexpression
###############################################################################

def compute_log_fold_change(df, threshold, genes_of_interest, threshold_type='fixed'):
    data = df.copy()
    if threshold_type=='fixed':
        data = data[~((data['Overall.Survival.Months']<=threshold) & (data['case.vital.status']=='alive'))]
        cutoff = threshold
    else:
        low, high = threshold
        data = data[~((data['Overall.Survival.Months']<=low) & (data['case.vital.status']=='alive'))]
        data = data[(data['Overall.Survival.Months']<=low) | (data['Overall.Survival.Months']>high)]
        cutoff = low

    data['group'] = np.where(data['Overall.Survival.Months']>cutoff, 'high','low')
    results = []
    for gene in genes_of_interest:
        if gene in data.columns:
            m_high = data.loc[data['group']=='high', gene].mean()
            m_low  = data.loc[data['group']=='low', gene].mean()
            log_fc = m_high - m_low
            results.append((gene, log_fc))
        else:
            results.append((gene, np.nan))

    return pd.DataFrame(results, columns=['Gene','logFC'])

from scipy.cluster.hierarchy import linkage, dendrogram
def coexpression_analysis(df, genes, dataset_name="Dataset", outdir="./results", suffix=""):
    os.makedirs(outdir, exist_ok=True)

    subset_genes = [g for g in genes if g in df.columns]
    if len(subset_genes) < 2:
        print(f"[{dataset_name}] Not enough genes to do correlation heatmap.")
        return

    corr_mat = df[subset_genes].corr()

    linkage_matrix = linkage(corr_mat, method='average')
    dendro = dendrogram(linkage_matrix, labels=subset_genes, leaf_rotation=90)

    clustered_genes = [subset_genes[i] for i in dendro['leaves']]
    clustered_corr_mat = corr_mat.loc[clustered_genes, clustered_genes]

    # Try turning off annot or reduce the figure size if the matrix is very large
    plt.figure(figsize=(max(8, len(subset_genes)), max(6, len(subset_genes))))
    sns.heatmap(clustered_corr_mat, annot=False, cmap='vlag', center=0)
    plt.title(f"Co-expression (Hierarchical Clustering): {dataset_name} {suffix}")
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()

    outname = os.path.join(outdir, f"coexpr_clustered_{dataset_name}_{suffix}.pdf")
    outname = str(outname)  # ensure plain string

    plt.savefig(outname, dpi=150)  # maybe lower DPI too if it's very large
    plt.close()


###############################################################################
# 5. Scatter + Bar Plots
###############################################################################

def plot_scatter_logfc_across_datasets(logfc_dict_1, logfc_dict_2, dataset_pair_label,
                                       threshold_label, outdir="./results"):
    """
    Scatter plot of logFC across two datasets.
    logfc_dict_1, logfc_dict_2 are dataframes with columns ['Gene','logFC'].
    """
    os.makedirs(outdir, exist_ok=True)
    df1 = logfc_dict_1.rename(columns={'logFC':'logFC_1'})
    df2 = logfc_dict_2.rename(columns={'logFC':'logFC_2'})

    merged = pd.merge(df1, df2, on='Gene', how='inner')

    plt.figure(figsize=(6,6))
    sns.scatterplot(data=merged, x='logFC_1', y='logFC_2')
    plt.axhline(0, color='black', linewidth=0.8)
    plt.axvline(0, color='black', linewidth=0.8)
    plt.title(f"{dataset_pair_label} logFC @ {threshold_label}")
    plt.xlabel(f"{dataset_pair_label.split('_vs_')[0]} logFC")
    plt.ylabel(f"{dataset_pair_label.split('_vs_')[1]} logFC")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, f"scatter_logFC_{dataset_pair_label}_{threshold_label}.png"), dpi=120)
    plt.close()


def plot_bar_logfc_across_datasets(logfc_tcga, logfc_cgga, logfc_glass, threshold_label, outdir="./results"):
    """
    Creates a grouped bar plot:
      x-axis: Genes
      groups: Datasets (TCGA, CGGA, GLASS)
      y-axis: logFC
    If many genes, we can do a facet or 2-row arrangement.
    """
    os.makedirs(outdir, exist_ok=True)

    merged = logfc_tcga.rename(columns={'logFC':'TCGA'})
    merged = pd.merge(merged, logfc_cgga.rename(columns={'logFC':'CGGA'}), on='Gene', how='outer')
    merged = pd.merge(merged, logfc_glass.rename(columns={'logFC':'GLASS'}), on='Gene', how='outer')

    melted = merged.melt(id_vars='Gene', value_vars=['TCGA','CGGA','GLASS'],
                         var_name='Dataset', value_name='logFC')

    num_genes = len(melted['Gene'].unique())
    # Determine the number of rows based on gene count
    if num_genes <= 6:
        rows = 1
        cols = num_genes
    else:
        rows = (num_genes + 5) // 6  # Ensure multiples of 6
        cols = 6
    fig, axes = plt.subplots(rows, 1, figsize=(max(10, 1.5 * cols), 6 * rows), constrained_layout=True)
    # If only one row, ensure axes is iterable
    if rows == 1:
        axes = [axes]
    for i, ax in enumerate(axes):
        start_idx = i * cols
        end_idx = start_idx + cols
        genes_subset = melted['Gene'].unique()[start_idx:end_idx]
        data_subset = melted[melted['Gene'].isin(genes_subset)]
        sns.barplot(data=data_subset, x='Gene', y='logFC', hue='Dataset', ax=ax)
        ax.set_title(f"LogFC Comparison at {threshold_label} (Genes {start_idx + 1}-{min(end_idx, num_genes)})")
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')

    plt.savefig(os.path.join(outdir, f"logFC_bar_{threshold_label}.png"), dpi=120)
    plt.close()


###############################################################################
# 6. Main data-loading functions + Outlier Removal + PCA
###############################################################################

def load_TCGA(tcga_expr_file, tcga_outcome_file, outdir):

    results_dir = os.path.join(outdir, "results")
    os.makedirs(results_dir, exist_ok=True)

    outcome_data = pd.read_excel(tcga_outcome_file, sheet_name='TCGA')
    expr = pd.read_csv(tcga_expr_file, sep=',')

    # Filter to intersection
    common_ids = set(expr['samples']).intersection(set(outcome_data['Sample ID']))
    expr = expr[expr['samples'].isin(common_ids)].copy()
    expr.rename(columns={'samples':'Sample.ID'}, inplace=True)

    # Keep only RNAseq + Agilent
    valid_platforms = ['RNAseq','agilent']
    expr = expr[expr['platform'].isin(valid_platforms)].copy()
    expr = remove_zero_genes(expr)   # presumably your custom function
    print("[TCGA] shape after platform filter:", expr.shape)

    sample_ids = expr['Sample.ID'].values
    platforms  = expr['platform'].values
    gene_cols  = [c for c in expr.columns if c not in ['Sample.ID','platform']]

    # Expression matrix before normalization
    df_before = expr[gene_cols].fillna(0)

    # ----- PCA before normalization -----
    pca = PCA(n_components=2)
    pca_scores_before = pca.fit_transform(df_before)
    pc_df_before = pd.DataFrame(pca_scores_before, columns=['PC1','PC2'])
    pc_df_before['Platform'] = platforms
    plt.figure()
    sns.scatterplot(data=pc_df_before, x='PC1', y='PC2', hue='Platform')
    plt.title("TCGA PCA Before Normalization")
    plt.tight_layout()
    plt.savefig(os.path.join(results_dir, "TCGA_pca_before_norm.png"), dpi=120)
    plt.close()

    # ----- Separate out the two platforms -----
    rna_mask = (expr['platform'] == 'RNAseq')
    agl_mask = (expr['platform'] == 'agilent')

    rna_df = df_before.loc[rna_mask]
    agl_df = df_before.loc[agl_mask]

    # ----- Normalize each platform separately -----
    scaler_rna = StandardScaler()
    rna_norm = scaler_rna.fit_transform(rna_df)

    scaler_agl = StandardScaler()
    agl_norm = scaler_agl.fit_transform(agl_df)

    # Reconstruct dataframes
    rna_norm_df = pd.DataFrame(rna_norm, index=rna_df.index, columns=rna_df.columns)
    agl_norm_df = pd.DataFrame(agl_norm, index=agl_df.index, columns=agl_df.columns)

    # Concatenate them back in the original sample order
    df_after = pd.concat([rna_norm_df, agl_norm_df], axis=0)
    # reorder to match original expr order
    df_after = df_after.loc[df_before.index]

    # ----- PCA after separate normalization -----
    pca2 = PCA(n_components=2)
    pca_scores_after = pca2.fit_transform(df_after)
    pc_df_after = pd.DataFrame(pca_scores_after, columns=['PC1','PC2'])
    pc_df_after['Platform'] = platforms
    plt.figure()
    sns.scatterplot(data=pc_df_after, x='PC1', y='PC2', hue='Platform')
    plt.title("TCGA PCA After Separate Normalization")
    plt.tight_layout()
    plt.savefig(os.path.join(results_dir, "TCGA_pca_after_norm.png"), dpi=120)
    plt.close()

    # ----- Build final DataFrame (after) -----
    df_norm = df_after.copy()
    df_norm['Sample.ID'] = expr['Sample.ID']
    df_norm.reset_index(drop=True, inplace=True)

    # Merge with outcome
    outcome_cols = ['Sample ID','Overall Survival (Months)','Diagnosis Age','Sex','Overall Survival Status']
    outcome_data = outcome_data[outcome_cols].copy()
    outcome_data.rename(columns={
        'Sample ID':'Sample.ID',
        'Overall Survival (Months)':'Overall.Survival.Months',
        'Diagnosis Age':'Diagnosis.Age',
        'Overall Survival Status':'case.vital.status'
    }, inplace=True)
    outcome_data['case.vital.status'] = outcome_data['case.vital.status'].replace({
        '1:DECEASED':'dead','0:LIVING':'alive'
    })

    merged = pd.merge(outcome_data, df_norm, on='Sample.ID')
    merged.drop_duplicates(subset=['Sample.ID'], inplace=True)
    print("[TCGA] final shape:", merged.shape)

    # Save
    #merged.to_csv(os.path.join(outdir, "TCGA_processed.csv"), index=False)

    # Return: final + separate (before, after) for expression distribution
    df_before_out = df_before.copy()
    df_before_out.index = expr['Sample.ID']
    return merged, df_before_out, df_after

def load_CGGA(outdir):
    """
    Example of loading CGGA 325 + 693, log transforming, scaling, combining.
    Returns combined_data_2 or whichever is your final integrated DF.
    Adjust file paths as needed.
    """
    #--- Load clinical
    results_dir = os.path.join(outdir, "results")
    os.makedirs(results_dir, exist_ok=True)

    # File paths
    clinical_325 = "/home/mkh062/Desktop/scratch/TCGA_project/CGGA-RAW DATA/CGGA.mRNAseq_325_clinical.20200506.txt"
    clinical_693 = "/home/mkh062/Desktop/scratch/TCGA_project/CGGA-RAW DATA/CGGA.mRNAseq_693_clinical.20200506.txt"

    expr_325 = "/home/mkh062/Desktop/scratch/TCGA_project/CGGA.mRNAseq_325.RSEM-genes.20200506.txt/CGGA.mRNAseq_325.RSEM-genes.20200506.txt"
    expr_693 = "/home/mkh062/Desktop/scratch/TCGA_project/CGGA.mRNAseq_693.RSEM-genes.20200506.txt/CGGA.mRNAseq_693.RSEM-genes.20200506.txt"

    # Load clinical
    c325 = pd.read_csv(clinical_325, sep='\t')
    c325 = c325[
        (c325['PRS_type'] == 'Primary') & (c325['Histology'] == 'GBM') & (c325['IDH_mutation_status'] == 'Wildtype')]
    c325 = c325[['CGGA_ID', 'Gender', 'Age', 'Censor (alive=0; dead=1)', 'OS']].dropna()
    c325['OS'] = c325['OS'] / 30.44
    c325.rename(columns={'CGGA_ID': 'Sample.ID'}, inplace=True)

    c693 = pd.read_csv(clinical_693, sep='\t')
    c693 = c693[
        (c693['PRS_type'] == 'Primary') & (c693['Histology'] == 'GBM') & (c693['IDH_mutation_status'] == 'Wildtype')]
    c693 = c693[['CGGA_ID', 'Gender', 'Age', 'Censor (alive=0; dead=1)', 'OS']].dropna()
    c693['OS'] = c693['OS'] / 30.44
    c693.rename(columns={'CGGA_ID': 'Sample.ID'}, inplace=True)

    # Load expression
    df_325 = pd.read_csv(expr_325, sep='\t')
    df_693 = pd.read_csv(expr_693, sep='\t')

    def reshape_expression(df):
        df_t = df.set_index('Gene_Name').T.reset_index()
        new_cols = ['Sample.ID'] + df['Gene_Name'].tolist()
        df_t.columns = new_cols
        return df_t

    df_325_t = reshape_expression(df_325)
    df_693_t = reshape_expression(df_693)

    # Filter to overlap
    common_325 = set(df_325_t['Sample.ID']).intersection(set(c325['Sample.ID']))
    df_325_t = df_325_t[df_325_t['Sample.ID'].isin(common_325)].copy()

    common_693 = set(df_693_t['Sample.ID']).intersection(set(c693['Sample.ID']))
    df_693_t = df_693_t[df_693_t['Sample.ID'].isin(common_693)].copy()
    # Remove zero variance genes
    df_325_t = remove_zero_genes(df_325_t)
    df_693_t = remove_zero_genes(df_693_t)

    common_genes = list(set(df_325_t.columns) & set(df_693_t.columns) - {'Sample.ID'})
    df_325_t = df_325_t[['Sample.ID'] + common_genes]
    df_693_t = df_693_t[['Sample.ID'] + common_genes]

    df_325_t[common_genes] = np.log1p(df_325_t[common_genes] + 1e-5)
    df_693_t[common_genes] = np.log1p(df_693_t[common_genes] + 1e-5)
    # Merge with clinical
    batch1 = pd.merge(c325, df_325_t, on='Sample.ID')
    batch2 = pd.merge(c693, df_693_t, on='Sample.ID')
    batch1['batch'] = 'batch325'
    batch2['batch'] = 'batch693'

    combined = pd.concat([batch1, batch2], ignore_index=True)

    #--- Outlier removal example (via PCA)
    gene_cols_for_pca = [g for g in combined.columns if
                         g not in ['Sample.ID', 'OS', 'Age', 'Gender', 'Censor (alive=0; dead=1)', 'batch']]
    X_before = combined[gene_cols_for_pca].values
    pca = PCA(n_components=2)
    pca_before = pca.fit_transform(X_before)
    combined['PC1_before'] = pca_before[:, 0]
    combined['PC2_before'] = pca_before[:, 1]

    plt.figure()
    sns.scatterplot(data=combined, x='PC1_before', y='PC2_before', hue='batch')
    plt.title("CGGA PCA Before Normalization")
    plt.tight_layout()
    plt.savefig(os.path.join(results_dir, "CGGA_pca_before_norm.png"), dpi=120)
    plt.close()

    outlier_indices = combined[(np.abs(combined['PC1_before'])>100) | (np.abs(combined['PC2_before'])>100)].index
    combined = combined.drop(outlier_indices).reset_index(drop=True)
    common_sample_ids_batch1 = set(combined['Sample.ID']).intersection(set(batch1['Sample.ID']))
    batch1 = batch1[batch1['Sample.ID'].isin(common_sample_ids_batch1)].reset_index(drop=True)
    # Take intersection of 'Sample ID' between batch2_data and combined_data
    common_sample_ids_batch2 = set(combined['Sample.ID']).intersection(set(batch2['Sample.ID']))
    batch2 = batch2[batch2['Sample.ID'].isin(common_sample_ids_batch2)].reset_index(drop=True)

    combined = pd.concat([batch1, batch2], ignore_index=True)
    df_before = combined[common_genes+['Sample.ID']].set_index('Sample.ID')

    from sklearn.preprocessing import StandardScaler

    scaler = StandardScaler()

    batch1[common_genes] = scaler.fit_transform(batch1[common_genes])
    batch2[common_genes] = scaler.fit_transform(batch2[common_genes])
    combined = pd.concat([batch1, batch2], ignore_index=True)
    #combined[common_genes] = scaler.fit_transform(combined[common_genes])
    df_after = combined[common_genes+['Sample.ID']].set_index('Sample.ID')
    # PCA after scale
    X_after = combined[common_genes].values
    pca2 = PCA(n_components=2)
    pca_after = pca2.fit_transform(X_after)
    combined = combined.copy()
    combined['PC1_after'] = pca_after[:, 0]
    combined['PC2_after'] = pca_after[:, 1]

    plt.figure()
    sns.scatterplot(data=combined, x='PC1_after', y='PC2_after', hue='batch')
    plt.title("CGGA PCA After Normalization and Outlier removal")
    plt.tight_layout()
    plt.savefig(os.path.join(results_dir, "CGGA_pca_after_norm.png"), dpi=120)
    plt.close()
    combined.rename(columns={
        'OS': 'Overall.Survival.Months',
        'Age': 'Diagnosis.Age',
        'Gender': 'Sex',
        'Censor (alive=0; dead=1)': 'case.vital.status'
    }, inplace=True)
    combined = combined.drop(columns=['batch', 'PC1_after',  'PC2_after'])
    print("[CGGA] final shape:", combined.shape)
    # Save
    #combined.to_csv(os.path.join(outdir, "CGGA_processed.csv"), index=False)
    return combined, df_before, df_after

def load_GLASS(outdir):
    """
    Loads GLASS data, merges with outcome, log1p transform, PCA before/after, saves final CSV.
    Also can do gene expression distribution plots if desired.
    """
    results_dir = os.path.join(outdir, "results")
    os.makedirs(results_dir, exist_ok=True)

    outcome_file = "/home/mkh062/Desktop/scratch/TCGA_project/Compiled data_TCGA_CGGA_GLASS_jh080424.xlsx"
    expr_file = "/home/mkh062/Desktop/scratch/TCGA_project/gene_tpm_matrix_all_samples.tsv"

    outdf = pd.read_excel(outcome_file, sheet_name='GLASS')
    expr = pd.read_csv(expr_file, sep='\t')

    # Transpose
    expr_t = expr.set_index('Gene_symbol').T.reset_index()
    new_cols = ['Sample ID'] + expr['Gene_symbol'].tolist()
    expr_t.columns = new_cols
    expr_t['Sample ID'] = expr_t['Sample ID'].apply(transform_sample_name_GLASS)

    common_ids = set(expr_t['Sample ID']).intersection(set(outdf['Sample ID']))
    expr_t = expr_t[expr_t['Sample ID'].isin(common_ids)].copy()
    expr_t = remove_zero_genes(expr_t)

    # merge
    outdf = outdf[['Sample ID','case_overall_survival_mo','case_age_diagnosis_years','case_sex','case_vital_status']]
    merged = pd.merge(outdf, expr_t, on='Sample ID')

    merged.rename(columns={
        'Sample ID': 'Sample.ID',
        'case_overall_survival_mo':'Overall.Survival.Months',
        'case_age_diagnosis_years':'Diagnosis.Age',
        'case_sex':'Sex',
        'case_vital_status':'case.vital.status'
    }, inplace=True)
    merged.drop_duplicates(subset=['Sample.ID'], inplace=True)
    # PCA before
    gene_cols = [c for c in merged.columns if c not in ['Sample.ID','Overall.Survival.Months','Diagnosis.Age','Sex','case.vital.status']]
    # log1p
    merged[gene_cols] = np.log1p(merged[gene_cols] + 1e-5)
    df_before = merged[gene_cols].copy()
    # scale
    scaler = StandardScaler()
    merged[gene_cols] = scaler.fit_transform(merged[gene_cols])
    df_after = merged[gene_cols].copy()
    print("[GLASS] final shape:", merged.shape)
    #merged.to_csv(os.path.join(outdir, "GLASS_processed.csv"), index=False)

    # Also return (before, after)
    df_before.index = merged['Sample.ID']
    df_after.index  = merged['Sample.ID']
    return merged, df_before, df_after

######################################
# 6. Bringing it all together
######################################
def gather_significant_genes(outdir, consortium, threshold_label):
    fname = os.path.join(outdir, f"significant_genes_{consortium}_{threshold_label}.csv")
    if not os.path.exists(fname):
        return set()
    df = pd.read_csv(fname, index_col=0)
    if df.empty:
        return set()
    return set(df.index)


def main():
    # define base directories
    BASE_DIR = "/home/mkh062/Desktop/scratch/TCGA_project"
    PROCESSED_DATA_DIR = os.path.join(BASE_DIR, "processed_data/try2")
    RESULTS_DIR = os.path.join(PROCESSED_DATA_DIR, "results")

    os.makedirs(PROCESSED_DATA_DIR, exist_ok=True)
    os.makedirs(RESULTS_DIR, exist_ok=True)

    # 1) Load datasets
    tcga_data, tcga_before, tcga_after = load_TCGA(
        tcga_expr_file=os.path.join(BASE_DIR, "GBMLGG_EB_RmDiffFullGenesRanRmDup.csv"),
        tcga_outcome_file=os.path.join(BASE_DIR, "Compiled data_TCGA_CGGA_GLASS_jh080424.xlsx"),
        outdir=PROCESSED_DATA_DIR
    )

    cgga_data, cgga_before, cgga_after = load_CGGA(
        outdir=PROCESSED_DATA_DIR
    )

    glass_data, glass_before, glass_after = load_GLASS(
        outdir=PROCESSED_DATA_DIR
    )

    # 2) Summaries
    print("\n--- Data Summaries ---")
    print("TCGA:", tcga_data.shape)
    print("CGGA:", cgga_data.shape)
    print("GLASS:", glass_data.shape)

    # 3) Make before/after expression subplots for the first 5 genes in each dataset
    #    We'll do a 2x5 subplot per dataset
    # plot_expression_subplots_before_after(tcga_before, tcga_after, genes_of_interest_5, "TCGA", RESULTS_DIR)
    # plot_expression_subplots_before_after(cgga_before, cgga_after, genes_of_interest_5, "CGGA", RESULTS_DIR)
    # plot_expression_subplots_before_after(glass_before, glass_after, genes_of_interest_5, "GLASS", RESULTS_DIR)

    # 4) Make 3x5 grid: each row = dataset, each column = one gene
    #    We'll do two big figures: one "before" for all 3 datasets, one "after"
    df_before_dict = {
        'TCGA': tcga_before,
        'CGGA': cgga_before,
        'GLASS': glass_before
    }
    df_after_dict = {
        'TCGA': tcga_after,
        'CGGA': cgga_after,
        'GLASS': glass_after
    }
    # all_genes_of_interest = [
    #     'TRIP4', 'EFEMP2', 'EFNB2', 'TRAF3', 'DYNLT3', 'MSN', 'C9orf64',
    #     'HSP90B1', 'VASN', 'PARVA', 'SWAP70', 'PALM2-AKAP2', 'LGALS8',
    #     'PDIA4', 'TBC1D1', 'RPAP3', 'DIRAS3', 'RAB37', 'AMIGO3',
    #     'TCHP',
    #     'EXOSC10', 'GRM5',
    #     'C1orf94', 'KCNV1', 'HCN1', 'STK36', 'MARCH8', 'HORMAD2', 'BARX1',
    #     'SDCBP', 'UBE2H', 'LAP3', 'TEX2', 'WWP1', 'SCRN1', 'USP32',
    #     'SNX10', 'GRB2', 'RAB21'
    # ]
    # genes_of_interest_5 = all_genes_of_interest[:5]
    # plot_expression_grid_three_datasets(df_before_dict, df_after_dict, genes_of_interest_5, RESULTS_DIR)

    # 5) Run limma for thresholds
    thresholds_fixed = [24, 60]
    thresholds_range = [(24,48), (24,60)]

    # fixed
    for threshold in thresholds_fixed:
        for ds_data, ds_name in [(tcga_data, "TCGA"), (cgga_data, "CGGA"), (glass_data, "GLASS")]:
            res = limma(ds_data.copy(),
                                       threshold=threshold,
                                       consortium=ds_name,
                                       threshold_type="fixed",
                                       outdir=RESULTS_DIR)

    # range
    for (low,high) in thresholds_range:
        for ds_data, ds_name in [(tcga_data, "TCGA"), (cgga_data, "CGGA"), (glass_data, "GLASS")]:
            res = limma(ds_data.copy(),
                                       threshold=(low,high),
                                       consortium=ds_name,
                                       threshold_type="range",
                                       outdir=RESULTS_DIR)

    # 6) For each threshold, gather union of sig genes, do coexpression & scatter
    # def process_threshold(th_label, threshold_type="fixed"):
    #     sig_tcga = gather_significant_genes(RESULTS_DIR, "TCGA", th_label)
    #     sig_cgga = gather_significant_genes(RESULTS_DIR, "CGGA", th_label)
    #     sig_glass= gather_significant_genes(RESULTS_DIR, "GLASS", th_label)
    #     union_genes = sig_tcga.union(sig_cgga).union(sig_glass)
    #     if not union_genes:
    #         print(f"No union sig genes at {th_label}.")
    #         return
    #     print(f"Union of sig genes at {th_label}: {len(union_genes)}")
    #
    #
    #     # Coexpr
    #     # coexpression_analysis(tcga_data, union_genes, dataset_name="TCGA", outdir=RESULTS_DIR, suffix=th_label)
    #     # coexpression_analysis(cgga_data, union_genes, dataset_name="CGGA", outdir=RESULTS_DIR, suffix=th_label)
    #     # coexpression_analysis(glass_data, union_genes, dataset_name="GLASS", outdir=RESULTS_DIR, suffix=th_label)
    #
    #     # scatter of logFC
    #     # parse threshold
    #     if threshold_type == "range":
    #         # range
    #         parts = th_label.split('-')
    #         low_ = int(parts[0])
    #         high_ = int(parts[1])
    #         df_tcga_logfc = compute_log_fold_change(tcga_data, (low_, high_), union_genes, threshold_type='range')
    #         df_cgga_logfc = compute_log_fold_change(cgga_data, (low_, high_), union_genes, threshold_type='range')
    #         df_glass_logfc= compute_log_fold_change(glass_data,(low_, high_), union_genes, threshold_type='range')
    #     else:
    #         # fixed
    #         thr_int = int(th_label.replace('-months',''))
    #         df_tcga_logfc = compute_log_fold_change(tcga_data, thr_int, union_genes, threshold_type='fixed')
    #         df_cgga_logfc = compute_log_fold_change(cgga_data, thr_int, union_genes, threshold_type='fixed')
    #         df_glass_logfc= compute_log_fold_change(glass_data, thr_int, union_genes, threshold_type='fixed')
    #
    #     plot_scatter_logfc_across_datasets(df_tcga_logfc, df_cgga_logfc, "TCGA_vs_CGGA", th_label, outdir=RESULTS_DIR)
    #     plot_scatter_logfc_across_datasets(df_tcga_logfc, df_glass_logfc, "TCGA_vs_GLASS",th_label, outdir=RESULTS_DIR)
    #     plot_scatter_logfc_across_datasets(df_cgga_logfc, df_glass_logfc,"CGGA_vs_GLASS",th_label, outdir=RESULTS_DIR)
    #
    #     # Also new bar plot for these union genes
    #     plot_bar_logfc_across_datasets(df_tcga_logfc, df_cgga_logfc, df_glass_logfc,
    #                                    threshold_label=th_label, outdir=RESULTS_DIR)
    # for (l,h) in thresholds_range:
    #     process_threshold(f"{l}-{h}-months", threshold_type="range")
    # # process thresholds
    # for thr in thresholds_fixed:
    #     process_threshold(f"{thr}-months", threshold_type="fixed")
    #
    # print("Done.")

if __name__=="__main__":
    main()


