import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
#matplotlib.use("Agg")  # force non-GUI backend
matplotlib.use("tkagg")
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from pipeline_CGGA import remove_zero_genes, log_transform, make_surv_groups_int
import helpers as hp
# R-related imports
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
pandas2ri.activate()
# File paths
clinical_325 = "/home/mkh062/Desktop/scratch/TCGA_project/CGGA-RAW DATA/CGGA.mRNAseq_325_clinical.20200506.txt"
clinical_693 = "/home/mkh062/Desktop/scratch/TCGA_project/CGGA-RAW DATA/CGGA.mRNAseq_693_clinical.20200506.txt"

expr_325 = "/home/mkh062/Desktop/scratch/TCGA_project/CGGA.mRNAseq_325.RSEM-genes.20200506.txt/CGGA.mRNAseq_325.RSEM-genes.20200506.txt"
expr_693 = "/home/mkh062/Desktop/scratch/TCGA_project/CGGA.mRNAseq_693.RSEM-genes.20200506.txt/CGGA.mRNAseq_693.RSEM-genes.20200506.txt"

OUTDIR           = "/home/mkh062/Desktop/scratch/TCGA_project/TCGA-combined-gene-exp-data/18052025"
BASE_DIR = "/home/mkh062/Desktop/scratch/TCGA_project"
tcga_expr_file = os.path.join(BASE_DIR, "GBMLGG_EB_RmDiffFullGenesRanRmDup.csv")
tcga_outcome_file=os.path.join(BASE_DIR, "Compiled data_TCGA_CGGA_GLASS_jh080424.xlsx")
# ════════════════════════════════════════════════════════════════════

SURV_THRESHOLD = 24      # months
PADJ_CUTOFF    = 0.1    # for “significant” genes

OUT_XLS     = os.path.join(OUTDIR, "limma_significant_genes2.xlsx")

def save_limma_sig_to_excel(results_dict: dict[str, pd.DataFrame],
                            xlsx_path: str,
                            padj_cutoff: float = 0.10) -> None:
    file_exists = os.path.isfile(xlsx_path)

    writer_opts = dict(engine="openpyxl")
    if file_exists:                 # append → allowed to use if_sheet_exists
        writer_opts.update(mode="a", if_sheet_exists="replace")

    with pd.ExcelWriter(xlsx_path, **writer_opts) as xls:
        for sheet, df in results_dict.items():
            sig = df[df["adj.P.Val"] < padj_cutoff].copy()
            sig.to_excel(xls, sheet_name=sheet, index=False)

    print(f"★ Significant limma hits written to: {xlsx_path}")
# -------------------------------------------------------------------
def run_limma(expr, meta, platform_label, padj=PADJ_CUTOFF):
    """Generic limma or voom‑limma depending on data type (RNA‑seq vs micro‑array)."""
    # expression (genes × samples) & metadata must be sample‑aligned
    r_expr = pandas2ri.py2rpy(expr)    # genes rows × samples cols
    r_meta = pandas2ri.py2rpy(meta)
    print("           platform_label    ", platform_label)
    print("########################################## hig low \n" , meta['group'].value_counts())
    # build design matrix ~ group + age + sex (+ batch if available)
    design_terms = ["group", "Diagnosis.Age", "Sex"]
    if 'batch' in meta.columns:
        design_terms.append("batch")
    design_formula = "~ " + " + ".join(design_terms)

    ro.globalenv["expr"] = r_expr
    ro.globalenv["meta"] = r_meta
    ro.r(f"""
        library(limma)
        library(edgeR)
        expr <- as.matrix(expr)
        meta <- as.data.frame(meta)
        meta$group <- factor(meta$group, levels=c('low','high'))
        design <- model.matrix({design_formula}, data=meta)
    """)
    #if platform_label == "CGGA-z" or platform_label == "CGGA-qn" :
    ro.r("""
           fit       <- lmFit(expr, design)
       """)
    ro.r("""
        fit  <- eBayes(fit)
        tt  <- topTable(fit, coef=2, number=Inf, adjust.method='fdr')
    """)
    tt = ro.r("tt")
    tt.index.name = "Gene"
    tt.reset_index(inplace=True)

    # # write to Excel (one sheet per dataset)
    # with pd.ExcelWriter(out_xls, engine="openpyxl", mode="a" if os.path.exists(out_xls) else "w") as xls:
    #     tt.to_excel(xls, sheet_name=f"{platform_label}_all", index=False)
    #     tt[tt['adj.P.Val'] < padj].to_excel(xls, sheet_name=f"{platform_label}_sig", index=False)

    print(f"[{platform_label}] limma done ‑ {(tt['adj.P.Val']<padj).sum()} significant genes (FDR<{padj})")
    return tt


def load_TCGA(care):
    outcome_data = pd.read_excel(tcga_outcome_file, sheet_name='TCGA')
    expr = pd.read_csv(tcga_expr_file, sep=',')
    valid_platforms = ['RNAseq','agilent']
    expr = expr[expr['platform'].isin(valid_platforms)].copy()
    if care:
        expr['samples'] = expr['samples'].str.replace(r'-01$', '', regex=True)
        outcome_data['Sample ID'] = outcome_data['Sample ID'].str.replace(r'-01$', '', regex=True)
        BASE_DIR = "/home/mkh062/Desktop/scratch/TCGA_project"
        id_file = os.path.join(BASE_DIR, "TCGA_filter_treatment1.xlsx")
        id_data = pd.read_excel(id_file, sheet_name='Sheet1')
        common_ids = set(expr['samples']) & set(id_data['Case ID'])
    else:
        common_ids = set(expr['samples']).intersection(set(outcome_data['Sample ID']))


    outcome_data = outcome_data[outcome_data['Sample ID'].isin(common_ids)].copy()
    expr = expr[expr['samples'].isin(common_ids)].copy()

    expr.rename(columns={'samples':'Sample.ID'}, inplace=True)
    outcome_data.rename(columns={'Sample ID':'Sample.ID'}, inplace=True)
    expr = remove_zero_genes(expr)
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
    if care:
        plt.title("TCGA PCA Before Normalization")
        plt.tight_layout()
        plt.show()
        # plt.savefig(os.path.join(results_dir, "care_TCGA_pca_before_norm.png"), dpi=300)
        # plt.close()
    else:
        plt.title("TCGA PCA Before Normalization")
        plt.tight_layout()
        plt.show()
        # plt.savefig(os.path.join(results_dir, "TCGA_pca_before_norm.png"), dpi=300)
        # plt.close()
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
    if care:
        plt.title("TCGA PCA After Separate Normalization")
        plt.tight_layout()
        plt.show()
        # plt.savefig(os.path.join(results_dir, "care_TCGA_pca_after_norm.png"), dpi=300)
        # plt.close()
    else:
        plt.title("TCGA PCA After Separate Normalization")
        plt.tight_layout()
        plt.show()
        # plt.savefig(os.path.join(results_dir, "TCGA_pca_after_norm.png"), dpi=300)
        # plt.close()

    # ----- Build final DataFrame (after) -----
    df_norm = df_after.copy()
    df_norm['Sample.ID'] = expr['Sample.ID']
    df_norm.reset_index(drop=True, inplace=True)
    # Merge with outcome
    outcome_cols = ['Sample.ID','Overall Survival (Months)','Diagnosis Age','Sex','Overall Survival Status']
    outcome_data = outcome_data[outcome_cols].copy()
    outcome_data.rename(columns={
        'Overall Survival (Months)':'Overall.Survival.Months',
        'Diagnosis Age':'Diagnosis.Age',
        'Overall Survival Status':'case.vital.status'
    }, inplace=True)
    outcome_data['case.vital.status'] = outcome_data['case.vital.status'].replace({
        '1:DECEASED':1,'0:LIVING':0
    })

    merged = pd.merge(outcome_data, df_norm, on='Sample.ID')
    merged.drop_duplicates(subset=['Sample.ID'], inplace=True)
    if 'platform' in merged.columns:
        merged.drop('platform', axis=1, inplace=True)
    if 'batch' in merged.columns:
        merged.drop('batch', axis=1, inplace=True)
    if 'Platform' in merged.columns:
        merged.drop('Platform', axis=1, inplace=True)


    if care:
        df_norm.set_index('Sample.ID', inplace=True)
        outcome_data.set_index('Sample.ID', inplace=True)
        outcome_data, df_norm = make_surv_groups_int(outcome_data, df_norm, SURV_THRESHOLD)
        outcome_data['Sex'] = outcome_data['Sex'].map({'Female': 0, 'Male': 1})
    else:
        df_norm.set_index('Sample.ID', inplace=True)
        outcome_data.set_index('Sample.ID', inplace=True)
        outcome_data, df_norm = make_surv_groups_int(outcome_data, df_norm, SURV_THRESHOLD)
        outcome_data['Sex'] = outcome_data['Sex'].map({'Female': 0, 'Male': 1})

    # Return: final + separate (before, after) for expression distribution
    df_before_out = df_before.copy()
    df_before_out.index = expr['Sample.ID']
    print("df_norm  ", df_norm)
    print("outcome_data  ", outcome_data)
    return outcome_data, df_norm

def load_CGGA(care):
    c325 = pd.read_csv(clinical_325, sep='\t')
    ############################################
    if care:
        # both treatments
        c325 = c325.loc[
            (c325['Radio_status (treated=1;un-treated=0)'] == 1) &
            (c325['Chemo_status (TMZ treated=1;un-treated=0)'] == 1)
            ].reset_index(drop=True)
    else:
        c325 = c325.loc[
            (c325['Radio_status (treated=1;un-treated=0)'] == 1) | (c325['Chemo_status (TMZ treated=1;un-treated=0)'] == 1)
            ].reset_index(drop=True)
    ############################################

    c325 = c325[
        (c325['PRS_type'] == 'Primary') & (c325['Histology'] == 'GBM') & (c325['IDH_mutation_status'] == 'Wildtype')]
    c325 = c325[['CGGA_ID', 'Gender', 'Age', 'Censor (alive=0; dead=1)', 'OS']].dropna()
    c325['OS'] = c325['OS'] / 30.44
    c325.rename(columns={'CGGA_ID': 'Sample.ID'}, inplace=True)
    ############################################
    c693 = pd.read_csv(clinical_693, sep='\t')
    if care:
        c693 = c693.loc[
            (c693['Radio_status (treated=1;un-treated=0)'] == 1) & (
                        c693['Chemo_status (TMZ treated=1;un-treated=0)'] == 1)
            ].reset_index(drop=True)
    else:
        c693 = c693.loc[
            (c693['Radio_status (treated=1;un-treated=0)'] == 1) | (c693['Chemo_status (TMZ treated=1;un-treated=0)'] == 1)
            ].reset_index(drop=True)
    ####################################################
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

    common_325 = set(df_325_t['Sample.ID']).intersection(set(c325['Sample.ID']))
    df_325_t = df_325_t[df_325_t['Sample.ID'].isin(common_325)].copy()
    c325 = c325[c325['Sample.ID'].isin(common_325)].copy()

    common_693 = set(df_693_t['Sample.ID']).intersection(set(c693['Sample.ID']))
    df_693_t = df_693_t[df_693_t['Sample.ID'].isin(common_693)].copy()
    c693 = c693[c693['Sample.ID'].isin(common_693)].copy()

    # Remove zero variance genes
    df_325_t = remove_zero_genes(df_325_t)
    df_693_t = remove_zero_genes(df_693_t)

    common_genes = list(set(df_325_t.columns) & set(df_693_t.columns) - {'Sample.ID'})
    df_325_t = df_325_t[['Sample.ID'] + common_genes]
    df_693_t = df_693_t[['Sample.ID'] + common_genes]

    return df_325_t, df_693_t, c325, c693

def PCA_cgga(batch14, batch24):
    combined = pd.concat([batch14, batch24], ignore_index=True)
    # --- Outlier removal example (via PCA)
    gene_cols_for_pca = [g for g in combined.columns if
                         g not in ['Sample.ID', 'OS', 'Age', 'Gender', 'Censor (alive=0; dead=1)', 'batch']]
    X_before = combined[gene_cols_for_pca].values
    pca = PCA(n_components=2)
    pca_before = pca.fit_transform(X_before)
    combined['PC1_before'] = pca_before[:, 0]
    combined['PC2_before'] = pca_before[:, 1]
    palette_map = {
        'batch325': 'red',
        'batch693': 'blue'
    }

    plt.figure(figsize=(6, 5))
    sns.scatterplot(data=combined, x='PC1_before', y='PC2_before', hue='batch', palette=palette_map, s=90,
                    edgecolor='w', alpha=0.8)
    plt.title("CGGA PCA")
    plt.tight_layout()
    plt.show()
    outlier_indices = combined[(np.abs(combined['PC1_before']) > 100) | (np.abs(combined['PC2_before']) > 100)].index
    combined = combined.drop(outlier_indices).reset_index(drop=True)
    common_sample_ids_batch1 = set(combined['Sample.ID']).intersection(set(batch14['Sample.ID']))
    batch14 = batch14[batch14['Sample.ID'].isin(common_sample_ids_batch1)].reset_index(drop=True)
    # Take intersection of 'Sample ID' between batch2_data and combined_data
    common_sample_ids_batch2 = set(combined['Sample.ID']).intersection(set(batch24['Sample.ID']))
    batch24 = batch24[batch24['Sample.ID'].isin(common_sample_ids_batch2)].reset_index(drop=True)
    return batch14, batch24


def CGGA_zscores(batch11, batch21):
    common_genes = list(set(batch11.columns) & set(batch21.columns) - {'Sample.ID', 'OS', 'Age', 'Gender', 'Censor (alive=0; dead=1)', 'batch'})

    from sklearn.preprocessing import StandardScaler

    scaler = StandardScaler()

    batch11[common_genes] = scaler.fit_transform(batch11[common_genes])
    batch21[common_genes] = scaler.fit_transform(batch21[common_genes])
    combined = pd.concat([batch11, batch21], ignore_index=True)


    X_after = combined[common_genes].values
    pca2 = PCA(n_components=2)
    pca_after = pca2.fit_transform(X_after)
    combined = combined.copy()
    combined['PC1_after'] = pca_after[:, 0]
    combined['PC2_after'] = pca_after[:, 1]
    palette_map = {
        'batch325': 'red',
        'batch693': 'blue'
    }

    plt.figure(figsize=(6, 5))
    sns.scatterplot(data=combined, x='PC1_after', y='PC2_after', hue='batch',  palette=palette_map, s=90,
                    edgecolor='w',alpha=0.8)
    plt.title("CGGA PCA After Normalization and Outlier removal")
    plt.tight_layout()
    plt.show()
    print(" z batch1 now ", batch11)
    print("z batch2 now  ", batch21)
    return batch11, batch21

def CGGA_QN_combat(batch13, batch23):
    print("QN batch1", batch13)
    print("QN batch2", batch23)
    concat= pd.concat([batch13, batch23], ignore_index=True)
    gene_cols = list(
        set(batch13.columns) & set(batch23.columns) - {'Sample.ID', 'OS', 'Age', 'Gender', 'Censor (alive=0; dead=1)',
                                                     'batch'})
    # --- 3) QUANTILE NORMALISE ---
    # transpose to genes×samples
    expr325 = (batch13.set_index('Sample.ID')[gene_cols]).T  # genes × samples
    expr693 = (batch23.set_index('Sample.ID')[gene_cols]).T   # genes × samples
    # combine them
    expr_combined = pd.concat([expr325, expr693], axis=1)
    # run single quantile-normalisation on the full matrix
    expr_qn = hp.quantile_normalize_combined(expr_combined)
    # recreate batch vector for downstream steps
    batch = pd.Series(
        ['batch325'] * expr325.shape[1] + ['batch693'] * expr693.shape[1],
        index=expr_qn.columns)
    print("  batch  C ", batch)
    # --- 5) COMBAT BATCH CORRECTION ---
    expr_cb = hp.combat(expr_qn, batch)
    # --- 6) PCA BEFORE & AFTER COMBAT ---
    hp.pca_plot(expr_qn, batch, title="PCA Before Combat", outfile="PCA_beforeComBat_CGGA.png")
    hp.pca_plot(expr_cb, batch, title="PCA After Combat", outfile="PCA_afterComBat_CGGA.png")

    # final expression (samples × genes):
    expr_df = expr_cb.T
    expr_df.reset_index(drop=True, inplace=True)
    print("  expr_df  ", expr_df)
    # 3. Sanity check: number of rows matches expr_df
    assert expr_df.shape[0] == concat.shape[0], \
        f"Row mismatch: expr_df has {expr_df.shape[0]} rows but combined meta has {concat.shape[0]}"
    concat.loc[:, gene_cols] = expr_df.loc[:, gene_cols].values
    grp = concat.groupby('batch')
    batch23 = grp.get_group('batch693')
    batch13 = grp.get_group('batch325')
    print("QN batch1  now ", batch13)
    print("QN batch2  now ", batch23)
    return batch13, batch23

def TCGA_QN(care):
    outcome_data = pd.read_excel(tcga_outcome_file, sheet_name='TCGA')
    expr = pd.read_csv(tcga_expr_file, sep=',')
    valid_platforms = ['RNAseq', 'agilent']
    expr = expr[expr['platform'].isin(valid_platforms)].copy()
    if care:
        expr['samples'] = expr['samples'].str.replace(r'-01$', '', regex=True)
        outcome_data['Sample ID'] = outcome_data['Sample ID'].str.replace(r'-01$', '', regex=True)
        BASE_DIR = "/home/mkh062/Desktop/scratch/TCGA_project"
        id_file = os.path.join(BASE_DIR, "TCGA_filter_treatment1.xlsx")
        id_data = pd.read_excel(id_file, sheet_name='Sheet1')
        common_ids = set(expr['samples']) & set(id_data['Case ID'])
    else:
        common_ids = set(expr['samples']).intersection(set(outcome_data['Sample ID']))

    outcome_data = outcome_data[outcome_data['Sample ID'].isin(common_ids)].copy()
    expr = expr[expr['samples'].isin(common_ids)].copy()

    expr.rename(columns={'samples': 'Sample.ID'}, inplace=True)
    outcome_data.rename(columns={'Sample ID': 'Sample.ID'}, inplace=True)
    expr = remove_zero_genes(expr)
    expr.set_index('Sample.ID', inplace=True)
    gene_cols = [c for c in expr.columns if c not in ['Sample.ID', 'platform']]

    # --- 3)  QUANTILE NORMALISATION -------------------------------------------
    expr_T = expr[gene_cols].T  # genes × samples
    expr_qn = hp.quantile_normalize_combined(expr_T)


    # --- 4)  ComBat (platform = batch) ----------------------------------------
    batch = expr['platform']
    assert (expr_qn.columns == batch.index).all(), "Sample order mismatch!"
    expr_cb = hp.apply_combat(expr_qn, batch)  # genes × samples



    # --- 5)  PCA before/after --------------------------------------------------
    hp.pca_TCGA(expr_qn, batch,
                title="PCA Before ComBat", outfile="PCA_beforeComBat_TCGA.png")
    hp.pca_TCGA(expr_cb, batch,
                title="PCA After  ComBat", outfile="PCA_afterComBat_TCGA.png")

    # --- 6)  Bring back to samples × genes ------------------------------------
    expr_final = expr_cb.T  # samples × genes  (index = sample IDs)

    # overwrite the expression columns inside `expr`
    expr.loc[expr_final.index, gene_cols] = expr_final[gene_cols]
    outcome_cols = ['Sample.ID', 'Overall Survival (Months)', 'Diagnosis Age', 'Sex', 'Overall Survival Status']
    outcome_data = outcome_data[outcome_cols].copy()
    outcome_data.rename(columns={
        'Overall Survival (Months)': 'Overall.Survival.Months',
        'Diagnosis Age': 'Diagnosis.Age',
        'Overall Survival Status': 'case.vital.status'
    }, inplace=True)
    outcome_data['case.vital.status'] = outcome_data['case.vital.status'].replace({
        '1:DECEASED': 1, '0:LIVING': 0
    })
    if 'platform' in expr.columns:
        expr.drop('platform', axis=1, inplace=True)

    outcome_data.set_index('Sample.ID', inplace=True)
    if 'platform' in outcome_data.columns:
        outcome_data.drop('platform', axis=1, inplace=True)
    outcome_data['Sex'] = outcome_data['Sex'].map({'Female': 0, 'Male': 1})
    outcome_data, expr = make_surv_groups_int(outcome_data, expr, SURV_THRESHOLD)

    print("   expr qnn  ", expr)
    print("   outcome_data qnn  ", outcome_data)
    return outcome_data, expr


def volcano_plot(tt: pd.DataFrame,
                 label: str,
                 outdir: str = OUTDIR,
                 padj_thresh: float = 0.1,
                 lfc_thresh: float = 1.0):
    # compute -log10 p
    x = tt['logFC']
    y = -np.log10(tt['P.Value'])

    # define significance
    sig = tt['adj.P.Val'] < padj_thresh

    plt.figure(figsize=(6,6))
    # non-significant in grey
    plt.scatter(x[~sig], y[~sig],
                color='lightgrey', s=10, alpha=0.6, label=f'n.s. (n={(~sig).sum()})')
    # significant in red
    plt.scatter(x[sig], y[sig],
                color='red',       s=10, alpha=0.8, label=f'signif (n={(sig).sum()})')

    # add threshold lines
    plt.axhline(-np.log10(padj_thresh), color='black', linestyle='--', linewidth=1)
    plt.axvline( lfc_thresh, color='black', linestyle='--', linewidth=1)
    plt.axvline(-lfc_thresh, color='black', linestyle='--', linewidth=1)

    plt.title(f"Volcano: {label}")
    plt.xlabel("log₂ fold change")
    plt.ylabel("-log₁₀(p-value)")
    plt.legend(frameon=False, fontsize=8, loc='upper left')

    plt.tight_layout()
    fn = f"{label}_volcano.png"
    #plt.show()
    # plt.savefig(os.path.join(outdir, fn), dpi=300)
    # plt.close()
    # print(f"Wrote {fn}")

def plot_pca_integration(expr_a   : pd.DataFrame,   # samples × genes
                         expr_b   : pd.DataFrame,   # samples × genes
                         label_a  : str,            # e.g. "TCGA"
                         label_b  : str,            # e.g. "CGGA"
                         outdir   : str,
                         tag      : str):           # "z" or "qn"
    # 1) keep only genes shared by *both* datasets --------------------------
    common_genes = expr_a.columns.intersection(expr_b.columns)
    expr_a = expr_a[common_genes]
    expr_b = expr_b[common_genes]

    # 2) build a simple batch label ----------------------------------------
    meta_a = pd.DataFrame(index=expr_a.index, data={'batch': label_a})
    meta_b = pd.DataFrame(index=expr_b.index, data={'batch': label_b})

    # 3) concatenate samples × genes ---------------------------------------
    expr_all = pd.concat([expr_a, expr_b], axis=0)
    meta_all = pd.concat([meta_a, meta_b], axis=0)

    # 4) (optional but recommended) centre / scale genes before PCA --------
    scaled = StandardScaler(with_mean=True, with_std=True).fit_transform(expr_all.values)

    # 5) run PCA ------------------------------------------------------------
    pcs = PCA(n_components=2, random_state=1).fit_transform(scaled)
    df_pcs = pd.DataFrame(pcs, index=expr_all.index, columns=['PC1', 'PC2'])
    df_pcs['batch'] = meta_all['batch'].values

    # 6) plot ---------------------------------------------------------------
    plt.figure(figsize=(6, 5))
    sns.scatterplot(
        data=df_pcs,
        x='PC1', y='PC2',
        hue='batch',
        palette={label_a: 'green', label_b: 'purple'},
        s=80, edgecolor='w', alpha=0.8
    )
    plt.title(f'PCA of Integrated {label_a} + {label_b} ({tag})')
    plt.tight_layout()
    plt.show()
    # 7) save ---------------------------------------------------------------
    # fname = os.path.join(outdir, f'PCA_{label_a}-{label_b}_{tag}.png')
    # plt.savefig(fname, dpi=300)
    # plt.close()
    # print(f'★ PCA plot saved → {fname}')

# -------------------------------------------------------------------
def main():
    #os.makedirs(OUTDIR, exist_ok=True)
    id_val = True
    df_325, df_693, c325, c693 = load_CGGA(id_val)
    df_325 = log_transform(df_325)
    df_693 = log_transform(df_693)
    batch1 = pd.merge(c325, df_325, on='Sample.ID')
    batch2 = pd.merge(c693, df_693, on='Sample.ID')
    batch1['batch'] = 'batch325'
    batch2['batch'] = 'batch693'
    batch1_p, batch2_p = PCA_cgga(batch1, batch2)
    batch1_qn, batch2_qn = CGGA_QN_combat(batch1_p, batch2_p)
    batch1_z, batch2_z = CGGA_zscores(batch1_p, batch2_p)
    rename_dict = {
        'OS': 'Overall.Survival.Months',
        'Age': 'Diagnosis.Age',
        'Gender': 'Sex',
        'Censor (alive=0; dead=1)': 'case.vital.status'
    }
    for df in [batch1_qn, batch2_qn, batch1_z, batch2_z]:
        df.rename(columns=rename_dict, inplace=True)
        df['Sex'] = df['Sex'].map({'Female': 0, 'Male': 1})
        df.drop('batch', axis=1, inplace=True)
    combined_qn = pd.concat([batch1_qn, batch2_qn], ignore_index=True)
    combined_z = pd.concat([batch1_z, batch2_z], ignore_index=True)
    clinical_cols = [
        'Overall.Survival.Months',
        'Diagnosis.Age',
        'Sex',
        'case.vital.status'
    ]
    meta_qn = combined_qn[['Sample.ID'] + clinical_cols]
    expr_qn = combined_qn.drop(columns=clinical_cols)
    meta_z = combined_z[['Sample.ID'] + clinical_cols]
    expr_z = combined_z.drop(columns=clinical_cols)

    for df in (meta_qn, expr_qn, meta_z, expr_z):
        df.set_index('Sample.ID', inplace=True)

    meta_qn, expr_qn =make_surv_groups_int(meta_qn, expr_qn, SURV_THRESHOLD)
    meta_z, expr_z =make_surv_groups_int(meta_z, expr_z, SURV_THRESHOLD)
    res_qn = run_limma(expr_qn.T, meta_qn, "CGGA-qn")
    res_z = run_limma(expr_z.T, meta_z, "CGGA-z")
    volcano_plot(res_qn, "CGGA-qn")
    volcano_plot(res_z, "CGGA-z")

    #################################################################################################################3
    id_val = False
    df_325_or, df_693_or, c325_or, c693_or = load_CGGA(id_val)
    df_325_or = log_transform(df_325_or)
    df_693_or = log_transform(df_693_or)
    batch1_or = pd.merge(c325_or, df_325_or, on='Sample.ID')
    batch2_or = pd.merge(c693_or, df_693_or, on='Sample.ID')
    batch1_or['batch'] = 'batch325'
    batch2_or['batch'] = 'batch693'
    batch1_p_or, batch2_p_or = PCA_cgga(batch1_or, batch2_or)
    batch1_qn_or, batch2_qn_or = CGGA_QN_combat(batch1_p_or, batch2_p_or)
    batch1_z_or, batch2_z_or = CGGA_zscores(batch1_p_or, batch2_p_or)
    rename_dict = {
        'OS': 'Overall.Survival.Months',
        'Age': 'Diagnosis.Age',
        'Gender': 'Sex',
        'Censor (alive=0; dead=1)': 'case.vital.status'
    }
    for df in [batch1_qn_or, batch2_qn_or, batch1_z_or, batch2_z_or]:
        df.rename(columns=rename_dict, inplace=True)
        df['Sex'] = df['Sex'].map({'Female': 0, 'Male': 1})
        df.drop('batch', axis=1, inplace=True)
    combined_qn_or = pd.concat([batch1_qn_or, batch2_qn_or], ignore_index=True)
    combined_z_or = pd.concat([batch1_z_or, batch2_z_or], ignore_index=True)
    clinical_cols = [
        'Overall.Survival.Months',
        'Diagnosis.Age',
        'Sex',
        'case.vital.status'
    ]
    meta_qn_or = combined_qn_or[['Sample.ID'] + clinical_cols]
    expr_qn_or = combined_qn_or.drop(columns=clinical_cols)
    meta_z_or = combined_z_or[['Sample.ID'] + clinical_cols]
    expr_z_or = combined_z_or.drop(columns=clinical_cols)

    for df in (meta_qn_or, expr_qn_or, meta_z_or, expr_z_or):
        df.set_index('Sample.ID', inplace=True)

    meta_qn_or, expr_qn_or = make_surv_groups_int(meta_qn_or, expr_qn_or, SURV_THRESHOLD)
    meta_z_or, expr_z_or = make_surv_groups_int(meta_z_or, expr_z_or, SURV_THRESHOLD)
    res_qn_or = run_limma(expr_qn_or.T, meta_qn_or, "CGGA-qn_or")
    res_z_or = run_limma(expr_z_or.T, meta_z_or, "CGGA_z_or")
    volcano_plot(res_qn_or, "CGGA-qn_or")
    volcano_plot(res_z_or, "CGGA-z_or")
    #################################################################################################################3
    id_val = True
    meta_z_TCGA, expr_z_TCGA = load_TCGA(id_val)
    meta_qn_TCGA, expr_qn_TCGA = TCGA_QN(id_val)
    res_z_TCGA = run_limma(expr_z_TCGA.T, meta_z_TCGA, "TCGA-z")
    res_qn_TCGA = run_limma(expr_qn_TCGA.T, meta_qn_TCGA, "TCGA-qn")
    volcano_plot(res_qn_TCGA , "TCGA-qn")
    volcano_plot(res_z_TCGA , "TCGA-z")
    #############################################################################################################
    id_val = False
    meta_z_TCGA_nc, expr_z_TCGA_nc = load_TCGA(id_val)
    res_z_TCGA_nc = run_limma(expr_z_TCGA_nc.T, meta_z_TCGA_nc, "TCGA-z_nc")
    volcano_plot(res_z_TCGA_nc , "TCGA-z_nc")
    ####################################################################################################
    # --------------------------------------------------------------------------
    # 2)  Build the *integrated* datasets and run limma on them
    # --------------------------------------------------------------------------
    # --- z-score integration ---------------------------------------------------
    common_genes_z = expr_z_or.columns.intersection(expr_z_TCGA.columns)

    expr_int_z = pd.concat(
        [expr_z_or[common_genes_z], expr_z_TCGA[common_genes_z]],
        axis=0
    )  # samples × genes

    meta_int_z = pd.concat(
        [meta_z_or.assign(batch="CGGA"), meta_z_TCGA.assign(batch="TCGA")],
        axis=0
    )

    res_int_z = run_limma(expr_int_z.T, meta_int_z, "CGGA+TCGA-z")
    volcano_plot(res_int_z , "CGGA+TCGA-z")

    # --- quantile-normalised integration --------------------------------------
    common_genes_qn = expr_qn_or.columns.intersection(expr_qn_TCGA.columns)

    expr_int_qn = pd.concat(
        [expr_qn_or[common_genes_qn], expr_qn_TCGA[common_genes_qn]],
        axis=0
    )

    meta_int_qn = pd.concat(
        [meta_qn_or.assign(batch="CGGA"), meta_qn_TCGA.assign(batch="TCGA")],
        axis=0
    )
    # 2) GLOBAL quantile normalisation  -----------------------------------------
    expr_int_qn_T = expr_int_qn.T  # genes × samples
    expr_int_qn_glob_T = hp.quantile_normalize_combined(expr_int_qn_T)
    expr_int_qn = expr_int_qn_glob_T.T

    res_int_qn = run_limma(expr_int_qn.T, meta_int_qn, "CGGA+TCGA-qn")
    volcano_plot(res_int_qn, "CGGA+TCGA-qn")
    # ----- PCA for the two integrated datasets ---------------------------------
    plot_pca_integration(expr_z_TCGA, expr_z, "TCGA", "CGGA", OUTDIR, tag="z")
    plot_pca_integration(expr_qn_TCGA, expr_qn, "TCGA", "CGGA", OUTDIR, tag="qn")
    # --------------------------------------------------------------------------
    # 3)  Collect everything and write the Excel workbook
    # --------------------------------------------------------------------------
    all_results = {
        "TCGA_z": res_z_TCGA,
        "TCGA_z_nc": res_z_TCGA_nc,
        "TCGA_qn": res_qn_TCGA,
        "CGGA_z": res_z,
        "CGGA_qn": res_qn,
        "CGGA_z_or": res_z_or,
        "CGGA_qn_or": res_qn_or,
        "Integrated_z": res_int_z,
        "Integrated_qn": res_int_qn,
    }

    save_limma_sig_to_excel(all_results, OUT_XLS)



# if __name__ == "__main__":
#     main()
#
#
def compare_gene_lists(xlsx_path: str,
                       sheet_a: str,
                       sheet_b: str,
                       gene_column: str = "Gene") -> pd.DataFrame:
    """
    Reads two sheets from an Excel workbook, extracts the gene column,
    counts totals and the intersection, and returns a tiny summary table.
    """
    s1 = set(pd.read_excel(xlsx_path, sheet_name=sheet_a)[gene_column].dropna())
    s2 = set(pd.read_excel(xlsx_path, sheet_name=sheet_b)[gene_column].dropna())

    summary = pd.DataFrame({
        "Sheet": [sheet_a, sheet_b, f"{sheet_a} ∩ {sheet_b}"],
        "n_genes": [len(s1), len(s2), len(s1 & s2)],
        "comment": ["", "", "genes in common"]
    })
    return summary


# ---- example --------------------------------------------------------------
xls = "/home/mkh062/Desktop/scratch/TCGA_project/TCGA-combined-gene-exp-data/18052025/limma_significant_genes2.xlsx"
#stats = compare_gene_lists(xls, sheet_a="TCGA_z", sheet_b="TCGA_qn")
stats = compare_gene_lists(xls, sheet_a="TCGA_z", sheet_b="Integrated_z")
print(stats.to_string(index=False))
