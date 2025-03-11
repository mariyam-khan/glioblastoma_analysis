import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from lifelines import CoxPHFitter, KaplanMeierFitter
from statsmodels.stats.multitest import multipletests
import os
# R-related imports
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
pandas2ri.activate()


def limma(data, threshold, consortium, threshold_type="fixed", outdir="./results"):
    """
    Runs a differential expression analysis in R using limma.
    """
    os.makedirs(outdir, exist_ok=True)
    data = data.copy()

    if threshold_type == "fixed":
        data = data[~((data['Overall.Survival.Months'] <= threshold) & (data['case.vital.status'] == 0))]
        print("less than thresh and alive", data[((data['Overall.Survival.Months'] <= threshold) & (data['case.vital.status'] == 0))])
        file_suffix = f"{consortium}_{threshold}-months"
        median_survival_r = threshold
    else:
        # range-based threshold
        low, high = threshold
        data = data[~((data['Overall.Survival.Months'] <= low) & (data['case.vital.status'] == 0))]
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

def create_event_column(df, status_col='case.vital.status'):
    """
    Transforms 'case.vital.status' text (e.g. 'alive'/'dead') into a binary event column.
    Returns the modified dataframe and the name of the created column.
    """
    df = df.copy()
    unique_values = df[status_col].unique()
    print(f"Unique values in '{status_col}':", unique_values)
    # 1 = event occurred (dead), 0 = censored (alive)
    if pd.api.types.is_numeric_dtype(df[status_col]):
        # rename that column to 'event'
        df.rename(columns={status_col: 'event'}, inplace=True)
        #df['event'] = df[status_col].map({0.0: 0, 1.0: 1}).astype(int)
        return df, 'event'
    else:
        # else we assume string 'alive'/'dead'.  Map to 0/1
        df['event'] = df[status_col].map({'alive': 0, 'dead': 1}).astype(int)
        return df, 'event'


def run_cox_for_all_genes(df, time_col='Overall.Survival.Months', event_col='event',
                          covariates=['Diagnosis.Age', 'Sex']):
    """
    Runs CoxPH for each gene in df, adjusting for the given covariates.

    Parameters:
    - df: DataFrame containing survival time, event, covariates, and gene expression columns
    - time_col: name of the column containing survival time
    - event_col: name of the column containing the event (1=dead, 0=alive)
    - covariates: list of columns (besides genes) to include as covariates in the Cox model
    - alpha: significance threshold for the *unadjusted* p-value (only used to speed up if we want).

    Returns:
    - results_df: DataFrame with columns [gene, coef, p, hazard_ratio, ...]
      Each row corresponds to one gene’s Cox regression summary.
    """
    # Identify gene columns.
    # We'll exclude ID columns, time/event columns, and known covariates from the regression.
    excluded_cols = set([time_col, event_col] + covariates)
    gene_cols = [c for c in df.columns if c not in excluded_cols]

    # (Optional) If you have other non-gene columns (e.g. 'Sample.ID'), exclude them as well:
    # for example:
    additional_exclude = ['Sample.ID', 'case.vital.status']  # adjust as needed
    excluded_cols.update(additional_exclude)

    gene_cols = [c for c in gene_cols if c not in excluded_cols]

    # Prepare a results container
    records = []

    cph = CoxPHFitter()
    print("df[[time_col, event_col, gene] + covariates]", df[[time_col, event_col] + covariates])
    df_copy = df[[time_col, event_col] + covariates].copy()
    # For each gene, fit a univariate or multivariate Cox model
    #   - If you want to adjust for the covariates, you must pass them in the fit dataframe.
    for gene in gene_cols:
        # Build a temporary dataframe with [time_col, event_col, gene] + covariates
        temp_df = df[[time_col, event_col, gene] + covariates].dropna()

        if temp_df.shape[0] < 5:
            # skip if not enough samples
            continue

        try:
            cph.fit(temp_df, duration_col=time_col, event_col=event_col)
            # We'll extract the row for this gene
            # The cox model summary has an index that includes the covariates. We only want the gene row.
            summary = cph.summary
            # The row name for the gene is exactly `gene` in cph.summary
            if gene in summary.index:
                row = summary.loc[gene]
                records.append({
                    'gene': gene,
                    'coef': row['coef'],
                    'p': row['p'],
                    'exp(coef)': row['exp(coef)'],  # hazard ratio
                    'coef lower 95%': row['coef lower 95%'],
                    'coef upper 95%': row['coef upper 95%'],
                    'z': row['z'],
                    'gene_n': temp_df.shape[0]  # how many patients included for that gene
                })
        except:
            # occasionally might fail to converge for certain genes
            continue

    results_df = pd.DataFrame(records)
    return results_df


def fdr_correction(results_df, p_col='p', alpha=0.10):
    """
    Apply multiple-testing correction (FDR / Benjamini-Hochberg).
    Returns the DataFrame with a new column 'p_adj' and a boolean 'significant' column.

    Parameters:
    - results_df: DataFrame with a column p_col for p-values
    - p_col: name of the p-value column in results_df
    - alpha: the threshold for adjusted p-value significance (e.g. 0.10)
    """
    df = results_df.copy()
    if df.empty:
        df['p_adj'] = np.nan
        df['significant'] = False
        return df
    print("alpha", alpha)
    corrected = multipletests(df[p_col].values, alpha=alpha, method='fdr_bh')
    df['p_adj'] = corrected[1]
    df['significant'] = corrected[0]
    return df


def plot_km_for_gene(df, gene, time_col='Overall.Survival.Months', event_col='event', outpath=None):
    """
    Plots Kaplan-Meier curve for 'High' vs. 'Low' expression groups,
    splitting by median expression for the specified gene.
    Saves plot to outpath if provided.
    """
    # Split by median
    median_expr = df[gene].median()
    low_group = df[df[gene] <= median_expr]
    high_group = df[df[gene] > median_expr]

    kmf_low = KaplanMeierFitter()
    kmf_high = KaplanMeierFitter()

    plt.figure(figsize=(5, 4))

    kmf_low.fit(durations=low_group[time_col], event_observed=low_group[event_col], label='Low ' + gene)
    ax = kmf_low.plot(ci_show=False)

    kmf_high.fit(durations=high_group[time_col], event_observed=high_group[event_col], label='High ' + gene)
    kmf_high.plot(ci_show=False, ax=ax)

    plt.title(f"KM curve for {gene}")
    plt.xlabel("Time (Months)")
    plt.ylabel("Survival Probability")
    plt.tight_layout()
    if outpath:
        plt.savefig(outpath, dpi=120)
        plt.close()
    else:
        plt.show()


def analyze_dataset(dataset_csv, dataset_name, outdir, alpha_fdr):
    """
    Loads the dataset CSV, runs Cox analysis on all genes,
    does FDR correction, and plots KM for significant genes.

    dataset_csv: e.g. "TCGA_processed.csv"
    dataset_name: e.g. "TCGA"
    outdir: folder to store results
    """
    os.makedirs(outdir, exist_ok=True)

    # 1. Load data
    df = dataset_csv
    print(f"\n--- Analyzing {dataset_name} ---")
    print(f"Initial shape: {df.shape}")

    # 2. Make sure we have the columns we need
    #    'Overall.Survival.Months', 'case.vital.status' in {alive, dead}
    #    'Diagnosis.Age', (optionally 'Sex'), ...
    required_cols = ['Overall.Survival.Months', 'case.vital.status', 'Diagnosis.Age']
    for rc in required_cols:
        if rc not in df.columns:
            raise ValueError(f"Missing required column '{rc}' in {dataset_name} dataset.")

    # 3. Create the event column
    df, event_col = create_event_column(df, status_col='case.vital.status')

    # Let's see if 'Sex' is in columns:
    covariates = ['Diagnosis.Age']
    if 'Sex' in df.columns:
        # Quick numeric code (male=1, female=0 for example)
        unique_values = df['Sex'].unique()
        print(f"Unique values in 'Sex':", unique_values)
        df['Sex'] = df['Sex'].map({'M': 1, 'Male': 1, 'male': 1, 'F': 0, 'Female': 0, 'female': 0}).fillna(0).astype(int)
        covariates.append('Sex')

    # raw_results = run_cox_for_all_genes(df,
    #                                     time_col='Overall.Survival.Months',
    #                                     event_col=event_col,
    #                                     covariates=covariates)
    #
    # if raw_results.empty:
    #     print(f"No valid genes for Cox model in {dataset_name}")
    #     return None
    #
    # # 5. Apply FDR correction
    # results_df = fdr_correction(raw_results, p_col='p', alpha=alpha_fdr)
    # print("alpha_fdr", alpha_fdr)
    # # 6. Filter for significant genes
    # sig_genes = results_df[(results_df['significant'] == True) & (results_df['p_adj'] < alpha_fdr)]
    # sig_df = results_df[(results_df['significant'] == True) & (results_df['p_adj'] < alpha_fdr)].copy()
    # outpath = os.path.join(outdir, f"{dataset_name}_cox_significant_genes.csv")
    # sig_df.to_csv(outpath, index=False)
    # print(f"Saved {len(sig_df)} significant Cox genes to {outpath}")
    # # Save the Cox results
    # df_sorted = results_df.sort_values(by="p_adj", ascending=True)
    # results_csv = os.path.join(outdir, f"{dataset_name}_cox_results.csv")
    # df_sorted.to_csv(results_csv, index=False)
    # print(f"Cox results saved to {results_csv}")
    #
    # if len(sig_df) == 0:
    #     print(f"No significant genes (adjusted p < {alpha_fdr}) for {dataset_name}.")
    #     return None
    #
    # print(f"Found {len(sig_genes)} significant genes with adjusted p < 0.1 for {dataset_name}.")
    # km_plots_dir = os.path.join(outdir, f"{dataset_name}_KM_plots")
    # os.makedirs(km_plots_dir, exist_ok=True)
    #
    # for i, row in sig_genes.iterrows():
    #     gene = row['gene']
    #     plot_km_for_gene(df, gene,
    #                      time_col='Overall.Survival.Months',
    #                      event_col=event_col,
    #                      outpath=os.path.join(km_plots_dir, f"{gene}_KM.png"))

    km_plots_dir = '/home/mkh062/Desktop/scratch/TCGA_project/processed_data/TCGA_GLASS_integrated/cox_results'
    genes = ['DKK3', 'COL22A1', 'RAB37', 'MSN']
    for i in genes:
        gene = i
        plot_km_for_gene(df, gene,
                         time_col='Overall.Survival.Months',
                         event_col=event_col,
                         outpath=os.path.join(km_plots_dir, f"TCGA_GLASS_{gene}_KM.png"))

def TCGA_GLASS_integrate():
    # Load the datasets
    glass = pd.read_csv('/home/mkh062/Desktop/scratch/TCGA_project/processed_data/GLASS_processed.csv')
    glass['Platform'] = 'GLASS'
    glass['case.vital.status'] = glass['case.vital.status'].map({'alive': 0, 'dead': 1}).astype(int)
    print("glass['case.vital.status']", glass['case.vital.status'].values)
    tcga = pd.read_csv('/home/mkh062/Desktop/scratch/TCGA_project/processed_data/TCGA_processed.csv')
    tcga['Platform'] = 'TCGA'
    tcga['case.vital.status'] = tcga['case.vital.status'].map({'alive': 0, 'dead': 1}).astype(int)
    # Identify the gene columns (excluding metadata columns)
    metadata_cols = ['Sample.ID', 'Overall.Survival.Months', 'Diagnosis.Age', 'Sex', 'case.vital.status', 'Platform']

    glass_gene_cols = set(glass.columns) - set(metadata_cols)
    tcga_gene_cols = set(tcga.columns) - set(metadata_cols)

    # Find common genes in both datasets
    common_genes = glass_gene_cols.intersection(tcga_gene_cols)

    # Keep only common genes and metadata columns
    glass_filtered = glass[list(metadata_cols) + list(common_genes)]
    tcga_filtered = tcga[list(metadata_cols) + list(common_genes)]

    # Combine the datasets
    combined_data = pd.concat([glass_filtered, tcga_filtered], ignore_index=True)
    sample_ids = combined_data['Sample.ID']
    platforms = combined_data['Platform']

    # The gene columns are everything except 'Sample ID' and 'platform'
    gene_cols = [c for c in combined_data.columns if c not in ['Sample.ID','Overall.Survival.Months','Diagnosis.Age','Sex','case.vital.status','Platform']]

    # Extract the expression matrix
    expr_matrix = combined_data[gene_cols].copy()

    # 3. PCA before normalization
    pca = PCA(n_components=2)
    pca_scores = pca.fit_transform(expr_matrix.fillna(0))  # handle missing data somehow
    pc_df = pd.DataFrame(pca_scores, columns=['PC1', 'PC2'])
    pc_df['Sample.ID'] = sample_ids.values
    pc_df['Platform'] = platforms.values

    # Plot the PCA before normalization
    plt.figure(figsize=(8, 6))
    sns.scatterplot(x='PC1', y='PC2', hue='Platform', data=pc_df, alpha=0.7)
    plt.title('PCA')
    plt.legend()
    output_dir ='/home/mkh062/Desktop/scratch/TCGA_project/processed_data/TCGA_GLASS_integrated/'
    pca_plot_path = os.path.join(output_dir, "PCA_plot.png")
    #plt.savefig(pca_plot_path, dpi=300, bbox_inches='tight')
    plt.show()
    combined_data = combined_data.drop(columns=['Platform'])
    return combined_data


def main():
    output_dir = "/home/mkh062/Desktop/scratch/TCGA_project/processed_data/TCGA_GLASS_integrated/"
    os.makedirs(output_dir, exist_ok=True)  # Ensure directory exists
    outdir = "/home/mkh062/Desktop/scratch/TCGA_project/processed_data/TCGA_GLASS_integrated/cox_results"
    os.makedirs(outdir, exist_ok=True)
    tcga_glass = TCGA_GLASS_integrate()
    print("tcga_glass", tcga_glass)
    #thresholds_fixed = [24, 60]
    thresholds_fixed = [36, 48]
    thresholds_range = [(24, 48), (24, 60)]

    for threshold in thresholds_fixed:
        for ds_data, ds_name in [(tcga_glass, "GLASS_TCGA")]:
            res = limma(ds_data.copy(),
                                       threshold=threshold,
                                       consortium=ds_name,
                                       threshold_type="fixed",
                                       outdir=output_dir)
    #
    # # range
    # for (low,high) in thresholds_range:
    #     for ds_data, ds_name in [(tcga_glass, "GLASS_TCGA")]:
    #         res = limma(ds_data.copy(),
    #                                    threshold=(low,high),
    #                                    consortium=ds_name,
    #                                    threshold_type="range",
    #                                    outdir=output_dir)
    #analyze_dataset(tcga_glass, "GLASS_TCGA", outdir, alpha_fdr=0.1)

if __name__=="__main__":
    main()

