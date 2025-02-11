import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from lifelines import CoxPHFitter, KaplanMeierFitter
from statsmodels.stats.multitest import multipletests


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
    #df_copy.to_csv("/home/mkh062/Desktop/scratch/TCGA_project/processed_data/results/GLASS_clinical.csv")
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
    df = pd.read_csv(dataset_csv)
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

    raw_results = run_cox_for_all_genes(df,
                                        time_col='Overall.Survival.Months',
                                        event_col=event_col,
                                        covariates=covariates)

    if raw_results.empty:
        print(f"No valid genes for Cox model in {dataset_name}")
        return None

    # 5. Apply FDR correction
    results_df = fdr_correction(raw_results, p_col='p', alpha=alpha_fdr)
    print("alpha_fdr", alpha_fdr)
    # 6. Filter for significant genes
    sig_genes = results_df[(results_df['significant'] == True) & (results_df['p_adj'] < alpha_fdr)]
    sig_df = results_df[(results_df['significant'] == True) & (results_df['p_adj'] < alpha_fdr)].copy()
    # Save the Cox results
    df_sorted = results_df.sort_values(by="p_adj", ascending=True)
    results_csv = os.path.join(outdir, f"{dataset_name}_cox_results.csv")
    df_sorted.to_csv(results_csv, index=False)
    print(f"Cox results saved to {results_csv}")

    if len(sig_df) == 0:
        print(f"No significant genes (adjusted p < {alpha_fdr}) for {dataset_name}.")
        return None

    outpath = os.path.join(outdir, f"{dataset_name}_cox_significant_genes.csv")
    sig_df.to_csv(outpath, index=False)
    print(f"Saved {len(sig_df)} significant Cox genes to {outpath}")

    print(f"Found {len(sig_genes)} significant genes with adjusted p < 0.1 for {dataset_name}.")
    km_plots_dir = os.path.join(outdir, f"{dataset_name}_KM_plots")
    os.makedirs(km_plots_dir, exist_ok=True)

    # genes_of_interest = [
    #   "RAB37"
    # ]
    # km_plots_dir = '/home/mkh062/Desktop/scratch/TCGA_project/processed_data/cox_results/km_for_sig_genes'
    # for i in genes_of_interest:
    #     gene = i
    #     plot_km_for_gene(df, gene,
    #                      time_col='Overall.Survival.Months',
    #                      event_col=event_col,
    #                      outpath=os.path.join(km_plots_dir, f"{dataset_name}_{gene}_KM.png"))

    for i, row in sig_genes.iterrows():
        gene = row['gene']
        plot_km_for_gene(df, gene,
                         time_col='Overall.Survival.Months',
                         event_col=event_col,
                         outpath=os.path.join(km_plots_dir, f"{gene}_KM.png"))

def summarize_threshold_overlaps(consortium, cox_sig_genes_file, outdir,
                                 thresholds_fixed=[24, 60],
                                 thresholds_range=[(24,48), (24,60)]):
    """
    Reads the Cox-significant genes for a given consortium,
    reads each limma-significant gene file for multiple thresholds,
    and produces a summary of how many genes are in each set and their overlap.

    Returns a dataframe with columns:
        ['consortium', 'threshold_label', 'n_limma', 'n_cox', 'n_common']

    Example threshold labels:
       "24-months" (single threshold)
       "60-months"
       "24-48-months" (range)
       "24-60-months"
    """
    results = []

    # 1. Load the Cox significant genes (if file exists)
    if not os.path.isfile(cox_sig_genes_file):
        print(f"[WARN] Cox file not found: {cox_sig_genes_file}")
        return pd.DataFrame()

    cox_df = pd.read_csv(cox_sig_genes_file)
    if 'gene' not in cox_df.columns:
        print(f"[WARN] 'gene' column not found in {cox_sig_genes_file}")
        return pd.DataFrame()

    cox_genes = set(cox_df['gene'].unique())
    n_cox = len(cox_genes)

    # 2. Build each limma file suffix and read them
    #    (a) fixed thresholds
    for t in thresholds_fixed:
        file_suffix = f"{consortium}_{t}-months"
        limma_file = os.path.join(outdir, f"significant_genes_{file_suffix}.csv")

        threshold_label = f"{t}-months"
        if not os.path.isfile(limma_file):
            # skip if file doesn't exist
            continue

        limma_df = pd.read_csv(limma_file)
        if 'gene' not in limma_df.columns:
            # adjust to your actual column name if needed
            continue

        limma_genes = set(limma_df['gene'].unique())
        n_limma = len(limma_genes)
        n_common = len(cox_genes.intersection(limma_genes))

        results.append({
            'consortium': consortium,
            'threshold_label': threshold_label,
            'n_limma': n_limma,
            'n_cox': n_cox,
            'n_common': n_common
        })

    #    (b) range thresholds
    for (low, high) in thresholds_range:
        file_suffix = f"{consortium}_{low}-{high}-months"
        limma_file = os.path.join(outdir, f"significant_genes_{file_suffix}.csv")

        threshold_label = f"{low}-{high}-months"
        if not os.path.isfile(limma_file):
            continue

        limma_df = pd.read_csv(limma_file)
        if 'gene' not in limma_df.columns:
            continue

        limma_genes = set(limma_df['gene'].unique())
        n_limma = len(limma_genes)
        n_common = len(cox_genes.intersection(limma_genes))

        results.append({
            'consortium': consortium,
            'threshold_label': threshold_label,
            'n_limma': n_limma,
            'n_cox': n_cox,
            'n_common': n_common
        })

    # Create a final DataFrame summary
    summary_df = pd.DataFrame(results)
    return summary_df


if __name__ == "__main__":
    # Paths to your processed CSVs
    tcga_csv = "/home/mkh062/Desktop/scratch/TCGA_project/processed_data/TCGA_processed.csv"
    cgga_csv = "/home/mkh062/Desktop/scratch/TCGA_project/processed_data/CGGA_processed.csv"
    glass_csv = "/home/mkh062/Desktop/scratch/TCGA_project/processed_data/GLASS_processed.csv"
    #
    # # Output folder
    outdir = "/home/mkh062/Desktop/scratch/TCGA_project/processed_data/cox_results"
    os.makedirs(outdir, exist_ok=True)
    #
    # # Analyze each dataset
    analyze_dataset(tcga_csv, "TCGA", outdir, alpha_fdr=0.10)
    analyze_dataset(glass_csv, "GLASS", outdir, alpha_fdr=0.4)
    analyze_dataset(cgga_csv, "CGGA", outdir, alpha_fdr=0.5)

    # # 2. Summaries of how many genes from LIMMA results overlap with the Cox results
    # thresholds_fixed = [24, 60]
    # thresholds_range = [(24, 48), (24, 60)]
    # results_dir1 = '/home/mkh062/Desktop/scratch/TCGA_project/processed_data/cox_results'
    # results_dir = '/home/mkh062/Desktop/scratch/TCGA_project/processed_data/results'
    # summary_frames = []
    # for consortium in ["TCGA", "CGGA", "GLASS"]:
    #     cox_sig_file = os.path.join(results_dir1, f"{consortium}_cox_significant_genes.csv")
    #     summary_df = summarize_threshold_overlaps(
    #         consortium=consortium,
    #         cox_sig_genes_file=cox_sig_file,
    #         outdir=results_dir,
    #         thresholds_fixed=thresholds_fixed,
    #         thresholds_range=thresholds_range
    #     )
    #     summary_frames.append(summary_df)
    #
    # # Combine all datasets' summaries
    # final_summary = pd.concat(summary_frames, ignore_index=True)
    # summary_outpath = os.path.join(results_dir, "cox_limma_summary.csv")
    # final_summary.to_csv(summary_outpath, index=False)
    #
    # print("\n=== Summary of limma & Cox overlaps ===")
    # print(final_summary)
    # print(f"\nSummary saved to: {summary_outpath}")
    #
    # print("Done!")
