import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


def main():
    # define base directories
    BASE_DIR = "/home/mkh062/Desktop/scratch/TCGA_project"
    PROCESSED_DATA_DIR = os.path.join(BASE_DIR, "processed_data")
    RESULTS_DIR = os.path.join(PROCESSED_DATA_DIR, "results")

    os.makedirs(PROCESSED_DATA_DIR, exist_ok=True)
    os.makedirs(RESULTS_DIR, exist_ok=True)
    # File paths

    cgga = pd.read_csv("/home/mkh062/Desktop/scratch/TCGA_project/processed_data/CGGA_processed.csv")
    tcga = pd.read_csv("/home/mkh062/Desktop/scratch/TCGA_project/processed_data/TCGA_processed.csv")
    glass = pd.read_csv("/home/mkh062/Desktop/scratch/TCGA_project/processed_data/GLASS_processed.csv")
    dfs = {
        "cgga": cgga,
        "tcga": tcga,
        "glass": glass
    }
    all_genes_of_interest = ['TRAF3', 'SDCBP', 'UBE2H', 'LAP3', 'TEX2', 'WWP1', 'TRIP4', 'SCRN1', 'USP32',
                             'SNX10', 'GRB2', 'RAB21']
    # ----------------------------------------------------
    threshold_scenarios = {
        'T24': 24,
        'T60': 60,
        'T24_60': (24, 60),  # special case
        'T24_48': (24, 48)  # special case
    }

    def gather_significant_genes(outdir, consortium, threshold_label):
        fname = os.path.join(outdir, f"significant_genes_{consortium}_{threshold_label}.csv")
        if not os.path.exists(fname):
            return set()
        df = pd.read_csv(fname, index_col=0)
        if df.empty:
            return set()
        return set(df.index)

    def compute_log2_fc_for_dataset(
            merged, genes, threshold_key='T24'
    ):


        # Filter out if threshold is the "range" scenario
        if threshold_key in ['T24_60', 'T24_48']:
            # Get lower and upper bounds
            lower_bound, upper_bound = threshold_scenarios[threshold_key]
            mask = (
                    ((merged['Overall.Survival.Months'] <= lower_bound) & (merged['case.vital.status'] != 'alive')) |
                    (merged['Overall.Survival.Months'] > upper_bound)
            )
            merged = merged[mask].copy()
            # Define "Group"
            merged['Group'] = np.where(merged['Overall.Survival.Months'] <= lower_bound, 'short', 'long')

        else:
            # For T24 or T60 thresholds
            thr_value = threshold_scenarios[threshold_key]
            mask = (
                    (merged['Overall.Survival.Months'] <= thr_value) & (merged['case.vital.status'] != 'alive') |
                    (merged['Overall.Survival.Months'] > thr_value)
            )
            merged = merged[mask].copy()
            merged['Group'] = np.where(merged['Overall.Survival.Months'] <= thr_value, 'short', 'long')

        # Compute log2 FC for each gene
        out = {}
        for g in genes:
            # mean expr in short group
            short_mean = merged.loc[merged['Group'] == 'short', g].mean()
            # mean expr in long group
            long_mean = merged.loc[merged['Group'] == 'long', g].mean()
            log2_fc = (long_mean - short_mean)
            out[g] = log2_fc

        return out

    BASE_DIR = "/home/mkh062/Desktop/scratch/TCGA_project"
    PROCESSED_DATA_DIR = os.path.join(BASE_DIR, "processed_data")
    RESULTS_DIR = os.path.join(PROCESSED_DATA_DIR, "results")


    logFCs = {key: {} for key in threshold_scenarios.keys()}
    for threshold_key in threshold_scenarios:
        for dataset_name in dfs:
            df = dfs[dataset_name]
            logFC_dict = compute_log2_fc_for_dataset(
                df, all_genes_of_interest, threshold_key
            )
            logFCs[threshold_key][dataset_name] = logFC_dict

    fig, axes = plt.subplots(nrows=4, ncols=12, figsize=(25, 15), sharey=True)
    threshold_order = ['T24', 'T60', 'T24_60', 'T24_48']


    for row_idx, tkey in enumerate(threshold_order):
        for col_idx, gene_name in enumerate(all_genes_of_interest):
            ax = axes[row_idx, col_idx]

            # Collect the logFC for each dataset
            values = []
            dsets = ["cgga", "tcga", "glass"]
            for d in dsets:
                val = logFCs[tkey][d][gene_name]
                #val = FCs[tkey][d][gene_name]
                values.append(val)

            # Make bar plot
            ax.bar(dsets, values, color='skyblue', edgecolor='k')
            ax.set_title(f"{gene_name} - {tkey}", fontsize=9)
            ax.axhline(0, color='black', linewidth=1)

            # Adjust x-axis for the bottom row
            if row_idx == 3:  # Last row
                ax.set_xticklabels(dsets, rotation=25)
            else:
                ax.set_xticklabels([])  # Remove x-ticks for other rows

            ax.set_ylabel("Log Fold‐Change" if col_idx == 0 else "")
            ax.set_xlabel("")  # Remove x-labels

    plt.tight_layout()
    plt.show()


    logFCs = {key: {} for key in threshold_scenarios.keys()}
    for threshold_key in threshold_scenarios:
        if threshold_key =='T24':
            t_label = '24-months'
        elif threshold_key =='T60':
            t_label = '60-months'
        elif threshold_key =='T24_60':
            t_label = '24-60-months'
        else:
            t_label = '24-48-months'

        sig_tcga = gather_significant_genes(RESULTS_DIR, "TCGA", t_label)
        sig_cgga = gather_significant_genes(RESULTS_DIR, "CGGA", t_label)
        sig_glass = gather_significant_genes(RESULTS_DIR, "GLASS", t_label)
        union_genes = sig_tcga.union(sig_cgga).union(sig_glass)

        for dataset_name in dfs:
            print("expression_dfs[dataset_name].", dfs[dataset_name].head())
            missing_genes = [gene for gene in union_genes if gene not in dfs[dataset_name].columns]
            if missing_genes:
                print(f"Missing genes in {dataset_name}: {missing_genes}")
                union_genes = union_genes - set(missing_genes)

        print("union_genes 1", len(union_genes))
        for dataset_name in dfs:
            df = dfs[dataset_name]
            print("dataset_name", dataset_name)
            logFC_dict = compute_log2_fc_for_dataset(
                df, union_genes, threshold_key
            )
            logFCs[threshold_key][dataset_name] = logFC_dict

    from itertools import combinations

    datasets = ["cgga", "tcga", "glass"]
    dataset_pairs = list(combinations(datasets, 2))

    nrows = len(threshold_scenarios)
    ncols = len(dataset_pairs)
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(ncols * 5, nrows *4), sharey=True)
    axes = axes.flatten()
    # Determine global min and max for symmetric axis limits
    global_min = float('inf')
    global_max = float('-inf')
    for idx, (tkey, (ds1, ds2)) in enumerate([(t, pair) for t in threshold_scenarios for pair in dataset_pairs]):
        ax = axes[idx]
        logFC_ds1 = pd.Series(logFCs[tkey][ds1])
        logFC_ds2 = pd.Series(logFCs[tkey][ds2])
        combined = pd.DataFrame({
            'logFC_ds1': logFC_ds1,
            'logFC_ds2': logFC_ds2
        }).dropna()
        if not combined.empty:
            current_min = min(combined['logFC_ds1'].min(), combined['logFC_ds2'].min())
            current_max = max(combined['logFC_ds1'].max(), combined['logFC_ds2'].max())
            global_min = min(global_min, current_min)
            global_max = max(global_max, current_max)
        global_limit = max(abs(global_min), abs(global_max))
        global_min, global_max = -global_limit, global_limit
        ax.scatter(combined['logFC_ds1'], combined['logFC_ds2'], alpha=0.7, color='blue', edgecolor='k')
        ax.axhline(0, color='gray', linestyle='--', linewidth=1)
        ax.axvline(0, color='gray', linestyle='--', linewidth=1)
        # ax.set_title(f"{tkey}: {ds1} vs {ds2}", fontsize=6)
        ax.set_title(f" ")
        ax.set_xlabel(f"log2 Fold Change ({ds1})")
        ax.set_ylabel(f"log2 Fold Change ({ds2})")
        ax.set_xlim(global_min, global_max)
        ax.set_ylim(global_min, global_max)
    # Remove unused axes
    for i in range(len(threshold_scenarios) * len(dataset_pairs), len(axes)):
        fig.delaxes(axes[i])

    plt.suptitle("Log2 Fold Change Comparisons Across Dataset Pairs", fontsize=12)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
