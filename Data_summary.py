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

    all_genes_of_interest = ['TRAF3', 'SDCBP', 'UBE2H', 'LAP3', 'TEX2', 'WWP1', 'TRIP4', 'SCRN1', 'USP32',
                             'SNX10', 'GRB2', 'RAB21']

    unique_genes = [
        'TRIP4', 'EFEMP2', 'EFNB2', 'TRAF3', 'DYNLT3', 'MSN', 'C9orf64',
        'HSP90B1', 'VASN', 'PARVA', 'SWAP70', 'PALM2-AKAP2', 'LGALS8',
        'PDIA4', 'TBC1D1', 'RPAP3', 'DIRAS3', 'RAB37', 'AMIGO3',
         'TCHP',
         'EXOSC10','GRM5',
        'C1orf94', 'KCNV1', 'HCN1',  'STK36', 'MARCH8', 'HORMAD2', 'BARX1',
         'SDCBP', 'UBE2H', 'LAP3', 'TEX2', 'WWP1', 'SCRN1', 'USP32',
                             'SNX10', 'GRB2', 'RAB21'
    ]

    clinical_325 = "/home/mkh062/Desktop/scratch/TCGA_project/CGGA-RAW DATA/CGGA.mRNAseq_325_clinical.20200506.txt"
    clinical_693 = "/home/mkh062/Desktop/scratch/TCGA_project/CGGA-RAW DATA/CGGA.mRNAseq_693_clinical.20200506.txt"

    expr_325 = "/home/mkh062/Desktop/scratch/TCGA_project/CGGA.mRNAseq_325.RSEM-genes.20200506.txt/CGGA.mRNAseq_325.RSEM-genes.20200506.txt"
    expr_693 = "/home/mkh062/Desktop/scratch/TCGA_project/CGGA.mRNAseq_693.RSEM-genes.20200506.txt/CGGA.mRNAseq_693.RSEM-genes.20200506.txt"

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

    c325.rename(columns={
        'Censor (alive=0; dead=1)': 'case.vital.status'
    }, inplace=True)
    c693.rename(columns={
        'Censor (alive=0; dead=1)': 'case.vital.status'
    }, inplace=True)

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

    expr_325 = df_325_t[['Sample.ID'] + all_genes_of_interest].copy()
    expr_693 = df_693_t[['Sample.ID'] + all_genes_of_interest].copy()
    expr_325_2 = df_325_t[['Sample.ID'] + unique_genes].copy()
    expr_693_2 = df_693_t[['Sample.ID'] + unique_genes].copy()


    clinical_tcga = pd.read_excel('/home/mkh062/Desktop/scratch/TCGA_project/Compiled data_TCGA_CGGA_GLASS_jh080424.xlsx', sheet_name='TCGA')
    clinical_tcga.rename(columns={
        'Overall Survival Status': 'case.vital.status'
    }, inplace=True)
    expr_tcga = pd.read_csv('/home/mkh062/Desktop/scratch/TCGA_project/GBMLGG_EB_RmDiffFullGenesRanRmDup.csv', sep=',')
    common_ids = set(expr_tcga['samples']).intersection(set(clinical_tcga['Sample ID']))
    expr = expr_tcga[expr_tcga['samples'].isin(common_ids)].copy()
    expr.rename(columns={'samples':'Sample.ID'}, inplace=True)

    # Keep only RNAseq + Agilent
    valid_platforms = ['RNAseq','agilent']
    expr_tcga1 = expr[expr['platform'].isin(valid_platforms)].copy()

    expr_tcga = expr_tcga1[['Sample.ID'] + all_genes_of_interest].copy()
    expr_tcga_2 = expr_tcga1[['Sample.ID'] + unique_genes].copy()

    clinical_glass = pd.read_excel('/home/mkh062/Desktop/scratch/TCGA_project/Compiled data_TCGA_CGGA_GLASS_jh080424.xlsx', sheet_name='GLASS')
    expr_glass = pd.read_csv( "/home/mkh062/Desktop/scratch/TCGA_project/gene_tpm_matrix_all_samples.tsv", sep='\t')

    expr_t = expr_glass.set_index('Gene_symbol').transpose().reset_index()
    new_columns = ['Sample ID'] + expr_glass['Gene_symbol'].tolist()
    expr_t.columns = new_columns
    def transform_sample_name_GLASS(name):
        if 'TP' in name:
            name = name.split('TP')[0] + 'TP'
        name = name.replace('.', '-')
        return name
    # Apply sample name transformation (assuming transform_sample_name is defined elsewhere)
    expr_t['Sample ID'] = expr_t['Sample ID'].apply(transform_sample_name_GLASS)
    common_ids = pd.Series(list(set(expr_t['Sample ID']).intersection(set(clinical_glass['Sample ID'])))).tolist()
    expr_t = expr_t[expr_t['Sample ID'].isin(common_ids)]
    expr_t.reset_index(drop=True, inplace=True)
    expr_t.rename(columns={'Sample ID': 'Sample.ID'}, inplace=True)
    expr_glass = expr_t[['Sample.ID'] + all_genes_of_interest].copy()
    expr_glass_2 = expr_t[['Sample.ID'] + unique_genes].copy()
    clinical_glass.rename(columns={
        'case_vital_status': 'case.vital.status'
    }, inplace=True)

    expression_dfs = {
        "325": expr_325,
        "693": expr_693,
        "tcga": expr_tcga,
        "glass": expr_glass
    }



    # Also put clinical info in dict:
    clinical_dfs = {
        "325": c325.rename(columns={'OS': 'OS_months'}),
        "693": c693.rename(columns={'OS': 'OS_months'}),
        "tcga": clinical_tcga.rename(columns={'Sample ID':'Sample.ID', 'Overall Survival (Months)': 'OS_months'}),
        "glass": clinical_glass.rename(columns={'Sample ID':'Sample.ID', 'case_overall_survival_mo': 'OS_months'})
    }

    # ----------------------------------------------------
    # 1) Plot raw expression distributions: 4x12
    # ----------------------------------------------------
    fig, axes = plt.subplots(nrows=4, ncols=12, figsize=(36, 12), sharex=False, sharey=False)

    for row_idx, dataset_name in enumerate(["325", "693", "tcga", "glass"]):
        df_expr = expression_dfs[dataset_name]
        for col_idx, gene_name in enumerate(all_genes_of_interest):
            ax = axes[row_idx, col_idx]
            sns.histplot(
                data=df_expr[gene_name],
                kde=True,  # overlay a density estimate
                ax=ax
            )
            ax.set_title(f"{gene_name} - {dataset_name}", fontsize=12)
            ax.set_xlabel('')

    plt.tight_layout()
    plt.show()

    # ----------------------------------------------------
    # 2) Compute log2 fold changes at 3 thresholds
    #    Thresholds:
    #    (A) 24 months
    #    (B) 60 months
    #    (C) Range: keep only OS<=24 or OS>60
    # ----------------------------------------------------
    threshold_scenarios = {
        'T24': 24,
        'T60': 60,
        'T24_60': (24, 60),  # special case
        'T24_48': (24, 48)  # special case
    }

    def compute_log2_fc_for_dataset(
            expr_df, clin_df, genes, threshold_key='T24'
    ):
        """
        For a single dataset, merges clinical & expression,
        splits into low/high groups for survival threshold,
        returns a dict { gene_name: log2FC }.
        """
        merged = pd.merge(expr_df, clin_df, how='inner', on='Sample.ID')

        # Filter out if threshold is the "range" scenario
        if threshold_key in ['T24_60', 'T24_48']:
            # Get lower and upper bounds
            lower_bound, upper_bound = threshold_scenarios[threshold_key]
            mask = (
                    ((merged['OS_months'] <= lower_bound) & (merged['case.vital.status'] != 'alive')) |
                    (merged['OS_months'] > upper_bound)
            )
            merged = merged[mask].copy()
            # Define "Group"
            merged['Group'] = np.where(merged['OS_months'] <= lower_bound, 'short', 'long')

        else:
            # For T24 or T60 thresholds
            thr_value = threshold_scenarios[threshold_key]
            mask = (
                    (merged['OS_months'] <= thr_value) & (merged['case.vital.status'] != 'alive') |
                    (merged['OS_months'] > thr_value)
            )
            merged = merged[mask].copy()
            merged['Group'] = np.where(merged['OS_months'] <= thr_value, 'short', 'long')


        # Compute log2 FC for each gene
        out = {}
        f_change = {}
        for g in genes:
            # mean expr in short group
            short_mean = merged.loc[merged['Group'] == 'short', g].mean()
            # mean expr in long group
            long_mean = merged.loc[merged['Group'] == 'long', g].mean()
            f = (long_mean / short_mean)
            if short_mean > 0 and long_mean > 0:
                fc = (long_mean / short_mean)
                log2_fc = np.log2(fc)
            else:

                log2_fc = np.nan
            out[g] = log2_fc
            f_change[g] = f
        return out

    logFCs = {key: {} for key in threshold_scenarios.keys()}
    for threshold_key in threshold_scenarios:
        for dataset_name in expression_dfs:
            df_expr = expression_dfs[dataset_name]
            df_clin = clinical_dfs[dataset_name]
            print("dataset_name", dataset_name)
            logFC_dict = compute_log2_fc_for_dataset(
                df_expr, df_clin, all_genes_of_interest, threshold_key
            )
            logFCs[threshold_key][dataset_name] = logFC_dict

    fig, axes = plt.subplots(nrows=4, ncols=12, figsize=(25, 15), sharey=True)
    threshold_order = ['T24', 'T60', 'T24_60', 'T24_48']


    for row_idx, tkey in enumerate(threshold_order):
        for col_idx, gene_name in enumerate(all_genes_of_interest):
            ax = axes[row_idx, col_idx]

            # Collect the logFC for each dataset
            values = []
            dsets = ["325", "693", "tcga", "glass"]
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

    threshold_scenarios = {
        'T24_60': (24, 60),
        'T24_48': (24, 48)  # Added this
    }
    # Define the list of unique genes

    expression_dfs = {
        "325": expr_325_2,
        "693": expr_693_2,
        "tcga": expr_tcga_2,
        "glass": expr_glass_2
    }
    for dataset_name in expression_dfs:
        print("expression_dfs[dataset_name].", expression_dfs[dataset_name].head())
        missing_genes = [gene for gene in unique_genes if gene not in expression_dfs[dataset_name].columns]
        if missing_genes:
            print(f"Missing genes in {dataset_name}: {missing_genes}")

    logFCs = {key: {} for key in threshold_scenarios.keys()}
    for threshold_key in threshold_scenarios:
        for dataset_name in expression_dfs:
            df_expr = expression_dfs[dataset_name]
            df_clin = clinical_dfs[dataset_name]
            print("dataset_name", dataset_name)
            logFC_dict = compute_log2_fc_for_dataset(
                df_expr, df_clin, unique_genes, threshold_key
            )
            logFCs[threshold_key][dataset_name] = logFC_dict

    from itertools import combinations
    datasets = ["325", "693", "tcga", "glass"]
    dataset_pairs = list(combinations(datasets, 2))

    nrows = len(threshold_scenarios)
    ncols = len(dataset_pairs)
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(ncols * 6, nrows * 5), sharey=True)
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
        ax.set_title(f"{tkey}: {ds1} vs {ds2}", fontsize=12)
        ax.set_xlabel(f"log2 Fold Change ({ds1})")
        ax.set_ylabel(f"log2 Fold Change ({ds2})")
        ax.set_xlim(global_min, global_max)
        ax.set_ylim(global_min, global_max)
    # Remove unused axes
    for i in range(len(threshold_scenarios) * len(dataset_pairs), len(axes)):
        fig.delaxes(axes[i])

    plt.suptitle("Log2 Fold Change Comparisons Across Dataset Pairs", fontsize=16)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()

# List of CSV filenames and corresponding sheet names
# csv_files = [
#     ("/home/mkh062/Desktop/scratch/TCGA_project/processed_data/results/significant_genes_CGGA_24-48-months.csv", "CGGA_24-48"),
#     ("/home/mkh062/Desktop/scratch/TCGA_project/processed_data/results/significant_genes_CGGA_24-60-months.csv", "CGGA_24-60"),
#     ("/home/mkh062/Desktop/scratch/TCGA_project/processed_data/results/significant_genes_GLASS_24-48-months.csv", "GLASS_24-48"),
#     ("/home/mkh062/Desktop/scratch/TCGA_project/processed_data/results/significant_genes_GLASS_24-60-months.csv", "GLASS_24-60"),
# ]
#
# # Create an Excel writer object
# output_excel = "/home/mkh062/Desktop/scratch/TCGA_project/processed_data/results/significant_genes.xlsx"
# with pd.ExcelWriter(output_excel, engine="xlsxwriter") as writer:
#     for csv_file, sheet_name in csv_files:
#         try:
#             # Read CSV file
#             df = pd.read_csv(csv_file)
#
#             # Write to Excel with a separate sheet for each CSV
#             df.to_excel(writer, sheet_name=sheet_name, index=False)
#             print(f"Converted {csv_file} to sheet '{sheet_name}'")
#         except FileNotFoundError:
#             print(f"Error: {csv_file} not found. Skipping...")
#
# print(f"✅ Excel file '{output_excel}' created successfully!")