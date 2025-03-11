import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random
from lifelines import CoxPHFitter, KaplanMeierFitter
from sklearn.decomposition import PCA
from statsmodels.stats.multitest import multipletests
import rpy2.robjects as ro
from rpy2.robjects import r
from rpy2.robjects import pandas2ri
from lifelines import KaplanMeierFitter


pandas2ri.activate()
r['source']('/home/mkh062/Desktop/scratch/PCL/lfdr_calc3.R')
lfdr_r = ro.globalenv['plotLFDR']
output_folder = '/home/mkh062/Desktop/scratch/PCL/lfdr_simulations_15Feb2025/brca_plots'

def get_LF(df, covariate):
    p_values = df['p'].values
    _, p_adjusted, _, _ = multipletests(p_values, method='fdr_bh')
    df['p_adjusted BH'] = p_adjusted
    df = df.sort_values(by='p', ascending=False)
    df.reset_index(drop=True, inplace=True)
    p_values_unique = ro.FloatVector(df['p'].tolist())
    plt.figure(figsize=(8, 6))
    plt.hist(df['p'], bins=30, edgecolor='black',
             density=True, label='Pvalues counts')  # Normalize the histogram
    # plt.axhline(y=0.8267012, color='r', linestyle='--',
    #             label='Estimate of pi0 (proportion of true null hypotheses): 0.8267012 ')
    plt.title(f'Normalized PDF of p_values (Covariate: {covariate})')
    plt.xlabel('P-values')
    plt.ylabel('Density')
    plt.legend()
    plt.savefig('/home/mkh062/Desktop/scratch/PCL/lfdr_simulations_15Feb2025/brca_plots/lfdr' + str(covariate) + '.png')
    # plt.show()
    lfdr_r_result_unique = lfdr_r(p_values_unique)
    lfdr_values_unique = list(lfdr_r_result_unique[0])
    cumsum_unique = np.cumsum(lfdr_values_unique)
    cumavg_lfdr_unique = cumsum_unique / np.arange(1, len(lfdr_values_unique) + 1)
    df['lfdr'] = lfdr_values_unique
    # df['cumavg_lfdr'] = cumavg_lfdr_unique
    print("df", df)
    # df = df.sort_values(by='cumavg_lfdr', ascending=False)
    return df
########################################
# Helper functions

def plot_kaplan_meier(data, gene, covariate_name, ax):
    """
    Plot Kaplan–Meier curves for a gene by dichotomizing expression at the median.
    """
    kmf_high = KaplanMeierFitter()
    kmf_low = KaplanMeierFitter()
    median_val = data[gene].median()
    high_expr = data[data[gene] >= median_val]
    low_expr = data[data[gene] < median_val]
    kmf_high.fit(high_expr['overall_survival'], event_observed=high_expr['Status'], label=f'High {gene}')
    kmf_low.fit(low_expr['overall_survival'], event_observed=low_expr['Status'], label=f'Low {gene}')
    kmf_high.plot(ax=ax)
    kmf_low.plot(ax=ax)
    ax.set_title(f"KM: {gene} ({covariate_name})")

def bootstrap_null_model(merged_data, covariate, gene_columns, n_iterations=100, significance_threshold=0.05):
    """
    For a given covariate, randomize gene expression values (thus the gene labels)
    and recompute the Cox regression for each gene over multiple iterations.
    Return a list with the number of significant genes (based on adjusted p-value)
    for each iteration.
    """
    significant_counts = []
    for i in range(n_iterations):
        shuffled_data = merged_data.copy()
        # Shuffle each gene column separately
        for col in gene_columns:
            shuffled_data[col] = np.random.permutation(shuffled_data[col].values)
        results_list = []
        for gene in gene_columns:
            cox_data = shuffled_data[['overall_survival', 'Status', covariate, gene]].dropna()
            cox_data.columns = ['overall_survival', 'Status', 'Covariate', 'Gene']
            try:
                cph = CoxPHFitter()
                cph.fit(cox_data, duration_col='overall_survival', event_col='Status')
                p_val = cph.summary.loc['Gene', 'p']
            except Exception as e:
                p_val = np.nan
            results_list.append(p_val)
        results_df = pd.DataFrame({'Gene': gene_columns, 'p': results_list}).dropna()
        if not results_df.empty:
            _, p_adjusted, _, _ = multipletests(results_df['p'].values, method='fdr_bh')
            significant_counts.append((p_adjusted < significance_threshold).sum())
        else:
            significant_counts.append(0)
    return significant_counts

########################################
# Main analysis function

def get_stats():
    # File paths
    fexpr = "/home/mkh062/Desktop/scratch/PCL/BRCA.exp.547.med.txt"
    fclin = "/home/mkh062/Desktop/scratch/PCL/Supplementary Tables 1-4.xlsx"

    # -------------------------------
    # Read gene expression data
    df = pd.read_csv(fexpr, sep='\t', header=0, index_col=0)
    df = df.dropna()  # Drop genes with missing values
    df = df.transpose()
    df.reset_index(inplace=True, drop=False)
    df.rename(columns={'index': 'NAME'}, inplace=True)
    # Standardize the NAME to first 12 characters
    df['NAME'] = df['NAME'].apply(lambda x: x[:12])
    print("Gene expression data (original):\n", df.head())

    # 1. Plot histograms for 10 randomly selected genes (2x5 grid)
    gene_cols = [col for col in df.columns if col != 'NAME']
    if len(gene_cols) >= 10:
        random_genes = random.sample(gene_cols, 10)
        fig, axes = plt.subplots(2, 5, figsize=(20, 8))
        axes = axes.flatten()
        for i, gene in enumerate(random_genes):
            axes[i].hist(df[gene].dropna(), bins=30, edgecolor='black')
            axes[i].set_title(gene)
        plt.tight_layout()
        plt.savefig('/home/mkh062/Desktop/scratch/PCL/lfdr_simulations_15Feb2025/brca_plots/gene_expression_histograms.png')
        plt.close()
    else:
        print("Not enough genes for histogram plotting.")

    # -------------------------------
    # Perform PCA on gene expression data to visualize clustering
    expression_data = df[gene_cols]
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(expression_data)
    plt.figure(figsize=(8, 6))
    plt.scatter(pca_result[:, 0], pca_result[:, 1], alpha=0.5)
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.title('PCA of Gene Expression Data')
    plt.savefig('/home/mkh062/Desktop/scratch/PCL/lfdr_simulations_15Feb2025/brca_plots/pca_gene_expression.png')
    plt.close()

    # -------------------------------
    # Read clinical data
    xf = pd.read_excel(fclin, sheet_name="SuppTable1", skiprows=1)
    xf = xf[xf['Gender'] != 'MALE']
    xf.reset_index(inplace=True, drop=True)

    sub_xf = xf[['Complete TCGA ID', 'Vital Status', 'Days to Date of Last Contact',
                 'Days to date of Death', 'Age at Initial Pathologic Diagnosis', 'AJCC Stage', 'PAM50 mRNA',
                 'ER Status', 'HER2 Final Status', 'PR Status', 'Converted Stage']]
    sub_xf = sub_xf.rename(columns={'Complete TCGA ID': 'NAME'})
    print("Clinical data (filtered):\n", sub_xf.head())

    # 2. Compute overall_survival
    def compute_overall_survival(row):
        death = row['Days to date of Death']
        last_contact = row['Days to Date of Last Contact']
        if row['Vital Status'] == 'DECEASED':
            if pd.notna(death) and pd.notna(last_contact) and (death != last_contact):
                print("Discrepancy in sample:", row['NAME'])
            return death
        else:
            return last_contact

    sub_xf['overall_survival'] = sub_xf.apply(compute_overall_survival, axis=1)
    # Drop missing or not available values in overall_survival
    sub_xf = sub_xf.dropna(subset=['overall_survival'])
    sub_xf = sub_xf[sub_xf['overall_survival'] != '[Not Available]']
    sub_xf.reset_index(drop=True, inplace=True)
    sub_xf['overall_survival'] = sub_xf['overall_survival'].astype(int)

    # Create Status column: 0 if LIVING, 1 if DECEASED
    sub_xf['Status'] = sub_xf['Vital Status'].apply(lambda x: 0 if x == 'LIVING' else 1)

    # Replace 'No_Conversion' or ambiguous stage with values from AJCC Stage
    sub_xf.loc[sub_xf['Converted Stage'] == 'No_Conversion', 'Converted Stage'] = sub_xf['AJCC Stage']
    sub_xf.loc[sub_xf['Converted Stage'] == 'Stage ', 'Converted Stage'] = sub_xf['AJCC Stage']
    mapping_dict = {'Stage I': 0, 'Stage IA': 0, 'Stage IB': 0, 'Stage II': 1, 'Stage IIA': 1, 'Stage IIB': 1,
                    'Stage III': 2, 'Stage IIIA': 2, 'Stage IIIB': 2, 'Stage IIIC': 2, 'Stage IV': 3}
    sub_xf['Converted Stage'] = sub_xf['Converted Stage'].map(mapping_dict)
    sub_xf.rename(columns={'Converted Stage': 'Stage'}, inplace=True)
    sub_xf = sub_xf.drop(['AJCC Stage', 'Days to date of Death', 'Days to Date of Last Contact', 'Vital Status'], axis=1)

    # -------------------------------
    # Define covariates for analysis
    covariates = {
        'Age at Initial Pathologic Diagnosis': 'Age',
        'Stage': 'Stage',
        'PAM50 mRNA': 'PAM50',
        'ER Status': 'ER',
        'PR Status': 'PR',
        'HER2 Final Status': 'HER2'
    }

    print("Clinical data after processing:\n", sub_xf.head())

    # Merge gene expression data with clinical data
    merged_data = pd.merge(sub_xf, df, on='NAME', how='inner')
    merged_data = merged_data.drop_duplicates(subset='NAME')
    print("Merged data:\n", merged_data.head())

    # -------------------------------
    # Function to fit Cox models for each gene using a specific covariate
    def fit_cox_models(data, covariate1):
        results_list = []
        # Assuming gene expression columns start from column index 9 (adjust as needed)
        gene_columns = data.columns[9:]
        for gene in gene_columns:
            cox_data = data[['overall_survival', 'Status', covariate1, gene]].dropna()
            cox_data.columns = ['overall_survival', 'Status', 'Covariate', 'Gene']
            try:
                cph = CoxPHFitter()
                cph.fit(cox_data, duration_col='overall_survival', event_col='Status')
                summary = cph.summary.loc['Gene']
                results_list.append({
                    'Gene': gene,
                    'coef': summary['coef'],
                    'exp(coef)': summary['exp(coef)'],
                    'se(coef)': summary['se(coef)'],
                    'z': summary['z'],
                    'p': summary['p'],
                    'lower 95%': summary['coef lower 95%'],
                    'upper 95%': summary['coef upper 95%']
                })
            except Exception as e:
                print(f"Error processing gene {gene} with covariate {covariate1}: {e}")
        return pd.DataFrame(results_list)

    all_results = []
    top_genes_dict = {}  # To store top 100 significant genes for each covariate

    # -------------------------------
    # Loop over each covariate and run Cox regression analyses
    for covariate, covariate_name in covariates.items():
        print(f"\nPerforming Cox regression using {covariate_name} as covariate...")
        data_cov = merged_data.copy()
        if covariate_name == 'Stage':
            data_cov = data_cov.dropna(subset=['Stage'])
            data_cov = data_cov[(data_cov['Stage'] != '[Not Available]') & (data_cov['Stage'] != 'Stage X')]
            data_cov.reset_index(drop=True, inplace=True)
        elif covariate_name == 'PAM50':
            data_cov = data_cov.dropna(subset=['PAM50 mRNA'])
            data_cov.reset_index(drop=True, inplace=True)
            pam50_mapping = {val: idx for idx, val in enumerate(data_cov['PAM50 mRNA'].unique())}
            data_cov['PAM50 mRNA'] = data_cov['PAM50 mRNA'].map(pam50_mapping)
            print("PAM50 mapping:", data_cov['PAM50 mRNA'].unique())
        elif covariate_name == 'ER':
            data_cov = data_cov.dropna(subset=['ER Status'])
            data_cov = data_cov[~data_cov['ER Status'].isin(['Performed but Not Available', 'Not Available'])]
            data_cov.reset_index(drop=True, inplace=True)
            status_mapping = {'Negative': -1, 'Positive': 1}
            data_cov['ER Status'] = data_cov['ER Status'].map(status_mapping)
        elif covariate_name == 'PR':
            data_cov = data_cov.dropna(subset=['PR Status'])
            data_cov = data_cov[~data_cov['PR Status'].isin(['Performed but Not Available', 'Not Available', 'Not Performed'])]
            data_cov.reset_index(drop=True, inplace=True)
            status_mapping = {'Negative': -1, 'Positive': 1}
            data_cov['PR Status'] = data_cov['PR Status'].map(status_mapping)
        elif covariate_name == 'HER2':
            data_cov = data_cov.dropna(subset=['HER2 Final Status'])
            data_cov = data_cov[~data_cov['HER2 Final Status'].isin(['Not Available', 'Equivocal'])]
            data_cov.reset_index(drop=True, inplace=True)
            status_mapping = {'Negative': -1, 'Positive': 1}
            data_cov['HER2 Final Status'] = data_cov['HER2 Final Status'].map(status_mapping)
        else:
            data_cov = data_cov.dropna(subset=['Age at Initial Pathologic Diagnosis'])
            data_cov.reset_index(drop=True, inplace=True)

        print("Data for covariate", covariate_name, ":\n", data_cov.head())
        results = fit_cox_models(data_cov, covariate)
        results = get_LF(results, covariate_name)
        results = results.sort_values(by=f'p_adjusted_BH_{covariate_name}', ascending=True)
        print(f"Cox model results for {covariate_name}:\n", results.head())
        results.to_csv(f'/home/mkh062/Desktop/scratch/PCL/lfdr_simulations_15Feb2025/brca_plots/all_genes_results_{covariate}.csv', index=False)

        # Save top 100 significant genes (by adjusted p-value)
        top_100 = results.head(100)
        top_genes_dict[covariate_name] = set(top_100['Gene'])
        max_adj_p = top_100[f'p_adjusted_BH_{covariate_name}'].max()
        print(f"Top 100 genes for {covariate_name}: max adjusted p-value = {max_adj_p}")
        all_results.append(results[['Gene', f'p_{covariate_name}', f'p_adjusted_BH_{covariate_name}', f'lfdr_{covariate_name}']])

        # 5. Plot Kaplan–Meier curves for top 4 significant genes
        top_genes = results.head(4)['Gene']
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        axes = axes.flatten()
        for i, gene in enumerate(top_genes):
            plot_kaplan_meier(data_cov, gene, covariate_name, axes[i])
        plt.tight_layout()
        plt.savefig(f'/home/mkh062/Desktop/scratch/PCL/lfdr_simulations_15Feb2025/brca_plots/km_curves_top4_{covariate_name}.png')
        plt.close()

    # -------------------------------
    # 4. Compare significant genes across covariates
    common_all = set.intersection(*top_genes_dict.values())
    print("Number of genes commonly significant across all covariates:", len(common_all))
    cov_names = list(top_genes_dict.keys())
    for i in range(len(cov_names)):
        for j in range(i+1, len(cov_names)):
            common_pair = top_genes_dict[cov_names[i]].intersection(top_genes_dict[cov_names[j]])
            print(f"Common genes between {cov_names[i]} and {cov_names[j]}: {len(common_pair)}")

    # Merge results from all covariates
    merged_results = all_results[0]
    for res in all_results[1:]:
        merged_results = merged_results.merge(res, on='Gene', how='outer')
    # Reorder columns
    column_order = ['Gene']
    for covariate_name in covariates.values():
        column_order.extend([f'p_{covariate_name}', f'p_adjusted_BH_{covariate_name}', f'lfdr_{covariate_name}'])
    merged_results = merged_results[column_order]
    merged_results.to_csv('/home/mkh062/Desktop/scratch/PCL/lfdr_simulations_15Feb2025/brca_plots/merged_results_all_genes.csv', index=False)
    merged_results.to_excel('/home/mkh062/Desktop/scratch/PCL/lfdr_simulations_15Feb2025/brca_plots/merged_result_all_genes.xlsx', index=False)
    print("Merged results:\n", merged_results.head())

    # -------------------------------
    # 6. Bootstrap null model to assess significance of random gene signatures
    # Using all gene columns (assumed to start at index 9)
    all_gene_columns = merged_data.columns[9:]
    for covariate, covariate_name in covariates.items():
        print(f"\nRunning bootstrap null model for covariate {covariate_name}...")
        significant_counts = bootstrap_null_model(merged_data, covariate, all_gene_columns, n_iterations=100)
        # (100 iterations used for speed; adjust as needed)
        print(f"Bootstrap null model for {covariate_name}:")
        print("Mean number of significant genes:", np.mean(significant_counts))
        print("Std of significant genes:", np.std(significant_counts))

# To run the analysis, simply call:
get_stats()
