import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random
from lifelines import CoxPHFitter, KaplanMeierFitter
from sklearn.decomposition import PCA
from statsmodels.stats.multitest import multipletests
import rpy2.robjects as ro
from rpy2.robjects import r, pandas2ri

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
    plt.hist(df['p'], bins=30, edgecolor='black', density=True, label='Pvalues counts')
    plt.title(f'Normalized PDF of p_values (Covariate: {covariate})')
    plt.xlabel('P-values')
    plt.ylabel('Density')
    plt.legend()
    plt.savefig(f'{output_folder}/lfdr_{covariate}.png')
    lfdr_r_result_unique = lfdr_r(p_values_unique)
    lfdr_values_unique = list(lfdr_r_result_unique[0])
    df['lfdr'] = lfdr_values_unique
    print("df", df.head())
    return df

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

###############################################################################
# New Helper Function: Multivariate Cox Regression
###############################################################################
def fit_cox_models_multivariate(data, covariates_list):
    """
    For each gene (assumed to be in the gene expression columns starting from index 9),
    fit a Cox proportional hazards model using the gene expression value along with
    the covariates in covariates_list. Returns a DataFrame with the summary for the gene coefficient.
    """
    results_list = []
    gene_columns = data.columns[9:]
    for gene in gene_columns:
        cols = ['overall_survival', 'Status'] + covariates_list + [gene]
        cox_data = data[cols].dropna()
        new_names = ['overall_survival', 'Status'] + covariates_list + ['Gene']
        cox_data.columns = new_names
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
            print(f"Error processing gene {gene}: {e}")
    return pd.DataFrame(results_list)

###############################################################################
# Extended Analysis Function: get_stats_extended
###############################################################################
def get_stats_extended():
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

    # Plot histograms for 10 randomly selected genes
    gene_cols = [col for col in df.columns if col != 'NAME']
    # -------------------------------
    # Read clinical data
    xf = pd.read_excel(fclin, sheet_name="SuppTable1", skiprows=1)
    xf = xf[xf['Gender'] != 'MALE']
    xf.reset_index(inplace=True, drop=True)
    # For extended analysis we keep the original PAM50 mRNA values (do not map them)
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
    sub_xf = sub_xf.dropna(subset=['overall_survival'])
    sub_xf = sub_xf[sub_xf['overall_survival'] != '[Not Available]']
    sub_xf.reset_index(drop=True, inplace=True)
    sub_xf['overall_survival'] = sub_xf['overall_survival'].astype(int)

    sub_xf['Status'] = sub_xf['Vital Status'].apply(lambda x: 0 if x == 'LIVING' else 1)
    sub_xf.loc[sub_xf['Converted Stage'] == 'No_Conversion', 'Converted Stage'] = sub_xf['AJCC Stage']
    sub_xf.loc[sub_xf['Converted Stage'] == 'Stage ', 'Converted Stage'] = sub_xf['AJCC Stage']
    mapping_dict = {'Stage I': 0, 'Stage IA': 0, 'Stage IB': 0, 'Stage II': 1, 'Stage IIA': 1, 'Stage IIB': 1,
                    'Stage III': 2, 'Stage IIIA': 2, 'Stage IIIB': 2, 'Stage IIIC': 2, 'Stage IV': 3}
    sub_xf['Converted Stage'] = sub_xf['Converted Stage'].map(mapping_dict)
    sub_xf.rename(columns={'Converted Stage': 'Stage'}, inplace=True)
    sub_xf = sub_xf.drop(['AJCC Stage', 'Days to date of Death', 'Days to Date of Last Contact', 'Vital Status'], axis=1)

    print("Clinical data after processing:\n", sub_xf.head())

    # Merge gene expression data with clinical data
    merged_data = pd.merge(sub_xf, df, on='NAME', how='inner')
    merged_data = merged_data.drop_duplicates(subset='NAME')
    print("Merged data:\n", merged_data.head())

    ##########################################################################
    # Analysis (a): Cox model with covariates "Stage" and "Age at Initial Pathologic Diagnosis"
    ##########################################################################
    covariates_a = ['Stage', 'Age at Initial Pathologic Diagnosis']
    results_a = fit_cox_models_multivariate(merged_data, covariates_a)
    results_a = get_LF(results_a, "Stage_Age")
    results_a = results_a.sort_values(by='p_adjusted BH', ascending=True).reset_index(drop=True)
    top_100_a = results_a.head(100)
    print("Analysis (a) - Top 100 genes (all samples):")
    print(top_100_a.head())
    results_a.to_excel(f'{output_folder}/cox_age+stage_results.xlsx', index=False)

    ##########################################################################
    # Analysis (b): Cox model within PAM50 mRNA subtypes "Luminal A" and "Luminal B"
    ##########################################################################
    # Print unique PAM50 mRNA subtypes and count samples for Luminal A and Luminal B.
    unique_subtypes = merged_data['PAM50 mRNA'].unique()
    print("Unique PAM50 mRNA subtypes:", unique_subtypes)
    for subtype in ['Luminal A', 'Luminal B']:
        subtype_data = merged_data[merged_data['PAM50 mRNA'] == subtype].copy()
        print(f"Number of samples for {subtype}: {len(subtype_data)}")
        if subtype_data.empty:
            continue
        results_b = fit_cox_models_multivariate(subtype_data, covariates_a)
        results_b = get_LF(results_b, f"Stage_Age_{subtype}")
        results_b = results_b.sort_values(by='p_adjusted BH', ascending=True).reset_index(drop=True)
        top_100_b = results_b.head(100)
        print(f"Analysis (b) - Top 100 genes for {subtype}:")
        print(top_100_b.head())
        results_b.to_excel(f'{output_folder}/cox_age+stage_results_{subtype.replace(" ", "_")}.xlsx', index=False)

# To run the extended analysis, call:
get_stats_extended()
