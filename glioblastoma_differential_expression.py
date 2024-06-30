import pandas as pd
import statsmodels.api as sm
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from scipy.stats import levene
from collections import Counter
from lifelines import CoxPHFitter
from lifelines import KaplanMeierFitter
from lifelines.exceptions import ConvergenceError
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri


def cox_reg(data):
    # columns_to_exclude = {'Sample.ID', 'Overall.Survival.Months', 'Diagnosis.Age', 'Sex', 'case.vital.status'}
    columns_to_exclude = {'Sample.ID', 'Overall.Survival.Months', 'Diagnosis.Age', 'Sex', 'case.vital.status',
                          'platform'}
    gene_data = data[[col for col in data.columns if col not in columns_to_exclude]]
    std_dev_genes = gene_data.std()
    sd_limit = np.percentile(std_dev_genes, 25)

    # Filtering genes based on standard deviation
    significant_genes_std = std_dev_genes[abs(std_dev_genes) >= sd_limit].index
    final_genes_list = list(significant_genes_std)

    merged_data = data[list(columns_to_exclude) + final_genes_list]
    merged_data = merged_data.drop_duplicates('Sample.ID')

    alive_count = merged_data['case.vital.status'].value_counts().get('alive', 0)
    dead_count = merged_data['case.vital.status'].value_counts().get('dead', 0)
    print(f'Alive: {alive_count}, Dead: {dead_count}')

    merged_data['Status'] = merged_data['case.vital.status'].apply(lambda x: 0 if x == 'alive' else 1)
    cols = list(columns_to_exclude) + ['Status'] + final_genes_list
    merged_data = merged_data[cols]
    merged_data = merged_data.drop(columns=['case.vital.status'])

    # merged_data['Sex'] = merged_data['Sex'].map({'male': 0, 'female': 1})
    merged_data['Sex'] = merged_data['Sex'].map({'Male': 0, 'Female': 1})
    unique_platforms = merged_data['platform'].unique()
    print("Unique values in 'platform' column:", unique_platforms)
    # Count samples in each group of 'platform'
    
    platform_counts = merged_data['platform'].value_counts()
    print("\nNumber of samples in each platform group:\n", platform_counts)

    merged_data = merged_data[merged_data['platform'] != 'x']
    merged_data.reset_index(drop=True, inplace=True)
    merged_data['platform'] = merged_data['platform'].map({'agilent': 0, 'RNAseq': 1})
    
    cph = CoxPHFitter()
    # Create a list to store the statistics dictionaries
    statistics_list = []
    """
    df (DataFrame) – a Pandas DataFrame with necessary columns duration_col and event_col (see below), 
    covariates columns, and special columns (weights, strata). duration_col refers to the lifetimes 
    of the subjects. event_col refers to whether the ‘death’ events was observed:
     1 if observed, 0 else (censored).
    """
    # Perform CoxPH regression for each gene
    for gene in merged_data.columns[6:]:  # Skip sample_id, overall_survival, Status
        print("gene", gene)
        # gene_data = merged_data[['Overall.Survival.Months', 'Status', 'Diagnosis.Age', 'Sex', gene]].dropna()
        gene_data = merged_data[['Overall.Survival.Months', 'Status', 'Diagnosis.Age', 'Sex', 'platform', gene]].dropna()
        try:
            cph.fit(gene_data, duration_col='Overall.Survival.Months', event_col='Status')
        except ConvergenceError as e:
            print(f"Convergence error for {gene}: {e}")
            continue

        # Extract statistics
        coef = cph.summary.loc[gene, 'coef']
        exp_coef = cph.summary.loc[gene, 'exp(coef)']
        se_coef = cph.summary.loc[gene, 'se(coef)']
        z = cph.summary.loc[gene, 'z']
        p = cph.summary.loc[gene, 'p']
        lower_95 = cph.summary.loc[gene, 'coef lower 95%']
        upper_95 = cph.summary.loc[gene, 'coef upper 95%']

        # Append statistics to the DataFrame
        # Append statistics to the list
        statistics_list.append({
            'gene': gene,
            'coef': coef,
            'exp(coef)': exp_coef,
            'se(coef)': se_coef,
            'z': z,
            'p': p,
            'coef lower 95': lower_95,
            'coef upper 95': upper_95
        })

    # Convert the list of statistics dictionaries to a DataFrame
    statistics = pd.DataFrame(statistics_list)
    # Adjust p-values for multiple testing using FDR method
    p_values = statistics['p'].values
    _, p_adjusted, _, _ = multipletests(p_values, method='fdr_bh')
    statistics['p_adjusted'] = p_adjusted

    # Filter significant genes (e.g., adjusted p < 0.05)
    significant_genes = statistics[statistics['p_adjusted'] < 0.1]
    print(significant_genes)

    # Save all statistics to a CSV file
    statistics.to_csv('/home/mkh062/Desktop/scratch/TCGA_project/survival_analysis/TCGA_gene_survival_statistics.csv',
                      index=False)

    # Save significant genes statistics to a CSV file
    significant_genes.to_csv(
        '/home/mkh062/Desktop/scratch/TCGA_project/survival_analysis/TCGA_significant_survival_statistics.csv',
        index=False)

    print("All statistics saved to 'gene_survival_statistics.csv'")
    print("Significant gene statistics saved to 'significant_gene_survival_statistics.csv'")
    # Plot the adjusted p-values
    plt.figure(figsize=(12, 6))
    sns.barplot(x='gene', y=-np.log10(statistics['p_adjusted']), data=statistics)
    plt.xticks(rotation=90)
    plt.title('Adjusted p-values for Differential Gene Expression Analysis TCGA data')
    plt.xlabel('Gene')
    plt.ylabel('-log10(adjusted p-value)')
    plt.show()


def check_variance_assumption(group1, group2):
    variance_results = {}
    for gene in group1.columns:
        stat, p_val = levene(group1[gene], group2[gene])
        variance_results[gene] = p_val

    # Convert results to DataFrame
    variance_results_df = pd.DataFrame(list(variance_results.items()), columns=['Gene', 'p-value'])

    # Determine if variances are equal
    variance_results_df['Equal Variance'] = variance_results_df['p-value'] > 0.05

    return variance_results_df


def t_test(D, outcome_d, num_iterations=500, threshold=0.3):
    D = D[D['platform'] == 'agilent']
    D = D.drop(columns=['platform'])
    D.reset_index(drop=True, inplace=True)
    # outcome_d = outcome_d[['Sample ID', 'Overall Survival (Months)']]
    outcome_d = outcome_d[
        ['Sample ID', 'Overall Survival (Months)', 'Overall Survival Status']]
    outcome_d = outcome_d[
        ~((outcome_d['Overall Survival (Months)'] <= 24.0) & (outcome_d['Overall Survival Status'] == '0:LIVING'))
    ]
    outcome_d.drop(columns=['Overall Survival Status'], inplace=True)
    data = pd.merge(outcome_d, D, on='Sample ID')
    median_survival = 24.0

    # Split the data into two groups based on median survival
    group1 = data[data['Overall Survival (Months)'] > median_survival]
    group2 = data[data['Overall Survival (Months)'] <= median_survival]

    # Reset indices of the groups
    group1.reset_index(drop=True, inplace=True)
    group2.reset_index(drop=True, inplace=True)

    # Print the number of samples in each group
    print("Number of samples in Group 1:", len(group1))
    print("Number of samples in Group 2:", len(group2))

    # Extract gene expression data (removing other columns)
    gene_expression_data = data.drop(columns=['Sample ID', 'Overall Survival (Months)'])
    print("gene_expression_data", gene_expression_data.head())

    # Dictionary to count significant occurrences of each gene
    significant_genes_count = Counter()

    for i in range(num_iterations):
        print(f"Running iteration {i + 1}/{num_iterations}")
        results = []
        for gene in gene_expression_data.columns:
            expr_group1 = group1[gene]
            expr_group2 = group2[gene].sample(n=len(group1), random_state=np.random.randint(0, 10000))
            # Check variance assumption
            stat, p_val_levene = levene(expr_group1, expr_group2)
            equal_var = p_val_levene > 0.05
            t_stat, p_val = ttest_ind(expr_group1, expr_group2, equal_var=equal_var)

            results.append((gene, t_stat, p_val))

        # Convert results to DataFrame
        results_df = pd.DataFrame(results, columns=['Gene', 't-statistic', 'p-value'])

        # Adjust p-values for multiple testing
        results_df['p-adjusted'] = multipletests(results_df['p-value'], method='fdr_bh')[1]

        # Find significant genes for this iteration
        iteration_significant_genes = results_df[results_df['p-adjusted'] < 0.2]['Gene']
        significant_genes_count.update(iteration_significant_genes)

    # Determine genes that are consistently significant
    consistent_significant_genes = [gene for gene, count in significant_genes_count.items() if
                                    count >= num_iterations * threshold]

    # Prepare final DataFrame
    final_results_df = pd.DataFrame(consistent_significant_genes, columns=['Gene'])
    final_results_df['Occurrences'] = final_results_df['Gene'].map(significant_genes_count)

    # Save the results to a CSV file
    final_results_df.to_csv('/home/mkh062/Desktop/scratch/TCGA_project/consistent_significant_genes_new.csv',
                            index=False)

    print(final_results_df)
    return final_results_df


def gene_plot(merged_data):
    # Select the first 10 genes
    selected_genes = merged_data.columns[4:14]
    print("selected_genes", selected_genes)
    selected_gene_expression_data = merged_data[selected_genes]
    # Plot the distribution of gene expression levels for the first 10 genes
    fig, axes = plt.subplots(2, 5, figsize=(20, 8))
    axes = axes.flatten()
    for i, gene in enumerate(selected_genes):
        sns.histplot(selected_gene_expression_data[gene], bins=30, kde=True, ax=axes[i])
        axes[i].set_title(gene)
        axes[i].set_xlabel('Expression Level')
        axes[i].set_ylabel('Frequency')

    plt.tight_layout()
    plt.savefig('/home/mkh062/Desktop/scratch/TCGA_project/gene_exp_dist_glass.png', dpi=300)
    plt.show()


def PCA_plot(gene_expression_data, platform):
    # Standardize the data
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(gene_expression_data)  # Transpose to have samples as rows

    # Perform PCA
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(scaled_data)

    # Plot PCA results
    plt.figure(figsize=(12, 6))
    platform_colors = {'agilent': 'red', 'RNAseq': 'blue'}

    for plat in platform.unique():
        indices = platform[platform == plat].index
        plt.scatter(pca_result[indices, 0], pca_result[indices, 1], label=plat, alpha=0.7)

    plt.title('PCA of Gene Expression Data')
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.legend(title='Platform')
    plt.savefig('/home/mkh062/Desktop/scratch/TCGA_project/PCA.png')
    plt.show()


def transform_sample_name(name):
    # Extract part up to and including 'TP'
    if 'TP' in name:
        name = name.split('TP')[0] + 'TP'
    # Replace dots with dashes
    name = name.replace('.', '-')
    return name


def limma(data, threshold):
    data = data[
        ~((data['Overall.Survival.Months'] <= threshold) & (data['case.vital.status'] == 'alive'))
    ]
    data.drop(columns=['case.vital.status'], inplace=True)
    data.reset_index(drop=True, inplace=True)
    print(f"Number of samples before filtering: ", len(data))
    less_than_threshold = data[data['Overall.Survival.Months'] <= threshold]
    more_than_threshold = data[data['Overall.Survival.Months'] > threshold]

    count_less_than_threshold = len(less_than_threshold)
    count_more_than_threshold = len(more_than_threshold)

    print(f"Number of samples with survival <= {threshold} months: {count_less_than_threshold}")
    print(f"Number of samples with survival > {threshold} months: {count_more_than_threshold}")

    # gene_data = data.drop(columns=['Sample.ID', 'Overall.Survival.Months', 'Diagnosis.Age', 'Sex'])
    gene_data = data.drop(columns=['Sample.ID', 'Overall.Survival.Months', 'Diagnosis.Age', 'Sex', 'platform'])
    non_zero_variance_genes = gene_data.loc[:, gene_data.var(axis=0) != 0]
    # data_filtered = pd.concat(
    #     [data[['Sample.ID', 'Overall.Survival.Months', 'Diagnosis.Age', 'Sex']], non_zero_variance_genes], axis=1)
    data_filtered = pd.concat(
        [data[['Sample.ID', 'Overall.Survival.Months', 'Diagnosis.Age', 'Sex', 'platform']], non_zero_variance_genes],
        axis=1)
    data = data_filtered
    print("data", data.head())
    pandas2ri.activate()
    r_data = pandas2ri.py2rpy(data)

    r_script = f"""
      if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')
      BiocManager::install('limma')
      library(limma)
      # gene_expression_data <- data[, !names(data) %in% c('Sample.ID', 'Overall.Survival.Months', 'Diagnosis.Age', 'Sex')]
      gene_expression_data <- data[, !names(data) %in% c('Sample.ID', 'Overall.Survival.Months', 'Diagnosis.Age', 'Sex', 'platform')]
      rownames(gene_expression_data) <- data$Sample.ID
      median_survival <- {threshold}
      data$group <- ifelse(data$Overall.Survival.Months > median_survival, 'high', 'low')
      data$group <- factor(data$group, levels = c('low', 'high'))
      # design <- model.matrix(~ group + Diagnosis.Age + Sex , data = data)
      design <- model.matrix(~ group + Diagnosis.Age + Sex + platform, data = data)
      rownames(design) <- data$Sample.ID
      gene_expression_matrix <- t(as.matrix(gene_expression_data))
      fit <- lmFit(gene_expression_matrix, design)
      fit <- eBayes(fit)
      results <- topTable(fit, coef = 'grouphigh', number = Inf, adjust.method = 'fdr')
      write.csv(results, file = '/home/mkh062/Desktop/scratch/TCGA_project/limma_analysis/results_limma_TCGA_{threshold}months.csv', row.names = TRUE)
      significant_genes <- results[results$adj.P.Val < 0.1, ]
      write.csv(significant_genes, file = '/home/mkh062/Desktop/scratch/TCGA_project/limma_analysis/significant_genes_limma_TCGA_{threshold}months.csv', row.names = TRUE)
    """

    # Execute the R script with the data
    ro.globalenv['data'] = r_data
    ro.r(r_script)


def TCGA_data():
    out = '/home/mkh062/Desktop/scratch/TCGA_project/Compiled data_TCGA_CGGA_GLASS_jh080424.xlsx'  # Update this path
    outcome_data = pd.read_excel(out, sheet_name='TCGA')
    gene_expression_file = '/home/mkh062/Desktop/scratch/TCGA_project/GBMLGG_EB_RmDiffFullGenesRanRmDup.csv'  # Path to the gene expression data CSV file
    D2 = pd.read_csv(gene_expression_file, sep=',')

    common_ids = pd.Series(list(set(D2['samples']).intersection(set(outcome_data['Sample ID'])))).tolist()
    D2_subset = D2[D2['samples'].isin(common_ids)]
    D2_subset.reset_index(drop=True, inplace=True)
    D2_subset.rename(columns={'samples': 'Sample ID'}, inplace=True)

    outcome_data = outcome_data[
        ['Sample ID', 'Overall Survival (Months)', 'Diagnosis Age', 'Sex', 'Overall Survival Status']]
    data = pd.merge(outcome_data, D2_subset, on='Sample ID')
    data.reset_index(drop=True, inplace=True)
    data['Overall Survival Status'] = data['Overall Survival Status'].replace({
        '1:DECEASED': 'dead',
        '0:LIVING': 'alive'
    })
    data.rename(columns={
        'Sample ID': 'Sample.ID',
        'Overall Survival (Months)': 'Overall.Survival.Months',
        'Diagnosis Age': 'Diagnosis.Age',
        'Overall Survival Status': 'case.vital.status'
    }, inplace=True)

    # Check for duplicated Sample IDs
    duplicated_ids = data['Sample.ID'][data['Sample.ID'].duplicated(keep=False)]

    # Print duplicated Sample IDs
    if not duplicated_ids.empty:
        print("Duplicated Sample IDs:")
        for id in duplicated_ids.unique():
            print(id)
    else:
        print("No duplicated Sample IDs found.")
    # Drop all duplicated Sample IDs
    data = data[~data['Sample.ID'].isin(duplicated_ids)]
    data.reset_index(drop=True, inplace=True)
    cox_reg(data)
    thresholds = [24, 36, 48, 60]
    for threshold in thresholds:
        limma(data, threshold)


def Glass_data():
    out = '/home/mkh062/Desktop/scratch/TCGA_project/Compiled data_TCGA_CGGA_GLASS_jh080424.xlsx'  # Update this path
    outcome_data = pd.read_excel(out, sheet_name='GLASS')

    gene_expression_file = '/home/mkh062/Desktop/scratch/TCGA_project/gene_tpm_matrix_all_samples.tsv'  # Path to the gene expression data CSV file
    D2 = pd.read_csv(gene_expression_file, sep='\t')

    D2_transposed = D2.set_index('Gene_symbol').transpose().reset_index()
    new_columns = ['Sample ID'] + D2['Gene_symbol'].tolist()
    D2_transposed.columns = new_columns

    D2_transposed['Sample ID'] = D2_transposed['Sample ID'].apply(transform_sample_name)
    common_ids = pd.Series(list(set(D2_transposed['Sample ID']).intersection(set(outcome_data['Sample ID'])))).tolist()

    D2_subset = D2_transposed[D2_transposed['Sample ID'].isin(common_ids)]
    D2_subset.reset_index(drop=True, inplace=True)

    outcome_data = outcome_data[
        ['Sample ID', 'case_overall_survival_mo', 'case_age_diagnosis_years', 'case_sex', 'case_vital_status']]
    data = pd.merge(outcome_data, D2_subset, on='Sample ID')
    data.reset_index(drop=True, inplace=True)

    data.rename(columns={
        'Sample ID': 'Sample.ID',
        'case_overall_survival_mo': 'Overall.Survival.Months',
        'case_age_diagnosis_years': 'Diagnosis.Age',
        'case_sex': 'Sex',
        'case_vital_status': 'case.vital.status'
    }, inplace=True)

    # Check for duplicated Sample IDs
    duplicated_ids = data['Sample.ID'][data['Sample.ID'].duplicated(keep=False)]

    # Print duplicated Sample IDs
    if not duplicated_ids.empty:
        print("Duplicated Sample IDs:")
        for id in duplicated_ids.unique():
            print(id)
    else:
        print("No duplicated Sample IDs found.")
    # Drop all duplicated Sample IDs
    data = data[~data['Sample.ID'].isin(duplicated_ids)]
    data.reset_index(drop=True, inplace=True)
    cox_reg(data)
    thresholds = [24, 36, 48, 60]
    for threshold in thresholds:
        limma(data, threshold)

    # random_test(D2_subset, outcome_data)
    # t_test2(D2_subset, outcome_data)
    # gene_exp_plot(D2_subset, outcome_data)
    # platform = D2_subset['platform']
    # D2_subset = D2_subset.drop(columns=['platform'])
    # D2_PCA = D2_subset.drop(columns=['samples'])
    # PCA_plot(D2_PCA, platform)
    # Find intersection of all significant genes sets, handle empty intersection


def plot_kaplan_meier_glass():
    # threshold = 60
    # data1 = pd.read_csv(
    #     f'/home/mkh062/Desktop/scratch/TCGA_project/limma_analysis/significant_genes_limma_GLASS_{threshold}months.csv')
    # print("data1", data1.head())
    data1 = pd.read_csv(
        f'/home/mkh062/Desktop/scratch/TCGA_project/survival_analysis/GLASS_significant_survival_statistics.csv')
    genes = data1['gene'].head(4).tolist()

    out = '/home/mkh062/Desktop/scratch/TCGA_project/Compiled data_TCGA_CGGA_GLASS_jh080424.xlsx'  # Update this path
    outcome_data = pd.read_excel(out, sheet_name='GLASS')

    gene_expression_file = '/home/mkh062/Desktop/scratch/TCGA_project/gene_tpm_matrix_all_samples.tsv'  # Path to the gene expression data CSV file
    D2 = pd.read_csv(gene_expression_file, sep='\t')

    D2_transposed = D2.set_index('Gene_symbol').transpose().reset_index()
    new_columns = ['Sample ID'] + D2['Gene_symbol'].tolist()
    D2_transposed.columns = new_columns

    D2_transposed['Sample ID'] = D2_transposed['Sample ID'].apply(transform_sample_name)
    common_ids = pd.Series(list(set(D2_transposed['Sample ID']).intersection(set(outcome_data['Sample ID'])))).tolist()

    D2_subset = D2_transposed[D2_transposed['Sample ID'].isin(common_ids)]
    D2_subset.reset_index(drop=True, inplace=True)

    outcome_data = outcome_data[
        ['Sample ID', 'case_overall_survival_mo', 'case_age_diagnosis_years', 'case_sex', 'case_vital_status']]
    data = pd.merge(outcome_data, D2_subset, on='Sample ID')
    data.reset_index(drop=True, inplace=True)

    data.rename(columns={
        'Sample ID': 'Sample.ID',
        'case_overall_survival_mo': 'Overall.Survival.Months',
        'case_age_diagnosis_years': 'Diagnosis.Age',
        'case_sex': 'Sex',
        'case_vital_status': 'case.vital.status'
    }, inplace=True)

    # Check for duplicated Sample IDs
    duplicated_ids = data['Sample.ID'][data['Sample.ID'].duplicated(keep=False)]

    # Print duplicated Sample IDs
    if not duplicated_ids.empty:
        print("Duplicated Sample IDs:")
        for id in duplicated_ids.unique():
            print(id)
    else:
        print("No duplicated Sample IDs found.")
    # Drop all duplicated Sample IDs
    data = data[~data['Sample.ID'].isin(duplicated_ids)]
    data.reset_index(drop=True, inplace=True)
    # data = data[
    #     ~((data['Overall.Survival.Months'] <= threshold) & (data['case.vital.status'] == 'alive'))
    # ]
    # data.reset_index(drop=True, inplace=True)
    # Set the font sizes
    plt.rcParams.update({'font.size': 14, 'axes.titlesize': 18, 'axes.labelsize': 16,
                         'xtick.labelsize': 14, 'ytick.labelsize': 14, 'legend.fontsize': 14})

    kmf = KaplanMeierFitter()
    num_genes = min(len(genes), 6)  # Ensure we plot no more than 6 genes
    fig, axs = plt.subplots(2, 2, figsize=(15, 15))

    for i, gene in enumerate(genes[:num_genes]):
        row = i // 2
        col = i % 2
        ax = axs[row, col]

        data['high_expression'] = data[gene] > data[gene].median()

        for name, grouped_df in data.groupby('high_expression'):
            kmf.fit(grouped_df['Overall.Survival.Months'], grouped_df['case.vital.status'] == 'dead',
                    label=f'{gene} {"High" if name else "Low"}')
            kmf.plot_survival_function(ax=ax)

        ax.set_title(f'Kaplan-Meier Curves for {gene}')
        ax.set_xlabel('Months')
        ax.set_ylabel('Survival Probability')
        ax.legend()

    # Hide any unused subplots
    for j in range(num_genes, 4):
        row = j // 2
        col = j % 2
        fig.delaxes(axs[row, col])

    plt.tight_layout()
    # Save the figure as a PNG file with DPI 300
    plt.savefig(f'/home/mkh062/Desktop/scratch/TCGA_project/survival_analysis/kmc_GLASS.png',
                dpi=300)
    plt.show()


def plot_kaplan_meier_tcga():
    # threshold = 60
    # data1 = pd.read_csv(
    #     f'/home/mkh062/Desktop/scratch/TCGA_project/limma_analysis/significant_genes_limma_TCGA_{threshold}months.csv')
    data1 = pd.read_csv(
        f'/home/mkh062/Desktop/scratch/TCGA_project/survival_analysis/TCGA_significant_survival_statistics.csv')
    print("data1", data1.head())
    genes = data1['gene'].head(4).tolist()

    out = '/home/mkh062/Desktop/scratch/TCGA_project/Compiled data_TCGA_CGGA_GLASS_jh080424.xlsx'  # Update this path
    outcome_data = pd.read_excel(out, sheet_name='TCGA')
    gene_expression_file = '/home/mkh062/Desktop/scratch/TCGA_project/GBMLGG_EB_RmDiffFullGenesRanRmDup.csv'  # Path to the gene expression data CSV file
    D2 = pd.read_csv(gene_expression_file, sep=',')

    common_ids = pd.Series(list(set(D2['samples']).intersection(set(outcome_data['Sample ID'])))).tolist()
    D2_subset = D2[D2['samples'].isin(common_ids)]
    D2_subset.reset_index(drop=True, inplace=True)
    D2_subset.rename(columns={'samples': 'Sample ID'}, inplace=True)

    outcome_data = outcome_data[
        ['Sample ID', 'Overall Survival (Months)', 'Diagnosis Age', 'Sex', 'Overall Survival Status']]
    data = pd.merge(outcome_data, D2_subset, on='Sample ID')
    data.reset_index(drop=True, inplace=True)
    data['Overall Survival Status'] = data['Overall Survival Status'].replace({
        '1:DECEASED': 'dead',
        '0:LIVING': 'alive'
    })
    data.rename(columns={
        'Sample ID': 'Sample.ID',
        'Overall Survival (Months)': 'Overall.Survival.Months',
        'Diagnosis Age': 'Diagnosis.Age',
        'Overall Survival Status': 'case.vital.status'
    }, inplace=True)

    # Check for duplicated Sample IDs
    duplicated_ids = data['Sample.ID'][data['Sample.ID'].duplicated(keep=False)]

    # Print duplicated Sample IDs
    if not duplicated_ids.empty:
        print("Duplicated Sample IDs:")
        for id in duplicated_ids.unique():
            print(id)
    else:
        print("No duplicated Sample IDs found.")
    # Drop all duplicated Sample IDs
    data = data[~data['Sample.ID'].isin(duplicated_ids)]
    data.reset_index(drop=True, inplace=True)
    # data = data[
    #     ~((data['Overall.Survival.Months'] <= threshold) & (data['case.vital.status'] == 'alive'))
    # ]
    # data.reset_index(drop=True, inplace=True)
    # Set the font sizes
    plt.rcParams.update({'font.size': 14, 'axes.titlesize': 18, 'axes.labelsize': 16,
                         'xtick.labelsize': 14, 'ytick.labelsize': 14, 'legend.fontsize': 14})

    kmf = KaplanMeierFitter()
    num_genes = min(len(genes), 4)  # Ensure we plot no more than 6 genes
    fig, axs = plt.subplots(2, 2, figsize=(15, 15))

    for i, gene in enumerate(genes[:num_genes]):
        row = i // 2
        col = i % 2
        ax = axs[row, col]

        data['high_expression'] = data[gene] > data[gene].median()

        for name, grouped_df in data.groupby('high_expression'):
            kmf.fit(grouped_df['Overall.Survival.Months'], grouped_df['case.vital.status'] == 'dead',
                    label=f'{gene} {"High" if name else "Low"}')
            kmf.plot_survival_function(ax=ax)

        ax.set_title(f'Kaplan-Meier Curves for {gene}')
        ax.set_xlabel('Months')
        ax.set_ylabel('Survival Probability')
        ax.legend()

    # Hide any unused subplots
    for j in range(num_genes, 4):
        row = j // 2
        col = j % 2
        fig.delaxes(axs[row, col])

    plt.tight_layout()
    # Save the figure as a PNG file with DPI 300
    plt.savefig(f'/home/mkh062/Desktop/scratch/TCGA_project/survival_analysis/kmc_TCGA.png',
                dpi=300)
    plt.show()


plot_kaplan_meier_glass()
