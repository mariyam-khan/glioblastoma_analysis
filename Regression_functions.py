import pandas as pd
import numpy as np
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from lifelines import CoxPHFitter
from lifelines import KaplanMeierFitter
from lifelines.exceptions import ConvergenceError
from statsmodels.stats.multitest import multipletests

def limma(data, threshold, consortium, threshold_type="fixed"):
    # Exclude alive patients with survival less than threshold
    if threshold_type == "fixed":
        print("original", data.shape[0])
        data = data[
            ~((data['Overall.Survival.Months'] <= threshold) & (data['case.vital.status'] == 'alive'))
        ]
        print("after drop", data.shape[0])
        file_suffix = f"{consortium}_{threshold}-months"
    elif threshold_type == "range":
        low, high = threshold
        print("original", data.shape[0])
        data = data[
            ~((data['Overall.Survival.Months'] <= low) & (data['case.vital.status'] == 'alive'))
        ]
        data = data[
            (data['Overall.Survival.Months'] <= low) | (data['Overall.Survival.Months'] > high)
            ]
        print("after drop", data.shape[0])
        file_suffix = f"{consortium}_{low}-{high}-months"

    data.drop(columns=['case.vital.status'], inplace=True)
    data.reset_index(drop=True, inplace=True)

    # Remove genes with zero variance
    if 'batch' in data.columns:
        print("batch in columns")
        clinical_cols = ['Sample.ID', 'Overall.Survival.Months', 'Diagnosis.Age', 'Sex', 'batch']
    else:
        clinical_cols = ['Sample.ID', 'Overall.Survival.Months', 'Diagnosis.Age', 'Sex']

    gene_cols = data.columns.difference(clinical_cols)
    gene_data = data[gene_cols]
    gene_data = gene_data.loc[:, gene_data.var(axis=0) != 0]
    data = pd.concat([data[clinical_cols], gene_data], axis=1)

    # R integration
    pandas2ri.activate()
    r_data = pandas2ri.py2rpy(data)
    design_formula = "~ group + Diagnosis.Age + Sex + batch" if 'batch' in data.columns else "~ group + Diagnosis.Age + Sex"
    exclude_cols_r = "c('Sample.ID', 'Overall.Survival.Months', 'Diagnosis.Age', 'Sex', 'batch')" if 'batch' in data.columns else "c('Sample.ID', 'Overall.Survival.Months', 'Diagnosis.Age', 'Sex')"

    r_script = f"""
    library(limma)
    gene_expression_data <- data[, !names(data) %in% {exclude_cols_r}]
    rownames(gene_expression_data) <- data$Sample.ID
    #data$group <- ifelse(data$Overall.Survival.Months > {threshold}, 'high', 'low')
    median_survival <- {threshold if threshold_type == "fixed" else threshold[0]}
    data$group <- ifelse(data$Overall.Survival.Months > median_survival, 'high', 'low')
    data$group <- factor(data$group, levels = c('low', 'high'))
    design <- model.matrix({design_formula}, data = data)
    rownames(design) <- data$Sample.ID
    fit <- lmFit(t(gene_expression_data), design)
    fit <- eBayes(fit)
    results <- topTable(fit, coef = 'grouphigh', number = Inf, adjust.method = 'fdr')
    #write.csv(results, file = '/home/mkh062/Desktop/scratch/TCGA_project/Dec20/results_limma_{file_suffix}.csv', row.names = TRUE)
    significant_genes <- results[results$adj.P.Val < 0.1, ]
    write.csv(significant_genes, file = '/home/mkh062/Desktop/scratch/TCGA_project/Jan20/significant_genes_{file_suffix}.csv', row.names = TRUE)
    """

    ro.globalenv['data'] = r_data
    ro.r(r_script)


def cox_reg(data):
    columns_to_drop = ['PC1', 'PC2', 'batch']

    if set(columns_to_drop).intersection(data.columns):
        data = data.drop(columns=columns_to_drop, errors='ignore')
    print("cgga", data)
    columns_to_exclude = {'Sample.ID', 'Overall.Survival.Months', 'Diagnosis.Age', 'Sex', 'case.vital.status'}
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
    print("merged_data['Sex']", merged_data['Sex'].values)

    merged_data['Sex'] = merged_data['Sex'].map({'Male': 0, 'Female': 1})
    cph = CoxPHFitter()
    # Create a list to store the statistics dictionaries
    statistics_list = []
    """
    df (DataFrame) – a Pandas DataFrame with necessary columns duration_col and event_col (see below), 
    covariates columns, and special columns (weights, strata). duration_col refers to the lifetimes 
    of the subjects. event_col refers to whether the ‘death’ events was observed:
     1 if observed, 0 else (censored).
    """
    print(merged_data)
    # Perform CoxPH regression for each gene
    for gene in merged_data.columns[5:]:  # Skip sample_id, overall_survival, Status
        print("gene", gene)
        print("gene", merged_data[gene].values)

        gene_data = merged_data[['Overall.Survival.Months', 'Status', 'Diagnosis.Age', 'Sex', gene]].dropna()
        print("gene_data", gene_data)
        # gene_data = merged_data[['Overall.Survival.Months', 'Status', 'Diagnosis.Age', 'Sex', 'platform', gene]].dropna()
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

    # Save significant genes statistics to a CSV file
    significant_genes.to_csv(
        '/home/mkh062/Desktop/scratch/TCGA_project/Dec20/survival_analysis/TCGA_CGGA_integrated_survival_statistics.csv',
        index=False)

    # print("All statistics saved to 'gene_survival_statistics.csv'")
    # print("Significant gene statistics saved to 'significant_gene_survival_statistics.csv'")
    # # Plot the adjusted p-values
    # plt.figure(figsize=(12, 6))
    # sns.barplot(x='gene', y=-np.log10(statistics['p_adjusted']), data=statistics)
    # plt.xticks(rotation=90)
    # plt.title('Adjusted p-values for Differential Gene Expression Analysis TCGA data')
    # plt.xlabel('Gene')

    # plt.ylabel('-log10(adjusted p-value)')
    # plt.show()