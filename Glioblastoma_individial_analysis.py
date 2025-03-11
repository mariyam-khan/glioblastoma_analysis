import pandas as pd
import numpy as np
from gliblastoma_random_functions import transform_sample_name_GLASS, remove_zero_genes, apply_log_transform, \
    apply_zscore_transform,log2_transform
import os
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


def TCGA_integrate_RNA_agilent():
    out = '/home/mkh062/Desktop/scratch/TCGA_project/Compiled data_TCGA_CGGA_GLASS_jh080424.xlsx'  # Update this path
    outcome_data = pd.read_excel(out, sheet_name='TCGA')
    gene_expression_file = '/home/mkh062/Desktop/scratch/TCGA_project/GBMLGG_EB_RmDiffFullGenesRanRmDup.csv'  # Path to the gene expression data CSV file
    D2 = pd.read_csv(gene_expression_file, sep=',')
    print("gene_expression_file TCGA", D2.head())
    common_ids = pd.Series(list(set(D2['samples']).intersection(set(outcome_data['Sample ID'])))).tolist()
    D2_subset = D2[D2['samples'].isin(common_ids)]
    D2_subset.reset_index(drop=True, inplace=True)
    D2_subset.rename(columns={'samples': 'Sample ID'}, inplace=True)
    print("DA", D2_subset.head())
    # 1. Subset to keep RNAseq and agilent only
    valid_platforms = ['RNAseq', 'agilent']
    D2_subset = D2_subset[D2_subset['platform'].isin(valid_platforms)].copy()
    gene_names = [col for col in D2_subset.columns if col not in ['Sample.ID', 'Platform']]
    #D2_subset = remove_zero_genes(D2_subset)
    # 2. Separate out the sample info and the gene expression matrix
    sample_ids = D2_subset['Sample ID']
    platforms = D2_subset['platform']

    # The gene columns are everything except 'Sample ID' and 'platform'
    gene_cols = [c for c in D2_subset.columns if c not in ['Sample ID', 'platform']]

    # Extract the expression matrix
    expr_matrix = D2_subset[gene_cols].copy()

    # 3. PCA before normalization
    pca = PCA(n_components=2)
    pca_scores = pca.fit_transform(expr_matrix.fillna(0))  # handle missing data somehow
    pc_df = pd.DataFrame(pca_scores, columns=['PC1', 'PC2'])
    pc_df['Sample ID'] = sample_ids.values
    pc_df['Platform'] = platforms.values

    # Plot the PCA before normalization
    plt.figure(figsize=(8, 6))
    sns.scatterplot(x='PC1', y='PC2', hue='Platform', data=pc_df, alpha=0.7)
    plt.title('PCA Before Normalization')
    plt.legend()
    plt.show()

    # 4. Normalize (z-score)
    # Note: For multi-platform data, you may need a more advanced batch correction method (e.g., ComBat)
    scaler = StandardScaler()
    expr_matrix_scaled = scaler.fit_transform(expr_matrix.fillna(0))

    # 5. PCA after normalization
    pca_norm = PCA(n_components=2)
    pca_scores_norm = pca_norm.fit_transform(expr_matrix_scaled)
    pc_norm_df = pd.DataFrame(pca_scores_norm, columns=['PC1', 'PC2'])
    pc_norm_df['Sample ID'] = sample_ids.values
    pc_norm_df['Platform'] = platforms.values

    # Plot the PCA after normalization
    plt.figure(figsize=(8, 6))
    sns.scatterplot(x='PC1', y='PC2', hue='Platform', data=pc_norm_df, alpha=0.7)
    plt.title('PCA After Normalization')
    plt.legend()
    plt.show()
    expr_matrix = pd.DataFrame(expr_matrix_scaled, index=sample_ids, columns=gene_cols)
    normalized_data = expr_matrix
    normalized_data.reset_index(drop=True, inplace=True)
    # Add sample and platform metadata
    normalized_data['Sample.ID'] = sample_ids.values
    print("normalized_data", normalized_data)
    # Merge normalized data with survival data
    outcome_data = outcome_data[
        ['Sample ID', 'Overall Survival (Months)', 'Diagnosis Age', 'Sex', 'Overall Survival Status']]
    outcome_data.rename(columns={
        'Sample ID': 'Sample.ID',
        'Overall Survival (Months)': 'Overall.Survival.Months',
        'Diagnosis Age': 'Diagnosis.Age',
        'Overall Survival Status': 'case.vital.status'
    }, inplace=True)

    # Replace survival status for consistency
    outcome_data['case.vital.status'] = outcome_data['case.vital.status'].replace({
        '1:DECEASED': 'dead',
        '0:LIVING': 'alive'
    })

    # Merge normalized gene data with outcome data
    merged_data = pd.merge(outcome_data, normalized_data, on='Sample.ID')
    merged_data.reset_index(drop=True, inplace=True)
    print("merged_data", merged_data)
    # Check for duplicated Sample IDs
    duplicated_ids = merged_data['Sample.ID'][merged_data['Sample.ID'].duplicated(keep=False)]

    # Print and drop duplicated Sample IDs
    if not duplicated_ids.empty:
        print("Duplicated Sample IDs:")
        for id in duplicated_ids.unique():
            print(id)
        # Drop all duplicated Sample IDs
        merged_data = merged_data[~merged_data['Sample.ID'].isin(duplicated_ids)]
    else:
        print("No duplicated Sample IDs found.")
    merged_data[gene_cols] = merged_data[gene_cols].loc[:, merged_data[gene_cols].var() != 0]
    # merged_data.to_csv('/home/mkh062/Desktop/scratch/TCGA_project/Dec20/TCGA_integrate/TCGA_integrated.csv', index=False)
    print("PCA plot saved and combined data written to file.")
    # thresholds = [24, 36, 48, 60]
    # datasets = [merged_data]
    # dataset_names = ["TCGA_integrated"]
    # for threshold in thresholds:
    #     for data, name in zip(datasets, dataset_names):
    #         print(f"########### {name} threshold {threshold}")
    #         limma(data, threshold, name, "fixed")


def Glass_data():
    out = '/home/mkh062/Desktop/scratch/TCGA_project/Compiled data_TCGA_CGGA_GLASS_jh080424.xlsx'  # Update this path
    outcome_data = pd.read_excel(out, sheet_name='GLASS')

    gene_expression_file = '/home/mkh062/Desktop/scratch/TCGA_project/gene_tpm_matrix_all_samples.tsv'  # Path to the gene expression data CSV file
    D2 = pd.read_csv(gene_expression_file, sep='\t')

    # Transpose the gene expression data and reset the index
    D2_transposed = D2.set_index('Gene_symbol').transpose().reset_index()
    new_columns = ['Sample ID'] + D2['Gene_symbol'].tolist()
    D2_transposed.columns = new_columns

    # Apply sample name transformation (assuming transform_sample_name is defined elsewhere)
    D2_transposed['Sample ID'] = D2_transposed['Sample ID'].apply(transform_sample_name_GLASS)
    common_ids = pd.Series(list(set(D2_transposed['Sample ID']).intersection(set(outcome_data['Sample ID'])))).tolist()

    D2_subset = D2_transposed[D2_transposed['Sample ID'].isin(common_ids)]
    D2_subset.reset_index(drop=True, inplace=True)
    print("GLASS gene expr", D2_subset.head())
    #####################################################################################3
    D2_subset = remove_zero_genes(D2_subset)
    # Continue with outcome data processing
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

    if not duplicated_ids.empty:
        print("Duplicated Sample IDs:")
        for id in duplicated_ids.unique():
            print(id)
    else:
        print("No duplicated Sample IDs found.")

    # Drop all duplicated Sample IDs
    data = data[~data['Sample.ID'].isin(duplicated_ids)]
    data.reset_index(drop=True, inplace=True)

    gene_cols = data.columns.difference(
        ['Sample.ID', 'Overall.Survival.Months', 'Diagnosis.Age', 'Sex', 'case.vital.status'])
    #data[gene_cols] = data[gene_cols].loc[:, data[gene_cols].var() != 0]
    data[gene_cols] = np.log1p(data[gene_cols] + 1e-5)
    scaler = StandardScaler()
    data[gene_cols] = scaler.fit_transform(data[gene_cols])

    #data.to_csv('/home/mkh062/Desktop/scratch/TCGA_project/Dec20/processed_glass_data.csv', index=False)
    return data

def CGGA_data():
    # Load clinical and filter
    clinical_CGGA_325 = pd.read_csv(
        '/home/mkh062/Desktop/scratch/TCGA_project/CGGA-RAW DATA/CGGA.mRNAseq_325_clinical.20200506.txt',
        sep='\t'
    )
    clinical_CGGA_325 = clinical_CGGA_325[(clinical_CGGA_325['PRS_type'] == 'Primary') &
                                          (clinical_CGGA_325['Histology'] == 'GBM') &
                                          (clinical_CGGA_325['IDH_mutation_status'] == 'Wildtype')]
    clinical_CGGA_325 = clinical_CGGA_325[['CGGA_ID', 'Gender', 'Age', 'Censor (alive=0; dead=1)', 'OS']].dropna()
    clinical_CGGA_325['OS'] = clinical_CGGA_325['OS'] / 30.44
    clinical_CGGA_325.rename(columns={'CGGA_ID':'Sample ID'}, inplace=True)

    # Load clinical data for CGGA 693
    clinical_CGGA_693 = pd.read_csv(
        '/home/mkh062/Desktop/scratch/TCGA_project/CGGA-RAW DATA/CGGA.mRNAseq_693_clinical.20200506.txt',
        sep='\t'
    )
    clinical_CGGA_693 = clinical_CGGA_693[(clinical_CGGA_693['PRS_type'] == 'Primary') &
                                          (clinical_CGGA_693['Histology'] == 'GBM') &
                                          (clinical_CGGA_693['IDH_mutation_status'] == 'Wildtype')]
    clinical_CGGA_693 = clinical_CGGA_693[['CGGA_ID', 'Gender', 'Age', 'Censor (alive=0; dead=1)', 'OS']].dropna()
    clinical_CGGA_693['OS'] = clinical_CGGA_693['OS'] / 30.44
    clinical_CGGA_693.rename(columns={'CGGA_ID':'Sample ID'}, inplace=True)

    # Load gene expression data (raw counts)
    gene_expression_CGGA_325 = pd.read_csv(
        '/home/mkh062/Desktop/scratch/TCGA_project/CGGA.mRNAseq_325.RSEM-genes.20200506.txt/CGGA.mRNAseq_325.RSEM-genes.20200506.txt',
        sep='\t'
    )
    # Load gene expression data for CGGA 693
    gene_expression_CGGA_693 = pd.read_csv(
        '/home/mkh062/Desktop/scratch/TCGA_project/CGGA.mRNAseq_693.RSEM-genes.20200506.txt/CGGA.mRNAseq_693.RSEM-genes.20200506.txt',
        sep='\t'
    )
    print("CGGA_325 gene expr", gene_expression_CGGA_325.head())
    print("CGGA_693 gene expr", gene_expression_CGGA_693.head())

    # Transpose so that rows are samples and columns are genes
    def reshape_expression(df):
        df_t = df.set_index('Gene_Name').T.reset_index()
        new_cols = ['Sample ID'] + df['Gene_Name'].tolist()
        df_t.columns = new_cols
        return df_t

    gene_expression_CGGA_325 = reshape_expression(gene_expression_CGGA_325)
    gene_expression_CGGA_693 = reshape_expression(gene_expression_CGGA_693)
    print("CGGA_325 gene expr res", gene_expression_CGGA_325)
    print("CGGA_693 gene expr res", gene_expression_CGGA_693)

    # Filter to common samples
    common_ids_325 = set(gene_expression_CGGA_325['Sample ID']).intersection(clinical_CGGA_325['Sample ID'])
    gene_expression_CGGA_325 = gene_expression_CGGA_325[gene_expression_CGGA_325['Sample ID'].isin(common_ids_325)]

    common_ids_693 = set(gene_expression_CGGA_693['Sample ID']).intersection(clinical_CGGA_693['Sample ID'])
    gene_expression_CGGA_693 = gene_expression_CGGA_693[gene_expression_CGGA_693['Sample ID'].isin(common_ids_693)]

    gene_expression_CGGA_325 = remove_zero_genes(gene_expression_CGGA_325)
    gene_expression_CGGA_693 = remove_zero_genes(gene_expression_CGGA_693)

    # Find common genes between batches
    common_genes = list(set(gene_expression_CGGA_325.columns) & set(gene_expression_CGGA_693.columns) - {'Sample ID'})
    batch1_expr = gene_expression_CGGA_325[['Sample ID'] + common_genes].copy()
    batch2_expr = gene_expression_CGGA_693[['Sample ID'] + common_genes].copy()


    batch1_expr = log2_transform(batch1_expr, common_genes)
    batch2_expr = log2_transform(batch2_expr, common_genes)

    # Merge with clinical data
    batch1_data = pd.merge(clinical_CGGA_325, batch1_expr, on='Sample ID')
    batch2_data = pd.merge(clinical_CGGA_693, batch2_expr, on='Sample ID')

    batch1_data['batch'] = 'batch325'
    batch2_data['batch'] = 'batch693'

    combined_data = pd.concat([batch1_data, batch2_data], ignore_index=True)

    # PCA before outlier removal
    gene_cols_for_pca = [g for g in combined_data.columns if g not in ['Sample ID','OS','Age','Gender','Censor (alive=0; dead=1)','batch']]
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(combined_data[gene_cols_for_pca])
    combined_data['PC1'] = pca_result[:,0]
    combined_data['PC2'] = pca_result[:,1]

    plt.figure(figsize=(10,8))
    sns.scatterplot(x='PC1', y='PC2', hue='batch', data=combined_data)
    plt.title('PCA before Outlier Removal')
    #plt.savefig('/home/mkh062/Desktop/scratch/TCGA_project/Dec20/PCA_before_outlier_removal.png')
    plt.show()

    # Identify outliers (arbitrary threshold shown as example)
    outlier_indices = combined_data[(np.abs(combined_data['PC1'])>200) | (np.abs(combined_data['PC2'])>200)].index
    combined_data = combined_data.drop(outlier_indices).reset_index(drop=True)
    print("combined_data", combined_data)
    common_sample_ids_batch1 = set(combined_data['Sample ID']).intersection(set(batch1_data['Sample ID']))
    batch1_data = batch1_data[batch1_data['Sample ID'].isin(common_sample_ids_batch1)].reset_index(drop=True)

    # Take intersection of 'Sample ID' between batch2_data and combined_data
    common_sample_ids_batch2 = set(combined_data['Sample ID']).intersection(set(batch2_data['Sample ID']))
    batch2_data = batch2_data[batch2_data['Sample ID'].isin(common_sample_ids_batch2)].reset_index(drop=True)

    print("batch1_data n", batch1_data)
    print("batch2_data n", batch2_data)
    # Re-run PCA after outlier removal
    pca_result = pca.fit_transform(combined_data[gene_cols_for_pca])
    combined_data['PC1'] = pca_result[:,0]
    combined_data['PC2'] = pca_result[:,1]

    plt.figure(figsize=(10,8))
    sns.scatterplot(x='PC1', y='PC2', hue='batch', data=combined_data)
    #plt.savefig('/home/mkh062/Desktop/scratch/TCGA_project/Dec20/PCA_after_outlier_removal.png')
    plt.title('PCA after Outlier Removal')
    plt.show()

    combined_data.rename(columns={
        'Sample ID': 'Sample.ID',
        'OS': 'Overall.Survival.Months',
        'Age': 'Diagnosis.Age',
        'Gender': 'Sex',
        'Censor (alive=0; dead=1)': 'case.vital.status'
    }, inplace=True)
    clinical_cols = ['Sample.ID','Overall.Survival.Months','Diagnosis.Age','Sex','case.vital.status','batch','PC1','PC2']
    gene_cols = [c for c in combined_data.columns if c not in clinical_cols]

    scaler = StandardScaler()
    combined_data_norm = combined_data
    combined_data_norm[gene_cols] = scaler.fit_transform(combined_data[gene_cols])
    print("combined_data_norm", combined_data_norm)

    pca_result = pca.fit_transform(combined_data_norm[gene_cols_for_pca])
    combined_data_norm['PC1'] = pca_result[:,0]
    combined_data_norm['PC2'] = pca_result[:,1]

    plt.figure(figsize=(10,8))
    sns.scatterplot(x='PC1', y='PC2', hue='batch', data=combined_data_norm)
    plt.title('PCA after Outlier Removal, concatenation and normalisation')
    #plt.savefig('/home/mkh062/Desktop/scratch/TCGA_project/Dec20/PCA_after_outlier_removal_combined_nor.png')
    plt.show()

    # Assuming gene_cols and other variables are already defined
    scaler = StandardScaler()

    # Scale the data and convert back to DataFrame
    batch1_data_norm = pd.DataFrame(
        scaler.fit_transform(batch1_data[gene_cols]),
        columns=gene_cols,
        index=batch1_data.index
    )

    batch2_data_norm = pd.DataFrame(
        scaler.fit_transform(batch2_data[gene_cols]),
        columns=gene_cols,
        index=batch2_data.index
    )
    batch1_data.rename(columns={
        'Sample ID': 'Sample.ID',
        'OS': 'Overall.Survival.Months',
        'Age': 'Diagnosis.Age',
        'Gender': 'Sex',
        'Censor (alive=0; dead=1)': 'case.vital.status'
    }, inplace=True)
    batch2_data.rename(columns={
        'Sample ID': 'Sample.ID',
        'OS': 'Overall.Survival.Months',
        'Age': 'Diagnosis.Age',
        'Gender': 'Sex',
        'Censor (alive=0; dead=1)': 'case.vital.status'
    }, inplace=True)
    # Add back the clinical columns to each batch
    batch1_data_norm = pd.concat([batch1_data[['Sample.ID', 'Overall.Survival.Months', 'Diagnosis.Age', 'Sex', 'case.vital.status', 'batch']], batch1_data_norm], axis=1)
    batch2_data_norm = pd.concat([batch2_data[['Sample.ID', 'Overall.Survival.Months', 'Diagnosis.Age', 'Sex', 'case.vital.status', 'batch']], batch2_data_norm], axis=1)

    # Combine the normalized data
    combined_data_2 = pd.concat([batch1_data_norm, batch2_data_norm], ignore_index=True)

    pca_result = pca.fit_transform(combined_data_2[gene_cols_for_pca])
    combined_data_2['PC1'] = pca_result[:,0]
    combined_data_2['PC2'] = pca_result[:,1]

    plt.figure(figsize=(10,8))
    sns.scatterplot(x='PC1', y='PC2', hue='batch', data=combined_data_2)
    plt.title('PCA after Outlier Removal, normalisation and concatenation')
    # plt.savefig('/home/mkh062/Desktop/scratch/TCGA_project/Dec20/PCA_after_outlier_removal_nor_combined.png')
    plt.show()
    #
    # print("combined_data_2", combined_data_2)
    # print("combined_data_norm", combined_data_norm)
    # # Save normalized data
    # combined_data_norm.to_csv('/home/mkh062/Desktop/scratch/TCGA_project/Dec20/combined_CGGA_combined_nor.csv', index=False)
    # combined_data_2.to_csv('/home/mkh062/Desktop/scratch/TCGA_project/Dec20/processed_cgga_nor_combined.csv', index=False)
    # batch1_data_norm.to_csv('/home/mkh062/Desktop/scratch/TCGA_project/Dec20/batch1_data_norm.csv', index=False)
    # batch2_data_norm.to_csv('/home/mkh062/Desktop/scratch/TCGA_project/Dec20/batch2_data_norm.csv', index=False)

TCGA_integrate_RNA_agilent()
CGGA_data()
Glass_data()