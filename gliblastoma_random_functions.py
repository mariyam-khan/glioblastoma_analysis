import pandas as pd
import numpy as np

def remove_zero_genes(df):
    # Example: remove genes that have zero (or near zero) expression across all samples
    # Assuming df columns after 'Sample ID' are genes
    gene_cols = [c for c in df.columns if c != 'Sample ID']

    mask = (df[gene_cols].sum(axis=0) > 0)
    return df[['Sample ID'] + mask[mask].index.tolist()]

def log2_transform(df, gene_cols):
    # If the data is already non-integers (RSEM/FPKM/TPM), just log2(x+1)
    df[gene_cols] = np.log2(df[gene_cols] + 1)
    return df

def transform_sample_name_GLASS(name):
    # Extract part up to and including 'TP'
    if 'TP' in name:
        name = name.split('TP')[0] + 'TP'
    # Replace dots with dashes
    name = name.replace('.', '-')
    return name


def transform_sample_name_TCGA(sample_name):
    return '-'.join(sample_name.split('-')[:3])


# Log and Z-Score Transformations
def apply_log_transform(data):
    return data.apply(lambda x: np.log1p(x) + 1e-5 if np.issubdtype(x.dtype, np.number) else x)


def apply_zscore_transform(data):
    return data.apply(lambda x: zscore(x) if np.issubdtype(x.dtype, np.number) else x)