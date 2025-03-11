import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans


def elbow_plot_kmeans(X, case, max_clusters=10):
    inertias = []
    cluster_range = range(1, max_clusters + 1)

    # Compute the inertia for each k.
    for k in cluster_range:
        kmeans = KMeans(n_clusters=k, random_state=42)
        kmeans.fit(X)
        inertias.append(kmeans.inertia_)
    reductions = [inertias[i - 1] - inertias[i] for i in range(1, len(inertias))]

    plt.figure(figsize=(8, 6))
    plt.plot(range(2, max_clusters + 1), reductions, marker='o')
    plt.xlabel("Number of Clusters (k)", fontsize=12)
    plt.ylabel("Reduction in WCSS", fontsize=12)
    plt.title("Reduction in Variation vs. Number of Clusters", fontsize=12)
    plt.xticks(range(2, max_clusters + 1))
    plt.show()
    #plt.savefig("/home/mkh062/Desktop/scratch/TCGA_project/processed_data/results/elbow_plot" + str(case) + ".png", dpi=300)

def load_glass_data(filepath):
    df = pd.read_csv(filepath)
    metadata_cols = ['Sample.ID', 'Overall.Survival.Months', 'Diagnosis.Age', 'Sex', 'case.vital.status']
    gene_cols = [col for col in df.columns if col not in metadata_cols]
    return df, gene_cols

def kmeans_cluster_samples(filepath, n_clusters=4):
    df, gene_cols = load_glass_data(filepath)
    print("sample k means", df[gene_cols].head())
    sample_data = df[gene_cols].values  # rows: samples, columns: genes
    case = "sample_kmeans"
    if n_clusters == 2:
        elbow_plot_kmeans(sample_data, case, max_clusters=10)
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    sample_clusters = kmeans.fit_predict(sample_data)
    df['Cluster'] = sample_clusters
    return df, kmeans

def kmeans_cluster_samples_pca(filepath, n_clusters=4, n_components=2):
    df, gene_cols = load_glass_data(filepath)
    sample_data = df[gene_cols].values
    pca = PCA(n_components=n_components)
    sample_pcs = pca.fit_transform(sample_data)
    case = "samplePCs_kmeans"
    # if n_clusters == 2:
    #     elbow_plot_kmeans(sample_pcs, case, max_clusters=10)
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    sample_clusters = kmeans.fit_predict(sample_pcs)
    plt.figure(figsize=(8, 6))
    scatter = plt.scatter(sample_pcs[:, 0], sample_pcs[:, 1], c=sample_clusters, cmap='viridis', alpha=0.7)
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.title("PCA of TCGA with K-means Clustering by the Samples, #clusters= " + str(n_clusters))
    plt.legend(*scatter.legend_elements(), title="Cluster")
    # plt.savefig("/home/mkh062/Desktop/scratch/TCGA_project/processed_data/results/PCA_plot" + str(case) + "_" + str(n_clusters) + ".png",
    #             dpi=300)
    plt.show()
    df['Cluster'] = sample_clusters
    return df, kmeans, pca


def load_de_genes(file1):
    # Reads the DE gene files and returns the union of genes (from the 'Gene' column)
    df1 = pd.read_csv(file1)
    genes1 = set(df1['Gene'].unique())
    de_genes = list(genes1)
    print(f"Total DE genes loaded: {len(de_genes)}")
    return de_genes

if __name__ == "__main__":
    filepath = "/home/mkh062/Desktop/scratch/TCGA_project/processed_data/TCGA_processed.csv"
    de_file1 = "/home/mkh062/Desktop/scratch/TCGA_project/processed_data/cox_results/TCGA_cox_significant_genes.csv"

    de_genes = load_de_genes(de_file1)
    no_clusters = [2,3,4,5,6]
    # K-means clustering for genes (columns)
    for k in no_clusters:
        df_clustered_samples, sample_kmeans = kmeans_cluster_samples(filepath, n_clusters=k)
        df_clustered_samples_pca, sample_kmeans_pca, sample_pca = kmeans_cluster_samples_pca(filepath, n_clusters=k,n_components=2)
        case = "sample_kmeans"
        case = "samplePCs_kmeans"
        print("################### ", k)
        #stat, p_val = kruskal_test_survival(df_clustered_samples, cluster_col='Cluster', survival_col='Overall.Survival.Months')