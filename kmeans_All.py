import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans


def load_de_genes(file1, file2):
    # Reads the DE gene files and returns the union of genes (from the 'Gene' column)
    df1 = pd.read_csv(file1)
    df2 = pd.read_csv(file2)
    genes1 = set(df1['Gene'].unique())
    genes2 = set(df2['Gene'].unique())
    de_genes = list(genes1.union(genes2))
    print(f"Total DE genes loaded: {len(de_genes)}")
    return de_genes


def plot_gene_pca_with_de_annotation(filepath, de_genes, n_clusters=4, n_components=2,output_file="pca_de_annotation.png"):
    df, gene_cols = load_glass_data(filepath)
    gene_data = df[gene_cols].transpose()  # rows: genes, columns: samples
    de_genes_in_data = [gene for gene in de_genes if gene in gene_data.index]
    if not de_genes_in_data:
        print("No DE genes found in the gene expression data.")
        return
    pca = PCA(n_components=n_components)
    gene_pcs = pca.fit_transform(gene_data)
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    gene_clusters = kmeans.fit_predict(gene_pcs)

    plt.figure(figsize=(8, 6))
    scatter = plt.scatter(gene_pcs[:, 0], gene_pcs[:, 1], c=gene_clusters, cmap='viridis', alpha=0.7)
    # Annotate only the DE genes in red
    for i, gene in enumerate(gene_data.index):
        if gene in de_genes_in_data:
            plt.annotate(gene, (gene_pcs[i, 0], gene_pcs[i, 1]), fontsize=8, color="red")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.title("PCA of Genes with DE Gene Annotation")
    plt.legend(*scatter.legend_elements(), title="Cluster")
    # plt.savefig(output_file, dpi=300)
    # plt.close()
# ------------------------------
# Data Loading Function
# ------------------------------
def load_glass_data(filepath):
    df = pd.read_csv(filepath)
    metadata_cols = ['Sample.ID', 'Overall.Survival.Months', 'Diagnosis.Age', 'Sex', 'case.vital.status']
    gene_cols = [col for col in df.columns if col not in metadata_cols]
    return df, gene_cols
# ------------------------------
# 1. K-means Clustering for Genes
# ------------------------------
def kmeans_cluster_genes(filepath, n_clusters=4):
    df, gene_cols = load_glass_data(filepath)
    # Transpose the gene expression matrix: rows become genes, columns are samples.
    gene_data = df[gene_cols].transpose()
    case = "gene_kmeans"
    # if n_clusters==2:
    #     elbow_plot_kmeans(gene_data, case, max_clusters=10)
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    gene_clusters = kmeans.fit_predict(gene_data)
    return df, kmeans

def kmeans_cluster_genes_pca(filepath, n_clusters=4, n_components=2):
    df, gene_cols = load_glass_data(filepath)
    gene_data = df[gene_cols].transpose()  # rows: genes, columns: samples
    pca = PCA(n_components=n_components)
    gene_pcs = pca.fit_transform(gene_data)
    case = "genePCs_kmeans"
    # if n_clusters == 2:
    #     elbow_plot_kmeans(gene_pcs, case, max_clusters=10)
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    gene_clusters = kmeans.fit_predict(gene_pcs)
    plt.figure(figsize=(8, 6))
    scatter = plt.scatter(gene_pcs[:, 0], gene_pcs[:, 1], c=gene_clusters, cmap='viridis', alpha=0.7)
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.title("PCA of GLASS with K-means Clustering by the genes, #clusters= " + str(n_clusters))
    plt.legend(*scatter.legend_elements(), title="Cluster")
    # plt.savefig("/home/mkh062/Desktop/scratch/TCGA_project/processed_data/results/PCA_plot" + str(case) + "_" + str(n_clusters) + ".png",
    #             dpi=300)
    #gene_cluster_df = pd.DataFrame({'Gene': gene_data.index, 'Cluster': gene_clusters})
    return df, kmeans, pca
# ------------------------------
# 2. K-means Clustering for Samples
# ------------------------------
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
    plt.title("PCA of GLASS with K-means Clustering by the Samples, #clusters= " + str(n_clusters))
    plt.legend(*scatter.legend_elements(), title="Cluster")
    # plt.savefig("/home/mkh062/Desktop/scratch/TCGA_project/processed_data/results/PCA_plot" + str(case) + "_" + str(n_clusters) + ".png",
    #             dpi=300)
    df['Cluster'] = sample_clusters
    return df, kmeans, pca

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
    plt.savefig("/home/mkh062/Desktop/scratch/TCGA_project/processed_data/results/elbow_plot" + str(case) + ".png", dpi=300)


def overlay_clinical_data(df, case_p, n_cluster, cluster_col='Cluster'):
    # Boxplot: Overall Survival Months by Cluster
    plt.figure(figsize=(10, 6))
    sns.boxplot(x=cluster_col, y='Overall.Survival.Months', data=df)
    sns.swarmplot(x=cluster_col, y='Overall.Survival.Months', data=df, color=".25")
    plt.title('Overall Survival Months by Cluster')
    #plt.savefig("/home/mkh062/Desktop/scratch/TCGA_project/processed_data/results/kmeans_survival" + str(case_p) + "_" + str(n_cluster) + ".png", dpi=300)

    # Bar plot: Case Vital Status distribution by Cluster
    plt.figure(figsize=(10, 6))
    status_counts = pd.crosstab(df[cluster_col], df['case.vital.status'])
    status_counts.plot(kind='bar', stacked=True)
    plt.title('Case Vital Status by Cluster')
    plt.xlabel('Cluster')
    plt.ylabel('Count')
    #plt.savefig("/home/mkh062/Desktop/scratch/TCGA_project/processed_data/results/kmeans_vital" + str(case_p) + "_" + str(n_cluster) + ".png", dpi=300)


from scipy.stats import kruskal  # For the non-parametric test
# ------------------------------
# 4. Kruskal-Wallis Test for Survival across Clusters
# ------------------------------
def kruskal_test_survival(df, cluster_col='Cluster', survival_col='Overall.Survival.Months'):
    """
    Groups samples by cluster and performs a Kruskal-Wallis test
    to see if the survival distributions differ significantly.
    """
    # Group survival times by cluster (dropping missing values if any)
    groups = [group[survival_col].dropna().values for _, group in df.groupby(cluster_col)]
    stat, p_val = kruskal(*groups)
    print("Kruskal-Wallis test result:")
    print(f"Statistic: {stat:.3f}, p-value: {p_val:.3f}")
    return stat, p_val

if __name__ == "__main__":
    filepath = "/home/mkh062/Desktop/scratch/TCGA_project/processed_data/GLASS_processed.csv"
    de_file1 = "/home/mkh062/Desktop/scratch/TCGA_project/processed_data/results/significant_genes_GLASS_24-48-months.csv"
    de_file2 = "/home/mkh062/Desktop/scratch/TCGA_project/processed_data/results/significant_genes_GLASS_24-60-months.csv"
    de_genes = load_de_genes(de_file1, de_file2)
    no_clusters = [2,3,4,5,6]
    # K-means clustering for genes (columns)
    for k in no_clusters:
        gene_clusters_df, gene_kmeans = kmeans_cluster_genes(filepath, n_clusters=k)
        gene_clusters_pca_df, gene_kmeans_pca, gene_pca = kmeans_cluster_genes_pca(filepath, n_clusters=k, n_components=2)
        df_clustered_samples, sample_kmeans = kmeans_cluster_samples(filepath, n_clusters=k)
        df_clustered_samples_pca, sample_kmeans_pca, sample_pca = kmeans_cluster_samples_pca(filepath, n_clusters=k,n_components=2)
        plot_gene_pca_with_de_annotation(filepath, de_genes, n_clusters=k, n_components=2,
                                         output_file=f"/home/mkh062/Desktop/scratch/TCGA_project/processed_data/results/PCA_de_annotation_k{k}.png")
        case = "sample_kmeans"
        overlay_clinical_data(df_clustered_samples, case, k, cluster_col='Cluster')
        case = "samplePCs_kmeans"
        overlay_clinical_data(df_clustered_samples_pca, case, k, cluster_col='Cluster')
        print("################### ", k)
        stat, p_val = kruskal_test_survival(df_clustered_samples, cluster_col='Cluster', survival_col='Overall.Survival.Months')