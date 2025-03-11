import os
import pandas as pd
import itertools

# Define the directory, platforms, and thresholds.
directory = '/home/mkh062/Desktop/scratch/TCGA_project/processed_data/results'
platforms = ['CGGA', 'TCGA', 'GLASS']
thresholds = ['24', '60', '36', '48', '24-48', '24-60']

# Loop over each threshold.
for threshold in thresholds:
    print(f"Threshold {threshold}:")

    # Dictionary to hold gene sets for each platform.
    platform_genes = {}

    # Read gene sets for each platform.
    for platform in platforms:
        filename = f"significant_genes_{platform}_{threshold}-months.csv"
        filepath = os.path.join(directory, filename)
        df = pd.read_csv(filepath)
        genes = set(df['Gene'].dropna().unique())
        platform_genes[platform] = genes

    # Compute common genes for every combination of two platforms.
    for pair in itertools.combinations(platforms, 2):
        common_genes = platform_genes[pair[0]].intersection(platform_genes[pair[1]])
        print(f"Common genes for {pair[0]} and {pair[1]}:")
        print(sorted(common_genes))
        print()
    print("---------------------------------------------------\n")

import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans


def load_glass_data(filepath):
    """
    Loads the GLASS gene expression dataset and identifies the gene expression columns.
    Assumes that the metadata columns are:
      ['Sample.ID','Overall.Survival.Months','Diagnosis.Age','Sex','case.vital.status']
    """
    df = pd.read_csv(filepath)
    metadata_cols = ['Sample.ID', 'Overall.Survival.Months', 'Diagnosis.Age', 'Sex', 'case.vital.status']
    # Gene expression columns are those not in metadata_cols.
    gene_cols = [col for col in df.columns if col not in metadata_cols]
    return df, gene_cols


def pca_by_survival(filepath):
    """
    Performs PCA on the gene expression data and plots the first two principal components,
    coloring each sample based on its Overall Survival Months.

    The survival groups are defined as:
      - 0-24 months:       red
      - >24 and ≤48 months: orange
      - >48 and ≤60 months: green
      - >60 months:         blue
    """
    df, gene_cols = load_glass_data(filepath)
    X = df[gene_cols].values

    # Perform PCA to reduce to 2 dimensions.
    pca = PCA(n_components=2)
    pcs = pca.fit_transform(X)


    # Define a function to assign survival groups.
    def assign_group(survival):
        if survival <= 24:
            return "0-24"
        elif 24 < survival <= 48:
            return "24-48"
        elif 48 < survival <= 60:
            return "48-60"
        else:
            return ">60"

    df['Survival_Group'] = df['Overall.Survival.Months'].apply(assign_group)
    print(df[['Survival_Group', 'Overall.Survival.Months']])


    # Define colors for the survival groups.
    group_colors = {"0-24": "red", "24-48": "orange", "48-60": "green", ">60": "blue"}
    colors = df['Survival_Group'].map(group_colors)

    plt.figure(figsize=(8, 6))
    plt.scatter(pcs[:, 0], pcs[:, 1], c=colors, alpha=0.7)
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.title("PCA of GLASS Gene Expression (Colored by Survival Group)")
    plt.savefig("/home/mkh062/Desktop/scratch/TCGA_project/processed_data/results/PCA_glass_survival.png", dpi=300)

    # Create a custom legend.
    import matplotlib.patches as mpatches
    legend_patches = [mpatches.Patch(color=group_colors[g], label=g) for g in group_colors]
    plt.legend(handles=legend_patches)
    plt.show()


def pca_by_age(filepath):
    """
    Performs PCA on the gene expression data and plots the first two principal components.
    Each point is annotated with its Diagnosis Age.
    """
    df, gene_cols = load_glass_data(filepath)
    X = df[gene_cols].values

    # Perform PCA.
    pca = PCA(n_components=2)
    pcs = pca.fit_transform(X)

    plt.figure(figsize=(10, 8))
    plt.scatter(pcs[:, 0], pcs[:, 1], alpha=0.7)
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.title("PCA of GLASS Gene Expression (Labeled by Diagnosis Age)")

    # Annotate each point with its Diagnosis Age.
    for i, age in enumerate(df['Diagnosis.Age']):
        plt.annotate(str(age), (pcs[i, 0], pcs[i, 1]), fontsize=10, alpha=0.8)
    plt.savefig("/home/mkh062/Desktop/scratch/TCGA_project/processed_data/results/PCA_glass_age.png", dpi=300)
    plt.show()


# ---------------------------------------------------------------------------
# 2. K-means Clustering and Elbow Plot
# ---------------------------------------------------------------------------


def elbow_plot_kmeans(filepath, max_clusters=10):
    """
    Computes the within-cluster sum of squares (WCSS) for a range of cluster numbers,
    calculates the reduction in variation when adding an extra cluster,
    and plots this reduction against the number of clusters.

    Parameters:
      filepath: Path to the GLASS_processed.csv file.
      max_clusters: Maximum number of clusters to evaluate (default is 10).
    """
    df, gene_cols = load_glass_data(filepath)
    X = df[gene_cols].values
    inertias = []
    cluster_range = range(1, max_clusters + 1)

    # Compute the inertia for each k.
    for k in cluster_range:
        kmeans = KMeans(n_clusters=k, random_state=42)
        kmeans.fit(X)
        inertias.append(kmeans.inertia_)

    # Calculate reduction differences between successive k's.
    # For k=2,...,max_clusters, reduction = inertia(k-1) - inertia(k)
    reductions = [inertias[i - 1] - inertias[i] for i in range(1, len(inertias))]

    plt.figure(figsize=(8, 6))
    plt.plot(range(2, max_clusters + 1), reductions, marker='o')
    plt.xlabel("Number of Clusters (k)", fontsize=12)
    plt.ylabel("Reduction in WCSS", fontsize=12)
    plt.title("Reduction in Variation vs. Number of Clusters", fontsize=12)
    plt.xticks(range(2, max_clusters + 1))
    plt.savefig("/home/mkh062/Desktop/scratch/TCGA_project/processed_data/results/elbow,plot.png", dpi=300)
    plt.show()

def kmeans_clustering(filepath, n_clusters=4):
    """
    Applies k-means clustering to the gene expression data and visualizes the clusters in a PCA plot.

    Parameters:
      filepath: Path to the GLASS_processed.csv file.
      n_clusters: The number of clusters for k-means (default is 4).

    Returns:
      df: DataFrame with an added 'Cluster' column.
      kmeans: The fitted KMeans model.
    """
    df, gene_cols = load_glass_data(filepath)
    X = df[gene_cols].values

    # Perform k-means clustering.
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    cluster_labels = kmeans.fit_predict(X)

    # For visualization, perform PCA on the gene expression data.
    pca = PCA(n_components=2)
    pcs = pca.fit_transform(X)

    plt.figure(figsize=(8, 6))
    scatter = plt.scatter(pcs[:, 0], pcs[:, 1], c=cluster_labels, cmap='viridis', alpha=0.7)
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.title(f"PCA of GLASS Gene Expression with K-means Clustering (k={n_clusters})")
    plt.legend(*scatter.legend_elements(), title="Cluster")
    plt.savefig("/home/mkh062/Desktop/scratch/TCGA_project/processed_data/results/PCA_kmeans_" + str({n_clusters}) + ".png", dpi=300)
    plt.show()

    # Add cluster labels to the dataframe.
    df['Cluster'] = cluster_labels
    return df, kmeans

#
# # Example usage:
# filepath = "/home/mkh062/Desktop/scratch/TCGA_project/processed_data/GLASS_processed.csv"
#
# # 1. PCA colored by Overall Survival groups.
# pca_by_survival(filepath)
#
# # 2. PCA with points labeled by Diagnosis Age.
# pca_by_age(filepath)
#
# # 3. K-means clustering (if PCA does not reveal clear clusters).
# # You can adjust n_clusters as needed.
# labels, centers = kmeans_clustering(filepath, n_clusters=4)

import seaborn as sns
from sklearn.decomposition import PCA
from scipy import stats


def overlay_clinical_data(df, cluster_col='Cluster'):
    """
    For each cluster, visualize clinical variable distributions.

    - Boxplots are drawn for numeric variables (Overall.Survival.Months and Diagnosis.Age).
    - Bar plots are drawn for categorical variables (Sex and case.vital.status).

    Parameters:
      df : pandas.DataFrame
           DataFrame containing both clinical data and a column with cluster assignments.
      cluster_col : str
           Column name for cluster labels.
    """
    # Boxplot: Overall Survival Months by Cluster
    plt.figure(figsize=(10, 6))
    sns.boxplot(x=cluster_col, y='Overall.Survival.Months', data=df)
    plt.title('Overall Survival Months by Cluster')
    plt.savefig("/home/mkh062/Desktop/scratch/TCGA_project/processed_data/results/kmeans_survival.png", dpi=300)
    plt.show()

    # Boxplot: Diagnosis Age by Cluster
    plt.figure(figsize=(10, 6))
    sns.boxplot(x=cluster_col, y='Diagnosis.Age', data=df)
    plt.title('Diagnosis Age by Cluster')
    plt.savefig("/home/mkh062/Desktop/scratch/TCGA_project/processed_data/results/kmeans_age.png", dpi=300)
    plt.show()

    # Bar plot: Sex distribution by Cluster
    plt.figure(figsize=(10, 6))
    sex_counts = pd.crosstab(df[cluster_col], df['Sex'])
    sex_counts.plot(kind='bar', stacked=True)
    plt.title('Sex Distribution by Cluster')
    plt.xlabel('Cluster')
    plt.ylabel('Count')
    plt.savefig("/home/mkh062/Desktop/scratch/TCGA_project/processed_data/results/kmeans_sex.png", dpi=300)
    plt.show()

    # Bar plot: Case Vital Status distribution by Cluster
    plt.figure(figsize=(10, 6))
    status_counts = pd.crosstab(df[cluster_col], df['case.vital.status'])
    status_counts.plot(kind='bar', stacked=True)
    plt.title('Case Vital Status by Cluster')
    plt.xlabel('Cluster')
    plt.ylabel('Count')
    plt.savefig("/home/mkh062/Desktop/scratch/TCGA_project/processed_data/results/kmeans_vital.png", dpi=300)
    plt.show()


def examine_cluster_centroids(kmeans_model, gene_columns):
    """
    Examine the centroids of each cluster (i.e., the average gene expression profile for each cluster).

    Parameters:
      kmeans_model : sklearn.cluster.KMeans
           A fitted KMeans model.
      gene_columns : list of str
           List of column names corresponding to gene expression data.

    Returns:
      centroids_df : pandas.DataFrame
           DataFrame with clusters as rows and gene expression values (centroid values) as columns.
      diff_genes : pandas.Series
           Genes sorted by their standard deviation across centroids, highlighting those
           that differ the most between clusters.
    """
    centroids = kmeans_model.cluster_centers_
    centroids_df = pd.DataFrame(centroids, columns=gene_columns)

    # Identify genes with highest variability across centroids.
    gene_std = centroids_df.std(axis=0)
    diff_genes = gene_std.sort_values(ascending=False)

    print("Top 10 genes with highest variation across centroids:")
    print(diff_genes.head(10))

    return centroids_df, diff_genes


def visualize_clusters(df, gene_columns, cluster_col='Cluster', top_n_genes=50):
    """
    Visualize clusters using PCA and a heatmap.

    - A PCA scatter plot is generated, with points colored by their cluster label.
    - A heatmap of gene expression is shown, with samples ordered by cluster assignment.

    Parameters:
      df : pandas.DataFrame
           DataFrame containing gene expression data and cluster labels.
      gene_columns : list of str
           List of column names for gene expression.
      cluster_col : str
           Column name that contains the cluster assignments.
    """
    X = df[gene_columns].values

    # Perform PCA for 2D visualization.
    pca = PCA(n_components=2)
    pcs = pca.fit_transform(X)

    # PCA scatter plot colored by cluster.
    plt.figure(figsize=(8, 6))
    scatter = plt.scatter(pcs[:, 0], pcs[:, 1], c=df[cluster_col], cmap='viridis', alpha=0.7)
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.title("PCA of Gene Expression Data Colored by Cluster")
    plt.legend(*scatter.legend_elements(), title="Cluster")
    plt.show()

    # # Heatmap: Order samples by cluster for visualizing gene expression patterns.
    # ordered_df = df.sort_values(by=cluster_col)
    # plt.figure(figsize=(12, 8))
    # sns.heatmap(ordered_df[gene_columns], cmap="viridis", yticklabels=False)
    # plt.title("Heatmap of Gene Expression (Ordered by Cluster)")
    # plt.xlabel("Genes")
    # plt.ylabel("Samples (ordered by cluster)")
    # plt.show()
    # Compute variance for each gene.
    gene_variances = df[gene_columns].var().sort_values(ascending=False)
    top_genes = gene_variances.head(top_n_genes).index.tolist()

    # Order samples by cluster for clarity.
    ordered_df = df.sort_values(by=cluster_col)

    plt.figure(figsize=(12, 8))
    sns.heatmap(ordered_df[top_genes], cmap="viridis", yticklabels=False)
    plt.title(f"Heatmap of Top {top_n_genes} Variable Genes (Ordered by Cluster)")
    plt.xlabel("Genes")
    plt.ylabel("Samples (ordered by cluster)")
    plt.show()


def statistical_validation(df, cluster_col='Cluster'):
    """
    Perform statistical tests to compare clinical variables across clusters.

    - For numeric variables (Overall.Survival.Months, Diagnosis.Age): One-way ANOVA is performed.
    - For categorical variables (Sex, case.vital.status): Chi-square tests are performed.

    Parameters:
      df : pandas.DataFrame
           DataFrame containing clinical variables and cluster assignments.
      cluster_col : str
           Column name for the cluster labels.

    Returns:
      results : dict
           Dictionary of p-values for each test.
    """
    results = {}

    # ANOVA for Overall.Survival.Months.
    survival_groups = [group['Overall.Survival.Months'].dropna().values for name, group in df.groupby(cluster_col)]
    f_val, p_val = stats.f_oneway(*survival_groups)
    results['Overall.Survival.Months_ANOVA'] = p_val
    print(f"ANOVA for Overall.Survival.Months: p-value = {p_val}")

    # ANOVA for Diagnosis.Age.
    age_groups = [group['Diagnosis.Age'].dropna().values for name, group in df.groupby(cluster_col)]
    f_val, p_val = stats.f_oneway(*age_groups)
    results['Diagnosis.Age_ANOVA'] = p_val
    print(f"ANOVA for Diagnosis.Age: p-value = {p_val}")

    # Chi-square test for Sex distribution.
    contingency_sex = pd.crosstab(df[cluster_col], df['Sex'])
    chi2, p_sex, dof, ex = stats.chi2_contingency(contingency_sex)
    results['Sex_Chi2'] = p_sex
    print(f"Chi-square test for Sex distribution: p-value = {p_sex}")

    # Chi-square test for case.vital.status distribution.
    contingency_status = pd.crosstab(df[cluster_col], df['case.vital.status'])
    chi2, p_status, dof, ex = stats.chi2_contingency(contingency_status)
    results['case.vital.status_Chi2'] = p_status
    print(f"Chi-square test for case.vital.status: p-value = {p_status}")

    return results


# ---------------------------------------------------------------------------
# 4. Main Pipeline: Integration and Usage
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    filepath = "/home/mkh062/Desktop/scratch/TCGA_project/processed_data/GLASS_processed.csv"

    # 1. Visualize PCA colored by Overall Survival groups.
    pca_by_survival(filepath)

    # 2. Visualize PCA with Diagnosis Age annotations.
    pca_by_age(filepath)

    # 3. Compute and display the Elbow Plot for determining optimal k.
    elbow_plot_kmeans(filepath, max_clusters=10)

    # 4. Apply k-means clustering with your chosen number of clusters (e.g., 4).
    df_clustered, kmeans_model = kmeans_clustering(filepath, n_clusters=5)

    # Retrieve gene expression columns.
    _, gene_cols = load_glass_data(filepath)

    # 5. Overlay clinical data to assess differences across clusters.
    overlay_clinical_data(df_clustered, cluster_col='Cluster')

    # 6. Examine the cluster centroids to identify driving genes.
    centroids_df, diff_genes = examine_cluster_centroids(kmeans_model, gene_cols)

    # 7. Visualize clusters using PCA and a heatmap of the top variable genes.
    visualize_clusters(df_clustered, gene_cols, cluster_col='Cluster')

    # 8. Statistically validate if clinical variables differ between clusters.
    results = statistical_validation(df_clustered, cluster_col='Cluster')