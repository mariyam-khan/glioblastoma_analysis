import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from lifelines import CoxPHFitter, KaplanMeierFitter
from statsmodels.stats.multitest import multipletests
import os
# R-related imports
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
pandas2ri.activate()


def limma(data, threshold, consortium, threshold_type="fixed", outdir="./results"):
    """
    Runs a differential expression analysis in R using limma.
    """
    os.makedirs(outdir, exist_ok=True)
    data = data.copy()

    if threshold_type == "fixed":
        data = data[~((data['Overall.Survival.Months'] <= threshold) & (data['case.vital.status'] == 0))]
        file_suffix = f"{consortium}_{threshold}-months"
        median_survival_r = threshold
    else:
        # range-based threshold
        low, high = threshold
        data = data[~((data['Overall.Survival.Months'] <= low) & (data['case.vital.status'] == 0))]
        data = data[(data['Overall.Survival.Months'] <= low) | (data['Overall.Survival.Months'] > high)]
        file_suffix = f"{consortium}_{low}-{high}-months"
        median_survival_r = low

    data['group'] = np.where(data['Overall.Survival.Months'] > median_survival_r, 'high', 'low')
    n_low = sum(data['group']=='low')
    n_high= sum(data['group']=='high')
    print(f"{consortium} {file_suffix} => #Low group = {n_low}, #High group = {n_high}")

    data.drop(columns=['case.vital.status'], inplace=True)

    if 'batch' in data.columns:
        design_formula = "~ group + Diagnosis.Age + Sex + batch"
        exclude_cols_r = "c('Sample.ID','Overall.Survival.Months','Diagnosis.Age','Sex','batch','group')"
    else:
        design_formula = "~ group + Diagnosis.Age + Sex"
        exclude_cols_r = "c('Sample.ID','Overall.Survival.Months','Diagnosis.Age','Sex','group')"

    r_data = pandas2ri.py2rpy(data)
    r_script = f"""
    library(limma)
    data <- as.data.frame(data)
    gene_expression_data <- data[, !names(data) %in% {exclude_cols_r}]
    rownames(gene_expression_data) <- data$Sample.ID
    data$group <- factor(data$group, levels = c('low','high'))

    design <- model.matrix({design_formula}, data = data)
    rownames(design) <- data$Sample.ID

    fit <- lmFit(t(gene_expression_data), design)
    fit <- eBayes(fit)
    results <- topTable(fit, coef='grouphigh', number=Inf, adjust.method='fdr')
    out_all <- file.path("{outdir}", paste0("all_genes_{file_suffix}.csv"))
    write.csv(results, file=out_all, row.names=TRUE)
    sig <- results[results$adj.P.Val<0.1, ]
    out_sig <- file.path("{outdir}", paste0("significant_genes_{file_suffix}.csv"))
    write.csv(sig, file=out_sig, row.names=TRUE)
    """

    ro.globalenv['data'] = r_data
    ro.r(r_script)

    results_df = pd.read_csv(os.path.join(outdir, f"significant_genes_{file_suffix}.csv"), index_col=0)
    return results_df
