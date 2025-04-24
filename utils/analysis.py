import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
import statsmodels.stats.multitest as smm
from sklearn.decomposition import PCA

def load_counts(path="data/mock_expression_data.csv"):
    return pd.read_csv(path, index_col=0)

def normalize_cpm(counts):
    counts = counts.div(counts.sum(axis=0), axis=1) * 1e6
    return np.log2(counts + 1)

def differential_expression(counts, groups):
    # groups: list of 'Control'/'Stress' matching columns
    grp = pd.Series(groups, index=counts.columns)
    ctrl = counts.loc[:, grp=='Control']
    stress = counts.loc[:, grp=='Stress']
    results = []
    for gene in counts.index:
        stat, p = ttest_ind(stress.loc[gene], ctrl.loc[gene], equal_var=False)
        # compute log2 fold change
        lfc = np.log2(stress.loc[gene].mean() + 1) - np.log2(ctrl.loc[gene].mean() + 1)
        results.append((gene, lfc, p))
    df = pd.DataFrame(results, columns=["gene","log2FoldChange","pvalue"]).set_index("gene")
    df["padj"] = smm.multipletests(df.pvalue, method="fdr_bh")[1]
    return df.sort_values("padj")

def run_pca(norm_counts, groups, n_components=2):
    pca = PCA(n_components=n_components)
    pcs = pca.fit_transform(norm_counts.T)
    df = pd.DataFrame(pcs, columns=[f"PC{i+1}" for i in range(n_components)], index=norm_counts.columns)
    df["group"] = groups
    return df, pca.explained_variance_ratio_

def top_genes_for_heatmap(de_df, top_n=20):
    return list(de_df.sort_values("padj").head(top_n).index)
