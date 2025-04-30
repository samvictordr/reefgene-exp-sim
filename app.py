import streamlit as st
import numpy as np
import plotly.express as px
import seaborn as sns
import matplotlib.pyplot as plt
from utils.generate_data import simulate_coral_expression
from utils.analysis import load_counts, normalize_cpm, differential_expression, run_pca, top_genes_for_heatmap

st.set_page_config(page_title="Coral Gene Expr. Analyzer", layout="wide")

st.title("🌊 Coral Reef Gene Expression Simulator")
data_path = "data/mock_expression_data.csv"
if not st.sidebar.checkbox("Use existing data", value=False):
    st.sidebar.write("🔄 Generating mock data…")
    data_path = simulate_coral_expression(num_genes=200, num_samples=6)

counts = load_counts(data_path)
st.sidebar.success(f"Loaded counts: {counts.shape[0]} genes × {counts.shape[1]} samples")

# Define groups
groups = []
for col in counts.columns:
    groups.append("Control" if col.startswith("Ctrl") else "Stress")

# Normalization
norm = normalize_cpm(counts)
st.sidebar.write("✅ Normalization: CPM & log2")

# Differential Expression
de = differential_expression(counts, groups)
st.sidebar.write(f"✅ DE done ({(de.padj<0.05).sum()} significant)")

# Layout
col1, col2 = st.columns(2)

with col1:
    st.subheader("📊 Volcano Plot")
    fig = px.scatter(
        de.reset_index(),
        x="log2FoldChange", y=-np.log10(de.pvalue),
        hover_name="gene",
        color=de.padj < 0.05,
        labels={'color':'padj<0.05','x':'log2 Fold Change','y':'-log10(p-value)'}
    )
    st.plotly_chart(fig, use_container_width=True)

with col2:
    st.subheader("🎯 PCA")
    pca_df, var_ratio = run_pca(norm, groups)
    fig2 = px.scatter(
        pca_df, x="PC1", y="PC2", color="group",
        labels={'PC1':f"PC1 ({var_ratio[0]*100:.1f}%)",
                'PC2':f"PC2 ({var_ratio[1]*100:.1f}%)"}
    )
    st.plotly_chart(fig2, use_container_width=True)

st.subheader("Heatmap of Top Differentially Expressed Genes")
top_genes = top_genes_for_heatmap(de, top_n=25)
heatmap_data = norm.loc[top_genes]
fig3, ax = plt.subplots(figsize=(10,8))
sns.heatmap(heatmap_data, cmap="vlag", xticklabels=True, yticklabels=True)
st.pyplot(fig3)

st.markdown("### 📝 Top 10 Genes")
st.dataframe(de.head(10).reset_index().rename(columns={
    "gene":"Gene","log2FoldChange":"log2FC","pvalue":"p-value","padj":"adj p-value"
}))
