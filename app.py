import streamlit as st
import numpy as np
import pandas as pd
import plotly.express as px
import seaborn as sns
import matplotlib.pyplot as plt
from pymongo import MongoClient
from utils.generate_data import simulate_coral_expression
from utils.analysis import (
    normalize_cpm,
    differential_expression,
    run_pca,
    top_genes_for_heatmap
)

# Constants
MONGO_URI = "mongodb+srv://shinythewitch:SVOAOs1a4PbUsD8h@cluster0.yz3qbn2.mongodb.net/reefDE?retryWrites=true&w=majority"
DB_NAME = "reefDE"
COLLECTION_NAME = "DE-uploads"

st.set_page_config(page_title="Coral Gene Expr. Analyzer", layout="wide")
st.title("🌊 Coral Reef Gene Expression Simulator")

# --- Data Input Section ---
st.sidebar.header("Data Input")
# DELETE/CLEAR: Clear session data and optionally rerun
if st.sidebar.button("🗑️ Clear Data"):
    for key in list(st.session_state.keys()):
        del st.session_state[key]
    if hasattr(st, 'experimental_rerun'):
        st.experimental_rerun()
    else:
        st.sidebar.info("Please refresh the page to complete data clear.")
        st.stop()

# CREATE & READ Operations
uploaded_file = st.sidebar.file_uploader("Upload counts CSV (genes × samples)", type="csv")
if uploaded_file is not None:
    try:
        counts = pd.read_csv(uploaded_file, index_col=0)
        st.sidebar.success(f"Loaded uploaded data: {counts.shape[0]} genes × {counts.shape[1]} samples")
    except Exception as e:
        st.sidebar.error(f"Error loading file: {e}")
        st.stop()
else:
    use_existing = st.sidebar.checkbox("Use existing mock data", value=False)
    if not use_existing:
        st.sidebar.write("🔄 Simulating mock counts and creating file…")
        data_path = simulate_coral_expression(num_genes=200, num_samples=6)
    else:
        data_path = "data/mock_expression_data.csv"
    counts = pd.read_csv(data_path, index_col=0)
    st.sidebar.success(f"Loaded counts: {counts.shape[0]} genes × {counts.shape[1]} samples")

# --- Analysis Pipeline (UPDATE) ---
groups = ["Control" if col.lower().startswith("ctrl") else "Stress" for col in counts.columns]
norm = normalize_cpm(counts)
st.sidebar.write("✅ Normalization: CPM & log2")
de = differential_expression(counts, groups)
st.sidebar.write(f"✅ DE done ({(de.padj < 0.05).sum()} significant genes)")

# DOWNLOAD (CREATE): Export DE results
de_csv = de.reset_index().to_csv(index=False).encode("utf-8")
st.download_button(
    label="📥 Download DE results as CSV",
    data=de_csv,
    file_name="differential_expression_results.csv",
    mime="text/csv"
)

# MONGO UPLOAD: Insert DE results into MongoDB
st.sidebar.header("Database Operations")
if st.sidebar.button("🚀 Upload DE to MongoDB"):
    try:
        client = MongoClient(MONGO_URI)
        db = client[DB_NAME]
        coll = db[COLLECTION_NAME]
        records = de.reset_index().to_dict("records")
        # add timestamp metadata
        timestamp = pd.Timestamp.now().isoformat()
        for rec in records:
            rec["uploaded_at"] = timestamp
        result = coll.insert_many(records)
        st.sidebar.success(f"Inserted {len(result.inserted_ids)} records into MongoDB.")
    except Exception as e:
        st.sidebar.error(f"Upload failed: {e}")

# --- Visualization Layout (READ) ---
col1, col2 = st.columns(2)
with col1:
    st.subheader("📊 Volcano Plot")
    fig_volcano = px.scatter(
        de.reset_index(),
        x="log2FoldChange",
        y=-np.log10(de.pvalue),
        hover_name="gene",
        color=de.padj < 0.05
    )
    st.plotly_chart(fig_volcano, use_container_width=True)
with col2:
    st.subheader("🎯 PCA")
    pca_df, var_ratio = run_pca(norm, groups)
    fig_pca = px.scatter(
        pca_df,
        x="PC1",
        y="PC2",
        color="group",
        labels={'PC1': f"PC1 ({var_ratio[0]*100:.1f}% var)", 'PC2': f"PC2 ({var_ratio[1]*100:.1f}% var)"}
    )
    st.plotly_chart(fig_pca, use_container_width=True)

st.subheader("🔥 Heatmap of Top Differentially Expressed Genes")
top_genes = top_genes_for_heatmap(de, top_n=25)
heatmap_data = norm.loc[top_genes]
fig_heatmap, ax = plt.subplots(figsize=(10, 8))
sns.heatmap(heatmap_data, cmap="vlag", xticklabels=True, yticklabels=True)
st.pyplot(fig_heatmap)

st.markdown("### 📝 Top 10 Genes")
st.dataframe(
    de.head(10)
      .reset_index()
      .rename(columns={"gene":"Gene", "log2FoldChange":"log2FC", "pvalue":"p-value", "padj":"adj p-value"})
)
