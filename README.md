# Coral Transcriptome Differential Expression Analysis Suite

A comprehensive computational pipeline for simulating and analyzing coral reef transcriptomic responses to environmental stressors, implementing standardized RNA-seq differential expression workflows with interactive visualization capabilities.

## Overview

A bioinformatics tool designed for marine genomics researchers studying cnidarian transcriptional responses to environmental perturbations. The application provides a complete analytical framework for coral gene expression profiling, incorporating established methodologies from the field of marine molecular ecology and transcriptomics.

### Key Features

- **Transcriptome Simulation Engine**: Generates realistic coral gene expression count matrices with configurable parameters mimicking RNA-seq data characteristics
- **Differential Expression Analysis**: Implements Welch's t-test with Benjamini-Hochberg FDR correction for identifying stress-responsive genes
- **CPM Normalization Pipeline**: Counts per million (CPM) transformation with log₂ scaling for variance stabilization
- **Principal Component Analysis**: Dimensionality reduction for sample clustering and batch effect detection
- **Interactive Volcano Plots**: Statistical significance visualization with log₂ fold-change thresholds
- **Expression Heatmaps**: Hierarchical clustering of top differentially expressed genes (DEGs)

## Scientific Background

Coral reef ecosystems face unprecedented environmental challenges including ocean acidification, thermal stress, and pollution exposure. Understanding the molecular mechanisms underlying coral stress responses requires comprehensive transcriptomic profiling to identify key regulatory pathways and biomarker genes.

This tool facilitates:
- **Stress Response Genomics**: Analysis of differential gene expression between control and stress-exposed coral samples
- **Biomarker Discovery**: Identification of candidate genes for environmental monitoring applications
- **Pathway Enrichment Preparation**: Generation of gene lists suitable for downstream functional annotation
- **Quality Control Assessment**: Sample clustering analysis for experimental validation

## Methodology

### Data Simulation Framework

The simulation engine generates synthetic coral transcriptome data using Poisson-distributed count models:

- **Baseline Expression**: Control samples drawn from Poisson(λ=50) distribution
- **Stress Response Modeling**: 
  - Upregulated genes (n=10): Additional Poisson(λ=40) counts
  - Downregulated genes (n=10): Reduced by Poisson(λ=20) counts
  - Background genes: Maintain baseline expression levels

### Statistical Analysis Pipeline

1. **Count Matrix Loading**: Import expression data in standard CSV format with genes as rows and samples as columns
2. **CPM Normalization**: 
   ```
   CPM = (raw_counts / library_size) × 10⁶
   log₂CPM = log₂(CPM + 1)
   ```
3. **Differential Expression Testing**:
   - Statistical method: Welch's two-sample t-test (unequal variances)
   - Multiple testing correction: Benjamini-Hochberg FDR (α = 0.05)
   - Effect size: log₂ fold-change calculation with pseudocount adjustment

4. **Dimensionality Reduction**:
   - Principal Component Analysis on log₂CPM-transformed data
   - Variance explained reporting for PC1 and PC2

### Visualization Components

- **Volcano Plot**: -log₁₀(p-value) vs log₂(fold-change) with significance thresholding
- **PCA Biplot**: Sample ordination colored by treatment group with variance contribution
- **Expression Heatmap**: Z-score normalized expression of top 25 DEGs with hierarchical clustering

## Installation & Dependencies

### System Requirements
- Python ≥ 3.8
- 4GB RAM minimum (recommended: 8GB for large datasets)

### Package Dependencies
```bash
pip install -r requirements.txt
```

**Core Libraries:**
- `streamlit`: Web application framework
- `pandas`: Data manipulation and analysis
- `numpy`: Numerical computing
- `scipy`: Statistical functions
- `statsmodels`: Advanced statistical modeling
- `scikit-learn`: Machine learning and PCA implementation
- `plotly`: Interactive plotting
- `seaborn`: Statistical data visualization
- `matplotlib`: Publication-quality figures

## Usage

### Web Application Launch
```bash
streamlit run app.py
```

The application will launch in your default web browser at `http://localhost:8501`.

### Data Input Options

1. **Simulated Data Generation** (Default):
   - Toggle "Use existing data" checkbox OFF
   - Generates fresh mock expression matrix (200 genes × 6 samples)
   - Sample naming convention: `Ctrl_1`, `Ctrl_2`, `Ctrl_3`, `Stress_1`, `Stress_2`, `Stress_3`

2. **Custom Data Upload**:
   - Toggle "Use existing data" checkbox ON
   - Requires CSV format with gene identifiers as row names
   - Column headers must follow `Ctrl_*` and `Stress_*` naming pattern

### Data Format Specifications

**Input CSV Structure:**
```
gene_id,Ctrl_1,Ctrl_2,Ctrl_3,Stress_1,Stress_2,Stress_3
Gene_001,45,52,48,89,76,84
Gene_002,67,61,59,45,52,48
...
```

**Requirements:**
- Gene identifiers in first column
- Sample names containing "Ctrl" for control samples
- Sample names containing "Stress" for treatment samples
- Raw count data (non-negative integers)

## Interpretation Guidelines

### Volcano Plot Analysis
- **X-axis**: log₂ fold-change (positive = upregulated in stress, negative = downregulated)
- **Y-axis**: -log₁₀(p-value) (higher = more statistically significant)
- **Color coding**: 
  - Red/highlighted points: FDR-adjusted p-value < 0.05
  - Gray points: Non-significant genes

### PCA Interpretation
- **Sample clustering**: Groups of samples with similar expression profiles
- **Treatment separation**: Distance between control and stress sample clusters
- **Variance explained**: Percentage of total expression variance captured by each PC
- **Quality assessment**: Tight within-group clustering indicates good biological replication

### Heatmap Analysis
- **Rows**: Top 25 most significantly differentially expressed genes
- **Columns**: Individual sample replicates
- **Color scale**: Z-score normalized expression (red = high, blue = low)
- **Clustering**: Hierarchical clustering reveals gene co-expression patterns

### Statistical Significance Thresholds
- **p-value**: Raw statistical significance from t-test
- **FDR-adjusted p-value**: Multiple testing corrected significance (α = 0.05)
- **Fold-change**: Biological effect size (typical threshold: |log₂FC| > 1)

## Output Data

The application generates several data products:

1. **Differential Expression Results Table**:
   - Gene identifiers
   - log₂ fold-change values
   - Raw p-values
   - FDR-adjusted p-values
   - Ranked by statistical significance

2. **Interactive Visualizations**:
   - Volcano plot (HTML/SVG export compatible)
   - PCA biplot with variance statistics
   - Expression heatmap with clustering dendrograms

## Biological Applications

### Marine Stress Biology
- **Thermal Stress Response**: Identification of heat shock proteins and chaperone genes
- **Ocean Acidification Studies**: Analysis of calcification and pH homeostasis genes
- **Pollution Exposure**: Detection of xenobiotic metabolism and detoxification pathways

### Comparative Transcriptomics
- **Species Comparisons**: Cross-species analysis of conserved stress responses
- **Developmental Studies**: Temporal expression profiling during coral development
- **Symbiosis Research**: Host-symbiont interaction transcriptomics

### Biomarker Development
- **Environmental Monitoring**: Selection of indicator genes for ecosystem health assessment
- **Diagnostic Applications**: Development of molecular diagnostic assays
- **Conservation Genomics**: Population-level transcriptomic monitoring

## Technical Specifications

### Algorithm Complexity
- **Time Complexity**: O(n×m) where n = genes, m = samples
- **Space Complexity**: O(n×m) for count matrix storage
- **Scalability**: Optimized for datasets up to 50K genes × 1K samples

### Statistical Assumptions
- **Count Distribution**: Assumes Poisson-distributed read counts
- **Independence**: Requires independent biological replicates
- **Normality**: log₂CPM transformation improves normality for parametric tests

## Limitations & Considerations

1. **Simplified Model**: Uses t-test instead of more sophisticated methods (DESeq2, edgeR)
2. **Single Factor Design**: Currently supports only two-group comparisons
3. **No Batch Correction**: Does not account for technical batch effects
4. **Pseudocount Addition**: May introduce bias for very low-count genes

## Future Enhancements

- Integration with established RNA-seq pipelines (DESeq2, limma-voom)
- Multi-factor experimental design support
- Gene set enrichment analysis (GSEA) functionality
- Pathway visualization and annotation
- Batch effect correction methods
- Support for single-cell RNA-seq data

## Citation

If you use ReefGene-Exp-Sim in your research, please cite:

```
ReefGene-Exp-Sim: A Streamlit-based tool for coral transcriptome differential expression analysis
GitHub: https://github.com/samvictordr/reefgene-exp-sim
```

## License

This project is released under the Creative Commons CC0 1.0 Universal License - see the [LICENSE](LICENSE) file for details.

