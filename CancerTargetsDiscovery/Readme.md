<h1>Cancer Target Discovery from Understudied (Tdark) Genes</h1>

<p>
Systematic computational exploration of poorly characterized human protein-coding genes to evaluate their potential relevance in cancer using bulk RNA-seq, survival analysis, and integrated single-cell RNA-seq datasets.
</p>

<hr>

<h2>Project Overview</h2>

<p>
Despite rapid advances in cancer genomics, many human protein-coding genes remain poorly characterized. Many are classified as <b>Tdark</b> targets under the 
<a href="https://commonfund.nih.gov/idg">Illuminating the Druggable Genome (IDG)</a> initiative.
</p>

<p>This project evaluates whether selected understudied genes show:</p>

<ul>
  <li>Tumor-specific expression changes</li>
  <li>Clinical survival associations</li>
  <li>RNA–protein concordance or discordance</li>
  <li>Reproducible expression across independent single-cell datasets</li>
  <li>Scalability across multiple scRNA-seq cohorts</li>
</ul>

<hr>

<h2>Data Sources</h2>

<ul>
  <li>TCGA + GTEx via <a href="http://gepia2.cancer-pku.cn">GEPIA2</a></li>
  <li><a href="https://www.proteinatlas.org">Human Protein Atlas (HPA)</a></li>
  <li><a href="https://pharos.nih.gov">Pharos (NIH IDG)</a></li>
  <li>Single-cell RNA-seq datasets:
    <ul>
      <li>GSE72056 (Melanoma, ~4,645 cells, 19 patients)</li>
      <li>GSE115978 (Melanoma, ~7,186 cells, 31 patients)</li>
    </ul>
  </li>
</ul>

<hr>

<h2>Candidate Gene Panel</h2>

<pre>
CASKIN2
EMC9
PDIK1L
DBNDD2
FAM171A2
C1orf174
LOC124903857
TMEM161B
ZNF808
</pre>

<p><b>Selection Criteria:</b></p>
<ul>
  <li>Tdark classification (minimal literature)</li>
  <li>Low UniProt annotation score</li>
  <li>Sparse Gene Ontology annotation</li>
  <li>Limited or unknown protein domains</li>
  <li>Detectable expression in TCGA/GTEx</li>
</ul>

<hr>

<h2>Methods</h2>

<h3>1. Pan-Cancer Bulk RNA Analysis</h3>

<ul>
  <li>Tumor vs normal comparisons (TCGA + GTEx)</li>
  <li>Boxplot visualization</li>
  <li>Kaplan–Meier survival analysis</li>
  <li>Hazard ratio and log-rank interpretation</li>
</ul>

<p><b>Key Observations:</b></p>
<ul>
  <li><b>EMC9</b>: RNA upregulated in DLBC and THYM but weak protein staining (RNA–protein mismatch).</li>
  <li><b>PDIK1L</b>: Significant survival association in kidney cancer (KIRC).</li>
  <li><b>DBNDD2</b>: Strong pan-cancer survival association; high expression linked to poorer outcomes.</li>
</ul>

<hr>

<h3>2. Single-Cell RNA-seq Integration (Seurat)</h3>

<p><b>Workflow:</b></p>
<ol>
  <li>Create Seurat objects</li>
  <li>Normalize and identify variable features</li>
  <li>Anchor-based dataset integration</li>
  <li>PCA + UMAP dimensionality reduction</li>
  <li>Clustering</li>
  <li>Marker gene identification</li>
  <li>Manual cluster-to-cell-type mapping</li>
</ol>

<p><b>Integration Methods Used:</b></p>
<ul>
  <li>FindIntegrationAnchors()</li>
  <li>IntegrateData()</li>
  <li>FindAllMarkers()</li>
  <li>AverageExpression()</li>
</ul>

<hr>

<h2>Reproducibility Assessment</h2>

<ul>
  <li>Average gene expression per dataset</li>
  <li>Heatmap visualization of candidate genes</li>
  <li>UMAP visualization per dataset</li>
  <li>Cluster-level marker comparison</li>
</ul>

<p><b>Findings:</b></p>
<ul>
  <li>Consistent directional trends across both datasets</li>
  <li>Expected magnitude differences due to sequencing depth</li>
  <li>Low-expression genes detectable but near threshold</li>
  <li>DBNDD2 shows strong neural/tumor retention</li>
</ul>

<hr>

<h2>Scalable Multi-Dataset Framework</h2>

<p>The workflow is designed to scale beyond two datasets.</p>

<ul>
  <li>Batch-wise preprocessing</li>
  <li>Memory-efficient anchor integration</li>
  <li>Automated marker detection</li>
  <li>Cluster-level summarization</li>
  <li>Pseudo-bulk averaging</li>
  <li>Automated heatmap generation</li>
</ul>

<p>This enables extension to 10+ scRNA-seq datasets.</p>

<hr>

<h2>Repository Structure</h2>

<pre>
data/
  ├── raw/
  └── processed/

results/
  ├── figures/
  ├── tables/
  └── objects/

scripts/
  ├── config.R
  ├── data_preprocess.R
  ├── Clustering_integration.R
  ├── Clustering_individual.R
  ├── Marker_Analysis.R
  └── Visualisation.R
</pre>

<hr>

<h2>How to Run</h2>

<h3>1. Install Dependencies</h3>

<pre>
install.packages(c("Seurat", "dplyr", "ggplot2", "pheatmap"))
</pre>

<h3>2. Place Raw Data</h3>

<pre>
data/raw/
</pre>

<h3>3. Run Pipeline</h3>

<pre>
source("scripts/config.R")
source("scripts/data_preprocess.R")
source("scripts/Clustering_integration.R")
source("scripts/Marker_Analysis.R")
source("scripts/Visualisation.R")
</pre>

<p>All outputs are automatically saved in:</p>

<pre>
results/
</pre>

<hr>

<h2>Key Takeaways</h2>

<ul>
  <li>Understudied genes can show clinically meaningful cancer associations.</li>
  <li>RNA–protein mismatches are common for low-annotation targets.</li>
  <li>DBNDD2 emerges as a potential pan-cancer risk marker.</li>
  <li>PDIK1L shows kidney cancer survival relevance.</li>
  <li>Anchor-based integration enables reproducible cross-dataset comparison.</li>
</ul>

<hr>

<h2>Limitations</h2>

<ul>
  <li>Protein validation limited by antibody availability.</li>
  <li>Survival analysis exploratory (no multivariate modeling).</li>
  <li>Manual cluster annotation introduces subjectivity.</li>
  <li>Functional mechanisms remain experimentally unvalidated.</li>
</ul>

<hr>

<h2>Future Directions</h2>

<ul>
  <li>Multi-cancer scRNA-seq expansion</li>
  <li>Spatial transcriptomics validation</li>
  <li>Pseudotime / trajectory analysis</li>
  <li>Functional CRISPR validation</li>
  <li>Proteomics integration</li>
</ul>

<hr>

<h2>Author</h2>

<p>
<b>Kedar Joshi</b><br>
M.S. Bioinformatics
</p>