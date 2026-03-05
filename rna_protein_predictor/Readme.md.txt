<h1>RNA–Protein Discordance in Cancer</h1>

<h2>Predicting Protein Detectability from RNA Expression and Biological Context</h2>

<p>
RNA expression is commonly used as a proxy for protein abundance in cancer studies.
However, RNA levels do not always translate into detectable protein. Understanding
when RNA fails to predict protein presence is important for biomarker discovery,
target prioritization, and interpretation of transcriptomic data.
</p>

<p>
This project builds a machine learning pipeline to predict protein detectability
from RNA expression and then investigates biological cases where RNA is high but
protein is not detected.
</p>

<hr>

<h1>Data Sources</h1>

<h2>RNA Expression</h2>

<p>
RNA expression data was obtained from the TCGA + TARGET + GTEx RNA-seq recompute
dataset hosted on the UCSC Xena platform.
</p>

<ul>
<li>Dataset: <b>TcgaTargetGtex_rsem_gene_tpm</b></li>
<li>Samples: ~19,000</li>
<li>Used for tumor RNA expression, normal RNA expression, and fold-change calculations</li>
</ul>

<h2>Protein Expression</h2>

<p>
Protein detection data was obtained from the <b>Human Protein Atlas (HPA)</b>
immunohistochemistry dataset.
</p>

<p>Input format:</p>

<pre>
Gene | Cancer | High | Medium | Low | Not detected
</pre>

<p>
These values were converted into a binary classification target:
</p>

<pre>
protein_detected = 1 if (High + Medium + Low) > 0
protein_detected = 0 if Not detected
</pre>

<h2>Subcellular Localization</h2>

<p>
Subcellular localization data was also obtained from the Human Protein Atlas
Subcellular Atlas.
</p>

<p>Example localizations include:</p>

<ul>
<li>Nucleoplasm</li>
<li>Cytosol</li>
<li>Mitochondria</li>
<li>Plasma membrane</li>
<li>Vesicles</li>
<li>Endoplasmic reticulum</li>
</ul>

<p>
These were converted into binary localization features for machine learning.
</p>

<hr>

<h1>Cancers Included</h1>

<p>
Seven cancer types were used where RNA expression and HPA protein data overlapped.
</p>

<table>
<tr>
<th>TCGA Code</th>
<th>Cancer</th>
</tr>
<tr><td>BRCA</td><td>Breast Cancer</td></tr>
<tr><td>COAD</td><td>Colon Cancer</td></tr>
<tr><td>GBM</td><td>Glioblastoma</td></tr>
<tr><td>LIHC</td><td>Liver Cancer</td></tr>
<tr><td>PAAD</td><td>Pancreatic Cancer</td></tr>
<tr><td>PRAD</td><td>Prostate Cancer</td></tr>
<tr><td>THCA</td><td>Thyroid Cancer</td></tr>
</table>

<hr>

<h1>Pipeline Overview</h1>

<pre>
RNA expression (TCGA/TARGET/GTEx)
        │
        ▼
Compute gene–cancer RNA statistics
        │
        ▼
Merge with HPA protein detectability
        │
        ▼
Build master dataset
        │
        ▼
Train machine learning models
        │
        ▼
Identify RNA–protein discordance
        │
        ▼
Analyze biological explanations
</pre>

<hr>

<h1>Feature Engineering</h1>

<h2>RNA Features</h2>

<ul>
<li>Mean tumor expression (log2 TPM)</li>
<li>Mean normal expression</li>
<li>RNA fold change (tumor vs normal)</li>
</ul>

<h2>Gene Features</h2>

<ul>
<li>Gene length</li>
<li>Protein coding indicator</li>
</ul>

<h2>Subcellular Localization Features</h2>

<p>Binary features representing protein localization:</p>

<pre>
loc_nucleoplasm
loc_cytosol
loc_mitochondria
loc_plasma_membrane
loc_vesicles
loc_golgi_apparatus
</pre>

<p>Total localization features: <b>27</b></p>

<hr>

<h1>Master Dataset</h1>

<p>The final dataset contains:</p>

<ul>
<li><b>Rows:</b> 106,932</li>
<li><b>Genes:</b> 15,209</li>
<li><b>Cancers:</b> 7</li>
</ul>

<p>Target distribution:</p>

<ul>
<li>Protein detected: 59,037</li>
<li>Protein not detected: 47,503</li>
</ul>

<p>
Each row represents a <b>gene–cancer pair</b>.
</p>

<hr>

<h1>Model Training</h1>

<h2>Train/Test Strategy</h2>

<p>
To evaluate generalization across cancer types, a group split was performed
by cancer.
</p>

<p>Held-out test cancers:</p>

<ul>
<li>GBM</li>
<li>PRAD</li>
</ul>

<hr>

<h1>Models Evaluated</h1>

<h2>Logistic Regression (RNA Only)</h2>

<p>Features:</p>

<ul>
<li>tumor_log2tpm_mean</li>
<li>rna_log2fc</li>
</ul>

<p>Performance:</p>

<ul>
<li>ROC-AUC: 0.70</li>
<li>PR-AUC: 0.66</li>
</ul>

<h2>Logistic Regression (RNA + Gene Features)</h2>

<p>Additional features:</p>

<ul>
<li>gene_length</li>
<li>protein_coding</li>
</ul>

<p>Performance:</p>

<ul>
<li>ROC-AUC: 0.704</li>
<li>PR-AUC: 0.661</li>
</ul>

<h2>Random Forest (RNA + Gene Features)</h2>

<p>Performance:</p>

<ul>
<li>ROC-AUC: 0.717</li>
<li>PR-AUC: 0.680</li>
</ul>

<h2>Random Forest + Subcellular Localization</h2>

<p>
Adding localization features significantly improved prediction.
</p>

<ul>
<li>ROC-AUC: <b>0.81</b></li>
<li>PR-AUC: <b>0.79</b></li>
</ul>

<hr>

<h1>Feature Importance</h1>

<p>Top predictors:</p>

<ul>
<li>tumor_log2tpm_mean</li>
<li>gene_length</li>
<li>rna_log2fc</li>
<li>loc_Nucleoplasm</li>
<li>loc_Mitochondria</li>
<li>loc_Cytosol</li>
<li>loc_Vesicles</li>
<li>loc_Plasma_membrane</li>
</ul>

<p>
Subcellular localization contributes significant predictive signal.
</p>

<hr>

<h1>Biological Pivot: RNA–Protein Discordance</h1>

<p>
After building predictive models, the project pivoted to analyze biological
cases where RNA is high but protein is not detected.
</p>

<p>Definition of discordance:</p>

<pre>
High RNA expression
AND
Protein detected = 0
</pre>

<hr>

<h1>Discordance Rate</h1>

<p>Overall discordance among high-RNA genes:</p>

<p><b>24.6%</b></p>

<p>
This suggests that approximately one in four highly expressed transcripts
does not produce detectable protein.
</p>

<hr>

<h1>Discordance by Cancer</h1>

<ul>
<li>GBM – 36.8%</li>
<li>LIHC – 28.7%</li>
<li>PAAD – 25.2%</li>
<li>PRAD – 23.8%</li>
<li>BRCA – 20.5%</li>
<li>THCA – 19.7%</li>
<li>COAD – 17.3%</li>
</ul>

<p>Glioblastoma showed the highest RNA–protein discordance.</p>

<hr>

<h1>Localization Enrichment</h1>

<p>Discordant genes were enriched in:</p>

<ul>
<li>Actin filaments</li>
<li>Primary cilium</li>
<li>Intermediate filaments</li>
<li>Focal adhesion sites</li>
<li>Plasma membrane</li>
<li>Vesicles</li>
</ul>

<p>Discordant genes were depleted in:</p>

<ul>
<li>Mitochondria</li>
<li>Nucleoplasm</li>
<li>Nuclear speckles</li>
</ul>

<hr>

<h1>Unexpected Missing Proteins</h1>

<p>
Using the trained model, genes predicted to produce detectable protein but
reported as absent in HPA were identified.
</p>

<p>Results:</p>

<ul>
<li>Total genes: <b>706</b></li>
<li>GBM: 497</li>
<li>PRAD: 209</li>
</ul>

<p>Example genes:</p>

<ul>
<li>ALDOA</li>
<li>COX6A1</li>
<li>TIMM17B</li>
<li>CDC42</li>
<li>BID</li>
<li>MRPL21</li>
</ul>

<p>
These genes may represent cases of post-transcriptional regulation,
rapid protein degradation, or experimental detection limitations.
</p>

<hr>

<h1>Key Findings</h1>

<ul>
<li>RNA expression moderately predicts protein detectability.</li>
<li>Gene structural features improve predictions slightly.</li>
<li>Subcellular localization significantly improves prediction accuracy.</li>
<li>~25% of highly expressed transcripts lack detectable protein.</li>
<li>Discordance varies across cancer types.</li>
</ul>

<hr>

<h1>Technologies Used</h1>

<ul>
<li>Python</li>
<li>Pandas</li>
<li>NumPy</li>
<li>Scikit-learn</li>
</ul>

<hr>

<h1>Future Work</h1>

<ul>
<li>Integrate proteomics datasets such as CPTAC</li>
<li>Add protein half-life features</li>
<li>Incorporate translation efficiency metrics</li>
<li>Analyze ribosome profiling datasets</li>
</ul>

<hr>

<h1>Conclusion</h1>

<p>
RNA expression alone is insufficient to explain protein presence in cancer.
</p>

<p>
By integrating RNA expression, gene features, and subcellular localization,
this project demonstrates how machine learning can both improve protein
detectability prediction and uncover biological drivers of RNA–protein
discordance.
</p>