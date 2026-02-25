<body>

<h1>Exploring Understudied Human Protein-Coding Genes in Cancer</h1>

<h2>Abstract</h2>
<p>In this project, I explored a group of understudied human protein-coding genes to evaluate their potential relevance in cancer. In addition to four provided genes, I selected five underrepresented genes: <strong>EMC9, FAM171A2, PDIK1L, DBNDD2, CASKIN2</strong>, based on limited functional annotation in UniProt, Pfam/InterPro, and the literature. For each gene, basic information on predicted domains, Gene Ontology (GO) terms, and expression patterns was compiled.</p>

<p>Using <strong>GEPIA2</strong>, tumor and normal expression across TCGA and GTEx datasets was compared. Notable trends include EMC9 overexpression in thymoma and DLBCL, and PDIK1L survival associations in kidney cancer. RNA-protein comparisons using the Human Protein Atlas (HPA) revealed some RNA-protein mismatches, which are expected for low-expression, poorly characterized genes. Single-cell RNA-seq datasets (GSE72056 and GSE115978) showed consistent cell-type expression patterns. A scalable Seurat workflow was developed to integrate multiple scRNA-seq datasets for cluster-level summaries. Overall, this study highlights how minimal public data can provide preliminary clues about uncharacterized genes, some of which may warrant further investigation in cancer biology.</p>

<h2>Cancer Target Discovery</h2>

<h3>Identifying and Characterizing Understudied Human Protein-Coding Genes</h3>

<h4>Rationale</h4>
<p>Genes were selected based on three criteria: (i) truly understudied, (ii) minimal known function, and (iii) detectable expression in human tissues or cancers suitable for analysis using TCGA and GTEx data.</p>

<h4>Seed Genes</h4>
<p>The exercise provided four initial genes: <strong>C1orf174, LOC124903857, TMEM161B, ZNF808</strong>. These served as anchors for defining criteria for poorly characterized protein-coding genes.</p>

<h4>Finding Tdark Genes</h4>
<p>Using the <strong>Pharos (NIH IDG)</strong> platform, genes were classified into Target Development Levels (TDLs):</p>
<ul>
    <li><strong>Tclin:</strong> established drug targets</li>
    <li><strong>Tchem:</strong> targets with small molecule interactions</li>
    <li><strong>Tbio:</strong> targets with substantial biological information</li>
    <li><strong>Tdark:</strong> targets with very little functional information, few publications, limited GO terms, and scarce protein domain data</li>
</ul>

<p>Selection filters for Tdark genes included:</p>
<ul>
    <li>2–10 PubMed publications</li>
    <li>Minimal Gene Ontology annotations</li>
    <li>Unreviewed UniProt (TrEMBL) entries with low annotation scores (≤2/5)</li>
    <li>Few or no recognized protein domains in InterPro/Pfam</li>
    <li>Limited proof of protein expression</li>
</ul>

<h4>Final Gene Panel</h4>
<p>Five Tdark genes were selected: <strong>EMC9, FAM171A2, PDIK1L, DBNDD2, CASKIN2</strong>.</p>

<table>
<tr>
<th>Gene</th><th>Gene ID</th><th>Description</th><th>Protein Family</th><th>Expression / Domains</th>
</tr>
<tr>
<td>C1orf174</td><td>339448</td><td>Chromosome 1 open reading frame 174</td><td>-</td><td>Prognostic marker in liver, lung, pancreatic cancers</td>
</tr>
<tr>
<td>LOC124903857</td><td>124903857</td><td>FAM231A/C-like protein</td><td>FAM231 family / Homeobox</td><td>-</td>
</tr>
<tr>
<td>TMEM161B</td><td>153396</td><td>Transmembrane protein 161B</td><td>Transmembrane protein family</td><td>Prognostic marker in kidney and rectal cancers</td>
</tr>
<tr>
<td>ZNF808</td><td>388558</td><td>Zinc Finger Protein 808</td><td>C2H2 zinc finger</td><td>Prognostic marker in kidney cancer</td>
</tr>
<tr>
<td>EMC9</td><td>51016</td><td>ER membrane protein complex subunit 9</td><td>ER membrane complex</td><td>Upregulated RNA in THYM & DLBC; weak/negative protein; EMC9 domain</td>
</tr>
<tr>
<td>FAM171A2</td><td>284069</td><td>Family with sequence similarity 171 member A2</td><td>FAM171 family</td><td>Prognostic marker in glioblastoma, kidney papillary carcinoma; Signal domain</td>
</tr>
<tr>
<td>PDIK1L</td><td>149420</td><td>PDLIM1-interacting kinase 1-like</td><td>PDLIM1-interacting kinase-like</td><td>Moderate expression; prognostic in kidney cancer; Phosphorylase kinase domain</td>
</tr>
<tr>
<td>DBNDD2</td><td>55861</td><td>Dysbindin domain containing 2</td><td>Dysbindin domain proteins</td><td>High RNA in brain, muscle, and adipose; Dysbindin domain</td>
</tr>
<tr>
<td>CASKIN2</td><td>57513</td><td>CASK interacting protein 2</td><td>CASK-interacting scaffold proteins</td><td>Low/moderate RNA; strong protein in ovarian/stomach cancers; Ankyrin repeats, SH3, SAM</td>
</tr>
</table>

<h3>Expression in Cancer</h3>
<p>Across all nine genes, functional annotation was sparse. HPA data and GEPIA2 analyses revealed potential cancer relevance:</p>

<ul>
<li><strong>EMC9:</strong> Low tissue specificity; modest RNA across normal tissues; protein moderate in colorectal, thyroid, prostate, and skin cancers. Discrepancies between RNA and protein suggest transcriptional activity with limited protein detection.</li>
<li><strong>FAM171A2:</strong> Highest in cerebral cortex and epididymis; mostly weak or negative staining in cancers, with moderate cytoplasmic positivity in liver and thyroid.</li>
<li><strong>PDIK1L:</strong> Low/moderate RNA across tissues; prognostic in kidney cancer (KIRC) with favorable survival association.</li>
<li><strong>DBNDD2:</strong> Broadly expressed; high in gliomas; high expression correlates with poor prognosis across multiple cancers.</li>
<li><strong>CASKIN2:</strong> Low/moderate RNA; strong protein in ovarian and stomach cancers; tumor-selective expression suggests potential tissue-specific roles.</li>
</ul>

<h3>Pan-Cancer Expression & Survival Analysis</h3>
<p>GEPIA2 boxplots showed RNA upregulation in certain cancers (e.g., EMC9 in DLBC and THYM). HPA protein data often showed weak or absent staining, consistent with low-expression Tdark genes. Survival analysis suggested:</p>
<ul>
<li>High DBNDD2 expression = worse outcomes in multiple cancers.</li>
<li>High EMC9 expression showed trends but did not reach statistical significance.</li>
</ul>

<h3>Single-Cell Profiling & Multi-Dataset Analysis</h3>
<p>Two melanoma scRNA-seq datasets (GSE72056: 19 patients, 4,645 cells; GSE115978: 31 patients, 7,186 cells) were analyzed. Observations:</p>
<ul>
<li>Consistent expression trends across datasets</li>
<li>Low-expressed genes (FAM171A2, CASKIN2) detectable but at low levels</li>
<li>Clustering shows reproducibility of gene expression across cell types</li>
<li>Seurat workflow allows integration of 10+ datasets, automated cluster-level visualization and pseudo-bulk summaries. Code available on <a href="https://github.com/">GitHub repository</a></li>
</ul>

<h3>Conclusion</h3>
<p>Systematic integration of public data reveals preliminary insights into understudied genes. Several candidates (DBNDD2, CASKIN2) show potential relevance in cancer biology and merit further experimental investigation.</p>

</body>
</html>
