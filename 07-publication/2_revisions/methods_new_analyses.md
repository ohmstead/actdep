# Draft Methods text — new analyses (revision)

For insertion into *Quantification and statistical analysis*. Subsection headings
follow the italicized style used in the submitted manuscript (Olmstead, King, &
Bloodgood 2025). Superscript citation numbers continue from ref. 71 in the
submitted reference list; see the reference block at the end.

---

### *Complementary pseudobulk differential expression (limma-voom)*

To confirm that DEG calls were not an artifact of a single statistical framework,
we repeated the pseudobulk analysis with limma-voom⁷² using the same unit of
analysis as DESeq2 (animal-level pseudobulk counts). For each Subclass, nuclei
were subset and grouped by sample of origin, genes expressed in fewer than 1% of
nuclei were excluded, and counts were aggregated with Seurat's
AggregateExpression(). Genes were further filtered with edgeR's filterByExpr()
and libraries normalized by the trimmed mean of M-values (TMM) method⁷³. A design
matrix without an intercept (expression ~ 0 + activity_condition) was fit with
voom precision weights and lmFit, and each contrast from the DESeq2 analysis was
tested with contrasts.fit() followed by empirical Bayes moderation
(eBayes(robust = TRUE)). Only contrasts meeting the same minimum per-group
nucleus threshold (n = 30 nuclei per condition) were tested. Genes with absolute
log2 fold change > 0.585 and fdr < 0.05 were classified as differentially
expressed, identical to the DESeq2 criteria. Ninety-five percent confidence
intervals on log2 fold change were retained for the forest plots. Concordance
between the DESeq2 and limma-voom DEG sets was quantified per Subclass and
contrast as the Jaccard index of the two DEG sets.

### *Negative-binomial mixed-effects models*

Pseudobulk aggregation collapses all nuclei from an animal into a single
observation. To confirm our results while modeling nuclei individually — and
without treating nuclei from the same animal as independent replicates — we fit
cell-level negative-binomial generalized linear mixed models (NB-GLMMs) with
glmmTMB⁷⁴, using the formula

    counts ~ activity_condition + offset(log(nCount_RNA)) + (1 | sample)

with a nbinom2 (quadratic) variance and one random intercept per animal. The
offset on total transcript count per nucleus accounts for differences in
sequencing depth. Models were fit separately for each Subclass, over an a priori
panel of the 15 canonical IEGs and 9 circadian genes (*Per1, Per2, Per3, Clock,
Bmal1, Cry1, Cry2, Nr1d1, Nr1d2*). Only conditions with at least 30 nuclei in a
given Subclass were included, matching the DESeq2 and limma-voom analyses. Fixed
effect estimates for each condition contrast against SE were converted to log2
fold change by dividing by log(2); 95% Wald confidence intervals were converted
identically. p-values were corrected with the Benjamini-Hochberg procedure within
each Subclass × contrast family.

### *Rank-rank hypergeometric overlap (RRHO)*

Agreement between statistical frameworks was assessed threshold-free by rank-rank
hypergeometric overlap⁷⁵. Within each gene list, genes were ranked by the signed
statistic sign(log2FC) × -log10(p), from most down-regulated to most
up-regulated. Over a 100 × 100 grid of rank cutoff pairs (i, j), the overlap
between the top i genes of the first list and the top j genes of the second was
tested against a hypergeometric null, and the resulting grid of -log10(p) values
plotted as a heatmap; concordant rankings produce a signal along the diagonal
running from down- to up-regulation. Hypergeometric tail probabilities were
evaluated on the log scale, as -log10(p) grows with list length and exceeds
double-precision range for genome-wide lists. As a negative control, each
comparison was repeated after randomly permuting the gene order of one list,
which abolishes the diagonal signal. Comparisons were performed for 016 CA1-ProS
Glut, EE30m vs. SE. DESeq2 and limma-voom results were compared genome-wide; for
the DESeq2 vs. NB-GLMM comparison, NB-GLMMs were fit to the 2,000 most highly
expressed genes passing the 1% expression filter, and both maps were restricted
to that common gene set so that the two panels share axis lengths and color
scale. Because map brightness scales with the number of genes tested and is
therefore not an effect size, quadrant statistics rather than map intensity are
reported: the minimum-p coordinate in each quadrant was located with RedRibbon's
evolutionary algorithm⁷⁶ with 96 permutations for the permutation-adjusted
p-value, and the p-value at that coordinate recomputed on the log scale.

### *Receiver-operating-characteristic and precision-recall analysis*

To confirm that the binary "active cell" definition (FIGURE 4) is not sensitive
to its thresholds, we compared it against a continuous IEG activity score. For
each cell, the log-normalized SCT expression of each of the 15 IEGs was z-scored
against the mean and standard deviation of that gene in SE cells of the CA1
Subclass, and the score taken as the equal-weight mean of the z-scored genes;
z-scoring per gene prevents IEGs with high baseline transcript abundance from
dominating the composite. Treating the binary call (90th percentile threshold,
≥ 3 of 15 IEGs induced) as the reference label and the continuous score as the
predictor, receiver-operating-characteristic (ROC) curves and the area under them
were computed with pROC⁷⁷, and precision-recall (PR) curves and the area under
them (AUPRC) with PRROC⁷⁸. Curves were computed separately for each activity
condition (SE, EE30m, KA30m) and, for ROC-AUC, separately for each CA1 Supertype;
Supertypes with fewer than 10 cells or fewer than 5 cells in either class were
not scored. Because a PR curve has no fixed null — a random classifier achieves
precision equal to the prevalence of active cells at every recall — the
prevalence baseline is drawn on each PR panel, and ROC-AUC, which is
prevalence-invariant, is the statistic compared across conditions. In addition,
the binary definition was swept across percentile cutoffs (80, 85, 90, 95%) and
minimum IEG requirements (2, 3, 4 genes) to confirm that the reported fraction of
active cells per Supertype and condition is stable within this parameter range.

### *Software*

Analyses were performed in R with limma v3.60.6, edgeR v4.2.2, glmmTMB v1.1.14,
broom.mixed v0.2.9.7, RedRibbon v1.4.1, pROC v1.19.0.1, and PRROC v1.4.

---

## New references (continuing from 71)

72. Law, C. W., Chen, Y., Shi, W. & Smyth, G. K. voom: precision weights unlock
    linear model analysis tools for RNA-seq read counts. *Genome Biology* 15, R29
    (2014).
73. Robinson, M. D. & Oshlack, A. A scaling normalization method for differential
    expression analysis of RNA-seq data. *Genome Biology* 11, R25 (2010).
74. Brooks, M. E. *et al.* glmmTMB balances speed and flexibility among packages
    for zero-inflated generalized linear mixed modeling. *The R Journal* 9,
    378-400 (2017).
75. Plaisier, S. B., Taschereau, R., Wong, J. A. & Graeber, T. G. Rank-rank
    hypergeometric overlap: identification of statistically significant overlap
    between gene-expression signatures. *Nucleic Acids Research* 38, e169 (2010).
76. Piron, A. *et al.* RedRibbon: A new rank-rank hypergeometric overlap for gene
    and transcript expression signatures. *Life Science Alliance* 7, e202302203
    (2024).
77. Robin, X. *et al.* pROC: an open-source package for R and S+ to analyze and
    compare ROC curves. *BMC Bioinformatics* 12, 77 (2011).
78. Grau, J., Grosse, I. & Keilwagen, J. PRROC: computing and visualizing
    precision-recall and receiver operating characteristic curves in R.
    *Bioinformatics* 31, 2595-2597 (2015).
