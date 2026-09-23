# CGS4144 Project
## Plant Root Development During Spaceflight

### Scientific Question

Which genes and biological pathways associated with root development show different expression patterns in *Arabidopsis thaliana* grown during spaceflight on the International Space Station compared with ground controls?

This analysis identifies statistical associations between spaceflight and gene expression. It does not establish that any gene or pathway directly causes a change in root development.

### 1. Dataset and Experimental Design

The data came from Refine.bio experiment SRP100712, which measured gene expression in *Arabidopsis thaliana* roots. The original expression table contained 32,309 genes and 594 SRR sequencing-run columns. These runs represented 32 biological samples defined by environment, age, genotype, and replicate:

- Environment: Ground or Spaceflight
- Age: 4 or 8 days
- Genotype: Col-0 or WS
- Replication: four biological replicates for each environment × age × genotype combination

Multiple SRR accessions mapped to the same biological replicate. Treating those runs as independent samples would have inflated the sample size and introduced pseudoreplication. The run-level columns were therefore averaged within each biological sample, producing a 32,309-gene × 32-sample matrix with 16 Ground and 16 Spaceflight samples.

Refine.bio supplied quantile-normalized continuous expression values rather than raw RNA-seq counts. The existing scale was preserved because it contained negative values and had already been normalized. No count-based model, integer rounding, or additional `log2(x + 1)` transformation was applied.

### 2. Expression Variability

Gene variability was summarized as the maximum minus minimum expression value across the 32 biological samples. The median gene expression range was 0.364844. The mean was 0.411838, the 25th percentile was 0.126166, and the 75th percentile was 0.557759. The largest observed range was 4.577637.

![Density distribution of per-gene expression ranges](../plots/expression_variability_density.png)

*Figure 1. Density distribution of expression ranges for 32,309 genes across 32 biological samples.*

Most genes had relatively small expression ranges, while a smaller right tail contained genes with much larger changes across samples. Half of the genes had a range at or below 0.365, and 75% had a range at or below 0.558. This pattern supports focusing exploratory analyses on genes that vary across the experimental conditions while retaining all genes for the adjusted regression analysis.

### 3. Dimensionality Reduction

#### PCA

PCA used all 32,309 non-zero-variance genes. The values were centered by PCA but were not standardized gene by gene. PC1 explained 29.346% of the expression variance, PC2 explained 19.720%, and together they explained 49.066%.

![PCA colored by environment](../plots/PCA_flight_vs_ground.png)

*Figure 2. PCA of the 32 biological samples, colored by Ground or Spaceflight environment.*

![PCA colored by environment and shaped by genotype](../plots/PCA_flight_vs_ground_genotype.png)

*Figure 3. PCA colored by environment, with marker shape indicating Col-0 or WS genotype.*

The PCA showed that genotype was the strongest global organizing feature, with age contributing additional structure. Ground and Spaceflight samples often occupied different local positions within genotype and age groups, but the two environments did not form two completely isolated global clusters.

#### t-SNE

t-SNE used the 2,000 genes with the highest variance, a perplexity of 5, and a fixed random seed of 4144.

![t-SNE colored by environment](../plots/tSNE_flight_vs_ground.png)

*Figure 4. t-SNE representation of the 32 samples using the 2,000 most variable genes.*

The t-SNE plot produced distinct local groups dominated by genotype-related structure. Ground and Spaceflight samples showed local separation within several groups, but environment did not define a single global division.

#### UMAP

UMAP used the same 2,000 genes, with 8 neighbors, a minimum distance of 0.2, and random seed 4144.

![UMAP colored by environment](../plots/UMAP_flight_vs_ground.png)

*Figure 5. UMAP representation of the 32 samples using the 2,000 most variable genes.*

UMAP also emphasized strong genotype-related organization. Age and environment contributed smaller within-group patterns. As with t-SNE, distances and axes do not represent percentages of explained variance.

Across all three methods, genotype was the clearest source of global structure, and age also affected sample arrangement. Environment-related differences appeared within that broader structure rather than as two fully isolated clusters. These plots are exploratory and do not provide tests of statistical significance.

### 4. Differential Expression

A separate ordinary least squares model was fitted for every gene:

`expression ~ environment + age + genotype`

Environment, age, and genotype were modeled as categorical variables. Ground was the reference environment, so `spaceflight_effect` represents Spaceflight minus Ground after adjustment for age and genotype. Positive values indicate higher normalized expression in Spaceflight, and negative values indicate lower normalized expression in Spaceflight. This value is not described as a log2 fold change because the exact transformation underlying the processed Refine.bio scale is not established.

Benjamini-Hochberg correction was applied across all 32,309 genes. At FDR < 0.05, 2,560 genes were significant. Of these, 1,117 had positive Spaceflight effects and 1,443 had negative effects.

The five most statistically significant genes were:

| Gene | Spaceflight effect | Standard error | t statistic | Raw p-value | Adjusted p-value |
|---|---:|---:|---:|---:|---:|
| AT4G23496 | -0.536245 | 0.040912 | -13.107388 | 1.806e-13 | 5.835e-09 |
| AT4G18422 | -0.247720 | 0.019482 | -12.715626 | 3.753e-13 | 6.063e-09 |
| AT3G52740 | -0.548170 | 0.044787 | -12.239446 | 9.334e-13 | 1.005e-08 |
| AT4G15248 | -0.300948 | 0.031533 | -9.543954 | 2.670e-10 | 1.791e-06 |
| AT3G24520 | -0.359490 | 0.037852 | -9.497294 | 2.968e-10 | 1.791e-06 |

![Volcano plot of adjusted Spaceflight effects](../plots/volcano_spaceflight_vs_ground.png)

*Figure 6. Volcano plot of adjusted Spaceflight effects against −log10 adjusted p-values. Red points passed FDR < 0.05; no effect-size threshold was imposed.*

The volcano plot shows significant genes on both sides of zero, although more significant genes had negative than positive Spaceflight effects. The strongest statistical signals were largely negative, while some genes had large positive effects. These results identify candidate environment-associated genes and do not show that the genes cause altered root growth.

### 5. Top 50 Differentially Expressed Genes

The following genes had the 50 smallest adjusted p-values. Effects are normalized-expression differences for Spaceflight minus Ground after adjustment for age and genotype.

| Gene | Spaceflight effect | Raw p-value | Adjusted p-value |
|---|---:|---:|---:|
| AT4G23496 | -0.536245 | 1.806e-13 | 5.835e-09 |
| AT4G18422 | -0.247720 | 3.753e-13 | 6.063e-09 |
| AT3G52740 | -0.548170 | 9.334e-13 | 1.005e-08 |
| AT4G15248 | -0.300948 | 2.670e-10 | 1.791e-06 |
| AT3G24520 | -0.359490 | 2.968e-10 | 1.791e-06 |
| AT4G25830 | -0.226281 | 3.617e-10 | 1.791e-06 |
| AT5G66170 | 0.596014 | 3.881e-10 | 1.791e-06 |
| AT1G10370 | -0.388380 | 6.896e-10 | 2.785e-06 |
| AT4G14690 | -0.597412 | 1.159e-09 | 3.396e-06 |
| AT5G54470 | -0.165231 | 1.189e-09 | 3.396e-06 |
| AT5G53980 | 0.303192 | 1.219e-09 | 3.396e-06 |
| AT5G67210 | -0.304810 | 1.286e-09 | 3.396e-06 |
| AT5G23280 | -0.272256 | 1.366e-09 | 3.396e-06 |
| AT2G44080 | 0.565546 | 1.668e-09 | 3.849e-06 |
| AT1G67030 | -0.365319 | 2.578e-09 | 5.553e-06 |
| AT2G27370 | -0.420603 | 2.823e-09 | 5.701e-06 |
| AT4G04750 | -0.213938 | 3.171e-09 | 5.877e-06 |
| AT1G54540 | -0.356947 | 3.274e-09 | 5.877e-06 |
| AT1G21320 | 0.209215 | 3.672e-09 | 6.040e-06 |
| AT5G54230 | -0.246023 | 3.739e-09 | 6.040e-06 |
| AT2G36870 | -0.339902 | 4.175e-09 | 6.424e-06 |
| AT1G30750 | -0.636911 | 4.831e-09 | 6.840e-06 |
| AT5G54490 | 0.277633 | 4.869e-09 | 6.840e-06 |
| AT2G36100 | -0.646233 | 7.556e-09 | 1.017e-05 |
| AT5G13900 | -0.508712 | 1.007e-08 | 1.301e-05 |
| AT5G57660 | -0.194684 | 1.115e-08 | 1.386e-05 |
| AT4G29930 | -0.276786 | 1.176e-08 | 1.408e-05 |
| AT4G21660 | 0.220779 | 1.304e-08 | 1.487e-05 |
| AT3G24020 | -0.519566 | 1.334e-08 | 1.487e-05 |
| AT3G17130 | -0.187572 | 1.458e-08 | 1.570e-05 |
| AT5G47450 | -0.670674 | 1.593e-08 | 1.655e-05 |
| AT5G41590 | -0.241570 | 1.678e-08 | 1.655e-05 |
| AT5G08050 | -0.345019 | 1.690e-08 | 1.655e-05 |
| AT4G13580 | -0.590292 | 1.807e-08 | 1.717e-05 |
| AT4G36610 | -0.331019 | 2.000e-08 | 1.846e-05 |
| AT3G55230 | -0.649711 | 2.367e-08 | 2.124e-05 |
| AT3G13650 | -0.300347 | 2.498e-08 | 2.182e-05 |
| AT3G21560 | -0.487060 | 2.698e-08 | 2.294e-05 |
| AT5G25350 | 0.451181 | 2.984e-08 | 2.460e-05 |
| AT1G73650 | -0.273971 | 3.045e-08 | 2.460e-05 |
| AT5G07475 | -0.355828 | 3.165e-08 | 2.494e-05 |
| AT2G39430 | -0.589221 | 3.701e-08 | 2.847e-05 |
| AT3G11550 | -0.488275 | 4.275e-08 | 3.212e-05 |
| AT3G22620 | -0.535947 | 4.489e-08 | 3.296e-05 |
| AT5G09520 | -0.310112 | 4.864e-08 | 3.484e-05 |
| AT3G59900 | 0.303662 | 5.148e-08 | 3.484e-05 |
| AT1G59760 | 0.232490 | 5.160e-08 | 3.484e-05 |
| AT2G48130 | -0.461489 | 5.175e-08 | 3.484e-05 |
| AT3G56290 | -0.488290 | 5.296e-08 | 3.492e-05 |
| AT3G25655 | -0.331438 | 5.864e-08 | 3.789e-05 |

### 6. Significant-Gene Heatmap

The full heatmap used all 2,560 genes that passed FDR < 0.05. A second heatmap used the top 100 genes ranked by adjusted p-value. Expression was centered and scaled separately for each gene across the 32 samples only for visualization. These z-scores did not replace or alter the processed expression matrix.

![Full heatmap of significant genes](../plots/heatmap_significant_genes_full.png)

*Figure 7. Clustered heatmap of all 2,560 significant genes. Rows and columns were hierarchically clustered; sample annotation bars show environment, genotype, and age. Gene labels were omitted because of the number of rows.*

![Top-100 significant-gene heatmap](../plots/heatmap_top100_significant_genes.png)

*Figure 8. Clustered heatmap of the 100 genes with the smallest adjusted p-values. Gene labels are displayed, and sample annotation bars indicate environment, genotype, and age.*

Within this selected significant-gene subset, environment was the strongest global grouping factor. The full heatmap placed Ground and Spaceflight samples into separate major branches, with genotype and age producing additional substructure. Ground and Spaceflight separation was also visible within genotype and age groups. Because the genes were selected using their adjusted environment association, this separation is expected and should not be treated as independent confirmation of an environment effect.

### 7. My Gene-Set Enrichment Analysis

The primary individual enrichment analysis used two-sided Mann-Whitney tests, also described as a Wilcoxon rank-sum approach, with Gene Ontology Biological Process gene sets. All 32,309 genes were evaluated using `spaceflight_effect`; the analysis was not restricted to significant genes. For each eligible GO term, effects for genes in the term were compared with effects for genes outside the term. Terms represented by 10–2,000 genes were tested, and Benjamini-Hochberg correction was applied across the 1,020 eligible terms.

Of the 1,020 terms tested, 236 passed FDR < 0.05. There were 155 significant terms with positive median effect differences and 81 with negative differences. A positive difference means that genes in the GO set tended to have more positive Spaceflight effects than genes outside it; a negative difference means that they tended to have more negative effects.

The ten strongest statistical results were:

| GO ID | Biological process | Set size | Effect difference | Adjusted p-value |
|---|---|---:|---:|---:|
| GO:0006412 | translation | 346 | 0.065873 | 1.149e-53 |
| GO:0000398 | mRNA splicing, via spliceosome | 170 | 0.058214 | 1.874e-32 |
| GO:0006364 | rRNA processing | 117 | 0.063519 | 6.256e-20 |
| GO:0006457 | protein folding | 211 | 0.049211 | 6.256e-20 |
| GO:0006338 | chromatin remodeling | 164 | 0.041491 | 1.181e-16 |
| GO:0008380 | RNA splicing | 88 | 0.040311 | 1.011e-14 |
| GO:0002181 | cytoplasmic translation | 69 | 0.070756 | 3.759e-14 |
| GO:0006413 | translational initiation | 75 | 0.054374 | 4.279e-13 |
| GO:0006396 | RNA processing | 182 | 0.021210 | 1.394e-12 |
| GO:0006397 | mRNA processing | 114 | 0.024686 | 1.960e-11 |

Several processes directly relevant to the scientific question were significant:

| GO ID | Biological process | Set size | Effect difference | Adjusted p-value | Significant |
|---|---|---:|---:|---:|---|
| GO:0009834 | plant-type secondary cell wall biogenesis | 80 | -0.047145 | 6.680e-09 | Yes |
| GO:0042546 | cell wall biogenesis | 65 | -0.042003 | 4.937e-06 | Yes |
| GO:0048364 | root development | 214 | 0.012414 | 5.001e-04 | Yes |
| GO:0009926 | auxin polar transport | 43 | 0.029051 | 5.868e-03 | Yes |
| GO:0010449 | root meristem growth | 16 | 0.048450 | 8.314e-03 | Yes |
| GO:0010082 | regulation of root meristem growth | 43 | 0.018425 | 3.040e-02 | Yes |
| GO:0080022 | primary root development | 45 | -0.008491 | 1.140e-01 | No |
| GO:0080147 | root hair cell development | 33 | -0.007292 | 2.459e-01 | No |
| GO:0009629 | response to gravity | 11 | -0.016859 | 3.748e-01 | No |
| GO:0009630 | gravitropism | 50 | 0.006597 | 3.757e-01 | No |
| GO:0048527 | lateral root development | 56 | 0.002370 | 9.355e-01 | No |

The positive shifts for root development, root meristem growth, and auxin polar transport are consistent with altered regulation of developmental processes in Spaceflight. The negative shifts for cell-wall biogenesis terms indicate that genes in these sets tended toward more negative adjusted Spaceflight effects. Primary-root, root-hair, gravity-response, gravitropism, and lateral-root terms were found in the exploratory search but were not significant after multiple-testing correction.

GO:BP annotations were available for 22,238 of the 32,309 genes in the MyGene.info mapping. The remaining 10,071 genes had no returned GO:BP annotation and remained in the outside group for each test.

### 8. Combined Team Enrichment

The team comparison combined my Wilcoxon/Mann-Whitney GO Biological Process results with the teammate's g:Profiler GO Biological Process results. Terms were matched by GO ID. Method-specific p-values were retained separately because the two methods use different inputs, null hypotheses, and statistical procedures; p-values were not averaged. `methods_tested` records how many result sets contained a term, and `methods_significant` records how many methods classified it as significant.

The Wilcoxon analysis contained 1,020 terms and the teammate g:Profiler analysis contained 966 terms. Their outer join produced 1,576 unique terms, including 410 present in both methods. Five terms were significant in both methods, 231 were significant only in Wilcoxon, and 61 were significant only in g:Profiler. In total, 297 terms were significant in at least one method.

The top 10 combined terms, ranked by number of significant methods, significant fraction, and then the strongest available method-specific significance value, were:

| GO ID | Biological process | Methods tested | Methods significant | Wilcoxon adjusted p-value | g:Profiler p-value |
|---|---|---:|---:|---:|---:|
| GO:0071456 | cellular response to hypoxia | 2 | 2 | 2.100e-02 | 1.370e-03 |
| GO:0009699 | phenylpropanoid biosynthetic process | 2 | 2 | 1.056e-02 | 2.372e-03 |
| GO:0009873 | ethylene-activated signaling pathway | 2 | 2 | 2.648e-03 | 1.254e-02 |
| GO:0035556 | intracellular signal transduction | 2 | 2 | 6.369e-03 | 1.498e-02 |
| GO:0006979 | response to oxidative stress | 2 | 2 | 9.761e-03 | 2.906e-02 |
| GO:0006364 | rRNA processing | 1 | 1 | 6.256e-20 | — |
| GO:0000463 | maturation of LSU-rRNA from tricistronic rRNA transcript (SSU-rRNA, 5.8S rRNA, LSU-rRNA) | 1 | 1 | 1.960e-11 | — |
| GO:0042254 | ribosome biogenesis | 1 | 1 | 5.366e-11 | — |
| GO:0042221 | response to chemical | 1 | 1 | — | 2.202e-10 |
| GO:0002183 | cytoplasmic translational initiation | 1 | 1 | 2.790e-10 | — |

The five terms significant in both methods were cellular response to hypoxia, phenylpropanoid biosynthetic process, ethylene-activated signaling pathway, intracellular signal transduction, and response to oxidative stress. These results show the strongest cross-method consistency in the team comparison.

Among terms selected for their relevance to root biology, 57 appeared in both result sets and 22 were significant in at least one method. None of these root-related terms were significant in both methods. GO:0048364, root development, appeared in both sets and was significant in Wilcoxon (adjusted p = 5.001268e-04), but it was not significant in the teammate g:Profiler analysis (p approximately 1.0). Root development therefore did not have cross-method significance in the team results.

### 9. Biological Interpretation

The adjusted regression identified substantial transcriptional differences associated with Spaceflight, with 2,560 genes passing FDR < 0.05. Genotype and age were important sources of global expression variation, which supports their inclusion as covariates and cautions against interpreting an unadjusted Ground-versus-Spaceflight comparison.

The strongest shared team signals were cellular response to hypoxia, phenylpropanoid biosynthetic process, ethylene-activated signaling pathway, intracellular signal transduction, and response to oxidative stress. Each was significant in both the Wilcoxon and teammate g:Profiler analyses. This agreement indicates that these processes were identified under two different enrichment frameworks, but it does not establish a causal role in the Spaceflight response.

Root-related findings were less consistent between methods. Root development, root meristem growth, regulation of root meristem growth, auxin-related regulation, and cell-wall processes were identified by the Wilcoxon analysis, but no root-related term was significant in both team methods. For root development specifically, the positive Wilcoxon effect difference indicated a shift toward more positive adjusted Spaceflight effects among genes in the set, while the teammate g:Profiler result was not significant. This should be described as support from one method rather than cross-method agreement.

The difference is consistent with what the methods test. Wilcoxon/Mann-Whitney compares the distribution of `spaceflight_effect` values for genes inside and outside a GO term and uses all tested genes. g:Profiler tests whether a GO term is over-represented in the list of significant genes. A distributed shift across a gene set can therefore be detected by Wilcoxon without producing significant over-representation in g:Profiler. Broader results involving stress, hormone signaling, secondary metabolism, translation, RNA processing, and protein folding suggest that the Spaceflight-associated expression pattern extends beyond root-specific processes.

The most significant locus identifiers and genes contributing to enriched terms are candidates for additional study. The present data support associations with Spaceflight but do not establish that these genes directly cause changes in root structure or development.

### 10. Limitations

- The analysis used processed, quantile-normalized continuous values rather than raw sequencing counts. This prevented count-based modeling and limited interpretation of the expression scale.
- `spaceflight_effect` is an adjusted normalized-expression difference, not a raw RNA-seq log2 fold change.
- GO annotation coverage was incomplete. MyGene.info returned GO:BP annotations for 22,238 of 32,309 genes.
- Enrichment results show associations between gene sets and Spaceflight-related expression effects; they do not demonstrate causality.
- Genotype and age strongly influenced global expression patterns. Although both were included in the regression, the single primary environment coefficient does not describe genotype- or age-specific interactions.
- The same dataset was used for differential expression, gene selection, heatmap visualization, and enrichment. Visual separation in the significant-gene heatmaps is therefore not independent validation.
- Wilcoxon/Mann-Whitney and g:Profiler evaluate different properties of a gene set, so disagreement between their significance results is expected and should not be interpreted as a direct contradiction.

### 11. Conclusion

Spaceflight was associated with widespread expression differences in *Arabidopsis thaliana* roots after adjustment for genotype and age. The primary analysis identified 2,560 significant genes, while dimensionality reduction showed that genotype and age remained major sources of overall expression structure. The team enrichment comparison found five processes that were significant in both methods: cellular response to hypoxia, phenylpropanoid biosynthetic process, ethylene-activated signaling pathway, intracellular signal transduction, and response to oxidative stress.

Root development, root meristem growth, auxin-related regulation, and cell-wall biology remained candidate processes associated with Spaceflight-related expression differences, but root-development-related processes were not consistently significant across both team methods. These findings support further study of the identified genes and pathways, but they do not establish direct causal mechanisms.
