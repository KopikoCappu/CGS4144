from gprofiler import GProfiler

gp = GProfiler(return_dataframe=True)

enrichment = gp.profile(
    organism="athaliana",
    query=list(sig_genes),
    sources=["GO:BP"],
    background=list(log_expr.index),
    domain_scope="custom",
    no_evidences=False
)

enrichment = enrichment.sort_values("p_value")
enrichment.to_csv("results/enrichment_gprofiler_GOBP.csv", index=False)

print("Significant terms:", (enrichment["p_value"] < 0.05).sum())
enrichment.head(20)
