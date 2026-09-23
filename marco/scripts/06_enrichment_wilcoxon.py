"""Mann-Whitney enrichment of GO Biological Process gene sets."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.analysis.enrichment import (  # noqa: E402
    annotation_summary,
    build_gene_sets,
    load_differential_expression,
    load_or_query_go_bp_annotations,
    run_mannwhitney_enrichment,
    select_root_related_terms,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Test GO Biological Process gene sets for shifts in adjusted "
            "Spaceflight effects using two-sided Mann-Whitney tests."
        )
    )
    parser.add_argument(
        "--differential-expression",
        type=Path,
        default=(
            PROJECT_ROOT
            / "results"
            / "tables"
            / "differential_expression_all.csv"
        ),
    )
    parser.add_argument(
        "--annotation-cache",
        type=Path,
        default=(
            PROJECT_ROOT
            / "data"
            / "processed"
            / "arabidopsis_go_bp_annotations.csv"
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=PROJECT_ROOT / "results" / "enrichment",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    warnings: list[str] = []

    differential_expression = load_differential_expression(
        args.differential_expression
    )
    annotations, queried_count = load_or_query_go_bp_annotations(
        differential_expression["Gene"], args.annotation_cache
    )
    mapping = annotation_summary(annotations, differential_expression["Gene"])
    gene_sets = build_gene_sets(
        annotations,
        pd.Index(differential_expression["Gene"]),
        minimum_size=10,
        maximum_size=2_000,
    )

    print(f"Testing {len(gene_sets):,} eligible GO:BP terms...")
    enrichment = run_mannwhitney_enrichment(
        differential_expression, gene_sets, fdr_threshold=0.05
    )
    significant = enrichment.loc[enrichment["significant"]].copy()
    top20 = enrichment.head(20).drop(columns=["statistic", "significant"])
    root_related = select_root_related_terms(enrichment)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    all_output = args.output_dir / "wilcoxon_GO_BP_all.csv"
    significant_output = args.output_dir / "wilcoxon_GO_BP_significant.csv"
    top20_output = args.output_dir / "wilcoxon_GO_BP_top20.csv"
    root_output = args.output_dir / "wilcoxon_GO_BP_root_related.csv"
    enrichment.to_csv(all_output, index=False)
    significant.to_csv(significant_output, index=False)
    top20.to_csv(top20_output, index=False)
    root_related.to_csv(root_output, index=False)

    significant_positive = int(
        (significant["effect_difference"] > 0).sum()
    )
    significant_negative = int(
        (significant["effect_difference"] < 0).sum()
    )
    significant_zero = int((significant["effect_difference"] == 0).sum())
    if significant_zero:
        warnings.append(
            f"{significant_zero} significant terms have zero median effect difference."
        )
    if mapping["unannotated_genes"]:
        warnings.append(
            f"{mapping['unannotated_genes']:,} genes had no GO:BP annotation in "
            "the cached MyGene.info response; they remained in each outside set."
        )

    print("\n=== ANNOTATION AND TESTING SUMMARY ===")
    print(f"Total genes in differential-expression results: {mapping['total_genes']:,}")
    print(f"Genes with at least one GO:BP annotation: {mapping['annotated_genes']:,}")
    print(f"Genes without GO:BP annotations: {mapping['unannotated_genes']:,}")
    print(f"Total unique GO:BP terms: {mapping['unique_terms']:,}")
    print(f"GO terms tested after size filtering: {len(enrichment):,}")
    print(f"Significant GO terms at FDR < 0.05: {len(significant):,}")
    print(f"Significant terms with positive effect difference: {significant_positive:,}")
    print(f"Significant terms with negative effect difference: {significant_negative:,}")
    print(f"Minimum raw p-value: {enrichment['p_value'].min():.6e}")
    print(
        "Minimum adjusted p-value: "
        f"{enrichment['adjusted_p_value'].min():.6e}"
    )
    print(
        "Annotation cache status: "
        + (
            f"queried {queried_count:,} genes and updated the cache"
            if queried_count
            else "used the existing cache; no annotation query was needed"
        )
    )

    display_columns = [
        "GO_ID",
        "GO_term",
        "gene_set_size",
        "effect_difference",
        "adjusted_p_value",
        "significant",
    ]
    print("\nTop 10 statistically significant GO:BP terms:")
    print(enrichment.head(10)[display_columns].to_string(index=False))

    print(f"\nRoot-development-related tested terms ({len(root_related):,}):")
    if root_related.empty:
        print("  None")
    else:
        print(root_related[display_columns].to_string(index=False))

    print("\n=== ASSIGNMENT-READY INTERPRETATION ===")
    print(
        f"A two-sided Wilcoxon rank-sum approach, implemented with Mann-Whitney "
        f"tests, evaluated {len(enrichment):,} Gene Ontology Biological Process "
        "sets. All 32,309 genes were evaluated using the adjusted Spaceflight-"
        "effect statistic from the linear models, rather than restricting the "
        "analysis to significant genes. Benjamini-Hochberg correction was applied "
        f"across all tested terms, identifying {len(significant):,} terms at FDR "
        "< 0.05. Positive effect differences indicate a shift toward more positive "
        "Spaceflight effects among genes in a term, while negative differences "
        "indicate a shift toward more negative effects. These enrichments describe "
        "associations and do not establish causal effects on root development."
    )

    print("\nWarnings:")
    if warnings:
        for warning in warnings:
            print(f"  - {warning}")
    else:
        print("  None")

    print("\nFiles written:")
    for path in (
        args.annotation_cache,
        all_output,
        significant_output,
        top20_output,
        root_output,
    ):
        print(f"  {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
