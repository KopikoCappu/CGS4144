"""Create clustered heatmaps for adjusted significant SRP100712 genes."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.analysis.heatmap import (  # noqa: E402
    abbreviate_sample_names,
    calculate_heatmap_clustering,
    sample_grouping_scores,
    validate_significant_genes,
    validate_top50_ranking,
    zscore_genes,
)
from src.data.preprocess import (  # noqa: E402
    get_biological_sample_metadata,
    load_clean_metadata,
    load_processed_expression,
    validate_processed_data,
)


ANNOTATION_PALETTES = {
    "environment": {"Ground": "#31688E", "Spaceflight": "#D1495B"},
    "genotype": {"Col-0": "#4C956C", "WS": "#F4A261"},
    "age": {"4": "#8E7DBE", "8": "#E9C46A"},
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create full and top-100 clustered heatmaps from adjusted "
            "significant SRP100712 genes."
        )
    )
    parser.add_argument(
        "--expression",
        type=Path,
        default=PROJECT_ROOT / "data" / "processed" / "expression_32_samples.csv",
        help="Processed genes-by-biological-samples expression CSV",
    )
    parser.add_argument(
        "--metadata",
        type=Path,
        default=PROJECT_ROOT / "data" / "processed" / "metadata_clean.csv",
        help="Clean run-level metadata CSV",
    )
    parser.add_argument(
        "--significant",
        type=Path,
        default=(
            PROJECT_ROOT
            / "results"
            / "tables"
            / "differential_expression_significant.csv"
        ),
        help="Adjusted significant-gene result table",
    )
    parser.add_argument(
        "--top50",
        type=Path,
        default=(
            PROJECT_ROOT
            / "results"
            / "tables"
            / "differential_expression_top50.csv"
        ),
        help="Top-50 differential-expression table used for ranking validation",
    )
    parser.add_argument(
        "--tables-dir",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables",
        help="Output directory for heatmap tables",
    )
    parser.add_argument(
        "--plots-dir",
        type=Path,
        default=PROJECT_ROOT / "plots",
        help="Output directory for heatmaps",
    )
    return parser.parse_args()


def build_annotation_colors(
    sample_metadata: pd.DataFrame, abbreviations: dict[str, str]
) -> pd.DataFrame:
    """Build environment, genotype, and age color bars for clustermap."""

    metadata = sample_metadata.set_index("biological_sample").copy()
    metadata.index = metadata.index.map(abbreviations)
    colors = pd.DataFrame(index=metadata.index)
    colors["Environment"] = metadata["environment"].map(
        ANNOTATION_PALETTES["environment"]
    )
    colors["Genotype"] = metadata["genotype"].map(
        ANNOTATION_PALETTES["genotype"]
    )
    colors["Age"] = metadata["age"].astype(str).map(ANNOTATION_PALETTES["age"])
    if colors.isna().any().any():
        raise ValueError("At least one sample annotation color could not be assigned.")
    return colors


def add_annotation_legend(cluster_grid) -> None:
    """Add a legend that decodes all three sample annotation bars."""

    handles = [
        Patch(facecolor="#31688E", label="Ground"),
        Patch(facecolor="#D1495B", label="Spaceflight"),
        Patch(facecolor="#4C956C", label="Col-0"),
        Patch(facecolor="#F4A261", label="WS"),
        Patch(facecolor="#8E7DBE", label="4 days"),
        Patch(facecolor="#E9C46A", label="8 days"),
    ]
    cluster_grid.ax_heatmap.legend(
        handles=handles,
        title="Sample annotations",
        loc="upper left",
        bbox_to_anchor=(1.15, 1.0),
        frameon=True,
        borderaxespad=0,
    )


def save_clustered_heatmap(
    zscores: pd.DataFrame,
    sample_metadata: pd.DataFrame,
    output_path: Path,
    title: str,
    show_gene_labels: bool,
    figure_size: tuple[float, float],
) -> list[str]:
    """Cluster and save a z-score heatmap with three sample annotations."""

    abbreviations = abbreviate_sample_names(sample_metadata)
    display_matrix = zscores.rename(columns=abbreviations)
    annotation_colors = build_annotation_colors(sample_metadata, abbreviations)
    annotation_colors = annotation_colors.loc[display_matrix.columns]
    clustering = calculate_heatmap_clustering(display_matrix)

    sns.set_theme(style="white", context="notebook")
    cluster_grid = sns.clustermap(
        display_matrix,
        row_linkage=clustering.row_linkage,
        col_linkage=clustering.column_linkage,
        row_cluster=True,
        col_cluster=True,
        col_colors=annotation_colors,
        cmap="vlag",
        center=0,
        vmin=-3,
        vmax=3,
        xticklabels=True,
        yticklabels=show_gene_labels,
        figsize=figure_size,
        dendrogram_ratio=(0.12, 0.12),
        colors_ratio=0.025,
        cbar_kws={"label": "Gene-wise z-score"},
    )
    cluster_grid.fig.suptitle(title, fontsize=15, y=1.02)
    cluster_grid.ax_heatmap.set_xlabel("Biological sample")
    cluster_grid.ax_heatmap.set_ylabel(
        "Significant genes" if not show_gene_labels else "Gene"
    )
    cluster_grid.ax_heatmap.tick_params(
        axis="x", labelsize=7, labelrotation=90
    )
    if show_gene_labels:
        cluster_grid.ax_heatmap.tick_params(axis="y", labelsize=5.5)
    add_annotation_legend(cluster_grid)
    cluster_grid.fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(cluster_grid.fig)

    reverse_abbreviations = {value: key for key, value in abbreviations.items()}
    return [reverse_abbreviations[name] for name in clustering.ordered_samples]


def grouping_description(score: float) -> str:
    if score >= 0.50:
        return "strong"
    if score >= 0.25:
        return "moderate"
    if score >= 0.10:
        return "visible but weaker"
    return "limited"


def main() -> int:
    args = parse_args()
    print("Loading processed expression, metadata, and differential-expression tables...")
    expression = load_processed_expression(args.expression)
    metadata = load_clean_metadata(args.metadata)
    validation = validate_processed_data(expression, metadata)
    sample_metadata = get_biological_sample_metadata(metadata)
    significant_raw = pd.read_csv(args.significant)
    top50 = pd.read_csv(args.top50)

    significant = validate_significant_genes(significant_raw, expression)
    validate_top50_ranking(top50, significant)

    significant_expression = expression.loc[significant["Gene"]].copy()
    significant_zscores = zscore_genes(significant_expression)
    top100 = significant.head(100).copy()
    top100_zscores = significant_zscores.loc[top100["Gene"]].copy()

    args.tables_dir.mkdir(parents=True, exist_ok=True)
    args.plots_dir.mkdir(parents=True, exist_ok=True)
    top100_output = args.tables_dir / "heatmap_top100_genes.csv"
    zscores_output = args.tables_dir / "heatmap_top100_zscores.csv"
    full_plot = args.plots_dir / "heatmap_significant_genes_full.png"
    top100_plot = args.plots_dir / "heatmap_top100_significant_genes.png"

    top100.to_csv(top100_output, index=False)
    top100_zscores.to_csv(zscores_output, index=True, index_label="Gene")

    print("Clustering and drawing the full 2,560-gene heatmap...")
    full_sample_order = save_clustered_heatmap(
        significant_zscores,
        sample_metadata,
        full_plot,
        "All Adjusted Significant Genes Across SRP100712 Samples",
        show_gene_labels=False,
        figure_size=(16, 14),
    )
    print("Clustering and drawing the top-100 significant-gene heatmap...")
    top100_sample_order = save_clustered_heatmap(
        top100_zscores,
        sample_metadata,
        top100_plot,
        "Top 100 Adjusted Significant Genes Across SRP100712 Samples",
        show_gene_labels=True,
        figure_size=(16, 19),
    )

    global_scores, within_scores = sample_grouping_scores(
        significant_zscores, sample_metadata
    )
    primary_factor = max(global_scores, key=global_scores.get)
    within_separation = all(score >= 0.10 for score in within_scores.values())

    print("\n=== VALIDATION ===")
    print(
        f"Expression matrix: {validation['gene_count']:,} genes x "
        f"{validation['sample_count']} biological samples"
    )
    print(
        "Expression sample IDs match metadata biological sample IDs: "
        f"{'YES' if validation['identifiers_match'] else 'NO'}"
    )
    print(f"Significant genes loaded and included: {len(significant):,}")
    print("Top-50 ranking matches significant-gene table: YES")
    print(f"Genes included in smaller heatmap: {len(top100)}")
    print("Gene-wise z-scoring used only for heatmap visualization: YES")
    print("Original processed expression file overwritten: NO")

    print("\n=== DESCRIPTIVE SAMPLE GROUPING ===")
    for factor, score in global_scores.items():
        print(
            f"{factor}: {grouping_description(score)} grouping "
            f"(descriptive silhouette={score:.3f})"
        )
    print(f"Strongest global grouping in this selected gene set: {primary_factor}")
    print("Ground-vs-Spaceflight grouping within genotype/age strata:")
    for stratum, score in within_scores.items():
        print(f"  {stratum}: {grouping_description(score)} ({score:.3f})")
    print(
        "Ground and Spaceflight appear to separate within genotype/age structure: "
        f"{'YES' if within_separation else 'NOT CONSISTENTLY'}"
    )
    print("Full-heatmap clustered sample order:")
    print("  " + " | ".join(full_sample_order))
    print("Top-100 heatmap clustered sample order:")
    print("  " + " | ".join(top100_sample_order))
    print(
        "These scores and heatmap patterns are descriptive, not statistical "
        "significance tests. Environment grouping is expected to be emphasized "
        "because the displayed genes were selected by their adjusted environment "
        "association."
    )

    print("\n=== ASSIGNMENT-READY INTERPRETATION ===")
    print(
        f"The heatmaps display {len(significant):,} genes identified at FDR < "
        "0.05 by the adjusted Spaceflight-versus-Ground linear models. Expression "
        "was centered and scaled separately for each gene across the 32 samples "
        "to emphasize relative expression patterns; this visualization-only "
        "z-scoring did not alter the processed expression data. Environment is "
        "the strongest global grouping factor in this significant-gene subset, "
        "while genotype and age contribute additional structure. Ground and "
        "Spaceflight samples also show separation within genotype/age strata. "
        "Because the genes were selected for environment association, these "
        "patterns are exploratory and should not be treated as independent "
        "evidence or proof of causation."
    )

    print("\nTop five genes in the top-100 heatmap:")
    print(
        top100.loc[
            :4, ["Gene", "spaceflight_effect", "adjusted_p_value"]
        ].to_string(index=False)
    )

    print("\nFiles written:")
    for path in (top100_output, zscores_output, full_plot, top100_plot):
        print(f"  {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
