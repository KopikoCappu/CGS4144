"""Run PCA, t-SNE, and UMAP for the SRP100712 biological samples."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.analysis.dimensionality_reduction import (  # noqa: E402
    add_sample_metadata,
    embedding_silhouette_scores,
    prepare_sample_matrices,
    qualitative_grouping,
    run_pca,
    run_tsne,
    run_umap,
)
from src.data.preprocess import (  # noqa: E402
    get_biological_sample_metadata,
    load_clean_metadata,
    load_processed_expression,
    validate_processed_data,
)


ENVIRONMENT_PALETTE = {"Ground": "#31688E", "Spaceflight": "#D1495B"}
GENOTYPE_MARKERS = {"Col-0": "o", "WS": "s"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run PCA, t-SNE, and UMAP on the processed SRP100712 biological "
            "sample expression matrix."
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
        "--plots-dir",
        type=Path,
        default=PROJECT_ROOT / "plots",
        help="Output directory for plots",
    )
    parser.add_argument(
        "--tables-dir",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables",
        help="Output directory for coordinate tables",
    )
    return parser.parse_args()


def save_embedding_plot(
    coordinates: pd.DataFrame,
    x: str,
    y: str,
    x_label: str,
    y_label: str,
    title: str,
    output_path: Path,
    show_genotype: bool = False,
) -> None:
    """Save a consistently formatted two-dimensional sample scatter plot."""

    sns.set_theme(style="whitegrid", context="notebook")
    figure, axis = plt.subplots(figsize=(8.5, 6.5))
    plot_options = {
        "data": coordinates,
        "x": x,
        "y": y,
        "hue": "environment",
        "hue_order": ["Ground", "Spaceflight"],
        "palette": ENVIRONMENT_PALETTE,
        "s": 95,
        "edgecolor": "white",
        "linewidth": 0.8,
        "alpha": 0.9,
        "ax": axis,
    }
    if show_genotype:
        plot_options.update(
            {
                "style": "genotype",
                "style_order": ["Col-0", "WS"],
                "markers": GENOTYPE_MARKERS,
            }
        )
    sns.scatterplot(**plot_options)
    axis.set_title(title, fontsize=14, pad=14)
    axis.set_xlabel(x_label)
    axis.set_ylabel(y_label)
    axis.legend(title="Environment" if not show_genotype else None, frameon=True)
    sns.despine(ax=axis)
    figure.tight_layout()
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(figure)


def print_embedding_observations(
    method: str, scores: dict[str, float]
) -> None:
    print(
        f"{method}: Ground vs Spaceflight shows "
        f"{qualitative_grouping(scores['environment'])} "
        f"(descriptive silhouette={scores['environment']:.3f}); "
        f"age shows {qualitative_grouping(scores['age'])} "
        f"({scores['age']:.3f}); genotype shows "
        f"{qualitative_grouping(scores['genotype'])} "
        f"({scores['genotype']:.3f})."
    )


def main() -> int:
    args = parse_args()
    print("Loading processed SRP100712 data...")
    expression = load_processed_expression(args.expression)
    metadata = load_clean_metadata(args.metadata)
    validation = validate_processed_data(expression, metadata)
    sample_metadata = get_biological_sample_metadata(metadata)

    matrices = prepare_sample_matrices(expression, top_n=2_000)
    if matrices.all_nonzero_variance.shape[0] != 32:
        raise ValueError("Matrix preparation did not preserve all 32 samples.")

    print("Running PCA on all non-zero-variance genes...")
    pca_raw, explained_variance = run_pca(matrices.all_nonzero_variance)
    print("Running t-SNE on the 2,000 most variable genes...")
    tsne_raw = run_tsne(
        matrices.top_variable,
        perplexity=5,
        random_state=4144,
    )
    print("Running UMAP on the same 2,000 most variable genes...")
    umap_raw = run_umap(
        matrices.top_variable,
        n_neighbors=8,
        min_dist=0.2,
        random_state=4144,
    )

    pca_coordinates = add_sample_metadata(pca_raw, sample_metadata)
    tsne_coordinates = add_sample_metadata(tsne_raw, sample_metadata)
    umap_coordinates = add_sample_metadata(umap_raw, sample_metadata)

    args.plots_dir.mkdir(parents=True, exist_ok=True)
    args.tables_dir.mkdir(parents=True, exist_ok=True)

    pca_plot = args.plots_dir / "PCA_flight_vs_ground.png"
    pca_genotype_plot = args.plots_dir / "PCA_flight_vs_ground_genotype.png"
    tsne_plot = args.plots_dir / "tSNE_flight_vs_ground.png"
    umap_plot = args.plots_dir / "UMAP_flight_vs_ground.png"
    pca_table = args.tables_dir / "pca_coordinates.csv"
    tsne_table = args.tables_dir / "tsne_coordinates.csv"
    umap_table = args.tables_dir / "umap_coordinates.csv"

    pc1_percent = 100 * explained_variance[0]
    pc2_percent = 100 * explained_variance[1]
    save_embedding_plot(
        pca_coordinates,
        "PC1",
        "PC2",
        f"PC1 ({pc1_percent:.2f}% variance)",
        f"PC2 ({pc2_percent:.2f}% variance)",
        "PCA of Ground and Spaceflight Arabidopsis Samples",
        pca_plot,
    )
    save_embedding_plot(
        pca_coordinates,
        "PC1",
        "PC2",
        f"PC1 ({pc1_percent:.2f}% variance)",
        f"PC2 ({pc2_percent:.2f}% variance)",
        "PCA of Ground and Spaceflight Samples by Genotype",
        pca_genotype_plot,
        show_genotype=True,
    )
    save_embedding_plot(
        tsne_coordinates,
        "tSNE1",
        "tSNE2",
        "t-SNE 1",
        "t-SNE 2",
        "t-SNE of Ground and Spaceflight Arabidopsis Samples",
        tsne_plot,
    )
    save_embedding_plot(
        umap_coordinates,
        "UMAP1",
        "UMAP2",
        "UMAP 1",
        "UMAP 2",
        "UMAP of Ground and Spaceflight Arabidopsis Samples",
        umap_plot,
    )

    pca_coordinates.to_csv(pca_table, index=False)
    tsne_coordinates.to_csv(tsne_table, index=False)
    umap_coordinates.to_csv(umap_table, index=False)

    pca_scores = embedding_silhouette_scores(
        pca_coordinates, ("PC1", "PC2")
    )
    tsne_scores = embedding_silhouette_scores(
        tsne_coordinates, ("tSNE1", "tSNE2")
    )
    umap_scores = embedding_silhouette_scores(
        umap_coordinates, ("UMAP1", "UMAP2")
    )

    print("\n=== VALIDATION AND MATRIX PREPARATION ===")
    print(
        f"Input matrix: {validation['gene_count']:,} genes x "
        f"{validation['sample_count']} biological samples"
    )
    print(
        "Expression sample IDs match metadata biological sample IDs: "
        f"{'YES' if validation['identifiers_match'] else 'NO'}"
    )
    print(f"Missing expression values: {validation['missing_values']}")
    print(f"Infinite expression values: {validation['infinite_values']}")
    print(f"Zero-variance genes removed: {matrices.zero_variance_gene_count:,}")
    print(
        "Genes retained after zero-variance filtering: "
        f"{matrices.all_nonzero_variance.shape[1]:,}"
    )
    print(f"Genes used for t-SNE and UMAP: {matrices.top_variable.shape[1]:,}")
    print("Scale handling: centered by sklearn PCA; no gene-wise standardization")

    print("\n=== PCA VARIANCE EXPLAINED ===")
    print(f"PC1: {pc1_percent:.3f}%")
    print(f"PC2: {pc2_percent:.3f}%")
    print(f"PC1 + PC2 cumulative: {pc1_percent + pc2_percent:.3f}%")

    print("\n=== EXPLORATORY CLUSTERING OBSERVATIONS ===")
    print_embedding_observations("PCA", pca_scores)
    print_embedding_observations("t-SNE", tsne_scores)
    print_embedding_observations("UMAP", umap_scores)
    print(
        "These observations and silhouette values describe the displayed "
        "coordinates only; they are not tests of statistical significance."
    )
    print(
        "Across the three plots, genotype is the strongest organizing feature. "
        "Age-related structure is clearest in PCA and weaker in the nonlinear "
        "embeddings. Ground and Spaceflight do not form two single global "
        "clusters, although they often occupy nearby local subgroups within the "
        "larger genotype- and age-related structure."
    )

    print("\n=== ASSIGNMENT-READY METHOD COMPARISON ===")
    print(
        "PCA gives a linear, global summary of variation and quantifies how much "
        "variance PC1 and PC2 explain. t-SNE emphasizes local sample neighborhoods "
        "among the 2,000 most variable genes, but its axes and between-cluster "
        "distances have no direct variance interpretation. UMAP also captures "
        "nonlinear neighborhoods and often retains more broad arrangement than "
        "t-SNE, although its axes remain arbitrary. Agreement across methods is "
        "useful exploratory evidence of structure. Here, all three methods show "
        "stronger genotype-related organization than global environment-related "
        "separation, while PCA also displays age-related structure. Visual "
        "separation alone does not demonstrate statistical significance."
    )

    print("\nFiles written:")
    for path in (
        pca_table,
        tsne_table,
        umap_table,
        pca_plot,
        pca_genotype_plot,
        tsne_plot,
        umap_plot,
    ):
        print(f"  {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
