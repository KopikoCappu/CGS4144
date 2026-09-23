"""Adjusted Spaceflight-vs-Ground differential-expression analysis."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402
from adjustText import adjust_text  # noqa: E402


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.analysis.differential_expression import (  # noqa: E402
    build_model_design,
    fit_gene_wise_ols,
    select_volcano_labels,
)
from src.data.preprocess import (  # noqa: E402
    get_biological_sample_metadata,
    load_clean_metadata,
    load_processed_expression,
    validate_processed_data,
)


FDR_THRESHOLD = 0.05


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Fit adjusted gene-wise OLS models for Spaceflight versus Ground "
            "in SRP100712."
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
        "--tables-dir",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables",
        help="Output directory for result tables",
    )
    parser.add_argument(
        "--plots-dir",
        type=Path,
        default=PROJECT_ROOT / "plots",
        help="Output directory for plots",
    )
    return parser.parse_args()


def create_volcano_plot(results: pd.DataFrame, output_path: Path) -> int:
    """Create the adjusted-p-value volcano plot and return label count."""

    plot_data = results.copy()
    safe_adjusted = np.clip(
        plot_data["adjusted_p_value"].to_numpy(dtype=float),
        np.finfo(float).tiny,
        1.0,
    )
    plot_data["negative_log10_fdr"] = -np.log10(safe_adjusted)

    sns.set_theme(style="whitegrid", context="notebook")
    figure, axis = plt.subplots(figsize=(10, 7))
    non_significant = plot_data.loc[~plot_data["significant"]]
    significant = plot_data.loc[plot_data["significant"]]

    axis.scatter(
        non_significant["spaceflight_effect"],
        non_significant["negative_log10_fdr"],
        s=18,
        color="#9E9E9E",
        alpha=0.55,
        edgecolors="none",
        label="Not significant",
        rasterized=True,
    )
    axis.scatter(
        significant["spaceflight_effect"],
        significant["negative_log10_fdr"],
        s=27,
        color="#C44E52",
        alpha=0.8,
        edgecolors="none",
        label="Significant (FDR < 0.05)",
        rasterized=True,
    )
    axis.axhline(
        -np.log10(FDR_THRESHOLD),
        color="#333333",
        linestyle="--",
        linewidth=1.2,
        label="FDR = 0.05",
    )
    axis.axvline(0, color="#666666", linestyle=":", linewidth=1.0)

    labels = select_volcano_labels(results, maximum_labels=10)
    texts = []
    plot_lookup = plot_data.set_index("Gene")
    for row in labels.itertuples(index=False):
        texts.append(
            axis.text(
                row.spaceflight_effect,
                plot_lookup.loc[row.Gene, "negative_log10_fdr"],
                row.Gene,
                fontsize=8,
                ha="center",
                va="bottom",
            )
        )
    if texts:
        adjust_text(
            texts,
            ax=axis,
            arrowprops={"arrowstyle": "-", "color": "#555555", "lw": 0.6},
        )

    axis.set_title(
        "Adjusted Spaceflight vs Ground Gene Expression Effects",
        fontsize=14,
        pad=14,
    )
    axis.set_xlabel("Spaceflight effect (adjusted normalized-expression difference)")
    axis.set_ylabel("−log10(adjusted p-value)")
    axis.legend(frameon=True)
    sns.despine(ax=axis)
    figure.tight_layout()
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(figure)
    return len(labels)


def main() -> int:
    args = parse_args()
    warnings: list[str] = []

    print("Loading processed SRP100712 data...")
    expression = load_processed_expression(args.expression)
    metadata = load_clean_metadata(args.metadata)
    validation = validate_processed_data(expression, metadata)
    sample_metadata = get_biological_sample_metadata(metadata)
    model_design = build_model_design(sample_metadata, expression.columns)

    print("Fitting 32,309 separate adjusted OLS models...")
    results = fit_gene_wise_ols(
        expression,
        model_design,
        fdr_threshold=FDR_THRESHOLD,
    )
    significant = results.loc[results["significant"]].copy()
    top_50 = results.head(50).copy()

    args.tables_dir.mkdir(parents=True, exist_ok=True)
    args.plots_dir.mkdir(parents=True, exist_ok=True)
    all_output = args.tables_dir / "differential_expression_all.csv"
    significant_output = (
        args.tables_dir / "differential_expression_significant.csv"
    )
    top50_output = args.tables_dir / "differential_expression_top50.csv"
    volcano_output = args.plots_dir / "volcano_spaceflight_vs_ground.png"

    results.to_csv(all_output, index=False)
    significant.to_csv(significant_output, index=False)
    top_50.to_csv(top50_output, index=False)
    label_count = create_volcano_plot(results, volcano_output)

    positive_count = int((results["spaceflight_effect"] > 0).sum())
    negative_count = int((results["spaceflight_effect"] < 0).sum())
    zero_count = int((results["spaceflight_effect"] == 0).sum())
    significant_positive = int(
        (significant["spaceflight_effect"] > 0).sum()
    )
    significant_negative = int(
        (significant["spaceflight_effect"] < 0).sum()
    )
    significant_zero = int((significant["spaceflight_effect"] == 0).sum())

    if zero_count:
        warnings.append(f"{zero_count} genes have an exactly zero Spaceflight effect.")
    if significant_zero:
        warnings.append(
            f"{significant_zero} significant genes have an exactly zero effect."
        )
    if (results["adjusted_p_value"] == 0).any():
        warnings.append(
            "Some adjusted p-values equal zero numerically; plotting clips them "
            "to the smallest positive float before -log10 conversion."
        )

    print("\n=== DATA AND MODEL VALIDATION ===")
    print(
        f"Expression matrix: {validation['gene_count']:,} genes x "
        f"{validation['sample_count']} biological samples"
    )
    print(
        "Expression sample IDs match metadata biological sample IDs: "
        f"{'YES' if validation['identifiers_match'] else 'NO'}"
    )
    print(f"Missing expression values: {validation['missing_values']}")
    print(f"Infinite expression values: {validation['infinite_values']}")
    print(f"Model formula: expression ~ {model_design.formula}")
    print("Design matrix columns:")
    for column in model_design.matrix.columns:
        print(f"  {column}")
    print(
        "Environment coefficient reported: "
        f"{model_design.environment_coefficient}"
    )
    print("Ground design value: 0; Spaceflight design value: 1")
    print("Confirmed contrast: Spaceflight - Ground (adjusted for age and genotype)")

    print("\n=== DIFFERENTIAL-EXPRESSION SUMMARY ===")
    print(f"Total genes tested: {len(results):,}")
    print(f"Significant at FDR < 0.05: {len(significant):,}")
    print(f"Positive Spaceflight effect: {positive_count:,}")
    print(f"Negative Spaceflight effect: {negative_count:,}")
    print(f"Significant with positive effect: {significant_positive:,}")
    print(f"Significant with negative effect: {significant_negative:,}")
    print(f"Minimum raw p-value: {results['p_value'].min():.6e}")
    print(
        "Minimum adjusted p-value: "
        f"{results['adjusted_p_value'].min():.6e}"
    )
    print(
        "Largest positive Spaceflight effect: "
        f"{results['spaceflight_effect'].max():.6f}"
    )
    print(
        "Largest negative Spaceflight effect: "
        f"{results['spaceflight_effect'].min():.6f}"
    )
    print(f"Genes labeled on volcano plot: {label_count}")

    print("\nFive most statistically significant genes:")
    print(
        results.loc[
            :4,
            ["Gene", "spaceflight_effect", "adjusted_p_value"],
        ].to_string(index=False)
    )

    print("\n=== ASSIGNMENT-READY INTERPRETATION ===")
    print(
        f"A separate ordinary least squares model was fitted for each of "
        f"{len(results):,} genes. The model estimated the normalized-expression "
        "difference between Spaceflight and Ground while including age and "
        "genotype as categorical covariates. Benjamini-Hochberg correction was "
        f"applied across all genes, identifying {len(significant):,} genes at "
        f"FDR < 0.05 ({significant_positive:,} with positive and "
        f"{significant_negative:,} with negative Spaceflight effects). Positive "
        "effects indicate higher adjusted normalized expression in Spaceflight; "
        "negative effects indicate lower adjusted normalized expression. The "
        "effect is not labeled a log2 fold change because the exact Refine.bio "
        "transformation is not established. Differential expression describes "
        "association with environment and does not prove a causal role in root "
        "development."
    )

    print("\nWarnings:")
    if warnings:
        for warning in warnings:
            print(f"  - {warning}")
    else:
        print("  None")

    print("\nFiles written:")
    for path in (
        all_output,
        significant_output,
        top50_output,
        volcano_output,
    ):
        print(f"  {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
