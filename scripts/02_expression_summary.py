"""Describe expression values and gene variability in SRP100712."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import seaborn as sns  # noqa: E402


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.data.preprocess import (  # noqa: E402
    assess_expression_scale,
    calculate_gene_expression_ranges,
    load_aggregated_metadata,
    load_clean_metadata,
    load_processed_expression,
    summarize_expression_distribution,
    summarize_gene_expression_ranges,
    validate_processed_data,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Validate the 32-sample SRP100712 matrix and summarize expression "
            "variability without changing the normalized scale."
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
        "--aggregated-metadata",
        type=Path,
        default=(
            PROJECT_ROOT
            / "data"
            / "raw"
            / "SRP100712"
            / "aggregated_metadata.json"
        ),
        help="Refine.bio aggregated metadata JSON",
    )
    parser.add_argument(
        "--tables-dir",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables",
        help="Output directory for summary tables",
    )
    parser.add_argument(
        "--plots-dir",
        type=Path,
        default=PROJECT_ROOT / "plots",
        help="Output directory for plots",
    )
    return parser.parse_args()


def create_density_plot(ranges, output_path: Path) -> None:
    """Save a readable KDE density plot of per-gene expression ranges."""

    sns.set_theme(style="whitegrid", context="notebook")
    figure, axis = plt.subplots(figsize=(9, 6))
    sns.kdeplot(
        data=ranges,
        x="expression_range",
        fill=True,
        color="#2A6F97",
        linewidth=1.5,
        cut=0,
        ax=axis,
    )
    axis.set_title(
        "Distribution of Gene Expression Variability Across Biological Samples",
        fontsize=13,
        pad=14,
    )
    axis.set_xlabel("Gene expression range")
    axis.set_ylabel("Density")
    axis.set_xlim(left=0)
    sns.despine(ax=axis)
    figure.tight_layout()
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(figure)


def format_statistic(name: str, value: float) -> str:
    if name.startswith("proportion_"):
        return f"{value:.6f} ({value:.2%})"
    return f"{value:.6f}"


def main() -> int:
    args = parse_args()
    print("Loading processed SRP100712 data...")
    expression = load_processed_expression(args.expression)
    metadata = load_clean_metadata(args.metadata)
    aggregated_metadata = load_aggregated_metadata(args.aggregated_metadata)

    validation = validate_processed_data(expression, metadata)
    distribution = summarize_expression_distribution(expression)
    scale_assessment = assess_expression_scale(distribution, aggregated_metadata)
    if not scale_assessment["preserve_existing_scale"]:
        raise ValueError(
            "Preprocessing scale could not be validated: "
            + scale_assessment["evidence"]
        )

    # No transformation and no sample aggregation are performed here.
    gene_ranges = calculate_gene_expression_ranges(expression)
    variability_summary = summarize_gene_expression_ranges(gene_ranges)

    args.tables_dir.mkdir(parents=True, exist_ok=True)
    args.plots_dir.mkdir(parents=True, exist_ok=True)
    ranges_output = args.tables_dir / "gene_expression_ranges.csv"
    summary_output = args.tables_dir / "expression_variability_summary.csv"
    plot_output = args.plots_dir / "expression_variability_density.png"

    gene_ranges.to_csv(ranges_output, index=False)
    variability_summary.to_csv(summary_output, index=False)
    create_density_plot(gene_ranges, plot_output)

    print("\n=== VALIDATION ===")
    print(
        f"Expression matrix: {validation['gene_count']:,} genes x "
        f"{validation['sample_count']} biological samples"
    )
    print(
        "Unique metadata biological samples: "
        f"{validation['metadata_biological_sample_count']}"
    )
    print(
        "Expression sample IDs match metadata biological sample IDs: "
        f"{'YES' if validation['identifiers_match'] else 'NO'}"
    )
    print(f"Missing expression values: {validation['missing_values']}")
    print(f"Infinite expression values: {validation['infinite_values']}")

    print("\n=== EXPRESSION-VALUE DISTRIBUTION ===")
    for statistic, value in distribution.items():
        print(f"{statistic}: {format_statistic(statistic, value)}")

    print("\n=== PREPROCESSING DECISION ===")
    print(scale_assessment["decision"])
    print(scale_assessment["evidence"])
    print("No log2(x + 1), scaling, filtering, or sample aggregation was applied.")

    print("\n=== PER-GENE EXPRESSION-RANGE SUMMARY ===")
    for row in variability_summary.itertuples(index=False):
        print(f"{row.statistic}: {row.expression_range:.6f}")

    median_range = float(
        variability_summary.loc[
            variability_summary["statistic"].eq("median"), "expression_range"
        ].iloc[0]
    )
    upper_quartile = float(
        variability_summary.loc[
            variability_summary["statistic"].eq("75th_percentile"),
            "expression_range",
        ].iloc[0]
    )
    interpretation = (
        "Most genes have relatively low variability across the 32 samples: half "
        f"have an expression range at or below {median_range:.3f}, and 75% are "
        f"at or below {upper_quartile:.3f}. A smaller upper tail contains the "
        "more variable genes."
    )

    print("\n=== SHORT ASSIGNMENT SUMMARY ===")
    print(
        f"The validated expression matrix contains {expression.shape[0]:,} genes "
        f"and {expression.shape[1]} biological samples. The median gene expression "
        f"range is {median_range:.6f}. {interpretation}"
    )

    print("\nFiles written:")
    print(f"  {ranges_output}")
    print(f"  {summary_output}")
    print(f"  {plot_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
