"""Validate and aggregate the Refine.bio SRP100712 expression dataset."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.data.load_data import (  # noqa: E402
    aggregate_runs_by_biological_sample,
    biological_sample_metadata,
    build_clean_metadata,
    load_aggregated_metadata,
    load_expression,
    load_metadata,
)


EXPECTED_GENES = 32_309
EXPECTED_RUNS = 594
EXPECTED_BIOLOGICAL_SAMPLES = 32


def parse_args() -> argparse.Namespace:
    raw_dir = PROJECT_ROOT / "data" / "raw" / "SRP100712"
    parser = argparse.ArgumentParser(
        description=(
            "Inspect SRP100712, parse sample conditions, and average SRR runs "
            "belonging to each biological sample."
        )
    )
    parser.add_argument(
        "--expression",
        type=Path,
        default=raw_dir / "SRP100712.tsv",
        help="Refine.bio expression TSV",
    )
    parser.add_argument(
        "--metadata",
        type=Path,
        default=raw_dir / "metadata_SRP100712.tsv",
        help="Refine.bio sample metadata TSV",
    )
    parser.add_argument(
        "--aggregated-metadata",
        type=Path,
        default=raw_dir / "aggregated_metadata.json",
        help="Refine.bio aggregated metadata JSON",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=PROJECT_ROOT / "data" / "processed",
        help="Directory for validated CSV outputs",
    )
    return parser.parse_args()


def yes_no(value: bool) -> str:
    return "YES" if value else "NO"


def print_series(series: pd.Series) -> None:
    print(series.to_string())


def main() -> int:
    args = parse_args()
    warnings: list[str] = []

    print("Loading SRP100712 files...")
    expression = load_expression(args.expression)
    metadata = load_metadata(args.metadata)
    aggregate_metadata = load_aggregated_metadata(args.aggregated_metadata)
    clean_metadata = build_clean_metadata(metadata)
    sample_metadata = biological_sample_metadata(clean_metadata)

    expression_ids = pd.Index(expression.columns.astype(str))
    metadata_ids = pd.Index(clean_metadata["run_id"].astype(str))
    expression_only = expression_ids.difference(metadata_ids)
    metadata_only = metadata_ids.difference(expression_ids)
    identifiers_align = not len(expression_only) and not len(metadata_only)
    identifier_order_aligns = expression_ids.equals(metadata_ids)

    missing_gene_ids = int(expression.index.isna().sum())
    duplicate_gene_ids = int(expression.index.duplicated().sum())
    missing_run_ids = int(clean_metadata["run_id"].isna().sum())
    duplicate_run_ids = int(clean_metadata["run_id"].duplicated().sum())
    missing_sample_ids = int(clean_metadata["biological_sample"].isna().sum())
    missing_parsed_fields = int(
        clean_metadata[["environment", "age", "genotype", "replicate"]]
        .isna()
        .any(axis=1)
        .sum()
    )
    missing_expression_values = int(expression.isna().to_numpy().sum())

    gene_count = len(expression)
    run_count = expression.shape[1]
    metadata_run_count = len(metadata)
    biological_sample_count = len(sample_metadata)

    expected_checks = [
        (gene_count == EXPECTED_GENES, f"Expected {EXPECTED_GENES:,} genes; found {gene_count:,}."),
        (run_count == EXPECTED_RUNS, f"Expected {EXPECTED_RUNS} expression runs; found {run_count}."),
        (
            metadata_run_count == EXPECTED_RUNS,
            f"Expected {EXPECTED_RUNS} metadata rows; found {metadata_run_count}.",
        ),
        (
            biological_sample_count == EXPECTED_BIOLOGICAL_SAMPLES,
            "Expected 32 biological samples; "
            f"found {biological_sample_count}.",
        ),
        (missing_gene_ids == 0, f"Found {missing_gene_ids} missing gene IDs."),
        (duplicate_gene_ids == 0, f"Found {duplicate_gene_ids} duplicate gene IDs."),
        (missing_run_ids == 0, f"Found {missing_run_ids} missing run IDs."),
        (duplicate_run_ids == 0, f"Found {duplicate_run_ids} duplicate run IDs."),
        (
            missing_sample_ids == 0,
            f"Found {missing_sample_ids} missing biological sample names.",
        ),
        (
            missing_parsed_fields == 0,
            f"Found {missing_parsed_fields} rows with missing parsed condition fields.",
        ),
        (
            missing_expression_values == 0,
            f"Found {missing_expression_values} missing expression values.",
        ),
        (
            identifiers_align,
            f"Run ID mismatch: {len(expression_only)} expression-only and "
            f"{len(metadata_only)} metadata-only.",
        ),
    ]
    warnings.extend(message for passed, message in expected_checks if not passed)

    json_sample_count = aggregate_metadata.get("num_samples")
    quantile_normalized = aggregate_metadata.get("quantile_normalized")
    scale_by = aggregate_metadata.get("scale_by")
    if json_sample_count != metadata_run_count:
        warnings.append(
            f"JSON num_samples is {json_sample_count!r}, but metadata has "
            f"{metadata_run_count} rows."
        )
    if quantile_normalized is not True:
        warnings.append(
            "JSON does not mark this dataset as quantile normalized "
            f"(value: {quantile_normalized!r})."
        )

    condition_counts = (
        sample_metadata.groupby(["environment", "age", "genotype"], sort=True)
        .size()
        .rename("biological_replicates")
    )
    expected_conditions = pd.MultiIndex.from_product(
        [["Ground", "Spaceflight"], [4, 8], ["Col-0", "WS"]],
        names=["environment", "age", "genotype"],
    )
    condition_counts = condition_counts.reindex(expected_conditions, fill_value=0)
    four_replicates_per_condition = bool((condition_counts == 4).all())
    if not four_replicates_per_condition:
        warnings.append(
            "At least one environment x age x genotype condition does not have "
            "exactly four biological replicates."
        )

    run_counts = (
        clean_metadata.groupby("biological_sample", sort=True)
        .size()
        .rename("SRR_runs")
    )

    # The JSON documents normalized run-level values and supplies no alternate
    # technical-run aggregation instruction, so use an equal-weight mean.
    aggregated = aggregate_runs_by_biological_sample(expression, clean_metadata)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    metadata_output = args.output_dir / "metadata_clean.csv"
    expression_output = args.output_dir / "expression_32_samples.csv"
    clean_metadata.to_csv(metadata_output, index=False)
    aggregated.to_csv(expression_output, index=True, index_label="Gene")

    print("\n=== SRP100712 VALIDATION REPORT ===")
    print(
        f"Original expression table: {expression.shape[0]:,} rows x "
        f"{expression.shape[1] + 1:,} columns (Gene + {run_count} SRR runs)"
    )
    print(
        f"Original numeric expression matrix: {expression.shape[0]:,} genes x "
        f"{run_count} SRR runs"
    )
    print(
        f"Original metadata table: {metadata.shape[0]:,} rows x "
        f"{metadata.shape[1]:,} columns"
    )
    print(f"Number of genes: {gene_count:,}")
    print(f"SRR runs (expression / metadata): {run_count} / {metadata_run_count}")
    print(f"Unique biological samples: {biological_sample_count}")
    print(f"Missing gene IDs: {missing_gene_ids}")
    print(f"Duplicate gene IDs: {duplicate_gene_ids}")
    print(f"Missing run IDs: {missing_run_ids}")
    print(f"Duplicate run IDs: {duplicate_run_ids}")
    print(f"Missing biological sample names: {missing_sample_ids}")
    print(f"Rows with missing parsed fields: {missing_parsed_fields}")
    print(f"Missing expression values: {missing_expression_values}")
    print(f"Expression and metadata ID sets align: {yes_no(identifiers_align)}")
    print(f"Expression and metadata ID order aligns: {yes_no(identifier_order_aligns)}")
    print(f"Expression-only IDs: {len(expression_only)}")
    print(f"Metadata-only IDs: {len(metadata_only)}")

    print("\nAggregated metadata JSON:")
    print(f"  num_samples: {json_sample_count}")
    print(f"  quantile_normalized: {quantile_normalized}")
    print(f"  scale_by: {scale_by}")
    print("  run aggregation used: arithmetic mean by refinebio_title")

    print("\nSRR runs per biological sample:")
    print_series(run_counts)

    print("\nFinal aggregated expression matrix:")
    print(f"  {aggregated.shape[0]:,} genes x {aggregated.shape[1]} biological samples")

    print("\nBiological sample counts by environment:")
    print_series(sample_metadata["environment"].value_counts().sort_index())
    print("\nBiological sample counts by age:")
    print_series(sample_metadata["age"].value_counts().sort_index())
    print("\nBiological sample counts by genotype:")
    print_series(sample_metadata["genotype"].value_counts().sort_index())
    print("\nBiological samples by environment x age x genotype:")
    print_series(condition_counts)
    print(
        "\nEvery environment x age x genotype condition has 4 replicates: "
        f"{yes_no(four_replicates_per_condition)}"
    )

    print("\nValidation warnings:")
    if warnings:
        for warning in warnings:
            print(f"  - {warning}")
    else:
        print("  None")

    print("\nFiles written:")
    print(f"  {metadata_output}")
    print(f"  {expression_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
