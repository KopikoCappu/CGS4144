"""Validation and descriptive preprocessing for the 32-sample dataset."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from pandas.api.types import is_numeric_dtype


EXPECTED_GENES = 32_309
EXPECTED_BIOLOGICAL_SAMPLES = 32


def load_processed_expression(path: str | Path) -> pd.DataFrame:
    """Load a processed genes-by-biological-samples expression CSV."""

    expression = pd.read_csv(Path(path), index_col="Gene")
    expression.index.name = "Gene"

    non_numeric = [
        column
        for column in expression.columns
        if not is_numeric_dtype(expression[column])
    ]
    if non_numeric:
        preview = ", ".join(map(str, non_numeric[:10]))
        raise ValueError(f"Non-numeric expression columns found: {preview}")
    return expression


def load_clean_metadata(path: str | Path) -> pd.DataFrame:
    """Load run-level clean metadata created by ``01_inspect_data.py``."""

    metadata = pd.read_csv(Path(path), dtype={"run_id": str, "age": int})
    required = {
        "run_id",
        "biological_sample",
        "environment",
        "age",
        "genotype",
        "replicate",
    }
    missing = sorted(required.difference(metadata.columns))
    if missing:
        raise ValueError(
            "Clean metadata is missing required columns: " + ", ".join(missing)
        )
    return metadata


def load_aggregated_metadata(path: str | Path) -> dict[str, Any]:
    """Load Refine.bio's aggregated metadata JSON."""

    with Path(path).open(encoding="utf-8") as handle:
        return json.load(handle)


def get_biological_sample_metadata(metadata: pd.DataFrame) -> pd.DataFrame:
    """Return one internally consistent row per biological sample."""

    condition_columns = ["environment", "age", "genotype", "replicate"]
    inconsistent = (
        metadata.groupby("biological_sample", dropna=False)[condition_columns]
        .nunique(dropna=False)
        .gt(1)
        .any(axis=1)
    )
    if inconsistent.any():
        names = inconsistent.index[inconsistent].astype(str).tolist()
        raise ValueError(
            "Inconsistent condition fields within biological samples: "
            + ", ".join(names[:10])
        )

    return (
        metadata[["biological_sample", *condition_columns]]
        .drop_duplicates()
        .reset_index(drop=True)
    )


def validate_processed_data(
    expression: pd.DataFrame, metadata: pd.DataFrame
) -> dict[str, Any]:
    """Validate dimensions, identifiers, and missing values.

    Raises ``ValueError`` when a validation needed for the descriptive analysis
    fails. Returns details used in the printed validation report.
    """

    problems: list[str] = []
    sample_metadata = get_biological_sample_metadata(metadata)

    if expression.shape[0] != EXPECTED_GENES:
        problems.append(
            f"expected {EXPECTED_GENES:,} genes, found {expression.shape[0]:,}"
        )
    if expression.shape[1] != EXPECTED_BIOLOGICAL_SAMPLES:
        problems.append(
            "expected 32 expression sample columns, "
            f"found {expression.shape[1]}"
        )
    if len(sample_metadata) != EXPECTED_BIOLOGICAL_SAMPLES:
        problems.append(
            "expected 32 unique metadata biological samples, "
            f"found {len(sample_metadata)}"
        )
    if expression.index.isna().any():
        problems.append(f"found {int(expression.index.isna().sum())} missing gene IDs")
    if expression.index.duplicated().any():
        problems.append(
            f"found {int(expression.index.duplicated().sum())} duplicate gene IDs"
        )
    if expression.columns.isna().any():
        problems.append("found missing expression sample IDs")
    if expression.columns.duplicated().any():
        problems.append(
            f"found {int(expression.columns.duplicated().sum())} duplicate sample IDs"
        )
    if metadata["biological_sample"].isna().any():
        problems.append(
            "found "
            f"{int(metadata['biological_sample'].isna().sum())} missing metadata "
            "biological sample IDs"
        )

    expression_ids = pd.Index(expression.columns.astype(str))
    metadata_ids = pd.Index(sample_metadata["biological_sample"].astype(str))
    expression_only = expression_ids.difference(metadata_ids)
    metadata_only = metadata_ids.difference(expression_ids)
    identifiers_match = not len(expression_only) and not len(metadata_only)
    if not identifiers_match:
        problems.append(
            f"sample ID mismatch ({len(expression_only)} expression-only; "
            f"{len(metadata_only)} metadata-only)"
        )

    missing_values = int(expression.isna().to_numpy().sum())
    infinite_values = int(np.isinf(expression.to_numpy(dtype=np.float64)).sum())
    if missing_values:
        problems.append(f"found {missing_values:,} missing expression values")
    if infinite_values:
        problems.append(f"found {infinite_values:,} infinite expression values")

    if problems:
        raise ValueError("Processed data validation failed: " + "; ".join(problems))

    return {
        "gene_count": expression.shape[0],
        "sample_count": expression.shape[1],
        "metadata_biological_sample_count": len(sample_metadata),
        "identifiers_match": identifiers_match,
        "expression_only_ids": expression_only.tolist(),
        "metadata_only_ids": metadata_only.tolist(),
        "missing_values": missing_values,
        "infinite_values": infinite_values,
    }


def summarize_expression_distribution(expression: pd.DataFrame) -> pd.Series:
    """Calculate global descriptive statistics across all matrix values."""

    values = expression.to_numpy(dtype=np.float64, copy=False)
    percentiles = np.percentile(values, [1, 5, 25, 50, 75, 95, 99])
    return pd.Series(
        {
            "minimum": float(values.min()),
            "1st_percentile": float(percentiles[0]),
            "5th_percentile": float(percentiles[1]),
            "25th_percentile": float(percentiles[2]),
            "median": float(percentiles[3]),
            "mean": float(values.mean()),
            "75th_percentile": float(percentiles[4]),
            "95th_percentile": float(percentiles[5]),
            "99th_percentile": float(percentiles[6]),
            "maximum": float(values.max()),
            "proportion_negative": float(np.mean(values < 0)),
            "proportion_zero": float(np.mean(values == 0)),
        },
        name="value",
    )


def assess_expression_scale(
    distribution: pd.Series, aggregated_metadata: dict[str, Any]
) -> dict[str, Any]:
    """Document whether another expression transformation is justified."""

    quantile_normalized = aggregated_metadata.get("quantile_normalized")
    scale_by = aggregated_metadata.get("scale_by")
    has_negative_values = distribution["proportion_negative"] > 0

    preserve_existing_scale = bool(quantile_normalized and has_negative_values)
    if preserve_existing_scale:
        decision = (
            "Use the existing Refine.bio normalized scale without an additional "
            "transformation."
        )
        evidence = (
            "aggregated_metadata.json reports quantile_normalized=true and "
            f"scale_by={scale_by!r}; the matrix contains negative values "
            f"({distribution['proportion_negative']:.2%} of observations). "
            "These are processed continuous values rather than raw counts. The "
            "metadata does not identify the exact upstream mathematical transform, "
            "so a second log2(x + 1) transform is not justified."
        )
    else:
        decision = "The available evidence does not establish the expected scale."
        evidence = (
            f"quantile_normalized={quantile_normalized!r}, scale_by={scale_by!r}, "
            f"proportion_negative={distribution['proportion_negative']:.6f}."
        )

    return {
        "preserve_existing_scale": preserve_existing_scale,
        "decision": decision,
        "evidence": evidence,
        "quantile_normalized": quantile_normalized,
        "scale_by": scale_by,
    }


def calculate_gene_expression_ranges(expression: pd.DataFrame) -> pd.DataFrame:
    """Calculate maximum minus minimum expression for every gene."""

    ranges = expression.max(axis=1) - expression.min(axis=1)
    return ranges.rename("expression_range").reset_index()


def summarize_gene_expression_ranges(ranges: pd.DataFrame) -> pd.DataFrame:
    """Return the requested descriptive statistics for per-gene ranges."""

    values = ranges["expression_range"].to_numpy(dtype=np.float64, copy=False)
    summary = [
        ("minimum", float(values.min())),
        ("25th_percentile", float(np.percentile(values, 25))),
        ("median", float(np.median(values))),
        ("mean", float(values.mean())),
        ("75th_percentile", float(np.percentile(values, 75))),
        ("maximum", float(values.max())),
    ]
    return pd.DataFrame(summary, columns=["statistic", "expression_range"])
