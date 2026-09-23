"""Loading and metadata preparation utilities for Refine.bio datasets."""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from pandas.api.types import is_numeric_dtype


TITLE_PATTERN = re.compile(
    r"^(?P<environment>Ground Control|Spaceflight), "
    r"(?P<age>4|8) days old, "
    r"(?P<genotype>Col-0|WS) "
    r"(?P<replicate>Rep[1-4])$"
)


def load_expression(path: str | Path) -> pd.DataFrame:
    """Load a genes-by-runs Refine.bio expression TSV.

    The returned dataframe is indexed by the ``Gene`` column. All remaining
    columns are required to contain numeric expression values.
    """

    path = Path(path)
    expression = pd.read_csv(path, sep="\t", index_col="Gene")
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


def load_metadata(path: str | Path) -> pd.DataFrame:
    """Load the Refine.bio sample metadata TSV as strings."""

    return pd.read_csv(Path(path), sep="\t", dtype=str)


def load_aggregated_metadata(path: str | Path) -> dict[str, Any]:
    """Load the Refine.bio aggregated metadata JSON."""

    with Path(path).open(encoding="utf-8") as handle:
        return json.load(handle)


def parse_refinebio_titles(titles: pd.Series) -> pd.DataFrame:
    """Parse SRP100712 biological condition fields from ``refinebio_title``."""

    parsed = titles.astype("string").str.strip().str.extract(TITLE_PATTERN)
    failed = parsed.isna().any(axis=1)
    if failed.any():
        bad_titles = titles.loc[failed].drop_duplicates().astype(str).tolist()
        preview = "; ".join(bad_titles[:10])
        raise ValueError(
            f"Could not parse {int(failed.sum())} refinebio_title rows: {preview}"
        )

    parsed["environment"] = parsed["environment"].replace(
        {"Ground Control": "Ground"}
    )
    parsed["age"] = parsed["age"].astype(int)
    return parsed


def build_clean_metadata(metadata: pd.DataFrame) -> pd.DataFrame:
    """Create run-level metadata with parsed biological sample information."""

    required = {"refinebio_accession_code", "refinebio_title"}
    missing = sorted(required.difference(metadata.columns))
    if missing:
        raise ValueError(f"Metadata is missing required columns: {', '.join(missing)}")

    parsed = parse_refinebio_titles(metadata["refinebio_title"])
    clean = pd.DataFrame(
        {
            "run_id": metadata["refinebio_accession_code"].astype("string").str.strip(),
            "biological_sample": metadata["refinebio_title"].astype("string").str.strip(),
            "environment": parsed["environment"],
            "age": parsed["age"],
            "genotype": parsed["genotype"],
            "replicate": parsed["replicate"],
        }
    )
    return clean


def biological_sample_metadata(clean_metadata: pd.DataFrame) -> pd.DataFrame:
    """Collapse run-level metadata to one row per biological sample."""

    sample_columns = [
        "biological_sample",
        "environment",
        "age",
        "genotype",
        "replicate",
    ]
    inconsistent = (
        clean_metadata.groupby("biological_sample", dropna=False)[sample_columns[1:]]
        .nunique(dropna=False)
        .gt(1)
        .any(axis=1)
    )
    if inconsistent.any():
        names = inconsistent.index[inconsistent].astype(str).tolist()
        raise ValueError(
            "Biological samples have inconsistent parsed metadata: "
            + ", ".join(names[:10])
        )

    return clean_metadata[sample_columns].drop_duplicates().reset_index(drop=True)


def aggregate_runs_by_biological_sample(
    expression: pd.DataFrame, clean_metadata: pd.DataFrame
) -> pd.DataFrame:
    """Average technical SRR run columns into biological sample columns.

    The biological sample order follows the first occurrence of each title in
    the expression file. The arithmetic mean gives every SRR run equal weight.
    """

    if clean_metadata["run_id"].isna().any():
        raise ValueError("Cannot aggregate because metadata contains missing run IDs.")
    if clean_metadata["run_id"].duplicated().any():
        raise ValueError("Cannot aggregate because metadata contains duplicate run IDs.")

    expression_ids = pd.Index(expression.columns.astype(str))
    metadata_by_run = clean_metadata.set_index("run_id")
    expression_only = expression_ids.difference(metadata_by_run.index)
    metadata_only = metadata_by_run.index.difference(expression_ids)
    if len(expression_only) or len(metadata_only):
        raise ValueError(
            "Expression and metadata run IDs do not align: "
            f"{len(expression_only)} expression-only and "
            f"{len(metadata_only)} metadata-only IDs."
        )

    ordered_samples = metadata_by_run.loc[expression_ids, "biological_sample"]
    sample_names = list(dict.fromkeys(ordered_samples.astype(str)))
    sample_codes = pd.Categorical(ordered_samples, categories=sample_names).codes
    if (sample_codes < 0).any():
        raise ValueError("Cannot aggregate because a biological sample name is missing.")

    run_counts = np.bincount(sample_codes, minlength=len(sample_names))
    weights = np.zeros((len(expression_ids), len(sample_names)), dtype=np.float64)
    weights[np.arange(len(expression_ids)), sample_codes] = 1.0 / run_counts[sample_codes]

    values = expression.to_numpy(dtype=np.float64, copy=False)
    aggregated_values = values @ weights
    return pd.DataFrame(
        aggregated_values,
        index=expression.index.copy(),
        columns=sample_names,
    )
