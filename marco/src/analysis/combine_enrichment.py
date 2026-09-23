"""Combine GO enrichment results while preserving method-specific statistics."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd


WILCOXON_REQUIRED = {
    "GO_ID",
    "GO_term",
    "gene_set_size",
    "effect_difference",
    "p_value",
    "adjusted_p_value",
    "significant",
}
GPROFILER_REQUIRED = {
    "source",
    "native",
    "name",
    "p_value",
    "significant",
    "term_size",
    "intersection_size",
}


def _validate_columns(data: pd.DataFrame, required: set[str], label: str) -> None:
    missing = sorted(required.difference(data.columns))
    if missing:
        raise ValueError(f"{label} is missing columns: {', '.join(missing)}")


def _coerce_boolean(values: pd.Series, label: str) -> pd.Series:
    if pd.api.types.is_bool_dtype(values):
        return values.astype("boolean")
    normalized = values.astype("string").str.strip().str.lower()
    mapped = normalized.map({"true": True, "false": False})
    if mapped.isna().any():
        unexpected = sorted(normalized.loc[mapped.isna()].dropna().unique())
        raise ValueError(f"{label} has invalid significance values: {unexpected}")
    return mapped.astype("boolean")


def _validate_term_ids(data: pd.DataFrame, column: str, label: str) -> None:
    if data[column].isna().any():
        raise ValueError(f"{label} contains missing GO IDs.")
    if data[column].duplicated().any():
        raise ValueError(f"{label} contains duplicate GO IDs.")
    if not data[column].astype(str).str.fullmatch(r"GO:\d{7}").all():
        raise ValueError(f"{label} contains invalid GO IDs.")


def load_wilcoxon_results(path: str | Path) -> pd.DataFrame:
    data = pd.read_csv(path)
    _validate_columns(data, WILCOXON_REQUIRED, "Wilcoxon results")
    _validate_term_ids(data, "GO_ID", "Wilcoxon results")
    data = data.copy()
    data["significant"] = _coerce_boolean(
        data["significant"], "Wilcoxon results"
    )
    values = data[
        ["p_value", "adjusted_p_value", "effect_difference"]
    ].to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("Wilcoxon results contain non-finite statistics.")
    return data.rename(
        columns={
            "GO_ID": "term_id",
            "GO_term": "wilcoxon_term_name",
            "gene_set_size": "wilcoxon_gene_set_size",
            "effect_difference": "wilcoxon_effect_difference",
            "p_value": "wilcoxon_p_value",
            "adjusted_p_value": "wilcoxon_adjusted_p_value",
            "significant": "wilcoxon_significant",
        }
    )[
        [
            "term_id",
            "wilcoxon_term_name",
            "wilcoxon_gene_set_size",
            "wilcoxon_effect_difference",
            "wilcoxon_p_value",
            "wilcoxon_adjusted_p_value",
            "wilcoxon_significant",
        ]
    ]


def load_team_gprofiler_results(path: str | Path) -> pd.DataFrame:
    data = pd.read_csv(path)
    _validate_columns(data, GPROFILER_REQUIRED, "Teammate g:Profiler results")
    if not data["source"].eq("GO:BP").all():
        sources = sorted(data["source"].dropna().astype(str).unique())
        raise ValueError(
            f"Teammate g:Profiler results contain non-GO:BP sources: {sources}"
        )
    _validate_term_ids(data, "native", "Teammate g:Profiler results")
    data = data.copy()
    data["significant"] = _coerce_boolean(
        data["significant"], "Teammate g:Profiler results"
    )
    values = data[["p_value", "term_size", "intersection_size"]].to_numpy(
        dtype=float
    )
    if not np.isfinite(values).all():
        raise ValueError("Teammate g:Profiler results contain non-finite statistics.")
    return data.rename(
        columns={
            "native": "term_id",
            "name": "gprofiler_term_name",
            "term_size": "gprofiler_term_size",
            "intersection_size": "gprofiler_intersection_size",
            "p_value": "gprofiler_p_value",
            "significant": "gprofiler_significant",
        }
    )[
        [
            "term_id",
            "gprofiler_term_name",
            "gprofiler_term_size",
            "gprofiler_intersection_size",
            "gprofiler_p_value",
            "gprofiler_significant",
        ]
    ]


def combine_results(
    wilcoxon: pd.DataFrame, gprofiler: pd.DataFrame
) -> tuple[pd.DataFrame, int]:
    combined = wilcoxon.merge(
        gprofiler,
        on="term_id",
        how="outer",
        validate="one_to_one",
        indicator=True,
    )
    combined["wilcoxon_tested"] = combined["_merge"].isin(["left_only", "both"])
    combined["gprofiler_tested"] = combined["_merge"].isin(["right_only", "both"])

    both_names = combined[
        combined["wilcoxon_term_name"].notna()
        & combined["gprofiler_term_name"].notna()
    ]
    name_conflicts = int(
        (
            both_names["wilcoxon_term_name"].str.casefold()
            != both_names["gprofiler_term_name"].str.casefold()
        ).sum()
    )
    combined["term_name"] = combined["wilcoxon_term_name"].combine_first(
        combined["gprofiler_term_name"]
    )

    combined["methods_tested"] = (
        combined["wilcoxon_tested"].astype(int)
        + combined["gprofiler_tested"].astype(int)
    )
    combined["methods_significant"] = (
        combined["wilcoxon_significant"].fillna(False).astype(int)
        + combined["gprofiler_significant"].fillna(False).astype(int)
    )
    combined["significant_fraction"] = (
        combined["methods_significant"] / combined["methods_tested"]
    )

    columns = [
        "term_id",
        "term_name",
        "methods_tested",
        "methods_significant",
        "significant_fraction",
        "wilcoxon_tested",
        "gprofiler_tested",
        "wilcoxon_gene_set_size",
        "wilcoxon_effect_difference",
        "wilcoxon_p_value",
        "wilcoxon_adjusted_p_value",
        "wilcoxon_significant",
        "gprofiler_term_size",
        "gprofiler_intersection_size",
        "gprofiler_p_value",
        "gprofiler_significant",
    ]
    return combined[columns], name_conflicts


def rank_combined_results(combined: pd.DataFrame) -> pd.DataFrame:
    ranked = combined.copy()
    ranked["_strongest_available_p"] = ranked[
        ["wilcoxon_adjusted_p_value", "gprofiler_p_value"]
    ].min(axis=1, skipna=True)
    ranked = ranked.sort_values(
        [
            "methods_significant",
            "significant_fraction",
            "_strongest_available_p",
            "term_id",
        ],
        ascending=[False, False, True, True],
        kind="stable",
    )
    return ranked.drop(columns="_strongest_available_p").reset_index(drop=True)


def select_root_related(combined: pd.DataFrame) -> pd.DataFrame:
    keywords = [
        "root",
        "root system",
        "root meristem",
        "lateral root",
        "root hair",
        "auxin",
        "gravitropism",
        "gravity",
        "cell wall",
        "development",
        "morphogenesis",
    ]
    pattern = "|".join(re.escape(keyword) for keyword in keywords)
    return combined.loc[
        combined["term_name"].str.contains(pattern, case=False, na=False, regex=True)
    ].copy()
