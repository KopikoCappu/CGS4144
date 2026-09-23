"""Heatmap preparation and clustering utilities for SRP100712."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import leaves_list, linkage
from sklearn.metrics import silhouette_score


EXPECTED_SIGNIFICANT_GENES = 2_560
REQUIRED_DE_COLUMNS = {
    "Gene",
    "spaceflight_effect",
    "standard_error",
    "t_statistic",
    "p_value",
    "adjusted_p_value",
    "significant",
}


@dataclass(frozen=True)
class HeatmapClustering:
    """Linkage matrices and ordered sample identifiers for one heatmap."""

    row_linkage: np.ndarray
    column_linkage: np.ndarray
    ordered_samples: list[str]


def validate_significant_genes(
    significant: pd.DataFrame,
    expression: pd.DataFrame,
    expected_count: int = EXPECTED_SIGNIFICANT_GENES,
) -> pd.DataFrame:
    """Validate and sort the adjusted differential-expression gene table."""

    missing_columns = sorted(REQUIRED_DE_COLUMNS.difference(significant.columns))
    if missing_columns:
        raise ValueError(
            "Significant-gene table is missing columns: "
            + ", ".join(missing_columns)
        )
    if len(significant) != expected_count:
        raise ValueError(
            f"Expected {expected_count:,} significant genes, found "
            f"{len(significant):,}."
        )
    if significant["Gene"].isna().any():
        raise ValueError("Significant-gene table contains missing gene IDs.")
    if significant["Gene"].duplicated().any():
        raise ValueError("Significant-gene table contains duplicate gene IDs.")
    if not significant["significant"].astype(bool).all():
        raise ValueError("Significant-gene table contains rows not marked significant.")
    if not (significant["adjusted_p_value"] < 0.05).all():
        raise ValueError("Significant-gene table contains adjusted p-values >= 0.05.")

    missing_genes = pd.Index(significant["Gene"]).difference(expression.index)
    if len(missing_genes):
        raise ValueError(
            f"{len(missing_genes)} significant genes are absent from expression: "
            + ", ".join(missing_genes.astype(str)[:10])
        )

    return significant.sort_values(
        ["adjusted_p_value", "p_value", "Gene"], kind="stable"
    ).reset_index(drop=True)


def validate_top50_ranking(
    top50: pd.DataFrame, ranked_significant: pd.DataFrame
) -> None:
    """Confirm the saved top-50 table agrees with the significant-gene ranking."""

    if "Gene" not in top50.columns:
        raise ValueError("Top-50 differential-expression table lacks a Gene column.")
    if len(top50) != 50:
        raise ValueError(f"Expected 50 rows in the top-50 table, found {len(top50)}.")
    expected = ranked_significant["Gene"].head(50).astype(str).tolist()
    observed = top50["Gene"].astype(str).tolist()
    if observed != expected:
        raise ValueError(
            "Top-50 gene order does not match the significant table ranked by "
            "adjusted p-value."
        )


def zscore_genes(expression: pd.DataFrame) -> pd.DataFrame:
    """Center and scale each gene across samples using population SD (ddof=0)."""

    means = expression.mean(axis=1)
    standard_deviations = expression.std(axis=1, ddof=0)
    zero_variance = standard_deviations.eq(0)
    if zero_variance.any():
        genes = standard_deviations.index[zero_variance].astype(str).tolist()
        raise ValueError(
            "Cannot z-score zero-variance significant genes without removing "
            "them: "
            + ", ".join(genes[:10])
        )

    zscores = expression.sub(means, axis=0).div(standard_deviations, axis=0)
    if not np.isfinite(zscores.to_numpy(dtype=float)).all():
        raise ValueError("Gene-wise z-scoring produced missing or infinite values.")
    return zscores


def calculate_heatmap_clustering(zscores: pd.DataFrame) -> HeatmapClustering:
    """Compute average-linkage Euclidean clustering for genes and samples."""

    row_linkage = linkage(
        zscores.to_numpy(dtype=float), method="average", metric="euclidean"
    )
    column_linkage = linkage(
        zscores.T.to_numpy(dtype=float), method="average", metric="euclidean"
    )
    sample_order = zscores.columns[leaves_list(column_linkage)].astype(str).tolist()
    return HeatmapClustering(
        row_linkage=row_linkage,
        column_linkage=column_linkage,
        ordered_samples=sample_order,
    )


def sample_grouping_scores(
    zscores: pd.DataFrame, sample_metadata: pd.DataFrame
) -> tuple[dict[str, float], dict[str, float]]:
    """Calculate descriptive silhouette scores globally and within strata."""

    metadata = sample_metadata.set_index("biological_sample").loc[zscores.columns]
    sample_values = zscores.T.to_numpy(dtype=float)
    global_scores = {
        field: float(silhouette_score(sample_values, metadata[field].astype(str)))
        for field in ("environment", "age", "genotype")
    }

    within_environment_scores: dict[str, float] = {}
    for (genotype, age), group in metadata.groupby(
        ["genotype", "age"], observed=True, sort=True
    ):
        group_values = zscores.loc[:, group.index].T.to_numpy(dtype=float)
        key = f"{genotype}, {age} days"
        within_environment_scores[key] = float(
            silhouette_score(group_values, group["environment"].astype(str))
        )
    return global_scores, within_environment_scores


def abbreviate_sample_names(sample_metadata: pd.DataFrame) -> dict[str, str]:
    """Create compact, unique column labels from the experimental design."""

    abbreviations: dict[str, str] = {}
    for row in sample_metadata.itertuples(index=False):
        environment = "G" if row.environment == "Ground" else "SF"
        genotype = "Col0" if row.genotype == "Col-0" else "WS"
        replicate = str(row.replicate).replace("Rep", "R")
        abbreviations[str(row.biological_sample)] = (
            f"{environment}_{row.age}d_{genotype}_{replicate}"
        )
    if len(set(abbreviations.values())) != len(abbreviations):
        raise ValueError("Abbreviated biological sample names are not unique.")
    return abbreviations
