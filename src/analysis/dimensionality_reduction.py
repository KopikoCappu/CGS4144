"""Dimensionality-reduction utilities for the SRP100712 expression matrix."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
import umap
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.metrics import silhouette_score


TOP_VARIABLE_GENES = 2_000


@dataclass(frozen=True)
class PreparedMatrices:
    """Matrices and gene variances used by the three embedding methods."""

    all_nonzero_variance: pd.DataFrame
    top_variable: pd.DataFrame
    gene_variances: pd.Series
    zero_variance_gene_count: int


def prepare_sample_matrices(
    expression: pd.DataFrame, top_n: int = TOP_VARIABLE_GENES
) -> PreparedMatrices:
    """Transpose to samples by genes, filter constants, and select variable genes.

    Variances are sample variances across the 32 biological samples. No scaling,
    transformation, or sample aggregation is performed.
    """

    sample_by_gene = expression.T.copy()
    sample_by_gene.index.name = "biological_sample"
    gene_variances = sample_by_gene.var(axis=0, ddof=1)
    nonzero_mask = gene_variances > 0
    filtered = sample_by_gene.loc[:, nonzero_mask]
    filtered_variances = gene_variances.loc[nonzero_mask]

    if filtered.shape[1] < top_n:
        raise ValueError(
            f"Requested {top_n:,} variable genes, but only "
            f"{filtered.shape[1]:,} non-zero-variance genes are available."
        )

    top_genes = filtered_variances.nlargest(top_n).index
    top_variable = filtered.loc[:, top_genes]
    return PreparedMatrices(
        all_nonzero_variance=filtered,
        top_variable=top_variable,
        gene_variances=gene_variances,
        zero_variance_gene_count=int((~nonzero_mask).sum()),
    )


def run_pca(matrix: pd.DataFrame) -> tuple[pd.DataFrame, np.ndarray]:
    """Run centered, unscaled PCA and return PC1/PC2 coordinates and variance."""

    model = PCA(n_components=2)
    coordinates = model.fit_transform(matrix)
    coordinate_frame = pd.DataFrame(
        coordinates,
        index=matrix.index.copy(),
        columns=["PC1", "PC2"],
    )
    return coordinate_frame, model.explained_variance_ratio_.copy()


def run_tsne(
    matrix: pd.DataFrame,
    perplexity: float = 5,
    random_state: int = 4144,
) -> pd.DataFrame:
    """Run two-dimensional t-SNE on the selected variable genes."""

    model = TSNE(
        n_components=2,
        perplexity=perplexity,
        random_state=random_state,
        init="pca",
    )
    coordinates = model.fit_transform(matrix)
    return pd.DataFrame(
        coordinates,
        index=matrix.index.copy(),
        columns=["tSNE1", "tSNE2"],
    )


def run_umap(
    matrix: pd.DataFrame,
    n_neighbors: int = 8,
    min_dist: float = 0.2,
    random_state: int = 4144,
) -> pd.DataFrame:
    """Run two-dimensional UMAP on the selected variable genes."""

    model = umap.UMAP(
        n_components=2,
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        random_state=random_state,
    )
    coordinates = model.fit_transform(matrix)
    return pd.DataFrame(
        coordinates,
        index=matrix.index.copy(),
        columns=["UMAP1", "UMAP2"],
    )


def add_sample_metadata(
    coordinates: pd.DataFrame, sample_metadata: pd.DataFrame
) -> pd.DataFrame:
    """Attach biological condition fields to embedding coordinates."""

    metadata = sample_metadata.set_index("biological_sample")
    missing = coordinates.index.difference(metadata.index)
    if len(missing):
        raise ValueError(
            "Cannot attach metadata; missing biological samples: "
            + ", ".join(missing.astype(str)[:10])
        )

    aligned_metadata = metadata.loc[coordinates.index]
    result = aligned_metadata.join(coordinates)
    return result.reset_index()


def embedding_silhouette_scores(
    coordinates: pd.DataFrame,
    coordinate_columns: tuple[str, str],
) -> dict[str, float]:
    """Calculate descriptive grouping scores in a two-dimensional embedding.

    These scores summarize visible grouping only and are not inferential tests.
    """

    values = coordinates.loc[:, list(coordinate_columns)].to_numpy()
    return {
        field: float(silhouette_score(values, coordinates[field].astype(str)))
        for field in ("environment", "age", "genotype")
    }


def qualitative_grouping(score: float) -> str:
    """Translate a silhouette score into cautious visual-description wording."""

    if score >= 0.50:
        return "strong visual separation"
    if score >= 0.25:
        return "moderate visual separation"
    if score >= 0.10:
        return "weak visual separation"
    return "no clear single global separation"
