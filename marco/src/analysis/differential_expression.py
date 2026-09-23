"""Adjusted gene-wise linear models for the SRP100712 experiment."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
import statsmodels.api as sm
from patsy import dmatrix
from statsmodels.stats.multitest import multipletests


DESIGN_FORMULA = (
    "1 + C(environment, Treatment(reference='Ground')) "
    "+ C(age, Treatment(reference='4')) "
    "+ C(genotype, Treatment(reference='Col-0'))"
)
ENVIRONMENT_COEFFICIENT = (
    "C(environment, Treatment(reference='Ground'))[T.Spaceflight]"
)


@dataclass(frozen=True)
class ModelDesign:
    """Validated OLS design matrix and the coefficient of primary interest."""

    matrix: pd.DataFrame
    environment_coefficient: str
    formula: str


def build_model_design(
    sample_metadata: pd.DataFrame, sample_order: pd.Index
) -> ModelDesign:
    """Construct and validate the categorical additive model design.

    Ground is explicitly the environment reference, 4 days is the age
    reference, and Col-0 is the genotype reference.
    """

    metadata = sample_metadata.set_index("biological_sample")
    missing = sample_order.difference(metadata.index)
    extra = metadata.index.difference(sample_order)
    if len(missing) or len(extra):
        raise ValueError(
            "Expression and metadata sample identifiers do not align: "
            f"{len(missing)} expression-only and {len(extra)} metadata-only."
        )

    metadata = metadata.loc[sample_order].copy()
    metadata.index.name = "biological_sample"
    metadata = metadata.reset_index()
    metadata["environment"] = pd.Categorical(
        metadata["environment"], categories=["Ground", "Spaceflight"]
    )
    metadata["age"] = pd.Categorical(
        metadata["age"].astype(str), categories=["4", "8"]
    )
    metadata["genotype"] = pd.Categorical(
        metadata["genotype"], categories=["Col-0", "WS"]
    )

    if metadata[["environment", "age", "genotype"]].isna().any().any():
        raise ValueError("Unexpected environment, age, or genotype category found.")

    design = dmatrix(DESIGN_FORMULA, metadata, return_type="dataframe")
    design.index = pd.Index(metadata["biological_sample"], name="biological_sample")

    if ENVIRONMENT_COEFFICIENT not in design.columns:
        raise ValueError(
            "Expected Spaceflight-vs-Ground coefficient is absent from the "
            f"design matrix. Columns: {design.columns.tolist()}"
        )

    environment_column = design[ENVIRONMENT_COEFFICIENT]
    ground_values = environment_column.loc[
        metadata.set_index("biological_sample")["environment"].eq("Ground")
    ]
    spaceflight_values = environment_column.loc[
        metadata.set_index("biological_sample")["environment"].eq("Spaceflight")
    ]
    reference_is_correct = bool(
        (ground_values == 0).all() and (spaceflight_values == 1).all()
    )
    if not reference_is_correct:
        raise ValueError(
            "Ground is not encoded as the reference for the environment "
            "coefficient; differential-expression analysis was stopped."
        )

    rank = np.linalg.matrix_rank(design.to_numpy())
    if rank != design.shape[1]:
        raise ValueError(
            f"Model design is rank deficient: rank {rank}, "
            f"{design.shape[1]} columns."
        )

    return ModelDesign(
        matrix=design,
        environment_coefficient=ENVIRONMENT_COEFFICIENT,
        formula=DESIGN_FORMULA,
    )


def fit_gene_wise_ols(
    expression: pd.DataFrame,
    design: ModelDesign,
    fdr_threshold: float = 0.05,
) -> pd.DataFrame:
    """Fit a separate statsmodels OLS model for every expression row."""

    if not expression.columns.equals(design.matrix.index):
        raise ValueError(
            "Expression columns must exactly match the ordered design-matrix index."
        )

    design_values = design.matrix.to_numpy(dtype=np.float64, copy=False)
    expression_values = expression.to_numpy(dtype=np.float64, copy=False)
    coefficient_index = design.matrix.columns.get_loc(
        design.environment_coefficient
    )

    gene_count = expression.shape[0]
    coefficients = np.empty(gene_count, dtype=np.float64)
    standard_errors = np.empty(gene_count, dtype=np.float64)
    t_statistics = np.empty(gene_count, dtype=np.float64)
    p_values = np.empty(gene_count, dtype=np.float64)

    for position, gene_values in enumerate(expression_values):
        fitted = sm.OLS(gene_values, design_values, hasconst=True).fit()
        coefficients[position] = fitted.params[coefficient_index]
        standard_errors[position] = fitted.bse[coefficient_index]
        t_statistics[position] = fitted.tvalues[coefficient_index]
        p_values[position] = fitted.pvalues[coefficient_index]

    arrays = (coefficients, standard_errors, t_statistics, p_values)
    if not all(np.isfinite(values).all() for values in arrays):
        raise ValueError(
            "At least one gene produced a non-finite coefficient, standard error, "
            "t statistic, or p-value; FDR correction was not performed."
        )

    rejected, adjusted_p_values, _, _ = multipletests(
        p_values,
        alpha=fdr_threshold,
        method="fdr_bh",
    )
    results = pd.DataFrame(
        {
            "Gene": expression.index.astype(str),
            "spaceflight_effect": coefficients,
            "standard_error": standard_errors,
            "t_statistic": t_statistics,
            "p_value": p_values,
            "adjusted_p_value": adjusted_p_values,
            "significant": rejected,
        }
    )
    return results.sort_values(
        ["adjusted_p_value", "p_value", "Gene"], kind="stable"
    ).reset_index(drop=True)


def select_volcano_labels(
    results: pd.DataFrame, maximum_labels: int = 10
) -> pd.DataFrame:
    """Select significant genes using both FDR and absolute effect magnitude."""

    candidates = results.loc[results["significant"]].copy()
    if candidates.empty:
        return candidates

    safe_adjusted = np.clip(
        candidates["adjusted_p_value"].to_numpy(dtype=float),
        np.finfo(float).tiny,
        1.0,
    )
    candidates["negative_log10_fdr"] = -np.log10(safe_adjusted)
    candidates["label_score"] = (
        candidates["negative_log10_fdr"].rank(pct=True)
        + candidates["spaceflight_effect"].abs().rank(pct=True)
    )
    return candidates.nlargest(maximum_labels, "label_score")
