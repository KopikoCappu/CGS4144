"""GO Biological Process rank-sum enrichment for SRP100712."""

from __future__ import annotations

import time
from pathlib import Path

import mygene
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests


ANNOTATION_COLUMNS = ["Gene", "GO_ID", "GO_term"]
RESULT_COLUMNS = [
    "GO_ID",
    "GO_term",
    "gene_set_size",
    "median_effect_in_set",
    "median_effect_outside_set",
    "effect_difference",
    "statistic",
    "p_value",
    "adjusted_p_value",
    "significant",
]


def load_differential_expression(path: str | Path) -> pd.DataFrame:
    results = pd.read_csv(path)
    required = {"Gene", "spaceflight_effect", "p_value", "adjusted_p_value"}
    missing = sorted(required.difference(results.columns))
    if missing:
        raise ValueError(
            "Differential-expression table is missing columns: "
            + ", ".join(missing)
        )
    if len(results) != 32_309:
        raise ValueError(f"Expected 32,309 genes, found {len(results):,}.")
    if results["Gene"].isna().any() or results["Gene"].duplicated().any():
        raise ValueError("Differential-expression gene identifiers are missing or duplicated.")
    if not np.isfinite(results["spaceflight_effect"].to_numpy(dtype=float)).all():
        raise ValueError("Spaceflight effects contain missing or infinite values.")
    results = results.copy()
    results["Gene"] = results["Gene"].astype(str).str.upper()
    return results


def _parse_bp_terms(record: dict) -> list[tuple[str, str]]:
    go = record.get("go") or {}
    bp = go.get("BP") if isinstance(go, dict) else None
    if not bp:
        return []
    entries = bp if isinstance(bp, list) else [bp]
    terms = {
        (str(entry["id"]), str(entry["term"]))
        for entry in entries
        if isinstance(entry, dict) and entry.get("id") and entry.get("term")
    }
    return sorted(terms)


def _query_annotation_batch(
    client: mygene.MyGeneInfo, genes: list[str], attempts: int = 3
) -> pd.DataFrame:
    last_error: Exception | None = None
    for attempt in range(attempts):
        try:
            records = client.querymany(
                genes,
                scopes="symbol,ensembl.gene",
                fields="go.BP",
                species=3702,
                as_dataframe=False,
                verbose=False,
            )
            break
        except Exception as error:
            last_error = error
            if attempt + 1 < attempts:
                time.sleep(2**attempt)
    else:
        raise RuntimeError(
            f"MyGene.info annotation query failed after {attempts} attempts."
        ) from last_error

    requested = set(genes)
    terms_by_gene: dict[str, set[tuple[str, str]]] = {
        gene: set() for gene in genes
    }
    for record in records:
        gene = str(record.get("query", "")).upper()
        if gene in requested and not record.get("notfound"):
            terms_by_gene[gene].update(_parse_bp_terms(record))

    rows: list[dict[str, object]] = []
    for gene in genes:
        terms = terms_by_gene[gene]
        if terms:
            rows.extend(
                {"Gene": gene, "GO_ID": go_id, "GO_term": term}
                for go_id, term in sorted(terms)
            )
        else:
            rows.append({"Gene": gene, "GO_ID": pd.NA, "GO_term": pd.NA})
    return pd.DataFrame(rows, columns=ANNOTATION_COLUMNS)


def _read_annotation_cache(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=ANNOTATION_COLUMNS)
    annotations = pd.read_csv(path, dtype="string")
    missing = sorted(set(ANNOTATION_COLUMNS).difference(annotations.columns))
    if missing:
        raise ValueError("Annotation cache is missing columns: " + ", ".join(missing))
    annotations = annotations[ANNOTATION_COLUMNS].copy()
    annotations["Gene"] = annotations["Gene"].str.upper()
    if annotations["Gene"].isna().any():
        raise ValueError("Annotation cache contains a missing gene identifier.")
    incomplete = annotations["GO_ID"].isna() ^ annotations["GO_term"].isna()
    if incomplete.any():
        raise ValueError("Annotation cache contains incomplete GO term rows.")
    return annotations.drop_duplicates().reset_index(drop=True)


def _write_annotation_cache(annotations: pd.DataFrame, path: Path) -> None:
    ordered = annotations.sort_values(
        ["Gene", "GO_ID", "GO_term"], na_position="last", kind="stable"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    ordered.to_csv(path, index=False)


def load_or_query_go_bp_annotations(
    genes: pd.Series | pd.Index | list[str],
    cache_path: str | Path,
    batch_size: int = 1_000,
) -> tuple[pd.DataFrame, int]:
    cache_path = Path(cache_path)
    requested_genes = list(dict.fromkeys(str(gene).upper() for gene in genes))
    annotations = _read_annotation_cache(cache_path)
    cached_genes = set(annotations["Gene"].astype(str))
    missing_genes = [gene for gene in requested_genes if gene not in cached_genes]
    queried_count = len(missing_genes)

    if missing_genes:
        client = mygene.MyGeneInfo()
        batch_total = (len(missing_genes) + batch_size - 1) // batch_size
        for batch_number, start in enumerate(
            range(0, len(missing_genes), batch_size), start=1
        ):
            batch = missing_genes[start : start + batch_size]
            print(
                f"Querying MyGene.info GO:BP annotations: batch "
                f"{batch_number}/{batch_total}"
            )
            new_rows = _query_annotation_batch(client, batch)
            annotations = pd.concat([annotations, new_rows], ignore_index=True)
            annotations = annotations.drop_duplicates().reset_index(drop=True)
            _write_annotation_cache(annotations, cache_path)

    requested_set = set(requested_genes)
    selected = annotations.loc[annotations["Gene"].isin(requested_set)].copy()
    represented = set(selected["Gene"].astype(str))
    absent = requested_set.difference(represented)
    if absent:
        raise ValueError(
            f"Annotation cache lacks query records for {len(absent):,} genes."
        )
    return selected.reset_index(drop=True), queried_count


def annotation_summary(
    annotations: pd.DataFrame, genes: pd.Series | pd.Index | list[str]
) -> dict[str, int]:
    total_genes = len(set(str(gene).upper() for gene in genes))
    mapped = annotations.loc[annotations["GO_ID"].notna()]
    annotated_genes = mapped["Gene"].nunique()
    return {
        "total_genes": total_genes,
        "annotated_genes": int(annotated_genes),
        "unannotated_genes": int(total_genes - annotated_genes),
        "unique_terms": int(mapped["GO_ID"].nunique()),
    }


def build_gene_sets(
    annotations: pd.DataFrame,
    universe: pd.Index,
    minimum_size: int = 10,
    maximum_size: int = 2_000,
) -> dict[str, tuple[str, set[str]]]:
    mapped = annotations.dropna(subset=["GO_ID", "GO_term"]).copy()
    mapped = mapped.loc[mapped["Gene"].isin(universe)].drop_duplicates(
        ["Gene", "GO_ID"]
    )

    gene_sets: dict[str, tuple[str, set[str]]] = {}
    for go_id, group in mapped.groupby("GO_ID", sort=True):
        genes = set(group["Gene"].astype(str))
        if minimum_size <= len(genes) <= maximum_size:
            term_counts = group["GO_term"].value_counts()
            term = sorted(term_counts[term_counts.eq(term_counts.max())].index)[0]
            gene_sets[str(go_id)] = (str(term), genes)
    return gene_sets


def run_mannwhitney_enrichment(
    differential_expression: pd.DataFrame,
    gene_sets: dict[str, tuple[str, set[str]]],
    fdr_threshold: float = 0.05,
) -> pd.DataFrame:
    effects = differential_expression.set_index("Gene")["spaceflight_effect"]
    universe = effects.index
    effect_values = effects.to_numpy(dtype=float)
    rows: list[dict[str, object]] = []

    for go_id, (term, genes) in gene_sets.items():
        in_set = universe.isin(genes)
        set_effects = effect_values[in_set]
        outside_effects = effect_values[~in_set]
        test = mannwhitneyu(
            set_effects,
            outside_effects,
            alternative="two-sided",
            method="asymptotic",
        )
        median_in = float(np.median(set_effects))
        median_out = float(np.median(outside_effects))
        rows.append(
            {
                "GO_ID": go_id,
                "GO_term": term,
                "gene_set_size": int(in_set.sum()),
                "median_effect_in_set": median_in,
                "median_effect_outside_set": median_out,
                "effect_difference": median_in - median_out,
                "statistic": float(test.statistic),
                "p_value": float(test.pvalue),
            }
        )

    results = pd.DataFrame(rows)
    if results.empty:
        return pd.DataFrame(columns=RESULT_COLUMNS)
    _, adjusted, _, _ = multipletests(
        results["p_value"].to_numpy(), alpha=fdr_threshold, method="fdr_bh"
    )
    results["adjusted_p_value"] = adjusted
    results["significant"] = results["adjusted_p_value"] < fdr_threshold
    return results.sort_values(
        ["adjusted_p_value", "p_value", "GO_ID"], kind="stable"
    ).reset_index(drop=True)


def select_root_related_terms(results: pd.DataFrame) -> pd.DataFrame:
    keywords = [
        "root",
        "root hair",
        "lateral root",
        "gravitropism",
        "gravity",
        "auxin",
        "meristem",
        "cell wall",
        "development",
        "morphogenesis",
    ]
    pattern = "|".join(keywords)
    return results.loc[
        results["GO_term"].str.contains(pattern, case=False, na=False, regex=True)
    ].copy()
