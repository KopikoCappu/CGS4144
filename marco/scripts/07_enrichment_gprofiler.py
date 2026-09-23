"""GO Biological Process enrichment with g:Profiler."""

from __future__ import annotations

import argparse
import re
import sys
import time
from pathlib import Path

import pandas as pd
from gprofiler import GProfiler


PROJECT_ROOT = Path(__file__).resolve().parents[1]
ORGANISM = "athaliana"
SOURCE = "GO:BP"
CORRECTION_METHOD = "g_SCS"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run Arabidopsis GO:BP enrichment with g:Profiler."
    )
    parser.add_argument(
        "--significant-genes",
        type=Path,
        default=(
            PROJECT_ROOT
            / "results"
            / "tables"
            / "differential_expression_significant.csv"
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=PROJECT_ROOT / "results" / "enrichment",
    )
    return parser.parse_args()


def load_significant_genes(path: Path) -> list[str]:
    results = pd.read_csv(path)
    if "Gene" not in results.columns:
        raise ValueError("Significant differential-expression table lacks a Gene column.")
    if len(results) != 2_560:
        raise ValueError(f"Expected 2,560 significant genes, found {len(results):,}.")
    if results["Gene"].isna().any():
        raise ValueError("Significant differential-expression table has missing gene IDs.")
    if results["Gene"].duplicated().any():
        raise ValueError("Significant differential-expression table has duplicate gene IDs.")
    if "significant" in results and not results["significant"].astype(bool).all():
        raise ValueError("Input table contains genes not marked significant.")
    return results["Gene"].astype(str).str.upper().tolist()


def query_gprofiler(genes: list[str], attempts: int = 3) -> list[dict]:
    client = GProfiler(return_dataframe=False)
    last_error: Exception | None = None
    for attempt in range(attempts):
        try:
            return client.profile(
                organism=ORGANISM,
                query=genes,
                sources=[SOURCE],
                user_threshold=0.05,
                all_results=True,
                ordered=False,
                no_evidences=False,
                domain_scope="annotated",
                significance_threshold_method=CORRECTION_METHOD,
            )
        except Exception as error:
            last_error = error
            if attempt + 1 < attempts:
                time.sleep(2**attempt)
    raise RuntimeError(
        f"g:Profiler query failed after {attempts} attempts."
    ) from last_error


def _join_values(values: object) -> str:
    if not isinstance(values, list):
        return ""
    return ";".join(str(value) for value in values)


def _format_evidence(genes: object, evidence: object) -> str:
    if not isinstance(genes, list) or not isinstance(evidence, list):
        return ""
    pairs = []
    for gene, codes in zip(genes, evidence):
        code_text = "|".join(str(code) for code in codes) if isinstance(codes, list) else str(codes)
        pairs.append(f"{gene}:{code_text}")
    return ";".join(pairs)


def format_results(records: list[dict]) -> pd.DataFrame:
    rows = []
    for record in records:
        intersections = record.get("intersections", [])
        api_p_value = float(record["p_value"])
        rows.append(
            {
                "GO_ID": record["native"],
                "GO_term": record["name"],
                "p_value": api_p_value,
                "adjusted_p_value": api_p_value,
                "significant": bool(record["significant"]),
                "term_size": int(record["term_size"]),
                "query_size": int(record["query_size"]),
                "intersection_size": int(record["intersection_size"]),
                "effective_domain_size": int(record["effective_domain_size"]),
                "precision": float(record["precision"]),
                "recall": float(record["recall"]),
                "contributing_genes": _join_values(intersections),
                "evidence_codes": _format_evidence(
                    intersections, record.get("evidences", [])
                ),
                "description": record.get("description", ""),
                "parents": _join_values(record.get("parents", [])),
                "source": record["source"],
                "organism": ORGANISM,
                "correction_method": CORRECTION_METHOD,
            }
        )

    results = pd.DataFrame(rows)
    if results.empty:
        raise ValueError("g:Profiler returned no GO:BP results.")
    if not results["source"].eq(SOURCE).all():
        raise ValueError("g:Profiler response contains a source other than GO:BP.")
    if not results["GO_ID"].str.startswith("GO:").all():
        raise ValueError("g:Profiler response contains a non-GO identifier.")
    if results["GO_ID"].duplicated().any():
        raise ValueError("g:Profiler returned duplicate GO:BP terms for one query.")
    return results.sort_values(
        ["adjusted_p_value", "GO_ID"], kind="stable"
    ).reset_index(drop=True)


def select_root_related(results: pd.DataFrame) -> pd.DataFrame:
    keywords = [
        "root",
        "root hair",
        "root meristem",
        "lateral root",
        "auxin",
        "gravitropism",
        "gravity",
        "cell wall",
        "development",
        "morphogenesis",
    ]
    pattern = "|".join(re.escape(keyword) for keyword in keywords)
    return results.loc[
        results["GO_term"].str.contains(pattern, case=False, na=False, regex=True)
    ].copy()


def main() -> int:
    args = parse_args()
    genes = load_significant_genes(args.significant_genes)
    print(
        f"Submitting {len(genes):,} Arabidopsis genes to g:Profiler "
        f"for {SOURCE} enrichment..."
    )
    results = format_results(query_gprofiler(genes))
    significant = results.loc[results["significant"]].copy()
    top20 = results.head(20).copy()
    root_related = select_root_related(results)
    root_related_significant = root_related.loc[root_related["significant"]]

    args.output_dir.mkdir(parents=True, exist_ok=True)
    all_output = args.output_dir / "gprofiler_GO_BP_all.csv"
    significant_output = args.output_dir / "gprofiler_GO_BP_significant.csv"
    top20_output = args.output_dir / "gprofiler_GO_BP_top20.csv"
    root_output = args.output_dir / "gprofiler_GO_BP_root_related.csv"
    results.to_csv(all_output, index=False)
    significant.to_csv(significant_output, index=False)
    top20.to_csv(top20_output, index=False)
    root_related.to_csv(root_output, index=False)

    warnings = []
    recognized_query_size = int(results["query_size"].max())
    if recognized_query_size < len(genes):
        warnings.append(
            f"g:Profiler used {recognized_query_size:,} of {len(genes):,} submitted "
            "genes in the GO:BP query domain."
        )

    display_columns = [
        "GO_ID",
        "GO_term",
        "term_size",
        "intersection_size",
        "adjusted_p_value",
    ]
    print("\n=== G:PROFILER GO:BP SUMMARY ===")
    print(f"Genes submitted: {len(genes):,}")
    print(f"Genes represented in the GO:BP query domain: {recognized_query_size:,}")
    print(f"GO:BP terms returned: {len(results):,}")
    print(f"Significant GO:BP terms: {len(significant):,}")
    print(f"Minimum g:SCS-adjusted p-value: {results['adjusted_p_value'].min():.6e}")
    print("The g:Profiler p_value field is g:SCS corrected in these outputs.")

    print("\nTop 10 enriched GO:BP terms:")
    print(results.head(10)[display_columns].to_string(index=False))

    print(
        f"\nRoot-related significant terms "
        f"({len(root_related_significant):,} of {len(root_related):,} keyword matches):"
    )
    if root_related_significant.empty:
        print("  None")
    else:
        print(root_related_significant[display_columns].to_string(index=False))

    print("\n=== ASSIGNMENT-READY INTERPRETATION ===")
    print(
        f"g:Profiler was used as a second enrichment method to test Gene Ontology "
        f"Biological Process terms for the {len(genes):,} genes significant in "
        "the adjusted Spaceflight-versus-Ground analysis. g:Profiler's g:SCS "
        f"correction identified {len(significant):,} significant GO:BP terms. "
        "These results describe functional enrichment associations among the "
        "significant genes and do not establish causal effects on root development."
    )

    print("\nWarnings:")
    if warnings:
        for warning in warnings:
            print(f"  - {warning}")
    else:
        print("  None")

    print("\nFiles written:")
    for path in (all_output, significant_output, top20_output, root_output):
        print(f"  {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
