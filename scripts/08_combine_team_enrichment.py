"""Combine Wilcoxon and teammate g:Profiler GO:BP results by GO ID."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.analysis.combine_enrichment import (  # noqa: E402
    combine_results,
    load_team_gprofiler_results,
    load_wilcoxon_results,
    rank_combined_results,
    select_root_related,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Combine team GO:BP enrichment results without merging p-values."
    )
    parser.add_argument(
        "--wilcoxon",
        type=Path,
        default=(
            PROJECT_ROOT
            / "results"
            / "enrichment"
            / "wilcoxon_GO_BP_all.csv"
        ),
    )
    parser.add_argument(
        "--team-gprofiler",
        type=Path,
        default=(
            PROJECT_ROOT
            / "results"
            / "enrichment"
            / "team"
            / "enrichment_gprofiler_GOBP.csv"
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=PROJECT_ROOT / "results" / "enrichment",
    )
    return parser.parse_args()


def yes_no(value: bool) -> str:
    return "YES" if value else "NO"


def main() -> int:
    args = parse_args()
    warnings: list[str] = []
    wilcoxon = load_wilcoxon_results(args.wilcoxon)
    gprofiler = load_team_gprofiler_results(args.team_gprofiler)
    combined, name_conflicts = combine_results(wilcoxon, gprofiler)
    ranked = rank_combined_results(combined)

    significant = ranked.loc[ranked["methods_significant"] >= 1].copy()
    top10 = ranked.head(10).copy()
    root_related = select_root_related(ranked)

    if name_conflicts:
        warnings.append(
            f"{name_conflicts} shared GO IDs had different term names across "
            "files; the Wilcoxon term name was retained."
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    all_output = args.output_dir / "combined_team_enrichment_all.csv"
    significant_output = (
        args.output_dir / "combined_team_enrichment_significant.csv"
    )
    top10_output = args.output_dir / "combined_team_enrichment_top10.csv"
    root_output = args.output_dir / "combined_team_enrichment_root_related.csv"
    ranked.to_csv(all_output, index=False)
    significant.to_csv(significant_output, index=False)
    top10.to_csv(top10_output, index=False)
    root_related.to_csv(root_output, index=False)

    shared = ranked["methods_tested"].eq(2)
    wilcoxon_significant = ranked["wilcoxon_significant"].fillna(False)
    gprofiler_significant = ranked["gprofiler_significant"].fillna(False)
    significant_both = wilcoxon_significant & gprofiler_significant
    significant_wilcoxon_only = wilcoxon_significant & ~gprofiler_significant
    significant_gprofiler_only = gprofiler_significant & ~wilcoxon_significant

    root_id = "GO:0048364"
    root_rows = ranked.loc[ranked["term_id"].eq(root_id)]
    if len(root_rows) != 1:
        raise ValueError(f"Expected one combined row for {root_id}, found {len(root_rows)}.")
    root_development = root_rows.iloc[0]

    display_columns = [
        "term_id",
        "term_name",
        "methods_tested",
        "methods_significant",
        "significant_fraction",
        "wilcoxon_adjusted_p_value",
        "gprofiler_p_value",
    ]
    print("\nCombined enrichment summary")
    print(f"Wilcoxon terms: {len(wilcoxon):,}")
    print(f"Teammate g:Profiler terms: {len(gprofiler):,}")
    print(f"Unique combined terms: {len(ranked):,}")
    print(f"Terms appearing in both methods: {int(shared.sum()):,}")
    print(f"Terms significant in both methods: {int(significant_both.sum()):,}")
    print(
        "Terms significant only in Wilcoxon: "
        f"{int(significant_wilcoxon_only.sum()):,}"
    )
    print(
        "Terms significant only in g:Profiler: "
        f"{int(significant_gprofiler_only.sum()):,}"
    )

    print("\nGO:0048364 root development check")
    print(
        "Present in Wilcoxon: "
        f"{yes_no(bool(root_development['wilcoxon_tested']))}"
    )
    print(
        "Significant in Wilcoxon: "
        f"{yes_no(bool(root_development['wilcoxon_significant']))}"
    )
    print(
        "Present in teammate g:Profiler: "
        f"{yes_no(bool(root_development['gprofiler_tested']))}"
    )
    print(
        "Significant in teammate g:Profiler: "
        f"{yes_no(bool(root_development['gprofiler_significant']))}"
    )
    print(
        "Wilcoxon adjusted p-value: "
        f"{root_development['wilcoxon_adjusted_p_value']:.6e}"
    )
    print(
        "Teammate g:Profiler p-value: "
        f"{root_development['gprofiler_p_value']:.6e}"
    )

    print("\nTop 10 terms with strongest cross-method support")
    print(top10[display_columns].to_string(index=False))

    root_shared = root_related.loc[root_related["methods_tested"].eq(2)]
    root_shared_significant = root_shared.loc[
        root_shared["methods_significant"] >= 1
    ]
    print(
        f"\nRoot-related shared terms: {len(root_shared):,}; significant in at "
        f"least one method: {len(root_shared_significant):,}"
    )
    if not root_shared_significant.empty:
        print(root_shared_significant[display_columns].to_string(index=False))

    print("\nAssignment-ready interpretation")
    print(
        "The individual analysis used Wilcoxon/Mann-Whitney tests, while the "
        "teammate analysis used g:Profiler. Results were matched by GO ID, and "
        "methods_significant records how many methods marked each term significant. "
        "Method-specific p-values were retained separately and were not averaged "
        "because the analyses use different inputs and statistical approaches. "
        "Agreement across methods supports consistency of an enrichment pattern, "
        "but it does not establish causation."
    )

    print("\nWarnings:")
    if warnings:
        for warning in warnings:
            print(f"  - {warning}")
    else:
        print("  None")

    print("\nFiles written:")
    for path in (all_output, significant_output, top10_output, root_output):
        print(f"  {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
