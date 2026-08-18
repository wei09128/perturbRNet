#!/usr/bin/env python3
"""Locate, classify, and optionally collect perturbRNet Supplementary Table sources.

Run:
  python find_perturbRNet_supplementary_tables.py --copy

The script writes:
  /mnt/c/Wei/RRN/PRNet/manuscript/supplement/table_source_inventory.csv
  /mnt/c/Wei/RRN/PRNet/manuscript/supplement/TABLE_SOURCE_REPORT.txt

With --copy, selected final/primary sources are copied to the supplement directory
and renamed table_s1_... through table_s6_.... Originals are never altered.
"""
from pathlib import Path
import argparse
import csv
import fnmatch
import hashlib
import re
import shutil


REPO = Path("/mnt/c/Wei/RRN/PRNet")
NASA = Path("/mnt/c/Wei/NASA/prn_analysis")
AZD = Path("/mnt/c/Wei/AZD/brain/prn_analysis")


def first_existing(*candidates):
    for candidate in candidates:
        path = Path(candidate)
        if path.exists():
            return path
    # Return the preferred new location so the report is informative.
    return Path(candidates[0])


NORMAN = first_existing(
    "/mnt/c/Wei/RRN/PrNet_realdata/norman_analysis",
    "/mnt/c/Wei/RRN/PRNet_realdata/norman_analysis",
    "/mnt/c/Wei/PrNet_realdata/norman_analysis",
)
OUT = REPO / "manuscript" / "supplement"


SPECS = {
    "S1": {
        "title": "Simulation scenarios and performance",
        "roots": [REPO],
        "patterns": [
            "*benchmark*summary*.csv", "*scenario*.csv", "*method*comparison*.csv",
            "*identifiability*.csv", "*uncertainty*.csv", "*end_to_end*.csv",
        ],
        "exclude": ["manuscript/supplement", "perturbRNet.Rcheck", ".git"],
    },
    "S2": {
        "title": "Norman threshold sensitivity",
        "roots": [NORMAN],
        "patterns": [
            "*threshold*sensitivity*.csv", "*threshold*grid*.csv", "*sensitivity*summary*.csv",
            "*threshold*.csv", "*expression*quintile*.csv", "*expression*filter*.csv",
        ],
    },
    "S3": {
        "title": "Dataset, sample, and availability summary",
        "roots": [Path("/mnt/c/Wei/NASA/pipeline_config"), AZD, NORMAN],
        "patterns": ["samples.csv", "*sample*metadata*.csv", "*manifest.json", "inventory_summary.json", "matrix_diagnostics.csv"],
    },
    "S4": {
        "title": "Alzheimer restoration candidates and calibration",
        "roots": [AZD],
        "patterns": [
            "restoration_classification_counts.csv", "*restoration*genes*.csv", "*candidate*.csv",
            "final_state_calibration.csv", "calibration_*.csv", "degree_association_summary.csv",
            "matched_set_test.csv", "threshold_sensitivity_summary.csv", "matched_set_diagnostics.csv",
        ],
    },
    "S5": {
        "title": "NASA comparison-level PRN and null results",
        "roots": [NASA],
        "patterns": [
            "nasa_prn_comparison_summary.csv", "nasa_prn_mapping_summary.csv",
            "*null*comparison*.csv", "*benchmark*summary*.csv", "BENCHMARK5_SUMMARY.txt",
            "prn_summary.csv", "degree_association_summary.csv",
        ],
    },
    "S6": {
        "title": "Norman validation and robustness",
        "roots": [NORMAN],
        "patterns": [
            "selected_benchmark_metrics.csv", "selected_pair_bootstrap.csv",
            "all_specification_balance.csv", "minimum_positive_sensitivity.csv",
            "matching_balance.csv", "benchmark_metrics.csv", "pair_bootstrap_summary.csv",
            "primary_and_order_sensitivity_metrics.csv", "optimal_pair_bootstrap_summary.csv",
            "crossfit_manifest.json", "manifest.json", "BALANCE_CONSTRAINED_MANIFEST.json",
        ],
    },
}


RAW_SUFFIXES = {".h5ad", ".h5", ".loom", ".mtx", ".rds", ".rda", ".rdata", ".bam", ".fastq", ".fq"}
RAW_NAME_FRAGMENTS = ("row_counts", "raw_counts", "count_matrix", "expression_matrix")


def is_raw(path):
    name = path.name.lower()
    return path.suffix.lower() in RAW_SUFFIXES or any(x in name for x in RAW_NAME_FRAGMENTS)


def matches(path, patterns):
    low = path.name.lower()
    return any(fnmatch.fnmatch(low, pat.lower()) for pat in patterns)


def is_simulation_output(path):
    """Include CSV summaries/results from top-level simulation benchmark outputs."""
    if path.suffix.lower() != ".csv":
        return False
    parts = [part.lower() for part in path.parts]
    benchmark_parts = [part for part in parts if part.startswith("benchmark_")]
    if not benchmark_parts:
        return False
    blocked = ("realdata", "nasa", "nulls")
    return not any(any(word in part for word in blocked) for part in benchmark_parts)


def is_threshold_output(path):
    return (
        path.suffix.lower() == ".csv"
        and "threshold_expression_sensitivity" in str(path).lower()
    )


def classify(table, path):
    text = str(path).lower()
    name = path.name.lower()
    if "perturbrnet.rcheck" in text:
        return "exclude duplicate package-check copy"
    if table == "S4":
        if name in {
            "restoration_classification_counts.csv",
            "azd_low_amplitude_candidates.csv", "azd_strong_candidates.csv",
            "shatm_low_amplitude_candidates.csv", "shatm_strong_candidates.csv",
            "degree_association_summary.csv", "matched_set_test.csv",
            "threshold_sensitivity_summary.csv", "matched_set_diagnostics.csv",
        }:
            return "primary/compact source"
        if path.stat().st_size > 1_000_000:
            return "machine-readable supporting table; do not embed in DOCX"
    if table == "S5":
        if "_quick" in text:
            return "exclude quick-run result"
        if name in {"nasa_prn_comparison_summary.csv", "nasa_prn_mapping_summary.csv", "benchmark5_null_comparison.csv"}:
            return "primary/final source"
        if name == "benchmark5_summary.txt":
            return "human-readable verification only"
        if name == "prn_summary.csv":
            return "redundant with comparison summary"
    return "inspect/assemble"


def selected_for_copy(table, path):
    """Choose manuscript sources while excluding superseded or redundant outputs."""
    text = str(path).lower()
    name = path.name.lower()

    if "perturbrnet.rcheck" in text or "/manuscript/supplement/" in text:
        return False

    if table == "S1":
        # Retain the final revision of each versioned simulation and all files
        # from the uncertainty and method-comparison benchmark outputs.
        if "benchmark_end_to_end_quick" in text and "v3.1" not in text:
            return False
        if "benchmark_identifiability_v4_quick" in text and "v4.1" not in text:
            return False
        allowed_dirs = (
            "benchmark_end_to_end_v3.1_quick",
            "benchmark_identifiability_v4.1_quick",
            "benchmark_method_comparison_quick",
            "benchmark_uncertainty_quick",
            "benchmark_counterdirection",
        )
        return any(part in text for part in allowed_dirs) and path.suffix.lower() == ".csv"

    if table == "S2":
        return "threshold_expression_sensitivity" in text and path.suffix.lower() == ".csv"

    if table == "S3":
        # These are inputs for constructing the dataset/sample table. Other
        # manifests remain discoverable in the inventory but are not copied.
        return (
            name == "samples.csv"
            or name == "inventory_summary.json"
            or name == "matrix_diagnostics.csv"
            or name == "atlas_manifest.json"
            or name == "crossfit_manifest.json"
        )

    if table == "S4":
        compact = {
            "restoration_classification_counts.csv",
            "azd_low_amplitude_candidates.csv", "azd_strong_candidates.csv",
            "shatm_low_amplitude_candidates.csv", "shatm_strong_candidates.csv",
            "degree_association_summary.csv", "matched_set_test.csv",
            "threshold_sensitivity_summary.csv", "matched_set_diagnostics.csv",
        }
        if name not in compact:
            return False
        if "restoration_classification_final" in text or "final_azd_prn_validation" in text:
            return True
        return "final_state_calibration_200" in text and name == "degree_association_summary.csv"

    if table == "S5":
        if "_quick" in text or name in {"prn_summary.csv", "benchmark5_summary.txt"}:
            return False
        if name in {"nasa_prn_comparison_summary.csv", "nasa_prn_mapping_summary.csv"}:
            return True
        return name == "benchmark5_null_comparison.csv" and "_nulls_final" in text

    if table == "S6":
        primary_names = {
            "selected_benchmark_metrics.csv", "selected_pair_bootstrap.csv",
            "all_specification_balance.csv", "minimum_positive_sensitivity.csv",
            "matching_balance.csv", "benchmark_metrics.csv", "pair_bootstrap_summary.csv",
            "primary_and_order_sensitivity_metrics.csv", "optimal_pair_bootstrap_summary.csv",
            "manifest.json", "crossfit_manifest.json", "balance_constrained_manifest.json",
        }
        if name not in primary_names:
            return False
        useful_dirs = (
            "balance_constrained_matching", "expression_aware_matching",
            "gemgroup_crossfit", "noise_standardized_crossfit",
            "prn_crossfit_benchmark", "prn_crossfit_robustness_v2",
            "crossfit_benchmark",
        )
        return any(part in text for part in useful_dirs)

    return False


def safe_token(value):
    value = re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_").lower()
    return value or "source"


def source_root(table, path):
    roots = SPECS[table]["roots"]
    for root in roots:
        try:
            path.relative_to(root)
            return root
        except ValueError:
            continue
    return path.parent


def destination_name(table, path):
    """Build deterministic names and preserve enough parent context for collisions."""
    root = source_root(table, path)
    try:
        rel = path.relative_to(root)
        parents = list(rel.parts[:-1])
    except ValueError:
        parents = [path.parent.name]
    # Preserve enough context to distinguish fixed-threshold gemgroup cross-fit
    # from noise-standardized cross-fit without opaque hash suffixes.
    context = "_".join(safe_token(x) for x in parents[-4:])
    stem = safe_token(path.stem)
    suffix = path.suffix.lower()
    bits = [f"table_{table.lower()}"]
    if context:
        bits.append(context)
    bits.append(stem)
    return "_".join(bits) + suffix


def same_file_contents(a, b):
    if not b.exists() or a.stat().st_size != b.stat().st_size:
        return False
    def digest(path):
        h = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                h.update(chunk)
        return h.digest()
    return digest(a) == digest(b)


def copy_selected(records):
    copied = []
    for record in records:
        table = record["supplementary_table"]
        source = Path(record["source_path"])
        if not selected_for_copy(table, source):
            record["copy_status"] = "not selected"
            record["copied_path"] = ""
            continue
        destination = OUT / destination_name(table, source)
        if destination.exists() and not same_file_contents(source, destination):
            # This should be rare because names contain parent context. Preserve
            # both versions rather than overwriting ambiguous content.
            short_hash = hashlib.sha256(str(source).encode()).hexdigest()[:8]
            destination = destination.with_name(
                f"{destination.stem}_{short_hash}{destination.suffix}"
            )
        if same_file_contents(source, destination):
            status = "already current"
        else:
            shutil.copy2(source, destination)
            status = "copied"
        record["copy_status"] = status
        record["copied_path"] = str(destination)
        copied.append((table, source, destination, status))
    return copied


def scan(root, spec, table):
    if not root.exists():
        return [], f"MISSING ROOT: {root}"
    hits = []
    for path in root.rglob("*"):
        if not path.is_file() or is_raw(path):
            continue
        if any(fragment in str(path) for fragment in spec.get("exclude", [])):
            continue
        if (
            matches(path, spec["patterns"])
            or (table == "S1" and is_simulation_output(path))
            or (table == "S2" and is_threshold_output(path))
        ):
            hits.append(path)
    return sorted(set(hits)), None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--copy", action="store_true",
        help="copy selected final/primary files into the supplement directory",
    )
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    records, report = [], []
    for table, spec in SPECS.items():
        report.append(f"\n===== Table {table}: {spec['title']} =====")
        all_hits = []
        for root in spec["roots"]:
            hits, error = scan(root, spec, table)
            if error:
                report.append(error)
            all_hits.extend(hits)
        all_hits = sorted(set(all_hits))
        if not all_hits:
            report.append("NO MATCHES FOUND")
        for path in all_hits:
            stat = path.stat()
            action = classify(table, path)
            records.append({
                "supplementary_table": table,
                "table_title": spec["title"],
                "source_path": str(path),
                "file_name": path.name,
                "size_bytes": stat.st_size,
                "raw_matrix_excluded": False,
                "recommended_action": action,
            })
            report.append(f"{stat.st_size:>12,d}  {path}")

    copied = copy_selected(records) if args.copy else []
    for record in records:
        record.setdefault("copy_status", "report only")
        record.setdefault("copied_path", "")

    inventory = OUT / "table_source_inventory.csv"
    with inventory.open("w", newline="", encoding="utf-8") as handle:
        fields = [
            "supplementary_table", "table_title", "source_path", "file_name",
            "size_bytes", "raw_matrix_excluded", "recommended_action",
            "copy_status", "copied_path",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)
    report_path = OUT / "TABLE_SOURCE_REPORT.txt"
    report_path.write_text("\n".join(report).lstrip() + "\n", encoding="utf-8")
    print(f"Found {len(records)} candidate source files")
    print(inventory)
    print(report_path)
    if args.copy:
        copy_manifest = OUT / "COPIED_TABLE_SOURCES.txt"
        lines = [
            f"{table}\t{status}\t{source}\t=>\t{destination}"
            for table, source, destination, status in copied
        ]
        copy_manifest.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")
        print(f"Copied/verified {len(copied)} selected files")
        print(copy_manifest)
    print("\nS4 and S5 are derived result summaries, never raw count matrices.")


if __name__ == "__main__":
    main()
