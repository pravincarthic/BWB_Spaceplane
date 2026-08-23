#!/usr/bin/env python3
"""
Gmsh Mesh Quality Metrics Report
=================================

Loads an existing mesh with the Gmsh Python API and computes every
element quality metric that Gmsh exposes, broken down per element
type present in the mesh. The mesh itself is never modified - it is
only read and queried.

In addition to summary statistics (min/max/mean/std) for all metrics,
this script also produces, for one chosen "primary" metric:
  - a histogram of the quality distribution
  - counts/percentages of elements below quality thresholds
  - a list of the worst N elements, with their tag and centroid
    coordinates, so they can be located in the Gmsh GUI

Usage:
    python mesh_quality_report.py mesh_file.msh
    python mesh_quality_report.py mesh_file.msh --summary-csv summary.csv
    python mesh_quality_report.py mesh_file.msh --bad-elements-csv bad.csv
    python mesh_quality_report.py mesh_file.msh --primary-metric gamma --worst-n 50
    python mesh_quality_report.py mesh_file.msh --no-histogram --no-worst

Requires:
    pip install gmsh numpy
"""

import argparse
import csv
import statistics
import sys

import numpy as np
import gmsh


# Every quality measure supported by gmsh.model.mesh.getElementQualities
QUALITY_METRICS = [
    "minDetJac",
    "maxDetJac",
    "minSICN",
    "minSIGE",
    "gamma",
    "innerRadius",
    "outerRadius",
    "minIsotropy",
    "angleShape",
    "minEdge",
    "maxEdge",
    "volume",
]

METRIC_DESCRIPTIONS = {
    "minDetJac": "Minimum of the Jacobian determinant",
    "maxDetJac": "Maximum of the Jacobian determinant",
    "minSICN": "Minimum Signed Inverse Condition Number (shape quality)",
    "minSIGE": "Minimum Signed Inverse Gradient Error",
    "gamma": "Inscribed/circumscribed radius ratio, normalized to [0,1]",
    "innerRadius": "Radius of the inscribed circle/sphere",
    "outerRadius": "Radius of the circumscribed circle/sphere",
    "minIsotropy": "Minimum isotropy measure",
    "angleShape": "Angle-based shape measure",
    "minEdge": "Minimum edge length",
    "maxEdge": "Maximum edge length",
    "volume": "Element volume (area in 2D, length in 1D)",
}

# Metrics that are bounded, quality-style measures where lower = worse.
# These are the only ones eligible for histogram / threshold / worst-element
# analysis, since the other metrics (volume, edge length, radii, jacobian)
# are physical quantities, not bounded quality scores.
NORMALIZED_METRICS = ["minSICN", "minSIGE", "gamma", "minIsotropy", "angleShape"]

ELEMENT_TYPE_NAMES = {
    1: "2-node line",
    2: "3-node triangle",
    3: "4-node quadrangle",
    4: "4-node tetrahedron",
    5: "8-node hexahedron",
    6: "6-node prism",
    7: "5-node pyramid",
    8: "3-node second order line",
    9: "6-node second order triangle",
    10: "9-node second order quadrangle",
    11: "10-node second order tetrahedron",
    12: "27-node second order hexahedron",
    15: "1-node point",
    16: "8-node second order quadrangle",
    17: "20-node second order hexahedron",
}

# Gmsh uses DBL_MAX (~1.797e308) as a sentinel for "not applicable" in some
# metric/element-type combinations (observed e.g. for minSIGE on lines).
# Any |value| at or above this is treated as invalid and excluded from
# statistics rather than silently corrupting min/max/mean.
INVALID_SENTINEL = 1.0e300


def compute_stats(values):
    """Return summary statistics for a numpy array of numeric values,
    filtering out gmsh's DBL_MAX 'not applicable' sentinel values."""
    if values is None or len(values) == 0:
        return None
    arr = np.asarray(values, dtype=float)
    valid_mask = np.abs(arr) < INVALID_SENTINEL
    valid = arr[valid_mask]
    n_excluded = arr.size - valid.size
    if valid.size == 0:
        return None
    return {
        "count": int(valid.size),
        "excluded": int(n_excluded),
        "min": float(valid.min()),
        "max": float(valid.max()),
        "mean": float(valid.mean()),
        "std": float(valid.std()),
    }


def print_stats_table(title, rows):
    print()
    print(title)
    print("-" * len(title))
    header = "{:<14}{:>10}{:>16}{:>16}{:>16}{:>16}".format(
        "Metric", "Count", "Min", "Max", "Mean", "StdDev"
    )
    print(header)
    print("-" * len(header))
    for name, stats in rows:
        if stats is None:
            print("{:<14}{:>10}".format(name, "N/A"))
            continue
        note = ""
        if stats["excluded"] > 0:
            note = "  ({} N/A)".format(stats["excluded"])
        print(
            "{:<14}{:>10}{:>16.6g}{:>16.6g}{:>16.6g}{:>16.6g}{}".format(
                name,
                stats["count"],
                stats["min"],
                stats["max"],
                stats["mean"],
                stats["std"],
                note,
            )
        )


def has_valid_data(values):
    """True if this metric array exists and contains at least one
    non-sentinel value."""
    if values is None or len(values) == 0:
        return False
    return bool(np.any(np.abs(values) < INVALID_SENTINEL))


def get_valid_quality_array(tags, metric):
    """Fetch a quality metric for a set of element tags and return it as
    a numpy array. Returns None if the metric is not defined for this
    element type (gmsh raises an exception in that case)."""
    try:
        raw = gmsh.model.mesh.getElementQualities(tags, metric)
    except Exception:
        return None
    return np.asarray(raw, dtype=float)


def print_histogram(values, metric_name, n_bins):
    valid_mask = np.abs(values) < INVALID_SENTINEL
    valid = values[valid_mask]
    if valid.size == 0:
        print("  (no valid values to histogram for {})".format(metric_name))
        return
    vmin = float(valid.min())
    vmax = float(valid.max())
    # Guard against floating-point near-equal ranges (e.g. a perfectly
    # regular mesh where all values are identical up to numerical noise),
    # which would otherwise make np.histogram raise on a degenerate range.
    scale = max(abs(vmin), abs(vmax), 1.0)
    if (vmax - vmin) <= scale * 1e-9:
        print("  All {} valid elements have {} = {:.6g} (no meaningful spread)".format(
            valid.size, metric_name, vmin
        ))
        return
    counts, edges = np.histogram(valid, bins=n_bins, range=(vmin, vmax))
    print()
    print("Histogram of {} ({} elements)".format(metric_name, valid.size))
    print("-" * 60)
    total = valid.size
    for i in range(n_bins):
        lo, hi = edges[i], edges[i + 1]
        c = int(counts[i])
        pct = 100.0 * c / total
        bar = "#" * int(round(pct / 2))
        print("  [{:>12.5g}, {:>12.5g}) : {:>10} ({:>6.2f}%) {}".format(
            lo, hi, c, pct, bar
        ))


def print_thresholds(values, metric_name, thresholds):
    valid_mask = np.abs(values) < INVALID_SENTINEL
    valid = values[valid_mask]
    if valid.size == 0:
        print("  (no valid values for threshold check on {})".format(metric_name))
        return
    total = valid.size
    print()
    print("Elements at or below quality thresholds for {} ({} valid elements)".format(
        metric_name, total
    ))
    print("-" * 60)
    for t in thresholds:
        c = int(np.sum(valid <= t))
        pct = 100.0 * c / total
        print("  <= {:<10g} : {:>10} elements ({:>6.2f}%)".format(t, c, pct))


def get_element_centroid(element_tag):
    """Look up a single element's node coordinates and return its centroid.
    Only called for a small number (worst-N) of elements, so per-element
    API calls are cheap even on very large meshes."""
    etype, node_tags, _dim, _entity_tag = gmsh.model.mesh.getElement(element_tag)
    coords_sum = np.zeros(3)
    n_nodes = len(node_tags)
    for nt in node_tags:
        coord, _param, _dim, _tag = gmsh.model.mesh.getNode(nt)
        coords_sum += np.asarray(coord[:3], dtype=float)
    return coords_sum / n_nodes


def find_worst_elements(tags, values, n):
    """Return indices of the n smallest (worst-quality) valid values,
    sorted ascending (worst first). Uses argpartition to stay efficient
    on very large element counts instead of a full sort."""
    valid_mask = np.abs(values) < INVALID_SENTINEL
    valid_idx = np.nonzero(valid_mask)[0]
    if valid_idx.size == 0:
        return np.array([], dtype=int)
    n = min(n, valid_idx.size)
    sub_values = values[valid_idx]
    if n < sub_values.size:
        part = np.argpartition(sub_values, n)[:n]
    else:
        part = np.arange(sub_values.size)
    order = part[np.argsort(sub_values[part])]
    return valid_idx[order]


def analyze_mesh(mesh_file, summary_csv=None, bad_elements_csv=None,
                  primary_metric="minSICN", worst_n=20, hist_bins=10,
                  thresholds=None, do_histogram=True, do_worst=True):
    if thresholds is None:
        thresholds = [0.01, 0.05, 0.1, 0.2, 0.3]

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)

    summary_rows = []
    bad_element_rows = []

    try:
        # Read-only: gmsh.open() loads the mesh file, it does not alter it.
        gmsh.open(mesh_file)

        print("=" * 70)
        print("GMSH MESH QUALITY REPORT")
        print("=" * 70)
        print("File: {}".format(mesh_file))

        element_types, element_tags, node_tags = gmsh.model.mesh.getElements()

        if len(element_types) == 0:
            print("No elements found in the mesh.")
            return

        total_elements = sum(len(tags) for tags in element_tags)
        print("Total elements: {}".format(total_elements))
        print("Element types present: {}".format(len(element_types)))

        for etype, tags in zip(element_types, element_tags):
            type_name = ELEMENT_TYPE_NAMES.get(etype, "type {}".format(etype))
            n_elems = len(tags)
            if n_elems == 0:
                continue

            tags_arr = np.asarray(tags, dtype=np.int64)

            print()
            print("=" * 70)
            print(
                "Element type: {} (gmsh type {}), {} elements".format(
                    type_name, etype, n_elems
                )
            )
            print("=" * 70)

            rows = []
            metric_arrays = {}
            for metric in QUALITY_METRICS:
                values = get_valid_quality_array(tags_arr, metric)
                metric_arrays[metric] = values
                stats = compute_stats(values) if values is not None else None
                rows.append((metric, stats))

                if summary_csv is not None:
                    if stats is None:
                        summary_rows.append(
                            [etype, type_name, n_elems, metric] + ["N/A"] * 6
                        )
                    else:
                        summary_rows.append(
                            [
                                etype, type_name, n_elems, metric,
                                stats["count"], stats["excluded"],
                                stats["min"], stats["max"],
                                stats["mean"], stats["std"],
                            ]
                        )

            print_stats_table("Quality metrics for {}".format(type_name), rows)

            # Pick the metric to use for histogram / threshold / worst-element
            # analysis on this element type: the requested primary metric if
            # it is defined here, otherwise fall back to the first normalized
            # metric that is actually available.
            chosen_metric = None
            chosen_values = None
            if has_valid_data(metric_arrays.get(primary_metric)):
                chosen_metric = primary_metric
                chosen_values = metric_arrays[primary_metric]
            else:
                for m in NORMALIZED_METRICS:
                    if has_valid_data(metric_arrays.get(m)):
                        chosen_metric = m
                        chosen_values = metric_arrays[m]
                        break

            if chosen_metric is None:
                print()
                print("No normalized quality metric with valid data for {} - "
                      "skipping histogram/threshold/worst-element analysis."
                      .format(type_name))
                continue

            if chosen_metric != primary_metric:
                print()
                print("Note: '{}' is not defined for {}; using '{}' instead "
                      "for detailed analysis below.".format(
                          primary_metric, type_name, chosen_metric))

            if do_histogram:
                print_histogram(chosen_values, chosen_metric, hist_bins)

            print_thresholds(chosen_values, chosen_metric, thresholds)

            if do_worst:
                worst_positions = find_worst_elements(tags_arr, chosen_values, worst_n)
                if worst_positions.size > 0:
                    print()
                    print("Worst {} elements by {} for {}".format(
                        worst_positions.size, chosen_metric, type_name))
                    print("-" * 78)
                    print("{:<14}{:>14}{:>16}{:>16}{:>16}".format(
                        "ElementTag", chosen_metric, "X", "Y", "Z"
                    ))
                    for pos in worst_positions:
                        etag = int(tags_arr[pos])
                        qval = float(chosen_values[pos])
                        cx, cy, cz = get_element_centroid(etag)
                        print("{:<14}{:>14.6g}{:>16.6g}{:>16.6g}{:>16.6g}".format(
                            etag, qval, cx, cy, cz
                        ))
                        if bad_elements_csv is not None:
                            bad_element_rows.append(
                                [etype, type_name, etag, chosen_metric, qval, cx, cy, cz]
                            )

        print()
        print("=" * 70)
        print("Metric descriptions")
        print("=" * 70)
        for metric in QUALITY_METRICS:
            print("{:<14}{}".format(metric, METRIC_DESCRIPTIONS[metric]))

        if summary_csv is not None:
            with open(summary_csv, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(
                    [
                        "elem_type_id", "elem_type_name", "n_elements", "metric",
                        "count", "excluded_na", "min", "max", "mean", "std",
                    ]
                )
                writer.writerows(summary_rows)
            print()
            print("Summary CSV written to: {}".format(summary_csv))

        if bad_elements_csv is not None:
            with open(bad_elements_csv, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(
                    [
                        "elem_type_id", "elem_type_name", "element_tag",
                        "metric", "quality_value", "x", "y", "z",
                    ]
                )
                writer.writerows(bad_element_rows)
            print("Worst-element CSV written to: {}".format(bad_elements_csv))

        print()
        print("Report complete. Mesh was not modified.")

    finally:
        gmsh.finalize()


def parse_thresholds(text):
    return [float(t.strip()) for t in text.split(",") if t.strip()]


def main():
    parser = argparse.ArgumentParser(
        description="Compute every Gmsh mesh quality metric without modifying the mesh."
    )
    parser.add_argument("mesh_file", help="Path to the mesh file (.msh, .geo, etc.)")
    parser.add_argument(
        "--summary-csv", dest="summary_csv", default=None,
        help="Optional path to write the full min/max/mean/std summary as CSV.",
    )
    parser.add_argument(
        "--bad-elements-csv", dest="bad_elements_csv", default=None,
        help="Optional path to write the worst-N elements (tag + centroid) as CSV.",
    )
    parser.add_argument(
        "--primary-metric", dest="primary_metric", default="minSICN",
        choices=NORMALIZED_METRICS,
        help="Quality metric used for histogram/threshold/worst-element analysis "
             "(default: minSICN).",
    )
    parser.add_argument(
        "--worst-n", dest="worst_n", type=int, default=20,
        help="Number of worst elements to list per element type (default: 20).",
    )
    parser.add_argument(
        "--hist-bins", dest="hist_bins", type=int, default=10,
        help="Number of histogram bins (default: 10).",
    )
    parser.add_argument(
        "--thresholds", dest="thresholds", default="0.01,0.05,0.1,0.2,0.3",
        help="Comma-separated quality thresholds to report element counts for "
             "(default: 0.01,0.05,0.1,0.2,0.3).",
    )
    parser.add_argument(
        "--no-histogram", dest="do_histogram", action="store_false",
        help="Skip histogram output.",
    )
    parser.add_argument(
        "--no-worst", dest="do_worst", action="store_false",
        help="Skip worst-element listing.",
    )
    args = parser.parse_args()

    try:
        thresholds = parse_thresholds(args.thresholds)
    except ValueError:
        print("Error: --thresholds must be a comma-separated list of numbers.")
        sys.exit(1)

    try:
        analyze_mesh(
            args.mesh_file,
            summary_csv=args.summary_csv,
            bad_elements_csv=args.bad_elements_csv,
            primary_metric=args.primary_metric,
            worst_n=args.worst_n,
            hist_bins=args.hist_bins,
            thresholds=thresholds,
            do_histogram=args.do_histogram,
            do_worst=args.do_worst,
        )
    except Exception as exc:
        print("Error: {}".format(exc))
        sys.exit(1)


if __name__ == "__main__":
    main()
