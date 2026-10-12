import argparse
import json
import sys
import os
import glob

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import numpy as np

# ── Tweak these two constants to adjust histogram appearance ──────────────────
HIST_EDGE_COLOR = "black"     # contour line color around each bar
HIST_BAR_GAP    = 0.02        # fractional gap between bars (rwidth = 1 - gap)
# ─────────────────────────────────────────────────────────────────────────────

# Fields plotted on a log scale, each with a fixed [min, max] axis range.
tenktet_LOG_FIELDS = {
    "runtime_s":                [1e-2,  1e5],
    "n_init_vertices":          [1e1,  1e4],
    "n_init_cells":             [1e2,  1e5],
    "init_min_worst_quality":   [1e-2, 1e0],
    "init_min_relative_volume": [1e-6, 1e1],
}
HC_LOG_FIELDS = {
    "runtime_s":                [8*1e1,  2*1e5],
    "n_init_vertices":          [1e2,  1e5],
    "n_init_cells":             [2000,  2*1e5],
    "init_min_worst_quality":   [1e-8, 1e-2],
    "init_min_relative_volume": [1e-14, 1e-5],
}
#2460, 87214

LOG_FIELDS = HC_LOG_FIELDS

FIELD_LABELS = {
    "runtime_s":                  "runtime (s)",
    "n_init_vertices":            "# vertices",
    "n_init_cells":               "# cells",
    "n_final_cells":              "# cells (final)",
    "init_min_worst_quality":     "min. worst quality",
    "init_min_relative_volume":   "min. relative volume",
    "refinement_rate":            "growth ratio",
    "runtime_per_n_init_cells_s": "runtime / cell (s)",
}

OPTI_RESULT_LABELS = {
    -1: "error",
    0: "success",
    1: "failure (stuck)",
    3: "timeout",
    6: "boundary conditions failure",
    9: "success (backup)",
}

TAB10 = plt.get_cmap("tab10")

def opti_color(val):
    if val == 1:
        return TAB10(3)
    if val == 3:
        return TAB10(1)
    if val == 6:
        return opti_color(1)
    else:
        return TAB10(val % 10)

def load_all_json(directory):
    records = []
    paths = sorted(glob.glob(os.path.join(directory, "*.json")))
    if not paths:
        print(f"No JSON files found in '{directory}'.")
    for p in paths:
        with open(p) as f:
            try:
                records.append(json.load(f))
            except json.JSONDecodeError as e:
                print(f"  Warning: could not parse {p}: {e}")
    return records

def log_bins(vmin, vmax):
    lo = np.log10(vmin)
    hi = np.log10(vmax)
    n  = int(round((hi - lo) * 10))
    return 10 ** np.linspace(lo, hi, n + 1)

def apply_log_axis(ax, field, axis="x"):
    vmin, vmax = LOG_FIELDS[field]
    set_scale = ax.set_xscale if axis == "x" else ax.set_yscale
    set_lim   = ax.set_xlim   if axis == "x" else ax.set_ylim
    set_ticks = ax.xaxis      if axis == "x" else ax.yaxis
    set_scale("log")
    set_lim(vmin, vmax)
    set_ticks.set_major_locator(ticker.LogLocator(base=10, numticks=20))
    set_ticks.set_major_formatter(ticker.LogFormatterSciNotation(base=10))

def plot_colored_hist(ax, field, data, opti_result, unique_results, title):
    if field in LOG_FIELDS:
        bins = log_bins(*LOG_FIELDS[field])
    else:
        n_bins = max(8, min(30, len(data) // 2))
        bins = np.histogram_bin_edges(data, bins=n_bins)
    for val in unique_results:
        mask = opti_result == val
        ax.hist(data[mask], bins=bins, color=opti_color(val), alpha=0.8,
                label=str(val), edgecolor=HIST_EDGE_COLOR, rwidth=1 - HIST_BAR_GAP)
    if field in LOG_FIELDS:
        apply_log_axis(ax, field, axis="x")
    ax.set_xlabel(FIELD_LABELS.get(field, field), fontsize=9, labelpad=0)
    ax.set_ylabel("count", fontsize=9)
    ax.set_yscale("log")

def plot_scatter(ax, field_x, field_y, data, opti_result, unique_results, title):
    fontsize=6
    for val in unique_results:
        mask = (opti_result == val) & (~np.isnan(data[field_x])) & (~np.isnan(data[field_y]))
        ax.scatter(data[field_x][mask], data[field_y][mask],
                   color=opti_color(val), label=str(val), alpha=0.7, s=10)
    if field_x in LOG_FIELDS: apply_log_axis(ax, field_x, axis="x")
    if field_y in LOG_FIELDS: apply_log_axis(ax, field_y, axis="y")
    ax.set_xlabel(FIELD_LABELS.get(field_x, field_x), fontsize=fontsize, labelpad=0)
    ax.set_ylabel(FIELD_LABELS.get(field_y, field_y), fontsize=fontsize)
    ax.tick_params(axis='both', which='major', labelsize=fontsize)

def add_table(ax, cell_text, col_labels, row_colors=None):
    ax.axis("off")
    tbl = ax.table(cellText=cell_text, colLabels=col_labels,
                   cellLoc="center", loc="center")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(7)
    tbl.scale(1, 1.1)
    if row_colors:
        for row_idx, color in enumerate(row_colors, start=1):
            for col_idx in range(len(col_labels)):
                tbl[row_idx, col_idx].set_facecolor(color)

def main():
    parser = argparse.ArgumentParser(description="Plot optimization result statistics.")
    parser.add_argument("--input-directory", default=".", help="Input DIR.")
    parser.add_argument("--plot-output-path", default="opti_results.png", help="Output PATH.")
    parser.add_argument("--show-plots", action="store_true")
    parser.add_argument("--skip-failure-table", action="store_true",
                        help="If set, the table of non-zero results will not be plotted.")
    parser.add_argument("--skip-runtime-tables", action="store_true",
                        help="If set, the extreme runtime/cell tables will not be plotted.")
    args = parser.parse_args()

    records = load_all_json(args.input_directory)
    if not records: sys.exit("No data loaded.")

    data = {
        "runtime_s":                np.array([r.get("runtime_s",                np.nan) for r in records], dtype=float),
        "n_init_vertices":          np.array([r.get("n_init_vertices",          np.nan) for r in records], dtype=float),
        "n_init_cells":             np.array([r.get("n_init_cells",             np.nan) for r in records], dtype=float),
        "n_final_cells":            np.array([r.get("n_final_cells",            np.nan) for r in records], dtype=float),
        "init_min_worst_quality":   np.array([r.get("init_min_worst_quality",   np.nan) for r in records], dtype=float),
        "init_min_relative_volume": np.array([r.get("init_min_relative_volume", np.nan) for r in records], dtype=float),
    }


    # Derived field: growth ratio
    with np.errstate(invalid="ignore", divide="ignore"):
        data["refinement_rate"] = np.where(
            data["n_init_cells"] > 0,
            data["n_final_cells"] / data["n_init_cells"],
            np.nan)
        data["runtime_per_n_init_cells_s"] = np.where(
            data["n_init_cells"] > 0,
            data["runtime_s"] / data["n_init_cells"],
            np.nan)

    print("average runtime per cell: ", np.average(data["runtime_per_n_init_cells_s"]))
    print("    max runtime per cell: ", np.max(data["runtime_per_n_init_cells_s"]))
    print("average refinement ratio: ", np.average(data["refinement_rate"]))
    print("    max refinement ratio: ", np.max(data["refinement_rate"]))

    opti_result = np.array([r["opti_result"] for r in records], dtype=int)

    n_successes = np.count_nonzero(opti_result == 0)
    print("success rate: ", n_successes, "/", len(opti_result), "=", n_successes/len(opti_result))

    backup_mode = np.array([r.get("backup_mode", False) for r in records], dtype=bool)
    # Recode success cases with backup_mode=True as synthetic code 10
    opti_result = np.where((opti_result == 0) & backup_mode, 9, opti_result)
    mesh_name   = [r.get("mesh_name", "?") for r in records]
    unique_results = sorted(set(opti_result.tolist()))

    # Row 0: opti_result hist | runtime hist | n_init_vertices hist | n_init_cells hist
    # Row 1: 4 scatter plots
    hist_specs = [
        ("runtime_s",       "Runtime (s)"),
        ("n_init_vertices", "#initial vertices"),
        ("n_init_cells",    "#initial cells"),
    ]
    scatter_specs = [
        ("init_min_worst_quality",   "runtime_s", "min. quality vs runtime (s)"),
        ("init_min_relative_volume", "runtime_s", "min. relative volume vs runtime (s)"),
        ("n_init_vertices",          "runtime_s", "#vertices vs runtime (s)"),
        ("n_init_cells",             "runtime_s", "#cells vs runtime (s)"),
    ]
    refinement_scatter_specs = [
        ("init_min_worst_quality",   "refinement_rate", "min. quality vs growth ratio"),
        ("init_min_relative_volume", "refinement_rate", "min. relative volume vs growth ratio"),
        ("n_init_vertices",          "refinement_rate", "#vertices vs growth ratio"),
        ("n_init_cells",             "refinement_rate", "cells vs growth ratio"),
    ]

    # ── Build table data upfront so layout knows what's needed ────────────────
    tbl_col_failure = ["mesh_name", "n_init_cells", "opti_result", "runtime_s", "runtime/cell"]
    tbl_col_runtime = ["mesh_name", "n_init_cells", "runtime_s", "runtime/cell"]

    def make_rt_row(i):
        n_cells = data["n_init_cells"][i]
        rt      = data["runtime_s"][i]
        ratio   = rt / n_cells if (not np.isnan(n_cells) and n_cells > 0) else np.nan
        return [mesh_name[i],
                str(int(n_cells)) if not np.isnan(n_cells) else "?",
                f"{rt:.3f}",
                f"{ratio:.2e}" if not np.isnan(ratio) else "?"]

    def make_failure_rt_row(i):
        row = make_rt_row(i)
        # insert opti_result after n_init_cells
        return [row[0], row[1], str(opti_result[i]), row[2], row[3]]

    failure_rows_by_code = {}
    if not args.skip_failure_table:
        for code in [-1, 1, 6]:
            idx = np.where(opti_result == code)[0]
            if len(idx) > 0:
                sorted_idx = idx[np.argsort(data["n_init_cells"][idx])]
                failure_rows_by_code[code] = [make_failure_rt_row(i) for i in sorted_idx[:3]]

    fast_rows = slow_rows = None
    if not args.skip_runtime_tables:
        success_idx = np.where((opti_result == 0) | (opti_result == 10))[0]
        if len(success_idx) > 0:
            ratios     = data["runtime_s"][success_idx] / data["n_init_cells"][success_idx]
            sorted_idx = success_idx[np.argsort(ratios)]
            fast_rows  = [make_rt_row(i) for i in sorted_idx[:3]]
            slow_rows  = [make_rt_row(i) for i in sorted_idx[-3:]]

    has_failure_table  = len(failure_rows_by_code) > 0
    has_runtime_tables = fast_rows is not None
    has_any_tables     = has_failure_table or has_runtime_tables

    # ── Layout ────────────────────────────────────────────────────────────────
    DPI        = 200
    fig_width  = 16
    n_rows     = 3 + (1 if has_any_tables else 0)
    height_ratios = [1, 1, 1] + ([8.0] if has_any_tables else [])
    fig_height = (7.50 + (2.60 if has_any_tables else 0))

    fig = plt.figure(figsize=(fig_width, fig_height), dpi=DPI)

    gs = gridspec.GridSpec(n_rows, 4, figure=fig, height_ratios=height_ratios,
                           hspace=0.2, wspace=0.3,
                           top=0.93, bottom=0.02, left=0.06, right=0.98)

    ax_orhist            = fig.add_subplot(gs[0, 0])
    hist_axes            = [fig.add_subplot(gs[0, c]) for c in range(1, 4)]
    scatter_axes         = [fig.add_subplot(gs[1, c]) for c in range(4)]
    refinement_axes      = [fig.add_subplot(gs[2, c]) for c in range(4)]

    # ── opti_result histogram ─────────────────────────────────────────────────
    all_possible = list(range(-2, 10))
    counts = [int(np.sum(opti_result == v)) for v in all_possible]
    ax_orhist.bar([str(v) for v in all_possible], counts,
                  color=[opti_color(v) for v in all_possible], edgecolor=HIST_EDGE_COLOR)
    ax_orhist.set_xlabel("opti_result", fontsize=9)
    ax_orhist.set_ylabel("count", fontsize=9)
    ax_orhist.set_yscale("log")

    # ── Field histograms ──────────────────────────────────────────────────────
    for ax, (field, title) in zip(hist_axes, hist_specs):
        mask = ~np.isnan(data[field])
        plot_colored_hist(ax, field, data[field][mask], opti_result[mask], unique_results, title)

    # ── Scatter plots ─────────────────────────────────────────────────────────
    for ax, (fx, fy, title) in zip(scatter_axes, scatter_specs):
        plot_scatter(ax, fx, fy, data, opti_result, unique_results, title)

    # ── growth ratio scatter plots ─────────────────────────────────────────
    for ax, (fx, fy, title) in zip(refinement_axes, refinement_scatter_specs):
        plot_scatter(ax, fx, fy, data, opti_result, unique_results, title)

    # ── Legend ────────────────────────────────────────────────────────────────
    handles = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=opti_color(v),
                   markersize=6, label=OPTI_RESULT_LABELS.get(v, str(v)))
        for v in unique_results
    ]
    fig.legend(handles=handles, loc="upper center", ncol=len(unique_results),
               bbox_to_anchor=(0.5, 0.995), fontsize=8)

    # ── Tables: failures on left half, best/worst ratio stacked on right half ─
    if has_any_tables:
        if has_failure_table:
            codes_present = [c for c in [-1, 1, 6] if c in failure_rows_by_code]
            gs_left = gridspec.GridSpecFromSubplotSpec(
                len(codes_present), 1, subplot_spec=gs[3, :2], hspace=0.3)
            for sub_idx, code in enumerate(codes_present):
                ax_f = fig.add_subplot(gs_left[sub_idx])
                rows = failure_rows_by_code[code]
                fail_colors = [opti_color(code)[:3] + (0.15,)] * len(rows)
                add_table(ax_f, rows, tbl_col_failure, row_colors=fail_colors)

        if has_runtime_tables:
            gs_right = gridspec.GridSpecFromSubplotSpec(
                2, 1, subplot_spec=gs[3, 2:], hspace=0.3)
            ax_worst = fig.add_subplot(gs_right[0])
            ax_best  = fig.add_subplot(gs_right[1])
            add_table(ax_worst, slow_rows, tbl_col_runtime,
                      row_colors=["#f8d7da"] * len(slow_rows))
            add_table(ax_best, fast_rows, tbl_col_runtime,
                      row_colors=["#d4edda"] * len(fast_rows))

    plt.savefig(args.plot_output_path, dpi=DPI, bbox_inches="tight")
    print(f"Saved figure to {args.plot_output_path}")
    if args.show_plots: plt.show()

if __name__ == "__main__":
    main()