#!/usr/bin/env python3
"""
Plot remeshing data: scatter plots vs runtime_s in a grid layout,
plus a ternary plot of swap/collapse/split counts.

Usage:
    python plot_remeshing.py --input-dirs-file dirs.json [--show-plots]

The input JSON has the form:
    [{"dir": "/path/to/dir", "label": "My Run"}]

Dependencies:
    pip install matplotlib python-ternary
"""

import argparse
import json
import math
import time
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import ternary
import matplotlib.patches as mpatches
import numpy as np

from matplotlib.colors import LinearSegmentedColormap as _LSCM
CMAP_HEATMAP = _LSCM.from_list("YlOrRd_trunc", plt.get_cmap("YlOrRd")(np.linspace(0.15, 1.0, 256)))

# --- Tuneable display parameters ---
MARKER_SIZE              = 30     # scatter marker size (pt^2)
MARKER_ALPHA             = 0.7    # alpha for all scatter markers
MARKER_EDGE_WIDTH        = 0.4    # contour line width for all scatter markers
TERNARY_MARKER_SIZE      = 15     # ternary plot marker size (bigger than scatter)
FIGURE_WIDTH             = 20     # inches  (scatter grid figure)
FIGURE_HEIGHT            = 10     # inches  (scatter grid figure)
TERNARY_SUBPLOT_WIDTH    = 5      # inches per ternary subplot (figure width = 2x this)
TERNARY_SUBPLOT_HEIGHT   = 4      # inches per ternary subplot
TERNARY_FONT_SIZE        = 12     # font size for ternary axis labels, ticks, and legend
TERNARY_TICK_OFFSET      = 0.02  # extra offset for ternary tick labels from axes
PERF_FIG_WIDTH           = 18     # inches  (performance plot)
PERF_FIG_HEIGHT          = 6      # inches  (performance plot)
DISTORTION_FIG_WIDTH     = 10     # inches  (distortion plot)
DISTORTION_FIG_HEIGHT    = 6      # inches  (distortion plot)
DISTORTION_HSPACE_INNER  = 0.35   # vertical gap within SDE / conf_dist pairs
DISTORTION_HSPACE_OUTER  = 0.45   # vertical gap between the SDE and conf_dist groups
DISTORTION_BAR_ALPHA     = 0.7    # alpha for "after" overlay bars in distortion plots
HISTOGRAM_N_BINS         = 10     # number of bins for distortion and quality histograms
HISTOGRAM_STATS_FONTSIZE = 10      # font size for the stats box in histograms
SCATTER_FONT_SIZE           = 16  # font size for scatter grid plot
PERFORMANCE_FONT_SIZE       = 19  # font size for performance plot
DISTORTION_FONT_SIZE        = 14  # font size for distortion histograms
QUALITY_FONT_SIZE           = 22  # font size for quality histograms
QUALITY_PROGRESS_FONT_SIZE  = 14  # font size for quality progress scatter
BOUNDARY_DIST_FONT_SIZE     = 14  # font size for boundary distortion scatter
QUALITY_SIZE_FONT_SIZE      = 14  # font size for quality-size scatter
REMESHING_PERF_FONT_SIZE    = 16  # font size for remeshing performance scatter
HISTOGRAM_ALL_COLOR      = "darkgrey"  # color for the merged-all-datasets column
HISTOGRAM_RUNTIME_GRADIENT = True      # color each bar segment by runtime/cell (plasma)
USE_INPUT_QUALITY_METRICS  = True      # use input_min_* instead of init_min_* fields
PLOT_POST_OPTI_HISTOGRAMS  = True      # overlay post-opti bars on distortion histograms
HISTOGRAM_TICKS_FORMULA  = lambda lo, hi: (hi - lo) // 3 + 2
QUALITY_FIG_WIDTH        = 10     # inches  (quality plot)
QUALITY_FIG_HEIGHT       = 10    # inches  (quality plot)
QUALITY_PROGRESS_FIG_WIDTH  = 12  # inches  (quality progress scatter)
QUALITY_PROGRESS_FIG_HEIGHT = 3   # inches  (quality progress scatter)
BOUNDARY_DISTORTION_FIG_WIDTH  = 10   # inches  (boundary distortion scatter, single plot)
BOUNDARY_DISTORTION_FIG_HEIGHT = 5    # inches  (boundary distortion scatter)
QUALITY_SIZE_FIG_WIDTH   = 22     # inches  (quality-size scatter)
QUALITY_SIZE_FIG_HEIGHT  = 3      # inches  (quality-size scatter, per row)

# Field caps (values above these are clipped before display)
FIELD_CAPS = {
    "init_min_relative_edge_length":  1e7,
    "input_min_relative_edge_length": 1e7,
    "final_min_relative_edge_length": 1e5,
    "max_boundary_SDE": 1e20,
    "avg_boundary_SDE":1e15,
}

# Dataset color map: maps dataset label to a tab10 index.
# Datasets not listed are assigned their natural order (0, 1, 2, ...).
DATASET_COLOR_MAP: dict[str, int] = {
    "Knotted Holes": 6,
    "Torsion Cubes": 9,
}

# Marker mapping: opti_result -> (marker, legend label)
OPTI_MARKERS = {
    0: ("o", "success"),
    1: ("^", "stuck"),       # triangle instead of cross
    3: ("s", "timeout (24h)"),
}

# Human-readable labels for data fields
FIELD_LABELS = {
    "n_init_vertices":               "Initial Vertex Count",
    "init_min_worst_quality":        "Init Min Worst Quality",
    "init_min_relative_volume":      "Init Min Relative Volume",
    "init_min_relative_edge_length": "Init Min Relative Edge Length",
    "input_min_worst_quality":       "input minimum quality",
    "input_min_relative_volume":     "input minimum relative volume",
    "input_min_relative_edge_length":"input minimum relative edge length",
    "n_init_cells":                  "#cells",
    "n_final_cells":                 "Final Cell Count",
    "growth_ratio":                  "growth ratio",
    "max_final_SDE":                 "maximum final SDE",
    "avg_final_SDE":                 "average final SDE",
    "max_post_opti_SDE":             "Max Post-opti SDE",
    "avg_post_opti_SDE":             "Avg Post-opti SDE",
    "max_final_conf_dist":           "maximum final conformal distortion",
    "avg_final_conf_dist":           "average final conformal distortion",
    "max_post_opti_conf_dist":       "Max Post-opti Conf. Dist.",
    "avg_post_opti_conf_dist":       "Avg Post-opti Conf. Dist.",
    "max_boundary_SDE":              "maximum boundary SDE",
    "avg_boundary_SDE":              "average boundary SDE",
    "final_min_worst_quality":       "final minimum quality",
    "final_min_relative_volume":     "final minimum relative volume",
    "final_min_relative_edge_length":"final minimum relative edge length",
    "runtime_s":                     "runtime (s)",
    "runtime_per_cell":              "runtime / cell (s)",
}

# Grid layout: each row is a list of x-axis fields; y-axis is always runtime/cell
PLOT_GRID = [
    [
        "init_min_worst_quality",
        "init_min_relative_volume",
        "init_min_relative_edge_length",
    ],
    [
        "final_min_worst_quality",
        "final_min_relative_volume",
        "final_min_relative_edge_length",
    ],
]

# Fields that should use a log scale
LOG_SCALE_FIELDS = set(FIELD_LABELS.keys())

# Mapping from init_ quality fields to their input_ replacements
INIT_TO_INPUT_FIELDS = {
    "init_min_worst_quality":        "input_min_worst_quality",
    "init_min_relative_volume":      "input_min_relative_volume",
    "init_min_relative_edge_length": "input_min_relative_edge_length",
}

# Fields and labels for the remeshing operations ternary plot
TERNARY_FIELDS = ("total_swap_count", "total_collapse_count", "total_split_count")
TERNARY_LABELS = ("Swap",             "Collapse",             "Split")

# Fields and labels for the runtime breakdown ternary plot
RUNTIME_TERNARY_FIELDS = ("total_opti_time_s", "total_remeshing_time_s", "pre_processing_time_s")
RUNTIME_TERNARY_LABELS = ("optimisation time",          "remeshing time",         "pre-remeshing time")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _runtime_per_cell(r):
    """Return runtime_s / n_init_cells, or None if not available."""
    rt = r.get("runtime_s")
    n  = r.get("n_init_cells")
    if rt is None or n is None or n <= 0:
        return None
    return rt / n


def _get_field(r, field):
    """Return field value, applying cap if defined, or None."""
    v = r.get(field)
    if v is None:
        return None
    cap = FIELD_CAPS.get(field)
    if cap is not None:
        v = min(v, cap)
    return v


def _scatter(ax, xs, ys, color, marker, **kwargs):
    """Scatter with black contour edge."""
    ax.scatter(xs, ys, color=color, marker=marker,
               edgecolors="black", linewidths=MARKER_EDGE_WIDTH,
               alpha=MARKER_ALPHA, **kwargs)


def _scatter_c(ax, xs, ys, cs, marker, **kwargs):
    """Colormap scatter sorted by value (highest last), with black contour edge."""
    if not xs:
        return
    # Sort so highest-colored values are on top
    order = sorted(range(len(xs)), key=lambda i: cs[i] if isinstance(cs[i], float) else sum(cs[i][:3]))
    xs_s = [xs[i] for i in order]
    ys_s = [ys[i] for i in order]
    cs_s = [cs[i] for i in order]
    ax.scatter(xs_s, ys_s, c=cs_s, marker=marker,
               edgecolors="black", linewidths=MARKER_EDGE_WIDTH,
               alpha=MARKER_ALPHA, **kwargs)


def _runtime_colormap(all_recs, use_per_cell=True):
    """Compute norm and cmap for runtime (or runtime/cell) coloring.
    Returns (norm, cmap, vmin, vmax) or None if no valid values."""
    from matplotlib.colors import Normalize
    if use_per_cell:
        vals = [math.log10(v) for r in all_recs
                if (v := _runtime_per_cell(r)) is not None and v > 0]
    else:
        vals = [math.log10(r["runtime_s"]) for r in all_recs
                if r.get("runtime_s") and r["runtime_s"] > 0]
    if not vals:
        return None
    vmin, vmax = min(vals), max(vals)
    norm = Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.get_cmap("plasma")
    return norm, cmap, vmin, vmax


def _add_colorbar(fig, cbar_ax, cmap, norm, vmin, vmax, label, fontsize=12):
    """Add a 10^n formatted colorbar with all integer powers of ten."""
    from matplotlib.cm import ScalarMappable
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.set_label(label, fontsize=fontsize)
    ticks = list(range(math.floor(vmin), math.ceil(vmax) + 1))
    cbar.ax.yaxis.set_major_locator(plt.FixedLocator(ticks))
    cbar.ax.yaxis.set_major_formatter(
        plt.FuncFormatter(lambda x, _: fr"$10^{{{int(x)}}}$")
    )
    cbar.ax.tick_params(labelsize=fontsize)
    return cbar


def _style_ax(ax, xlabel, ylabel, xlog=True, ylog=True, fontsize=12):
    if xlog:
        ax.set_xscale("log")
    if ylog:
        ax.set_yscale("log")
    ax.xaxis.set_minor_locator(plt.NullLocator())
    ax.yaxis.set_minor_locator(plt.NullLocator())
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    ax.grid(True, which="major", linestyle="--", alpha=0.4)
    ax.tick_params(axis="both", labelsize=fontsize)


def _legend_handles(dir_labels, colors):
    return [mpatches.Patch(color=colors[i % len(colors)], label=l)
            for i, l in enumerate(dir_labels)]


def _marker_handles():
    return [plt.Line2D([0], [0], marker=m, color="black", linestyle="None",
                       markersize=7, label=l)
            for _, (m, l) in OPTI_MARKERS.items()]


def _print_field_caps(all_records, dir_labels):
    """Print max values and cap info for capped fields."""
    for field, cap in FIELD_CAPS.items():
        print(f"\n=== Field cap: {field} (cap={cap:.0e}) ===")
        overall_max = None
        for dir_label in dir_labels:
            vals = [r[field] for r in all_records[dir_label]
                    if r.get(field) is not None]
            if not vals:
                continue
            dmax = max(vals)
            n_clipped = sum(1 for v in vals if v > cap)
            print(f"  {dir_label}: max={dmax:.3e}  clipped={n_clipped}")
            overall_max = max(overall_max, dmax) if overall_max else dmax
        if overall_max:
            print(f"  Overall max: {overall_max:.3e}")


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_json_files(directory: Path, trim: int | None = None) -> list[tuple[Path, dict]]:
    data_dir = directory / "remeshing-maps-data"
    if not data_dir.exists():
        print(f"  [warning] Directory not found, skipping: {data_dir}", file=sys.stderr)
        return []
    results = []
    json_files = sorted(data_dir.glob("*.json"))
    if trim is not None:
        json_files = json_files[:trim]
    for json_file in json_files:
        try:
            with open(json_file) as f:
                data = json.load(f)
            results.append((json_file, data))
        except (json.JSONDecodeError, OSError) as e:
            print(f"  [warning] Could not read {json_file}: {e}", file=sys.stderr)
    return results


# ---------------------------------------------------------------------------
# Scatter grid
# ---------------------------------------------------------------------------

def plot_scatter_grid(all_records, dir_labels, colors, output_path):
    n_rows = len(PLOT_GRID)
    n_cols = max(len(row) for row in PLOT_GRID)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(FIGURE_WIDTH, FIGURE_HEIGHT), squeeze=False)

    for row_idx, row_fields in enumerate(PLOT_GRID):
        for col_idx, x_field in enumerate(row_fields):
            ax = axes[row_idx][col_idx]
            for dir_idx, dir_label in enumerate(dir_labels):
                color = colors[dir_idx % len(colors)]
                for opti_result, (marker, _) in OPTI_MARKERS.items():
                    xs, ys = [], []
                    for r in all_records[dir_label]:
                        if r["opti_result"] != opti_result:
                            continue
                        x = _get_field(r, x_field)
                        y = _runtime_per_cell(r)
                        if x is None or y is None:
                            continue
                        xs.append(x)
                        ys.append(y)
                    if xs:
                        _scatter(ax, xs, ys, color, marker, s=MARKER_SIZE)

            if x_field in LOG_SCALE_FIELDS:
                ax.set_xscale("log")
            ax.set_yscale("log")
            ax.xaxis.set_minor_locator(plt.NullLocator())
            ax.yaxis.set_minor_locator(plt.NullLocator())
            ax.set_xlabel(FIELD_LABELS.get(x_field, x_field), fontsize=SCATTER_FONT_SIZE)
            ax.set_ylabel(FIELD_LABELS["runtime_per_cell"], fontsize=SCATTER_FONT_SIZE)
            ax.grid(True, which="major", linestyle="--", alpha=0.4)
            ax.tick_params(axis="both", labelsize=SCATTER_FONT_SIZE)

    fig.legend(handles=_legend_handles(dir_labels, colors) + _marker_handles(),
               title="Datasets / Result",
               loc="upper right", bbox_to_anchor=(1.0, 1.0),
               borderaxespad=0.5, fontsize=8, title_fontsize=9)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"Scatter grid saved to: {output_path}")


# ---------------------------------------------------------------------------
# Ternary
# ---------------------------------------------------------------------------

def _setup_tax(fig, pos, labels):
    sub_ax = fig.add_subplot(pos)
    _, tax = ternary.figure(ax=sub_ax, scale=1.0)
    tax.boundary(linewidth=1.5)
    tax.gridlines(multiple=0.2, color="gray", linewidth=0.5, alpha=0.5)
    tax.bottom_axis_label(labels[0], fontsize=TERNARY_FONT_SIZE, offset=0.06)
    tax.right_axis_label(labels[1], fontsize=TERNARY_FONT_SIZE, offset=0.14)
    tax.left_axis_label(labels[2], fontsize=TERNARY_FONT_SIZE, offset=0.14)
    tax.ticks(axis="lbr", multiple=0.2, linewidth=1.5, tick_formats="%.1f",
              fontsize=TERNARY_FONT_SIZE - 2, offset=TERNARY_TICK_OFFSET)
    tax.get_axes().axis("off")
    return tax


def _collect_ternary_points(all_records, dir_labels, fields):
    points, runtimes, dir_idxs = [], [], []
    no_remeshing = {label: 0 for label in dir_labels}
    for dir_idx, dir_label in enumerate(dir_labels):
        for r in all_records[dir_label]:
            if r["opti_result"] != 0:
                continue
            a, b, c = r.get(fields[0]), r.get(fields[1]), r.get(fields[2])
            if a is None or b is None or c is None:
                continue
            total = a + b + c
            if total == 0:
                no_remeshing[dir_label] += 1
                continue
            points.append((a / total, b / total, c / total))
            runtimes.append(r["runtime_s"])
            dir_idxs.append(dir_idx)
    return points, runtimes, dir_idxs, no_remeshing


def plot_ternary(all_records, dir_labels, colors, output_path,
                 fields=TERNARY_FIELDS, labels=TERNARY_LABELS,
                 zero_label="cases without remeshing"):
    from matplotlib.colors import Normalize
    from matplotlib.cm import ScalarMappable

    marker, _ = OPTI_MARKERS[0]
    _, _, _, no_remeshing = _collect_ternary_points(all_records, dir_labels, fields)

    print(f"\n=== Ternary plot: {zero_label} (all counts zero) ===")
    for dir_label in dir_labels:
        print(f"  {dir_label}: {no_remeshing[dir_label]} {zero_label}")

    # Left: colored by dataset; Middle: n_init_cells (cividis); Right: runtime/cell (plasma)
    pts_by_dir = {di: [] for di in range(len(dir_labels))}
    pts_mid, ncells_list = [], []
    pts_right, rtpc_list = [], []

    for dir_idx, dir_label in enumerate(dir_labels):
        for r in all_records[dir_label]:
            if r["opti_result"] != 0:
                continue
            a, b, c = r.get(fields[0]), r.get(fields[1]), r.get(fields[2])
            if a is None or b is None or c is None:
                continue
            total = a + b + c
            if total == 0:
                continue
            pt = (a / total, b / total, c / total)
            pts_by_dir[dir_idx].append(pt)
            nc = r.get("n_init_cells")
            if nc is not None and nc > 0:
                pts_mid.append(pt)
                ncells_list.append(nc)
            rtpc = _runtime_per_cell(r)
            if rtpc is not None and rtpc > 0:
                pts_right.append(pt)
                rtpc_list.append(rtpc)

    from matplotlib.colors import Normalize

    def _sort_color(pts, vals, cmap_name):
        if not pts:
            return [], [], None, None, None, None
        log_vals = [math.log10(v) for v in vals]
        vmin, vmax = min(log_vals), max(log_vals)
        norm = Normalize(vmin=vmin, vmax=vmax)
        cmap = plt.get_cmap(cmap_name)
        order = sorted(range(len(pts)), key=lambda i: log_vals[i])
        return ([pts[i] for i in order], [cmap(norm(log_vals[i])) for i in order],
                norm, cmap, vmin, vmax)

    pts_ms, cms, norm_m, cmap_m, vmin_m, vmax_m = _sort_color(pts_mid,   ncells_list, "hot")
    pts_rs, crs, norm_r, cmap_r, vmin_r, vmax_r = _sort_color(pts_right, rtpc_list,   "plasma")

    TERNARY_CBAR_SPACER = 0.12
    fig = plt.figure(figsize=(TERNARY_SUBPLOT_WIDTH * 3, TERNARY_SUBPLOT_HEIGHT))
    gs = fig.add_gridspec(1, 8, width_ratios=[1, 0.06, TERNARY_CBAR_SPACER,
                                              1, 0.06, TERNARY_CBAR_SPACER,
                                              1, 0.06], wspace=0.25)
    tax_left  = _setup_tax(fig, gs[0], labels)
    _sp1 = fig.add_subplot(gs[1]); _sp1.set_visible(False)
    _sp2 = fig.add_subplot(gs[2]); _sp2.set_visible(False)
    tax_mid   = _setup_tax(fig, gs[3], labels)
    cbar_mid  = fig.add_subplot(gs[4])
    _sp3 = fig.add_subplot(gs[5]); _sp3.set_visible(False)
    tax_right = _setup_tax(fig, gs[6], labels)
    cbar_right = fig.add_subplot(gs[7])

    for dir_idx in range(len(dir_labels)):
        for pt in pts_by_dir[dir_idx]:
            tax_left.scatter([pt], marker=marker, color=colors[dir_idx % len(colors)],
                             s=TERNARY_MARKER_SIZE, alpha=MARKER_ALPHA, zorder=2,
                             edgecolors="black", linewidths=MARKER_EDGE_WIDTH)
    tax_left.get_axes().legend(handles=_legend_handles(dir_labels, colors),
                               loc="upper right",
                               bbox_to_anchor=(1.5, 1.1),
                               fontsize=TERNARY_FONT_SIZE - 2)

    for pt, c in zip(pts_ms, cms):
        tax_mid.scatter([pt], marker=marker, color=c,
                        s=TERNARY_MARKER_SIZE, alpha=MARKER_ALPHA, zorder=2,
                        edgecolors="black", linewidths=MARKER_EDGE_WIDTH)
    if pts_ms:
        _add_colorbar(fig, cbar_mid, cmap_m, norm_m, vmin_m, vmax_m,
                      FIELD_LABELS["n_init_cells"])

    for pt, c in zip(pts_rs, crs):
        tax_right.scatter([pt], marker=marker, color=c,
                          s=TERNARY_MARKER_SIZE, alpha=MARKER_ALPHA, zorder=2,
                          edgecolors="black", linewidths=MARKER_EDGE_WIDTH)
    if pts_rs:
        _add_colorbar(fig, cbar_right, cmap_r, norm_r, vmin_r, vmax_r,
                      FIELD_LABELS["runtime_per_cell"])

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"Ternary plot saved to: {output_path}")


def plot_remeshing_ternary(all_records, dir_labels, output_path,
                           fields=TERNARY_FIELDS, labels=TERNARY_LABELS,
                           zero_label="cases without remeshing"):
    """Ternary plot pair specific to remeshing operations:
      left  — colored by total_remeshing_time_s (viridis)
      right — colored by runtime/cell (plasma)
    """
    from matplotlib.colors import Normalize

    marker, _ = OPTI_MARKERS[0]
    _, _, _, no_remeshing = _collect_ternary_points(all_records, dir_labels, fields)

    print(f"\n=== Remeshing ternary plot: {zero_label} (all counts zero) ===")
    for dir_label in dir_labels:
        print(f"  {dir_label}: {no_remeshing[dir_label]} {zero_label}")

    pts_left, rmt_list = [], []
    pts_right, rtpc_list = [], []

    for dir_label in dir_labels:
        for r in all_records[dir_label]:
            if r["opti_result"] != 0:
                continue
            a, b, c = r.get(fields[0]), r.get(fields[1]), r.get(fields[2])
            if a is None or b is None or c is None:
                continue
            total = a + b + c
            if total == 0:
                continue
            pt = (a / total, b / total, c / total)
            rmt = r.get("total_remeshing_time_s")
            if rmt is not None and rmt > 0:
                pts_left.append(pt)
                rmt_list.append(rmt)
            rtpc = _runtime_per_cell(r)
            if rtpc is not None and rtpc > 0:
                pts_right.append(pt)
                rtpc_list.append(rtpc)

    def _sorted_colored(pts, vals, cmap_name):
        if not pts:
            return [], [], None, None, None, None
        log_vals = [math.log10(v) for v in vals]
        vmin, vmax = min(log_vals), max(log_vals)
        norm = Normalize(vmin=vmin, vmax=vmax)
        cmap = plt.get_cmap(cmap_name)
        order = sorted(range(len(pts)), key=lambda i: log_vals[i])
        return ([pts[i] for i in order],
                [cmap(norm(log_vals[i])) for i in order],
                norm, cmap, vmin, vmax)

    pts_l, cl, norm_l, cmap_l, vmin_l, vmax_l = _sorted_colored(pts_left,  rmt_list,  "viridis")
    pts_r, cr, norm_r, cmap_r, vmin_r, vmax_r = _sorted_colored(pts_right, rtpc_list, "plasma")

    TERNARY_CBAR_SPACER = 0.12
    fig = plt.figure(figsize=(TERNARY_SUBPLOT_WIDTH * 2 + 1.0, TERNARY_SUBPLOT_HEIGHT))
    gs = fig.add_gridspec(1, 5, width_ratios=[1, 0.06, TERNARY_CBAR_SPACER,
                                              1, 0.06], wspace=0.25)
    tax_left  = _setup_tax(fig, gs[0], labels)
    cbar_left = fig.add_subplot(gs[1])
    _sp1 = fig.add_subplot(gs[2]); _sp1.set_visible(False)
    tax_right = _setup_tax(fig, gs[3], labels)
    cbar_right = fig.add_subplot(gs[4])

    for pt, c in zip(pts_l, cl):
        tax_left.scatter([pt], marker=marker, color=c,
                         s=TERNARY_MARKER_SIZE, alpha=MARKER_ALPHA, zorder=2,
                         edgecolors="black", linewidths=MARKER_EDGE_WIDTH)
    if pts_l:
        _add_colorbar(fig, cbar_left, cmap_l, norm_l, vmin_l, vmax_l, "remeshing time (s)")

    for pt, c in zip(pts_r, cr):
        tax_right.scatter([pt], marker=marker, color=c,
                          s=TERNARY_MARKER_SIZE, alpha=MARKER_ALPHA, zorder=2,
                          edgecolors="black", linewidths=MARKER_EDGE_WIDTH)
    if pts_r:
        _add_colorbar(fig, cbar_right, cmap_r, norm_r, vmin_r, vmax_r,
                      FIELD_LABELS["runtime_per_cell"])

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"Remeshing ternary plot saved to: {output_path}")


# ---------------------------------------------------------------------------
# Performance
# ---------------------------------------------------------------------------

def _scatter_on_ax(ax, all_records, dir_labels, colors, get_xy):
    for dir_idx, dir_label in enumerate(dir_labels):
        color = colors[dir_idx % len(colors)]
        for opti_result, (marker, _) in OPTI_MARKERS.items():
            xs, ys = [], []
            for r in all_records[dir_label]:
                if r["opti_result"] != opti_result:
                    continue
                xy = get_xy(r)
                if xy is None:
                    continue
                xs.append(xy[0])
                ys.append(xy[1])
            if xs:
                _scatter(ax, xs, ys, color, marker, s=MARKER_SIZE)


def plot_performance(all_records, dir_labels, colors, output_path):
    fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(PERF_FIG_WIDTH, PERF_FIG_HEIGHT))

    def xy_runtime(r):
        x = r.get("n_init_cells")
        y = r.get("runtime_s")
        return (x, y) if x is not None and y is not None else None

    def xy_growth(r):
        n_init  = r.get("n_init_cells")
        n_final = r.get("n_final_cells")
        if n_init is None or n_final is None or n_init == 0:
            return None
        return (n_init, n_final / n_init)

    _scatter_on_ax(ax_left,  all_records, dir_labels, colors, xy_runtime)
    _scatter_on_ax(ax_right, all_records, dir_labels, colors, xy_growth)
    _style_ax(ax_left,  FIELD_LABELS["n_init_cells"], FIELD_LABELS["runtime_s"], fontsize=PERFORMANCE_FONT_SIZE)
    _style_ax(ax_right, FIELD_LABELS["n_init_cells"], FIELD_LABELS["growth_ratio"], fontsize=PERFORMANCE_FONT_SIZE)

    # Stats box for growth ratio
    all_growth = [n_final / n_init
                  for dl in dir_labels for r in all_records[dl]
                  if (n_init := r.get("n_init_cells")) and (n_final := r.get("n_final_cells"))
                  and n_init > 0]
    if all_growth:
        g_arr = np.array(all_growth)
        stats_text = (r"$\mu$" + f"={g_arr.mean():.2f}\n"
                      + r"$\sigma$" + f"={g_arr.std():.2f}")
        ax_right.text(0.98, 0.98, stats_text,
                      transform=ax_right.transAxes,
                      fontsize=PERFORMANCE_FONT_SIZE, family="monospace",
                      verticalalignment="top", horizontalalignment="right",
                      bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                                edgecolor="gray", alpha=0.8))

    fig.legend(handles=_legend_handles(dir_labels, colors) + _marker_handles(),
               title="Datasets / Result",
               loc="upper left", bbox_to_anchor=(1., 0.9),
               borderaxespad=0, fontsize=PERFORMANCE_FONT_SIZE, title_fontsize=PERFORMANCE_FONT_SIZE)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"Performance plot saved to: {output_path}")


def _xy_runtime(r):
    x = r.get("n_init_cells"); y = r.get("runtime_s")
    return (x, y) if x is not None and y is not None else None

def _xy_growth(r):
    n_init  = r.get("n_init_cells"); n_final = r.get("n_final_cells")
    if n_init is None or n_final is None or n_init == 0: return None
    return (n_init, n_final / n_init)

def _growth_stats_box(ax, all_records, dir_labels):
    all_growth = [n_final / n_init
                  for dl in dir_labels for r in all_records[dl]
                  if (n_init := r.get("n_init_cells")) and (n_final := r.get("n_final_cells"))
                  and n_init > 0]
    if not all_growth:
        return
    g_arr = np.array(all_growth)
    stats_text = (r"$\mu$" + f"={g_arr.mean():.2f}\n" + r"$\sigma$" + f"={g_arr.std():.2f}")
    ax.text(0.98, 0.98, stats_text,
            transform=ax.transAxes, fontsize=PERFORMANCE_FONT_SIZE, family="monospace",
            verticalalignment="top", horizontalalignment="right",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="gray", alpha=0.8))


def plot_performance_heatmap(all_records, dir_labels, colors, output_path):
    fig, axes = plt.subplots(2, 2, figsize=(PERF_FIG_WIDTH, PERF_FIG_HEIGHT * 2))
    ax_left, ax_right = axes[0]
    ax_left_h, ax_right_h = axes[1]

    # Row 1: scatter
    _scatter_on_ax(ax_left,  all_records, dir_labels, colors, _xy_runtime)
    _scatter_on_ax(ax_right, all_records, dir_labels, colors, _xy_growth)
    _style_ax(ax_left,  FIELD_LABELS["n_init_cells"], FIELD_LABELS["runtime_s"], fontsize=PERFORMANCE_FONT_SIZE)
    _style_ax(ax_right, FIELD_LABELS["n_init_cells"], FIELD_LABELS["growth_ratio"], fontsize=PERFORMANCE_FONT_SIZE)
    _growth_stats_box(ax_right, all_records, dir_labels)

    # Row 2: heatmaps
    all_recs = [r for dl in dir_labels for r in all_records[dl]]

    # Pre-compute global h_max across both heatmaps so colors are comparable
    _perf_h_max = 1
    for _get_xy in [_xy_runtime, _xy_growth]:
        _pts = [_get_xy(r) for r in all_recs]
        _pts = [(x, y) for xy in _pts if xy is not None for x, y in [xy] if x > 0 and y > 0]
        if _pts:
            _xl = np.array([math.log10(x) for x, y in _pts])
            _yl = np.array([math.log10(y) for x, y in _pts])
            _h, _, _ = np.histogram2d(_xl, _yl, bins=80,
                                      range=[[_xl.min(), _xl.max()], [_yl.min(), _yl.max()]])
            _perf_h_max = max(_perf_h_max, int(_h.max()))

    def _perf_heatmap(ax, get_xy, xlabel, ylabel):
        pts = [get_xy(r) for r in all_recs]
        pts = [(x, y) for xy in pts if xy is not None for x, y in [xy] if x > 0 and y > 0]
        if not pts:
            ax.set_visible(False); return
        xs_log = np.array([math.log10(x) for x, y in pts])
        ys_log = np.array([math.log10(y) for x, y in pts])
        n_bins = 80
        lo_x, hi_x = xs_log.min(), xs_log.max()
        lo_y, hi_y = ys_log.min(), ys_log.max()
        h, xedges, yedges = np.histogram2d(xs_log, ys_log, bins=n_bins,
                                           range=[[lo_x, hi_x], [lo_y, hi_y]])
        h_log_vals = np.where(h > 0, np.log10(h + 1), 0)
        h_norm = np.where(h > 0, h_log_vals / math.log10(_perf_h_max + 1), 0)
        h_masked = np.ma.masked_where(h == 0, h_norm)
        ax.pcolormesh(xedges, yedges, h_masked.T, cmap=CMAP_HEATMAP,
                      shading="auto", zorder=1, vmin=0, vmax=1)
        ax.set_xlim(lo_x, hi_x); ax.set_ylim(lo_y, hi_y)
        for axis, lo, hi in [(ax.xaxis, lo_x, hi_x), (ax.yaxis, lo_y, hi_y)]:
            ticks = list(range(math.floor(lo), math.ceil(hi) + 1))
            axis.set_major_locator(plt.FixedLocator(ticks))
            axis.set_major_formatter(plt.FuncFormatter(lambda v, _: fr"$10^{{{int(v)}}}$"))
            axis.set_minor_locator(plt.NullLocator())
        ax.set_xlabel(xlabel, fontsize=PERFORMANCE_FONT_SIZE)
        ax.set_ylabel(ylabel, fontsize=PERFORMANCE_FONT_SIZE)
        ax.tick_params(axis="both", labelsize=PERFORMANCE_FONT_SIZE)
        ax.grid(True, which="major", linestyle="--", alpha=0.4)

    _perf_heatmap(ax_left_h,  _xy_runtime, FIELD_LABELS["n_init_cells"], FIELD_LABELS["runtime_s"])
    _perf_heatmap(ax_right_h, _xy_growth,  FIELD_LABELS["n_init_cells"], FIELD_LABELS["growth_ratio"])
    _growth_stats_box(ax_right_h, all_records, dir_labels)

    # Colorbar for heatmap row
    from matplotlib.colors import Normalize as _Norm
    from matplotlib.cm import ScalarMappable as _SM
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    _sm = _SM(cmap=CMAP_HEATMAP, norm=_Norm(vmin=0, vmax=1))
    _sm.set_array([])
    _divider = make_axes_locatable(ax_right_h)
    _cax = _divider.append_axes("right", size="4%", pad=0.08)
    _cbar = fig.colorbar(_sm, cax=_cax)
    _cbar.set_label("count (log scale)", fontsize=PERFORMANCE_FONT_SIZE)
    _cbar.ax.yaxis.set_major_locator(plt.FixedLocator([0, 1.0]))
    _cbar.ax.yaxis.set_major_formatter(plt.FuncFormatter(
        lambda v, _: "1" if v <= 0 else str(_perf_h_max)
    ))
    _cbar.ax.tick_params(labelsize=PERFORMANCE_FONT_SIZE)

    fig.legend(handles=_legend_handles(dir_labels, colors) + _marker_handles(),
               title="Datasets / Result",
               loc="upper left", bbox_to_anchor=(1., 0.9),
               borderaxespad=0, fontsize=PERFORMANCE_FONT_SIZE, title_fontsize=PERFORMANCE_FONT_SIZE)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"Performance heatmap plot saved to: {output_path}")


# ---------------------------------------------------------------------------
# Distortion histograms
# ---------------------------------------------------------------------------

DISTORTION_ROWS = [
    ("max_final_SDE",       "max_post_opti_SDE"),
    ("avg_final_SDE",       "avg_post_opti_SDE"),
    ("max_final_conf_dist", "max_post_opti_conf_dist"),
    ("avg_final_conf_dist", "avg_post_opti_conf_dist"),
]
DISTORTION_TIGHT_PAIRS = [(0, 1), (2, 3)]
DISTORTION_CLAMP_LOG   = -3   # values <= 10^this are grouped into the leftmost bin


def _distortion_histogram(ax, values, color, n_total=None, overlay_values=None, overlay_color=None, rtpc_values=None, overlay_rtpc_values=None, fontsize=12, clamp_lo=None, clamp_hi=None):
    """Draw a ratio-of-total histogram on ax.
    X ticks: 3 if the integer log range is even, 4 if odd.
    Y ticks: all integer powers of ten from min ratio up to 10^0.
    """
    if not values:
        ax.set_visible(False)
        return

    lo = math.floor(min(values))
    hi = math.ceil(max(values))
    if clamp_hi is not None:
        hi = max(hi, clamp_hi)
    counts, bin_edges = np.histogram(values, bins=HISTOGRAM_N_BINS, range=(lo, hi))

    N = n_total if n_total and n_total > 0 else len(values)
    ratios = np.where(counts > 0, counts / N, np.nan)

    bin_step    = bin_edges[1] - bin_edges[0]
    has_overlay = bool(overlay_values)
    # Split bins only when there is an overlay; otherwise use full width
    if has_overlay:
        half = bin_step * 0.85 / 2
        before_width  = half
        before_offset = 0.0
        after_width   = half
        after_offset  = half
    else:
        before_width  = bin_step * 0.85
        before_offset = 0.0
        after_width   = 0.0
        after_offset  = 0.0

    def _draw_bars(positions, bar_ratios, bw, color_, rtpc_vals_=None, alpha=1.0, source_values=None):
        src_vals = source_values if source_values is not None else values
        if HISTOGRAM_RUNTIME_GRADIENT and rtpc_vals_ is not None and len(rtpc_vals_) == len(src_vals):
            from matplotlib.colors import Normalize as _Norm
            _log_rtpc = np.array([math.log10(v) for v in rtpc_vals_ if v is not None and v > 0])
            if len(_log_rtpc) > 0:
                _vmin, _vmax = _log_rtpc.min(), _log_rtpc.max()
                _norm = _Norm(vmin=_vmin, vmax=_vmax)
                _cmap = plt.get_cmap("plasma")
                bin_indices = np.clip(
                    np.digitize(src_vals, bin_edges[:-1], right=False) - 1,
                    0, HISTOGRAM_N_BINS - 1)
                for b_idx, x0 in enumerate(positions):
                    mask = bin_indices == b_idx
                    if not mask.any():
                        continue
                    bin_rtpc = np.sort([math.log10(rtpc_vals_[i])
                                        for i in range(len(src_vals))
                                        if mask[i] and rtpc_vals_[i] is not None and rtpc_vals_[i] > 0])
                    if len(bin_rtpc) == 0:
                        continue
                    seg_h = 1.0 / N
                    bottom_ = 0.0
                    for rv in bin_rtpc:
                        ax.bar(x0, seg_h, width=bw, align="edge",
                               color=_cmap(_norm(rv)), edgecolor="none",
                               linewidth=0, bottom=bottom_, alpha=alpha)
                        bottom_ += seg_h
                    ax.bar(x0, bottom_, width=bw, align="edge",
                           color="none", edgecolor="black", linewidth=0.5, bottom=0)
                return
        # Plain colored bars
        ax.bar(positions, bar_ratios, width=bw, align="edge",
               color=color_, edgecolor="black", linewidth=0.5, alpha=alpha)

    before_positions = bin_edges[:-1] + before_offset
    _draw_bars(before_positions, ratios, before_width, color, rtpc_values)

    if has_overlay:
        overlay_counts, _ = np.histogram(overlay_values, bins=bin_edges)
        overlay_ratios = np.where(overlay_counts > 0, overlay_counts / N, np.nan)
        after_positions = bin_edges[:-1] + after_offset
        ov_color = overlay_color if overlay_color else color
        overlay_vals_arr = np.array(overlay_values) if overlay_values else None
        _draw_bars(after_positions, overlay_ratios, after_width, ov_color,
                   rtpc_vals_=overlay_rtpc_values, alpha=0.5,
                   source_values=overlay_vals_arr)

    ax.set_yscale("log")
    ax.yaxis.set_minor_locator(plt.NullLocator())

    # Y ticks: all integer powers from floor(log10(min ratio)) up to 0 (inclusive)
    # Force y upper limit to 1 so 10^0 is always visible
    ymin_val = np.nanmin(ratios[ratios > 0]) if np.any(ratios > 0) else 1e-4
    y_lo = math.floor(math.log10(ymin_val))
    y_ticks = list(range(y_lo, 1))  # range(y_lo, 1) gives [..., -1, 0]
    if not y_ticks:
        y_ticks = [y_lo, 0]
    if 0 not in y_ticks:
        y_ticks.append(0)
    ax.set_ylim(bottom=10**y_lo * 0.5, top=2.0)  # top > 1 ensures 10^0 tick is visible
    ax.yaxis.set_major_locator(plt.FixedLocator([10**t for t in y_ticks]))
    ax.yaxis.set_major_formatter(
        plt.FuncFormatter(lambda y, _: f"{y * 100:.4g}%")
    )

    # X axis: always lo as first tick, then every x_step, then one tick at or past hi.
    # clamp_hi is always forced as a tick if set.
    x_step = 2 if (hi - lo) <= 8 else 3
    labels = [lo]
    nxt = lo + x_step
    while nxt <= hi:
        labels.append(nxt)
        nxt += x_step
    # If last label didn't reach hi, add one more step past hi
    if labels[-1] < hi:
        labels.append(labels[-1] + x_step)
    # Guarantee at least two ticks
    if len(labels) < 2:
        labels.append(labels[0] + x_step)
    ax.xaxis.set_major_locator(plt.FixedLocator(labels))
    def _xfmt(x, _):
        xi = round(x)
        if clamp_lo is not None and xi == clamp_lo:
            return fr"$\leq10^{{{xi}}}$"
        if clamp_hi is not None and xi == clamp_hi:
            return fr"$\geq10^{{{xi}}}$"
        return fr"$10^{{{xi}}}$"
    ax.xaxis.set_major_formatter(plt.FuncFormatter(_xfmt))
    ax.tick_params(axis="both", labelsize=fontsize)
    ax.grid(True, which="major", axis="y", linestyle="--", alpha=0.4)

    # Stats box: avg and std for baseline (and overlay if present)
    vals_arr = np.array(values)
    lines = [r"$\mu_{\text{log}}$" + f"={vals_arr.mean():.2f}  " + r"$\sigma_{\text{log}}$" + f"={vals_arr.std():.2f}"]
    if overlay_values:
        ov = np.array(overlay_values)
        lines.append(r"$\mu^{\prime}_{\text{log}}$" + f"={ov.mean():.2f}  " + r"$\sigma^{\prime}_{\text{log}}$" + f"={ov.std():.2f}")
    ax.text(0.98, 0.98, "\n".join(lines),
            transform=ax.transAxes, fontsize=HISTOGRAM_STATS_FONTSIZE,
            verticalalignment="top", horizontalalignment="right",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                      edgecolor="gray", alpha=0.8))


def _add_histogram_runtime_colorbar(fig, all_recs, anchor_ax, fontsize=12):
    """Add a plasma colorbar attached to the last of anchor_ax."""
    from matplotlib.colors import Normalize
    from matplotlib.cm import ScalarMappable
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    rtpc_vals = [math.log10(v) for r in all_recs
                 if (v := _runtime_per_cell(r)) is not None and v > 0]
    if not rtpc_vals:
        return
    vmin, vmax = min(rtpc_vals), max(rtpc_vals)
    norm = Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.get_cmap("plasma")
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    last_ax = anchor_ax[-1] if isinstance(anchor_ax, list) else anchor_ax
    divider = make_axes_locatable(last_ax)
    cax = divider.append_axes("right", size="8%", pad=0.05)
    cbar = fig.colorbar(sm, cax=cax)
    cbar.set_label(FIELD_LABELS["runtime_per_cell"], fontsize=fontsize)
    ticks = list(range(math.floor(vmin), math.ceil(vmax) + 1))
    cbar.ax.yaxis.set_major_locator(plt.FixedLocator(ticks))
    cbar.ax.yaxis.set_major_formatter(
        plt.FuncFormatter(lambda x, _: fr"$10^{{{int(x)}}}$")
    )
    cbar.ax.tick_params(labelsize=fontsize)


def plot_distortion(all_records, dir_labels, colors, output_path, show_post_opti=False, only_overall=False):
    all_recs_merged = [r for dl in dir_labels for r in all_records[dl]]
    n_rows = len(DISTORTION_ROWS)

    fontsize = DISTORTION_FONT_SIZE
    if only_overall:
        # 2x2 grid: [avg_SDE, max_SDE / avg_conf, max_conf]
        # Reorder rows for 2x2: avg_SDE, max_SDE, avg_conf, max_conf
        rows_2x2 = [
            ("avg_final_SDE",       "avg_post_opti_SDE"),
            ("max_final_SDE",       "max_post_opti_SDE"),
            ("avg_final_conf_dist", "avg_post_opti_conf_dist"),
            ("max_final_conf_dist", "max_post_opti_conf_dist"),
        ]
        fig, ax_grid = plt.subplots(2, 2,
                                    figsize=(DISTORTION_FIG_WIDTH, DISTORTION_FIG_HEIGHT),
                                    squeeze=False)
        axes_all = [ax_grid[0][0], ax_grid[0][1], ax_grid[1][0], ax_grid[1][1]]
    else:
        n_cols = len(dir_labels) + 1  # +1 for merged-all column
        fig = plt.figure(figsize=(DISTORTION_FIG_WIDTH * n_cols / max(len(dir_labels), 1), DISTORTION_FIG_HEIGHT))
        outer = fig.add_gridspec(2, 1, hspace=DISTORTION_HSPACE_OUTER)
        axes = [[None] * n_cols for _ in range(n_rows)]
        for group_idx, row_pair in enumerate([(0, 1), (2, 3)]):
            inner = outer[group_idx].subgridspec(2, n_cols, hspace=DISTORTION_HSPACE_INNER)
            for inner_row, row_idx in enumerate(row_pair):
                for col_idx in range(n_cols):
                    axes[row_idx][col_idx] = fig.add_subplot(inner[inner_row, col_idx])

    _rows = rows_2x2 if only_overall else DISTORTION_ROWS
    for row_idx, (field_before, field_after) in enumerate(_rows):
        if not only_overall:
            for col_idx, dir_label in enumerate(dir_labels):
                ax    = axes[row_idx][col_idx]
                color = colors[col_idx % len(colors)]
                recs  = all_records[dir_label]
                n_total = len(recs)
                vals_before = [max(math.log10(r[field_before]), DISTORTION_CLAMP_LOG)
                               for r in recs if r.get(field_before) and r[field_before] > 0]
                vals_after  = [max(math.log10(r[field_after]), DISTORTION_CLAMP_LOG)
                               for r in recs if r.get(field_after) and r[field_after] > 0]
                _cap = FIELD_CAPS.get(field_before)
                _cap_hi = int(math.log10(_cap)) if _cap else None
                _distortion_histogram(ax, vals_before, color, n_total=n_total,
                                      overlay_values=vals_after if show_post_opti else None,
                                      overlay_color=color, fontsize=fontsize,
                                      clamp_lo=DISTORTION_CLAMP_LOG, clamp_hi=_cap_hi)
                if row_idx == 0:
                    ax.set_title(dir_label, fontsize=fontsize)

        # Merged-all column (or only column when only_overall)
        ax_all = axes_all[row_idx] if only_overall else axes[row_idx][len(dir_labels)]
        n_total_all = len(all_recs_merged)
        tmp_all = [(max(math.log10(r[field_before]), DISTORTION_CLAMP_LOG), _runtime_per_cell(r))
                   for r in all_recs_merged if r.get(field_before) and r[field_before] > 0]
        vals_all_before = [v[0] for v in tmp_all]
        rtpc_all_dist   = [v[1] for v in tmp_all]
        vals_all_after  = [max(math.log10(r[field_after]), DISTORTION_CLAMP_LOG)
                           for r in all_recs_merged if r.get(field_after) and r[field_after] > 0]
        overlay_rtpc_all = [_runtime_per_cell(r)
                            for r in all_recs_merged if r.get(field_after) and r[field_after] > 0]
        _cap = FIELD_CAPS.get(field_before)
        _cap_hi = int(math.log10(_cap)) if _cap else None
        _distortion_histogram(ax_all, vals_all_before, HISTOGRAM_ALL_COLOR,
                              n_total=n_total_all,
                              overlay_values=vals_all_after if show_post_opti else None,
                              overlay_color=HISTOGRAM_ALL_COLOR,
                              rtpc_values=rtpc_all_dist,
                              overlay_rtpc_values=overlay_rtpc_all if show_post_opti else None,
                              fontsize=fontsize, clamp_lo=DISTORTION_CLAMP_LOG, clamp_hi=_cap_hi)
        ax_all.set_title(FIELD_LABELS.get(field_before, field_before) if only_overall else "All",
                         fontsize=fontsize)

        if not only_overall:
            axes[row_idx][0].set_ylabel(
                FIELD_LABELS.get(field_before, field_before), fontsize=fontsize, rotation=90, labelpad=6)

    all_col_ax_list = axes_all if only_overall else [axes[i][len(dir_labels)] for i in range(n_rows)]
    if HISTOGRAM_RUNTIME_GRADIENT:
        _add_histogram_runtime_colorbar(fig, all_recs_merged, all_col_ax_list, fontsize=fontsize)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"Distortion plot saved to: {output_path}")


# ---------------------------------------------------------------------------
# Quality histograms
# ---------------------------------------------------------------------------

QUALITY_ROWS = [
    "final_min_worst_quality",
    "final_min_relative_volume",
    "final_min_relative_edge_length",
]


def _make_quality_axes(fig, n_rows, n_data_cols, with_cbar):
    """Build a gridspec for quality plots: n_data_cols histogram cols + optional cbar col."""
    if with_cbar:
        gs = fig.add_gridspec(n_rows, n_data_cols + 1,
                              width_ratios=[1] * n_data_cols + [0.04],
                              wspace=0.35, hspace=0.4)
        axes = [[fig.add_subplot(gs[r, c]) for c in range(n_data_cols)]
                for r in range(n_rows)]
        cbar_axes = [fig.add_subplot(gs[r, n_data_cols]) for r in range(n_rows)]
        # Merge cbar axes visually by hiding all but the first, then use fig.add_axes later
        return axes, cbar_axes
    else:
        gs = fig.add_gridspec(n_rows, n_data_cols, wspace=0.35, hspace=0.4)
        axes = [[fig.add_subplot(gs[r, c]) for c in range(n_data_cols)]
                for r in range(n_rows)]
        return axes, None


def _add_histogram_runtime_colorbar_gs(fig, all_recs, cbar_axes, fontsize=12):
    """Add a plasma colorbar using a dedicated gridspec column spanning all rows."""
    from matplotlib.colors import Normalize
    from matplotlib.cm import ScalarMappable
    from matplotlib.gridspec import SubplotSpec
    rtpc_vals = [math.log10(v) for r in all_recs
                 if (v := _runtime_per_cell(r)) is not None and v > 0]
    if not rtpc_vals:
        return
    vmin, vmax = min(rtpc_vals), max(rtpc_vals)
    norm = Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.get_cmap("plasma")
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    # Hide all but the first placeholder axes, then span them with a single colorbar
    for ax in cbar_axes:
        ax.set_visible(False)
    # Get the gridspec rows spanned by cbar_axes
    ss_top = cbar_axes[0].get_subplotspec()
    ss_bot = cbar_axes[-1].get_subplotspec()
    from matplotlib.gridspec import GridSpecFromSubplotSpec
    gs_parent = ss_top.get_gridspec()
    row_top = ss_top.rowspan.start
    row_bot = ss_bot.rowspan.stop
    col     = ss_top.colspan.start
    merged_ss = gs_parent[row_top:row_bot, col]
    cax = fig.add_subplot(merged_ss)
    cbar = fig.colorbar(sm, cax=cax)
    cbar.set_label(FIELD_LABELS["runtime_per_cell"], fontsize=fontsize)
    ticks = list(range(math.floor(vmin), math.ceil(vmax) + 1))
    cbar.ax.yaxis.set_major_locator(plt.FixedLocator(ticks))
    cbar.ax.yaxis.set_major_formatter(
        plt.FuncFormatter(lambda x, _: fr"$10^{{{int(x)}}}$")
    )
    cbar.ax.tick_params(labelsize=fontsize)


def plot_quality(all_records, dir_labels, colors, output_path, only_overall=False):
    n_rows = len(QUALITY_ROWS)
    all_recs_merged = [r for dl in dir_labels for r in all_records[dl]]
    recs_success    = [r for r in all_recs_merged if r.get("opti_result") == 0]
    recs_fail       = [r for r in all_recs_merged if r.get("opti_result") != 0]
    fontsize = QUALITY_FONT_SIZE

    if only_overall:
        # 2 data cols (success, failures) + cbar
        n_data = 2
        fig = plt.figure(figsize=(QUALITY_FIG_WIDTH * n_data, QUALITY_FIG_HEIGHT))
        with_cbar = HISTOGRAM_RUNTIME_GRADIENT
        # Outer gridspec: data area + cbar area; inner for data columns with own wspace
        if with_cbar:
            outer_gs = fig.add_gridspec(1, 2, width_ratios=[n_data, 0.05], wspace=0.04)
            inner_gs = outer_gs[0].subgridspec(n_rows, n_data, wspace=0.3, hspace=0.95)
            cbar_gs  = outer_gs[1].subgridspec(n_rows, 1, hspace=0.75)
        else:
            inner_gs = fig.add_gridspec(n_rows, n_data, wspace=0.45, hspace=0.95)
        axes = [[fig.add_subplot(inner_gs[r, c]) for c in range(n_data)] for r in range(n_rows)]

        for row_idx, field in enumerate(QUALITY_ROWS):
            for c_idx, (recs_sub, title) in enumerate([
                (recs_success, "successes"), (recs_fail, "failures")]):
                ax = axes[row_idx][c_idx]
                sub_vals_rtpc = [(_get_field(r, field), _runtime_per_cell(r))
                                 for r in recs_sub
                                 if _get_field(r, field) is not None and _get_field(r, field) > 0]
                sub_vals   = [math.log10(v[0]) for v in sub_vals_rtpc]
                rtpc_sub_q = [v[1] for v in sub_vals_rtpc]
                _oo_cap = FIELD_CAPS.get(field)
                _oo_cap_hi = int(math.log10(_oo_cap)) if _oo_cap else None
                _distortion_histogram(ax, sub_vals, HISTOGRAM_ALL_COLOR,
                                      n_total=len(recs_sub),
                                      rtpc_values=rtpc_sub_q, fontsize=fontsize,
                                      clamp_hi=_oo_cap_hi)
                ax.set_xlabel(FIELD_LABELS.get(field, field), fontsize=fontsize)
                if row_idx == 0:
                    ax.set_title(title, fontsize=fontsize, pad=30)

        if with_cbar:
            cbar_col_axes = [fig.add_subplot(cbar_gs[r, 0]) for r in range(n_rows)]
            _add_histogram_runtime_colorbar_gs(fig, all_recs_merged, cbar_col_axes, fontsize=fontsize)

    else:
        n_data = len(dir_labels) + 2  # datasets + success-all + failures-all
        with_cbar = HISTOGRAM_RUNTIME_GRADIENT
        fig_w = QUALITY_FIG_WIDTH * n_data / max(len(dir_labels), 1)
        fig = plt.figure(figsize=(fig_w, QUALITY_FIG_HEIGHT))
        # Outer gridspec: data area + cbar area; inner for data columns with own wspace
        if with_cbar:
            outer_gs = fig.add_gridspec(1, 2, width_ratios=[n_data, 0.12], wspace=0.04)
            inner_gs = outer_gs[0].subgridspec(n_rows, n_data, wspace=0.55, hspace=0.95)
            cbar_gs  = outer_gs[1].subgridspec(n_rows, 1, hspace=0.75)
        else:
            inner_gs = fig.add_gridspec(n_rows, n_data, wspace=0.55, hspace=0.95)
        axes = [[fig.add_subplot(inner_gs[r, c]) for c in range(n_data)] for r in range(n_rows)]

        for row_idx, field in enumerate(QUALITY_ROWS):
            for col_idx, dir_label in enumerate(dir_labels):
                ax    = axes[row_idx][col_idx]
                color = colors[col_idx % len(colors)]
                recs  = all_records[dir_label]
                values = [math.log10(_get_field(r, field))
                          for r in recs
                          if _get_field(r, field) is not None and _get_field(r, field) > 0]
                _qcap = FIELD_CAPS.get(field)
                _qcap_hi = int(math.log10(_qcap)) if _qcap else None
                _distortion_histogram(ax, values, color, n_total=len(recs), fontsize=fontsize,
                                      clamp_hi=_qcap_hi)
                if row_idx == 0:
                    ax.set_title(dir_label, fontsize=fontsize)
                if col_idx == 0:
                    ax.set_ylabel(FIELD_LABELS.get(field, field), fontsize=fontsize,
                                  rotation=90, labelpad=6)
            # All (success)
            ax_all = axes[row_idx][len(dir_labels)]
            all_vals_rtpc = [(_get_field(r, field), _runtime_per_cell(r))
                             for r in recs_success
                             if _get_field(r, field) is not None and _get_field(r, field) > 0]
            all_vals   = [math.log10(v[0]) for v in all_vals_rtpc]
            rtpc_all_q = [v[1] for v in all_vals_rtpc]
            _qcap = FIELD_CAPS.get(field)
            _qcap_hi = int(math.log10(_qcap)) if _qcap else None
            _distortion_histogram(ax_all, all_vals, HISTOGRAM_ALL_COLOR,
                                  n_total=len(recs_success), rtpc_values=rtpc_all_q,
                                  fontsize=fontsize, clamp_hi=_qcap_hi)
            if row_idx == 0:
                ax_all.set_title("All (success)", fontsize=fontsize)
            # All (failures)
            ax_fail = axes[row_idx][len(dir_labels) + 1]
            fail_vals_rtpc = [(_get_field(r, field), _runtime_per_cell(r))
                              for r in recs_fail
                              if _get_field(r, field) is not None and _get_field(r, field) > 0]
            fail_vals   = [math.log10(v[0]) for v in fail_vals_rtpc]
            rtpc_fail_q = [v[1] for v in fail_vals_rtpc]
            _distortion_histogram(ax_fail, fail_vals, HISTOGRAM_ALL_COLOR,
                                  n_total=len(recs_fail), rtpc_values=rtpc_fail_q,
                                  fontsize=fontsize, clamp_hi=_qcap_hi)
            if row_idx == 0:
                ax_fail.set_title("All (failures)", fontsize=fontsize)

        if with_cbar:
            cbar_col_axes = [fig.add_subplot(cbar_gs[r, 0]) for r in range(n_rows)]
            _add_histogram_runtime_colorbar_gs(fig, all_recs_merged, cbar_col_axes, fontsize=fontsize)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"Quality plot saved to: {output_path}")


# ---------------------------------------------------------------------------
# Quality progress scatter
# ---------------------------------------------------------------------------

def plot_quality_progress_scatter(all_records, dir_labels, output_path):
    PAIRS = [
        ("init_min_worst_quality",       "final_min_worst_quality"),
        ("init_min_relative_volume",     "final_min_relative_volume"),
        ("init_min_relative_edge_length","final_min_relative_edge_length"),
    ]
    fig = plt.figure(figsize=(QUALITY_PROGRESS_FIG_WIDTH, QUALITY_PROGRESS_FIG_HEIGHT))
    gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 1, 0.05], wspace=0.35)
    axes = [fig.add_subplot(gs[i]) for i in range(3)]
    cbar_ax = fig.add_subplot(gs[3])

    all_recs = [r for dl in dir_labels for r in all_records[dl]]
    rc = _runtime_colormap(all_recs, use_per_cell=True)
    if not rc:
        print("[warning] No valid runtime values for quality progress scatter.")
        return
    norm, cmap, vmin, vmax = rc

    for ax, (x_field, y_field) in zip(axes, PAIRS):
        for opti_result, (marker, _) in OPTI_MARKERS.items():
            xs, ys, cs = [], [], []
            for r in all_recs:
                if r.get("opti_result") != opti_result:
                    continue
                x  = _get_field(r, x_field)
                y  = _get_field(r, y_field)
                rtpc = _runtime_per_cell(r)
                if x is None or y is None or rtpc is None or rtpc <= 0:
                    continue
                xs.append(x)
                ys.append(y)
                cs.append(math.log10(rtpc))
            if xs:
                order = sorted(range(len(xs)), key=lambda i: cs[i])
                _scatter_c(ax, [xs[i] for i in order], [ys[i] for i in order],
                           [cmap(norm(cs[i])) for i in order], marker, s=MARKER_SIZE)

        _style_ax(ax, FIELD_LABELS.get(x_field, x_field), FIELD_LABELS.get(y_field, y_field),
                  xlog=x_field in LOG_SCALE_FIELDS, ylog=y_field in LOG_SCALE_FIELDS,
                  fontsize=QUALITY_PROGRESS_FONT_SIZE)

    _add_colorbar(fig, cbar_ax, cmap, norm, vmin, vmax, FIELD_LABELS["runtime_per_cell"],
                  fontsize=QUALITY_PROGRESS_FONT_SIZE)
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"Quality progress scatter saved to: {output_path}")


# ---------------------------------------------------------------------------
# Boundary distortion
# ---------------------------------------------------------------------------

def plot_boundary_distortion(all_records, dir_labels, colors, output_path):
    """Single scatter plot: max_boundary_SDE (x) vs runtime/cell (y)."""
    fig, ax = plt.subplots(1, 1, figsize=(BOUNDARY_DISTORTION_FIG_WIDTH,
                                          BOUNDARY_DISTORTION_FIG_HEIGHT))
    x_field = "max_boundary_SDE"
    for dir_idx, dir_label in enumerate(dir_labels):
        color = colors[dir_idx % len(colors)]
        for opti_result, (marker, _) in OPTI_MARKERS.items():
            xs, ys = [], []
            for r in all_records[dir_label]:
                if r.get("opti_result") != opti_result:
                    continue
                x    = _get_field(r, x_field)
                y    = _runtime_per_cell(r)
                if x is None or y is None:
                    continue
                xs.append(x)
                ys.append(y)
            if xs:
                _scatter(ax, xs, ys, color, marker, s=MARKER_SIZE)

    _style_ax(ax, FIELD_LABELS.get(x_field, x_field), FIELD_LABELS["runtime_per_cell"],
              fontsize=BOUNDARY_DIST_FONT_SIZE)

    fig.legend(handles=_legend_handles(dir_labels, colors) + _marker_handles(),
               title="Datasets / Result",
               loc="upper left", bbox_to_anchor=(0.8, 0.9),
               borderaxespad=0, fontsize=BOUNDARY_DIST_FONT_SIZE, title_fontsize=BOUNDARY_DIST_FONT_SIZE)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"Boundary distortion plot saved to: {output_path}")


# ---------------------------------------------------------------------------
# Quality-size scatter
# ---------------------------------------------------------------------------

def plot_quality_size(all_records, dir_labels, colors, output_path):
    import matplotlib.colors as mcolors

    ROW1_Y_FIELDS = [
        "max_boundary_SDE",
        "init_min_worst_quality",
        "init_min_relative_volume",
        "init_min_relative_edge_length",
    ]
    # Extra composite plot appended at end of row 1
    ROW2_PAIRS = [
        ("max_boundary_SDE",            "init_min_worst_quality"),
        ("max_boundary_SDE",            "init_min_relative_edge_length"),
        ("max_boundary_SDE",            "init_min_relative_volume"),
        ("init_min_relative_volume",    "init_min_relative_edge_length"),
        ("init_min_relative_volume",    "init_min_worst_quality"),
    ]

    n_row1 = len(ROW1_Y_FIELDS) + 1   # +1 for composite
    n_row2 = len(ROW2_PAIRS)
    n_cols  = max(n_row1, n_row2)

    all_recs = [r for dl in dir_labels for r in all_records[dl]]
    rc = _runtime_colormap(all_recs, use_per_cell=True)
    if not rc:
        print("[warning] No valid runtime values for quality-size scatter.")
        return
    norm, cmap, vmin, vmax = rc

    fig = plt.figure(figsize=(QUALITY_SIZE_FIG_WIDTH, QUALITY_SIZE_FIG_HEIGHT * 2))
    gs = fig.add_gridspec(2, n_cols + 1,
                          width_ratios=[1] * n_cols + [0.05],
                          wspace=0.35, hspace=0.4)
    axes_row1 = [fig.add_subplot(gs[0, i]) for i in range(n_row1)]
    axes_row2 = [fig.add_subplot(gs[1, i]) for i in range(n_row2)]
    # Hide unused row2 slots
    for i in range(n_row2, n_cols):
        fig.add_subplot(gs[1, i]).set_visible(False)
    cbar_ax = fig.add_subplot(gs[0, n_cols])

    # --- Row 1: n_init_cells vs quality field, colored by runtime/cell ---
    for ax, y_field in zip(axes_row1[:len(ROW1_Y_FIELDS)], ROW1_Y_FIELDS):
        for opti_result, (marker, _) in OPTI_MARKERS.items():
            xs, ys, cs = [], [], []
            for r in all_recs:
                if r.get("opti_result") != opti_result:
                    continue
                x    = _get_field(r, "n_init_cells")
                y    = _get_field(r, y_field)
                rtpc = _runtime_per_cell(r)
                if x is None or y is None or rtpc is None or rtpc <= 0:
                    continue
                xs.append(x)
                ys.append(y)
                cs.append(math.log10(rtpc))
            if xs:
                order = sorted(range(len(xs)), key=lambda i: cs[i])
                _scatter_c(ax, [xs[i] for i in order], [ys[i] for i in order],
                           [cmap(norm(cs[i])) for i in order], marker, s=MARKER_SIZE)
        _style_ax(ax, FIELD_LABELS.get("n_init_cells"), FIELD_LABELS.get(y_field, y_field),
                  xlog=True, ylog=y_field in LOG_SCALE_FIELDS, fontsize=QUALITY_SIZE_FONT_SIZE)

    # --- Row 1 extra: composite metric vs n_init_cells ---
    ax_comp = axes_row1[-1]
    for opti_result, (marker, _) in OPTI_MARKERS.items():
        xs, ys, cs = [], [], []
        for r in all_recs:
            if r.get("opti_result") != opti_result:
                continue
            sde  = _get_field(r, "max_boundary_SDE")
            wq   = _get_field(r, "init_min_worst_quality")
            el   = _get_field(r, "init_min_relative_edge_length")
            vol  = _get_field(r, "init_min_relative_volume")
            n    = _get_field(r, "n_init_cells")
            rtpc = _runtime_per_cell(r)
            if any(v is None or v <= 0 for v in [sde, wq, el, vol, n]) or rtpc is None or rtpc <= 0:
                continue
            x_val = (math.log10(n)
                     * (1 + math.log10(sde))
                     * (-math.log10(wq))
                     * (abs(math.log10(el)) + 1)
                     * (abs(math.log10(vol)) + 1))
            if x_val <= 0:
                continue
            xs.append(x_val)
            ys.append(n)
            cs.append(math.log10(rtpc))
        if xs:
            order = sorted(range(len(xs)), key=lambda i: cs[i])
            _scatter_c(ax_comp, [xs[i] for i in order], [ys[i] for i in order],
                       [cmap(norm(cs[i])) for i in order], marker, s=MARKER_SIZE)
    _style_ax(ax_comp,
              "log(n)·(1+log SDE)·(-log WQ)·(|log EL|+1)·(|log Vol|+1)",
              FIELD_LABELS.get("n_init_cells"), xlog=True, ylog=True,
              fontsize=QUALITY_SIZE_FONT_SIZE)
    ax_comp.set_xlabel("Composite metric", fontsize=QUALITY_SIZE_FONT_SIZE)

    _add_colorbar(fig, cbar_ax, cmap, norm, vmin, vmax, FIELD_LABELS["runtime_per_cell"],
                  fontsize=QUALITY_SIZE_FONT_SIZE)

    # --- Row 2: quality pairs, colored by runtime/cell ---
    for ax, (x_field, y_field) in zip(axes_row2, ROW2_PAIRS):
        for opti_result, (marker, _) in OPTI_MARKERS.items():
            xs, ys, cs = [], [], []
            for r in all_recs:
                if r.get("opti_result") != opti_result:
                    continue
                x    = _get_field(r, x_field)
                y    = _get_field(r, y_field)
                rtpc = _runtime_per_cell(r)
                if x is None or y is None or rtpc is None or rtpc <= 0:
                    continue
                xs.append(x)
                ys.append(y)
                cs.append(math.log10(rtpc))
            if xs:
                order = sorted(range(len(xs)), key=lambda i: cs[i])
                _scatter_c(ax, [xs[i] for i in order], [ys[i] for i in order],
                           [cmap(norm(cs[i])) for i in order], marker, s=MARKER_SIZE)
        _style_ax(ax, FIELD_LABELS.get(x_field, x_field), FIELD_LABELS.get(y_field, y_field),
                  xlog=x_field in LOG_SCALE_FIELDS, ylog=y_field in LOG_SCALE_FIELDS,
                  fontsize=QUALITY_SIZE_FONT_SIZE)

    # Row 2 shared colorbar reuses row 1's — no second cbar needed
    # Add dataset+marker legend outside
    legend_ax = fig.add_subplot(gs[1, n_cols])
    legend_ax.axis("off")
    legend_ax.legend(handles=_legend_handles(dir_labels, colors) + _marker_handles(),
                     title="Datasets / Result", loc="center left", fontsize=QUALITY_SIZE_FONT_SIZE, title_fontsize=QUALITY_SIZE_FONT_SIZE)

    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"Quality-size scatter saved to: {output_path}")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Remeshing performance plot
# ---------------------------------------------------------------------------

REMESHING_PERF_FIG_WIDTH  = 12   # inches
REMESHING_PERF_FIG_HEIGHT = 5    # inches


def plot_remeshing_performance(all_records, dir_labels, output_path):
    """Two scatter plots: min/avg quality before vs after each remeshing step.

    Uses min_quality_before_remeshing_history and min_quality_after_remeshing_history,
    each being an array of values (one per remeshing step).
    """
    avg_before, avg_after = [], []
    min_before, min_after = [], []
    n_cases_with_data = 0

    for dir_label in dir_labels:
        for r in all_records[dir_label]:
            min_bh = r.get("min_quality_before_remeshing_history")
            min_ah = r.get("min_quality_after_remeshing_history")
            avg_bh = r.get("avg_quality_before_remeshing_history")
            avg_ah = r.get("avg_quality_after_remeshing_history")
            if not min_bh or not min_ah or not avg_bh or not avg_ah:
                continue
            n_cases_with_data += 1
            for b, a in zip(min_bh, min_ah):
                if b is not None and a is not None and b > 0 and a > 0:
                    min_before.append(b); min_after.append(a)
            for b, a in zip(avg_bh, avg_ah):
                if b is not None and a is not None and b > 0 and a > 0:
                    avg_before.append(b); avg_after.append(a)

    print("\n=== Remeshing performance ===")
    print(f"  Cases with remeshing history data: {n_cases_with_data}")
    print(f"  Total before/after pairs: {len(min_before)}")

    fig = plt.figure(figsize=(REMESHING_PERF_FIG_WIDTH, REMESHING_PERF_FIG_HEIGHT))
    gs_perf = fig.add_gridspec(1, 4, width_ratios=[1, 1, -0.7, 0.06], wspace=1.1)
    ax_avg   = fig.add_subplot(gs_perf[0])
    ax_max   = fig.add_subplot(gs_perf[1])
    _sp      = fig.add_subplot(gs_perf[2]); _sp.set_visible(False)
    cbar_ax  = fig.add_subplot(gs_perf[3])
    color = "mediumseagreen"

    print("  Improvement stats (log10(after/before)):")
    # Pre-compute global max bin count for colorbar label
    h_max_global = 1
    for _, xs, ys, _ in [
        (None, avg_before, avg_after, None),
        (None, min_before, min_after, None),
    ]:
        if not xs:
            continue
        log_xs = np.array([math.log10(v) for v in xs])
        log_ys = np.array([math.log10(v) for v in ys])
        lo = min(log_xs.min(), log_ys.min())
        hi = max(log_xs.max(), log_ys.max())
        h, _, _ = np.histogram2d(log_xs, log_ys, bins=80, range=[[lo, hi], [lo, hi]])
        h_max_global = max(h_max_global, int(h.max()))
    for ax, xs, ys, label in [
        (ax_avg, avg_before, avg_after, "average quality"),
        (ax_max, min_before, min_after, "minimum quality"),
    ]:
        fontsize = REMESHING_PERF_FONT_SIZE
        if not xs:
            ax.set_visible(False)
            continue
        print(f"    {label}: {len(xs)} data points plotted")

        # Log-space 2D histogram heatmap
        log_xs = np.array([math.log10(v) for v in xs])
        log_ys = np.array([math.log10(v) for v in ys])
        lo = min(log_xs.min(), log_ys.min())
        hi = max(log_xs.max(), log_ys.max())
        n_bins = 80
        h, xedges, yedges = np.histogram2d(log_xs, log_ys, bins=n_bins,
                                           range=[[lo, hi], [lo, hi]])
        # Normalise against global max so both plots share the same scale
        h_log = np.where(h > 0, np.log10(h + 1), 0)
        h_log_max = math.log10(h_max_global + 1)
        h_norm = np.where(h > 0, 0. + 1.0 * h_log / h_log_max, 0)
        h_masked = np.ma.masked_where(h == 0, h_norm)
        ax.pcolormesh(xedges, yedges, h_masked.T, cmap="YlOrRd",
                      shading="auto", zorder=1, vmin=0, vmax=1)

        # Identity line (no change)
        ax.plot([lo, hi], [lo, hi], color="black", linewidth=0.8,
                linestyle="--", alpha=0.6, zorder=2)
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        # Format ticks as 10^n
        ticks = list(range(math.floor(lo), math.ceil(hi) + 1))
        ax.xaxis.set_major_locator(plt.FixedLocator(ticks))
        ax.yaxis.set_major_locator(plt.FixedLocator(ticks))
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: fr"$10^{{{int(x)}}}$"))
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: fr"$10^{{{int(x)}}}$"))
        ax.xaxis.set_minor_locator(plt.NullLocator())
        ax.yaxis.set_minor_locator(plt.NullLocator())
        ax.set_xlabel(f"{label} before remeshing", fontsize=fontsize)
        ax.set_ylabel(f"{label} after remeshing",  fontsize=fontsize)
        ax.grid(True, which="major", linestyle="--", alpha=0.4)
        ax.tick_params(axis="both", labelsize=fontsize)

        # Compute improvement stats: log10(y) - log10(x)
        deltas = np.array([math.log10(y) - math.log10(x)
                           for x, y in zip(xs, ys) if x > 0 and y > 0])
        if len(deltas) == 0:
            continue
        d_min  = deltas.min()
        d_max  = deltas.max()
        d_avg  = deltas.mean()
        d_std  = deltas.std()

        # Print to terminal
        print(f"    {label}:")
        print(f"      min={d_min:+.3f}  avg={d_avg:+.3f}  max={d_max:+.3f}  std={d_std:.3f}")

        # Legend box in lower right corner
        stats_text = "\n".join([r"$min_{\text{log}}$"+f" {d_min:+.3f}", r"$max_{\text{log}}$"+f" {d_max:+.3f}", r"$\mu_{\text{log}}$"+f" {d_avg:+.3f}",  r"$\sigma_{\text{log}}$" + f" {d_std:+.3f}"])
        ax.text(0.98, 0.02, stats_text,
                transform=ax.transAxes,
                fontsize=fontsize, family="monospace",
                verticalalignment="bottom", horizontalalignment="right",
                bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                          edgecolor="gray", alpha=0.8))

    # Colorbar spanning full height
    from matplotlib.colors import Normalize as _Norm
    from matplotlib.cm import ScalarMappable as _SM
    from matplotlib.colors import LinearSegmentedColormap
    _base = plt.get_cmap("YlOrRd")
    _cmap_trunc = LinearSegmentedColormap.from_list(
        "YlOrRd_trunc", _base(np.linspace(0.15, 1.0, 256)))
    _sm = _SM(cmap=_cmap_trunc, norm=_Norm(vmin=0, vmax=1))
    _sm.set_array([])
    cbar = fig.colorbar(_sm, cax=cbar_ax)
    cbar.set_label("count (log scale)", fontsize=REMESHING_PERF_FONT_SIZE)
    cbar.ax.yaxis.set_major_locator(plt.FixedLocator([0, 1.0]))
    cbar.ax.yaxis.set_major_formatter(plt.FuncFormatter(
        lambda v, _: "1" if v <= 0 else str(h_max_global)
    ))
    cbar.ax.tick_params(labelsize=REMESHING_PERF_FONT_SIZE)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"Remeshing performance plot saved to: {output_path}")


def _print_cell_count_ranges(all_records, dir_labels):
    """Print the min/max n_init_cells for each dataset."""
    print("\n=== n_init_cells range per dataset ===")
    for dir_label in dir_labels:
        vals = [r["n_init_cells"] for r in all_records[dir_label]
                if r.get("n_init_cells") is not None]
        if not vals:
            print(f"  {dir_label}: no data")
        else:
            print(f"  {dir_label}: {min(vals)} \u2013 {max(vals)}")


def _print_success_rate(all_records, dir_labels, runtime_cap_s=6 * 3600):
    """Print success rate (opti_result==0 and runtime_s <= cap) per dataset."""
    print(f"\n=== Success rate (runtime <= {runtime_cap_s/3600:.0f}h) ===")
    for dir_label in dir_labels:
        records = all_records[dir_label]
        total = len(records)
        if total == 0:
            print(f"  {dir_label}: no data")
            continue
        n_success = sum(1 for r in records
                        if r.get("opti_result") == 0
                        and r.get("runtime_s", float("inf")) <= runtime_cap_s)
        pct = 100 * n_success / total
        print(f"  {dir_label}: {n_success}/{total}  ({pct:.1f}%)")


def _print_low_runtime_shrunk_cases(all_records, dir_labels):
    """Print the 10 cases with lowest runtime_s where volume shrank significantly."""
    FINAL_VOL_CAP  = 1e-9
    INIT_VOL_FLOOR = 1e-6
    candidates = []
    for dir_label in dir_labels:
        for r in all_records[dir_label]:
            fv = r.get("final_min_relative_volume")
            iv = r.get("init_min_relative_volume")
            rt = r.get("runtime_s")
            mn  = r.get("mesh_name", "<no mesh_name>")
            fn  = r.get("_filename", "<unknown>")
            if fv is None or iv is None or rt is None:
                continue
            if fv < FINAL_VOL_CAP and iv > INIT_VOL_FLOOR:
                candidates.append((rt, dir_label, mn, fn, iv, fv))
    candidates.sort(key=lambda x: x[0])
    print(f"\n=== Top 10 lowest runtime_s: final_vol < {FINAL_VOL_CAP:.0e} & init_vol > {INIT_VOL_FLOOR:.0e} ===")
    if not candidates:
        print("  (none found)")
        return
    for rt, dir_label, mesh_name, filename, iv, fv in candidates[:10]:
        print(f"  runtime={rt:.2f}s  init_vol={iv:.2e}  final_vol={fv:.2e}  "
              f"file={filename}  dataset={dir_label}  mesh={mesh_name}")


QUALITY_RESULT_FIELDS = [
    "max_boundary_SDE",
    "avg_boundary_SDE",
    #"input_min_worst_quality",
    #"input_min_relative_volume",
    #"input_min_relative_edge_length",
    "n_init_cells",
]
QUALITY_RESULT_FIG_WIDTH  = 10   # inches
QUALITY_RESULT_FIG_HEIGHT  = 3    # inches
QUALITY_RESULT_FONT_SIZE  = 16   # font size for quality-result plot


def plot_quality_result(all_records, dir_labels, colors, output_path):
    """One histogram per field, bars split by opti_result==0 vs !=0."""
    all_recs = [r for dl in dir_labels for r in all_records[dl]]
    n_fields = len(QUALITY_RESULT_FIELDS)
    fontsize  = QUALITY_RESULT_FONT_SIZE

    fig, axes = plt.subplots(n_fields, 1,
                             figsize=(QUALITY_RESULT_FIG_WIDTH, QUALITY_RESULT_FIG_HEIGHT * n_fields),
                             squeeze=False)

    # Two colors: success (opti_result==0) and non-success
    color_success = "steelblue"
    color_fail    = "tomato"

    for col_idx, field in enumerate(QUALITY_RESULT_FIELDS):
        ax = axes[col_idx][0]

        vals_success = [math.log10(v) for r in all_recs
                        if r.get("opti_result") == 0
                        and (v := _get_field(r, field)) is not None and v > 0]
        vals_fail    = [math.log10(v) for r in all_recs
                        if r.get("opti_result") != 0
                        and (v := _get_field(r, field)) is not None and v > 0]

        if not vals_success and not vals_fail:
            ax.set_visible(False)
            continue

        all_vals = vals_success + vals_fail
        lo = math.floor(min(all_vals))
        hi = math.ceil(max(all_vals))
        cap = FIELD_CAPS.get(field)
        cap_hi = int(math.log10(cap)) if cap else None
        if cap_hi is not None:
            hi = max(hi, cap_hi)
        N  = len(all_recs)

        counts_s, bin_edges = np.histogram(vals_success, bins=HISTOGRAM_N_BINS, range=(lo, hi))
        counts_f, _         = np.histogram(vals_fail,    bins=HISTOGRAM_N_BINS, range=(lo, hi))

        N_s = len(vals_success) if vals_success else 1
        N_f = len(vals_fail)    if vals_fail    else 1
        ratios_s = np.where(counts_s > 0, counts_s / N_s, np.nan)
        ratios_f = np.where(counts_f > 0, counts_f / N_f, np.nan)

        bin_step = bin_edges[1] - bin_edges[0]
        half     = bin_step * 0.85 / 2

        ax.bar(bin_edges[:-1],        ratios_s, width=half, align="edge",
               color=color_success, edgecolor="black", linewidth=0.5, label="success")
        ax.bar(bin_edges[:-1] + half, ratios_f, width=half, align="edge",
               color=color_fail,    edgecolor="black", linewidth=0.5, alpha=0.7, label="stuck/timeout")

        ax.set_yscale("log")
        ax.yaxis.set_minor_locator(plt.NullLocator())

        ymin_val = min((v for v in [*ratios_s, *ratios_f] if v > 0), default=1e-4)
        y_lo = math.floor(math.log10(ymin_val))
        y_ticks = list(range(y_lo, 1))
        if not y_ticks:
            y_ticks = [y_lo, 0]
        if 0 not in y_ticks:
            y_ticks.append(0)
        ax.set_ylim(bottom=10**y_lo * 0.5, top=2.0)
        ax.yaxis.set_major_locator(plt.FixedLocator([10**t for t in y_ticks]))
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f"{y * 100:.4g}%"))

        # At most 5 evenly-spaced tick labels, cap_hi always included
        log_range = hi - lo
        for x_step in range(1, log_range + 2):
            tick_labels = list(range(lo, hi + 1, x_step))
            if len(tick_labels) <= 5:
                break
        if len(tick_labels) < 2:
            tick_labels = [lo, lo + x_step]
        if cap_hi is not None and cap_hi not in tick_labels:
            tick_labels = [l for l in tick_labels if abs(l - cap_hi) > x_step / 2]
            tick_labels.append(cap_hi)
            tick_labels.sort()
        def _qr_xfmt(x, _, _cap_hi=cap_hi):
            xi = round(x)
            if _cap_hi is not None and xi == _cap_hi:
                return fr"$\geq10^{{{xi}}}$"
            return fr"$10^{{{xi}}}$"
        ax.xaxis.set_major_locator(plt.FixedLocator(tick_labels))
        ax.xaxis.set_major_formatter(plt.FuncFormatter(_qr_xfmt))
        ax.tick_params(axis="both", labelsize=fontsize)
        ax.grid(True, which="major", axis="y", linestyle="--", alpha=0.4)
        if col_idx == 0:
            ax.legend(fontsize=fontsize - 2)
        ax.set_xlabel(FIELD_LABELS.get(field, field), fontsize=fontsize)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"Quality-result plot saved to: {output_path}")


def main():
    t_script_start = time.perf_counter()
    parser = argparse.ArgumentParser(description="Plot remeshing scatter data.")
    parser.add_argument("--input-dirs-file", required=True, metavar="FILE",
                        help='Path to a JSON file: [{"dir": "...", "label": "..."}]')
    parser.add_argument("--show-plots", action="store_true")
    parser.add_argument("--trim", metavar="N", type=int, default=None)
    parser.add_argument("--only-plot-overall-quality", action="store_true",
                        help="Use vertical/2x2 layout showing only the merged-all histograms.")
    parser.add_argument("--scatter-plot-path",                metavar="PATH", default=None)
    parser.add_argument("--performance-plot-path",            metavar="PATH", default=None)
    parser.add_argument("--performance-heatmap-plot-path",    metavar="PATH", default=None)
    parser.add_argument("--quality-size-plot-path",           metavar="PATH", default=None)
    parser.add_argument("--boundary-distortion-plot-path",    metavar="PATH", default=None)
    parser.add_argument("--quality-progress-scatter-plot-path", metavar="PATH", default=None)
    parser.add_argument("--quality-plot-path",                metavar="PATH", default=None)
    parser.add_argument("--distortion-plot-path",             metavar="PATH", default=None)
    parser.add_argument("--runtime-ternary-plot-path",        metavar="PATH", default=None)
    parser.add_argument("--remeshing-ternary-plot-path",      metavar="PATH", default=None)
    parser.add_argument("--remeshing-performance-plot-path",  metavar="PATH", default=None)
    parser.add_argument("--quality-result-plot-path",         metavar="PATH", default=None)
    args = parser.parse_args()

    plot_args = [args.scatter_plot_path, args.remeshing_ternary_plot_path,
                 args.performance_plot_path, args.performance_heatmap_plot_path,
                 args.distortion_plot_path,
                 args.quality_plot_path, args.runtime_ternary_plot_path,
                 args.quality_progress_scatter_plot_path,
                 args.boundary_distortion_plot_path, args.quality_size_plot_path,
                 args.remeshing_performance_plot_path,
                 args.quality_result_plot_path]
    if not any(plot_args):
        parser.error("At least one plot output argument must be set.")

    if USE_INPUT_QUALITY_METRICS:
        for init_f, input_f in INIT_TO_INPUT_FIELDS.items():
            if input_f in FIELD_LABELS:
                FIELD_LABELS[init_f] = FIELD_LABELS[input_f]

    input_file = Path(args.input_dirs_file)
    if not input_file.exists():
        print(f"Error: input file not found: {input_file}", file=sys.stderr)
        sys.exit(1)
    try:
        entries = json.loads(input_file.read_text())
    except json.JSONDecodeError as e:
        print(f"Error: could not parse JSON input file: {e}", file=sys.stderr)
        sys.exit(1)
    if not entries:
        print("Error: no entries found in input file.", file=sys.stderr)
        sys.exit(1)

    directories: list[tuple[Path, str]] = [
        (Path(entry["dir"]), entry["label"]) for entry in entries
    ]

    # Collect data
    all_records: dict[str, list[dict]] = {}
    unknown_opti: list[tuple[str, str, object]] = []

    x_fields = [field for row in PLOT_GRID for field in row]
    extra_fields = list(TERNARY_FIELDS) + list(RUNTIME_TERNARY_FIELDS) + [
        "n_init_cells",     "n_init_vertices",
        "input_min_worst_quality", "input_min_relative_volume",
        "input_min_relative_edge_length",
        "n_final_cells",
        "max_final_SDE",        "avg_final_SDE",
        "max_post_opti_SDE",    "avg_post_opti_SDE",
        "max_final_conf_dist",  "avg_final_conf_dist",
        "max_post_opti_conf_dist", "avg_post_opti_conf_dist",
        "final_min_worst_quality", "final_min_relative_volume",
        "final_min_relative_edge_length",
        "max_boundary_SDE",     "avg_boundary_SDE",
        "min_quality_before_remeshing_history", "min_quality_after_remeshing_history",
        "avg_quality_before_remeshing_history", "avg_quality_after_remeshing_history",
    ]
    all_fields = list(dict.fromkeys(x_fields + extra_fields))  # deduplicated

    for directory, dir_label in directories:
        all_records[dir_label] = []
        print(f"\nLooking into directory: {directory}  [{dir_label}]")
        t_start = time.perf_counter()
        json_entries = load_json_files(directory, trim=args.trim)
        if not json_entries:
            print("  Loaded 0 files in 0.00s")
            continue
        loaded = 0
        for json_file, data in json_entries:
            runtime_s   = data.get("runtime_s")
            opti_result = data.get("opti_result")
            mesh_name   = data.get("mesh_name", "<no mesh_name>")
            if runtime_s is None or opti_result is None:
                print(f"\n  [warning] Missing runtime_s or opti_result in {json_file}, skipping.",
                      file=sys.stderr)
                continue
            if opti_result not in OPTI_MARKERS:
                unknown_opti.append((dir_label, mesh_name, opti_result))
                continue
            record = {"runtime_s": runtime_s, "opti_result": opti_result,
                      "_filename": json_file.stem}
            for field in all_fields:
                record[field] = data.get(field)
            if USE_INPUT_QUALITY_METRICS:
                for init_f, input_f in INIT_TO_INPUT_FIELDS.items():
                    if data.get(input_f) is not None:
                        record[init_f] = data.get(input_f)
            all_records[dir_label].append(record)
            loaded += 1
            print(f"  Loading... {loaded}", end="\r", flush=True)
        elapsed = time.perf_counter() - t_start
        print(f"  Loaded {loaded} files in {elapsed:.2f}s")

    # Summary
    print("\n=== Data summary per directory ===")
    dir_labels = list(all_records.keys())
    for dir_label, records in all_records.items():
        total = len(records)
        print(f"\n  {dir_label}  ({total} cases)")
        if total == 0:
            print("    (no valid entries)")
            continue
        for opti_result, (marker, _) in OPTI_MARKERS.items():
            count = sum(1 for r in records if r["opti_result"] == opti_result)
            filled = round(20 * count / total) if total else 0
            bar = "█" * filled + "░" * (20 - filled)
            pct = f"{100 * count / total:5.1f}%" if total else "    -"
            print(f"    opti_result={opti_result} [{marker}]  {count:>4}  ({pct})  {bar}")
        unknown_in_dir = sum(1 for d, _, _ in unknown_opti if d == dir_label)
        if unknown_in_dir:
            print(f"    opti_result=?  (unexpected)  {unknown_in_dir:>4}")

    if unknown_opti:
        print("\n=== JSON entries with unexpected opti_result ===")
        for dir_label, mesh_name, opti_result in unknown_opti:
            print(f"  directory  : {dir_label}")
            print(f"  mesh_name  : {mesh_name}")
            print(f"  opti_result: {opti_result}")
            print()
    else:
        print("\nAll opti_result values are among {0, 1, 3}.")

    _print_field_caps(all_records, dir_labels)
    _print_cell_count_ranges(all_records, dir_labels)
    _print_success_rate(all_records, dir_labels)
    _print_low_runtime_shrunk_cases(all_records, dir_labels)

    t_plots_start = time.perf_counter()

    # Render
    colors = plt.cm.tab10.colors

    colors = [
        plt.cm.tab10.colors[DATASET_COLOR_MAP.get(label, i) % 10]
        for i, label in enumerate(dir_labels)
    ]
    print()

    if args.scatter_plot_path:
        p = Path(args.scatter_plot_path); p.parent.mkdir(parents=True, exist_ok=True)
        plot_scatter_grid(all_records, dir_labels, colors, p)

    if args.performance_plot_path:
        p = Path(args.performance_plot_path); p.parent.mkdir(parents=True, exist_ok=True)
        plot_performance(all_records, dir_labels, colors, p)

    if args.performance_heatmap_plot_path:
        p = Path(args.performance_heatmap_plot_path); p.parent.mkdir(parents=True, exist_ok=True)
        plot_performance_heatmap(all_records, dir_labels, colors, p)

    if args.quality_size_plot_path:
        p = Path(args.quality_size_plot_path); p.parent.mkdir(parents=True, exist_ok=True)
        plot_quality_size(all_records, dir_labels, colors, p)

    if args.boundary_distortion_plot_path:
        p = Path(args.boundary_distortion_plot_path); p.parent.mkdir(parents=True, exist_ok=True)
        plot_boundary_distortion(all_records, dir_labels, colors, p)

    if args.quality_progress_scatter_plot_path:
        p = Path(args.quality_progress_scatter_plot_path); p.parent.mkdir(parents=True, exist_ok=True)
        plot_quality_progress_scatter(all_records, dir_labels, p)

    if args.quality_plot_path:
        p = Path(args.quality_plot_path); p.parent.mkdir(parents=True, exist_ok=True)
        plot_quality(all_records, dir_labels, colors, p,
                     only_overall=args.only_plot_overall_quality)

    if args.distortion_plot_path:
        p = Path(args.distortion_plot_path); p.parent.mkdir(parents=True, exist_ok=True)
        plot_distortion(all_records, dir_labels, colors, p,
                        show_post_opti=PLOT_POST_OPTI_HISTOGRAMS,
                        only_overall=args.only_plot_overall_quality)

    if args.runtime_ternary_plot_path:
        p = Path(args.runtime_ternary_plot_path); p.parent.mkdir(parents=True, exist_ok=True)
        plot_ternary(all_records, dir_labels, colors, p,
                     fields=RUNTIME_TERNARY_FIELDS, labels=RUNTIME_TERNARY_LABELS,
                     zero_label="cases without runtime data")

    if args.remeshing_ternary_plot_path:
        p = Path(args.remeshing_ternary_plot_path); p.parent.mkdir(parents=True, exist_ok=True)
        plot_remeshing_ternary(all_records, dir_labels, p)

    if args.remeshing_performance_plot_path:
        p = Path(args.remeshing_performance_plot_path); p.parent.mkdir(parents=True, exist_ok=True)
        plot_remeshing_performance(all_records, dir_labels, p)

    if args.quality_result_plot_path:
        p = Path(args.quality_result_plot_path); p.parent.mkdir(parents=True, exist_ok=True)
        plot_quality_result(all_records, dir_labels, colors, p)

    t_end = time.perf_counter()
    print(f"\nPlot generation time : {t_end - t_plots_start:.2f}s")
    print(f"Total script time    : {t_end - t_script_start:.2f}s")

    if args.show_plots:
        plt.show()


if __name__ == "__main__":
    main()