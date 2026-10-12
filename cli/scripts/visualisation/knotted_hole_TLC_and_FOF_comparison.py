#!/usr/bin/env python3
"""
Summarize algorithm results for knotted_hole_i cases.

Rows (methods): VE, TLC, FOF, TLC with VE output, FOF with VE output
Columns: cases 1..10

VE row shows timing (seconds); all others show a green (success) or red (failure) square.
"""

import argparse
import json
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np


# ── Parsers ──────────────────────────────────────────────────────────────────

def parse_stats_txt(path: Path) -> dict:
    """
    Parse a _stats.txt file.
    Expected format (8 lines):
        name
        n_vertices
        n_edges
        n_faces
        n_cells
        n_degenerate
        n_flipped
        0
    Returns dict with keys: n_degenerate, n_flipped, success
    """
    lines = path.read_text().strip().splitlines()
    if len(lines) < 8:
        raise ValueError(f"Unexpected format in {path}: only {len(lines)} lines")
    n_degenerate = int(lines[5])
    n_flipped    = int(lines[6])
    return {
        "n_degenerate": n_degenerate,
        "n_flipped":    n_flipped,
        "success":      int(n_degenerate == 0 and n_flipped == 0),
    }


def parse_ve_json(path: Path) -> dict:
    """
    Parse a boundary_deformation_map_output.json file.
    Returns dict with keys: success, timing (seconds).
    """
    data = json.loads(path.read_text())
    opti_result = data.get("opti_result", -1)
    timing      = data.get("runtime_s", None)
    return {
        "opti_result": opti_result,
        "success":     int(opti_result == 0),
        "timing":      timing,
    }


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    # Set non-interactive backend early if we will not be showing the plot
    import sys
    if "--show-plots" not in sys.argv:
        matplotlib.use("Agg")

    parser = argparse.ArgumentParser(description="Summarize knotted_hole results.")
    parser.add_argument("--input-directory", required=True, type=Path,
                        help="Directory containing the result files.")
    parser.add_argument("--plot-output-path", type=Path, default=None,
                        help="Path to save the summary image as a .png file. "
                             "Defaults to results_summary.png inside --input-directory.")
    parser.add_argument("--show-plots", action="store_true",
                        help="Display the plot interactively instead of (or in addition to) saving.")
    args = parser.parse_args()

    base: Path = args.input_directory
    if not base.is_dir():
        raise SystemExit(f"ERROR: {base} is not a directory.")

    n_cases = 10
    cases   = range(1, n_cases + 1)

    # Storage: row → list of cell values (str)
    rows = {
        "REVE":              [],
        "TLC":             [],
        "FOF":             [],
        "TLC\'": [],
        "FOF\'": [],
    }

    all_ok = True

    for i in cases:
        case = f"knotted_hole_{i}"

        # ── VE ────────────────────────────────────────────────────────────
        ve_json = base / f"{case}_boundary_deformation_map_output.json"
        if ve_json.exists():
            ve = parse_ve_json(ve_json)
            if ve["success"] != 1:
                print(f"WARNING: VE did NOT succeed for case {i} "
                      f"(opti_result={ve['opti_result']})")
                all_ok = False
            timing = ve["timing"]
            rows["REVE"].append(f"{timing:.0f}s" if timing is not None else "?s")
        else:
            print(f"WARNING: missing VE json for case {i}: {ve_json}")
            rows["REVE"].append("N/A")

        # ── TLC ───────────────────────────────────────────────────────────
        tlc_txt = base / f"{case}_TLC_stats.txt"
        if tlc_txt.exists():
            rows["TLC"].append(str(parse_stats_txt(tlc_txt)["success"]))
        else:
            print(f"WARNING: missing TLC stats for case {i}: {tlc_txt}")
            rows["TLC"].append("N/A")

        # ── FOF ───────────────────────────────────────────────────────────
        fof_txt = base / f"{case}_FOF_stats.txt"
        if fof_txt.exists():
            rows["FOF"].append(str(parse_stats_txt(fof_txt)["success"]))
        else:
            print(f"WARNING: missing FOF stats for case {i}: {fof_txt}")
            rows["FOF"].append("N/A")

        # ── TLC on VE output ──────────────────────────────────────────────
        tlc_ve_txt = base / f"{case}_boundary_deformation_map_output_TLC_stats.txt"
        if tlc_ve_txt.exists():
            rows["TLC\'"].append(str(parse_stats_txt(tlc_ve_txt)["success"]))
        else:
            print(f"WARNING: missing TLC-on-VE stats for case {i}: {tlc_ve_txt}")
            rows["TLC\'"].append("N/A")

        # ── FOF on VE output ──────────────────────────────────────────────
        fof_ve_txt = base / f"{case}_boundary_deformation_map_output_FOF_stats.txt"
        if fof_ve_txt.exists():
            rows["FOF\'"].append(str(parse_stats_txt(fof_ve_txt)["success"]))
        else:
            print(f"WARNING: missing FOF-on-VE stats for case {i}: {fof_ve_txt}")
            rows["FOF\'"].append("N/A")

    # ── Build image ───────────────────────────────────────────────────────────
    method_names = list(rows.keys())
    n_rows = len(method_names)
    n_cols = n_cases

    GREEN  = "#81c784"
    RED    = "#e57373"

    # ── Visual parameters (tweak to taste) ───────────────────────────────
    cell_size         = 0.9   # inches per cell (width)
    cell_height_scale = 0.5   # row height relative to cell_size (e.g. 0.5 = half as tall)
    left_pad          = 0.9   # inches reserved for row labels
    top_pad           = 0.7   # inches reserved for column headers
    bottom_pad        = 0.1
    right_pad         = 0.1
    fontsize_cell     = 14    # font size for VE timing labels inside squares
    fontsize_headers  = 15    # font size for row/column labels
    fontsize_titles   = 16    # font size for "digging index" and "method" titles
    gap               = 0.08  # half-gap between adjacent squares (in data units)
    method_title_x    = -left_pad  # x position of the rotated "method" title (more negative = further left)
    sq_pad            = 0.01  # extra inset on each side to shrink the square
    # ─────────────────────────────────────────────────────────────────────

    cell_h = cell_size * cell_height_scale
    fig_w  = left_pad + n_cols * cell_size + right_pad
    fig_h  = top_pad  + n_rows * cell_h    + bottom_pad

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_xlim(0, n_cols)
    ax.set_ylim(0, n_rows * cell_height_scale)
    ax.axis("off")

    # Coordinate helpers (row 0 = top)
    def cx(col): return col + 0.5
    def cy(row): return (n_rows - row - 0.5) * cell_height_scale

    sq_offset_x = gap + sq_pad
    sq_offset_y = (gap + sq_pad) * cell_height_scale
    square_w    = 1 - 2 * sq_offset_x
    square_h    = cell_height_scale - 2 * sq_offset_y

    for r, method in enumerate(method_names):
        vals = rows[method]
        for c, val in enumerate(vals):
            if method == "VE":
                color = GREEN
                label = val          # e.g. "12s"
            elif val == "N/A":
                color = "#bdbdbd"
                label = "N/A"
            elif val == "1":
                color = GREEN
                label = ""
            else:
                color = RED
                label = ""

            rect = mpatches.FancyBboxPatch(
                (c + sq_offset_x, cy(r) - 0.5 * cell_height_scale + sq_offset_y),
                square_w, square_h,
                boxstyle="round,pad=0.04",
                linewidth=0,
                color=color,
                transform=ax.transData,
                clip_on=False,
            )
            ax.add_patch(rect)

            text_color = "black" if method == "VE" else "white"
            ax.text(cx(c), cy(r), label,
                    ha="center", va="center",
                    fontsize=fontsize_cell, fontweight="bold",
                    color=text_color, transform=ax.transData)

    # ── Column headers ("digging index") ─────────────────────────────────────
    for c in range(n_cols):
        ax.text(cx(c), n_rows * cell_height_scale + 0.15, str(c + 1),
                ha="center", va="bottom",
                fontsize=fontsize_headers, color="#333333",
                transform=ax.transData)

    ax.text(n_cols / 2, n_rows * cell_height_scale + 0.52, "digging index",
            ha="center", va="bottom",
            fontsize=fontsize_titles, fontweight="bold", color="#333333",
            transform=ax.transData)

    # ── Row headers ───────────────────────────────────────────────────────────
    for r, method in enumerate(method_names):
        ax.text(-0.15, cy(r), method,
                ha="right", va="center",
                fontsize=fontsize_headers, color="#333333",
                transform=ax.transData)

    # Rotated "method" axis title
    ax.text(method_title_x, n_rows * cell_height_scale / 2, "method",
            ha="center", va="center",
            fontsize=fontsize_titles, fontweight="bold", color="#333333",
            rotation=90, transform=ax.transData)

    # ── Layout & save ─────────────────────────────────────────────────────────
    fig.subplots_adjust(
        left   = left_pad  / fig_w,
        right  = 1 - right_pad / fig_w,
        bottom = bottom_pad / fig_h,
        top    = 1 - top_pad / fig_h,
    )

    if args.plot_output_path or not args.show_plots:
        out_path = args.plot_output_path if args.plot_output_path else Path(args.input_directory) / "results_summary.png"
        fig.savefig(out_path, dpi=150, bbox_inches="tight",
                    facecolor="white", edgecolor="none")
        print(f"Image saved to: {out_path}")
    if args.show_plots:
        plt.show()
    plt.close(fig)

    if all_ok:
        print("All VE runs succeeded ✓")


if __name__ == "__main__":
    main()