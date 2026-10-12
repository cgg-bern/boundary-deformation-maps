#!/usr/bin/env python3
import argparse
import json
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec
import numpy as np

def main():
    parser = argparse.ArgumentParser(description="Plot optimization history from JSON file.")
    parser.add_argument("--json-input-path", required=True, help="Path to the input JSON file.")
    parser.add_argument("--plot-output-path", required=True, help="Path to save the output PNG plot.")
    parser.add_argument("--show-plots", action="store_true", help="Also display the plot on screen.")
    args = parser.parse_args()

    with open(args.json_input_path, "r") as f:
        data = json.load(f)

    # --- Fields ---
    # Plot 1: top (log scale)
    top_fields = {
        "unweighted_initial_primary_virtual_energy_history": "unweighted initial primary penalty",
        "unweighted_final_primary_virtual_energy_history": "unweighted final primary penalty",
        "unweighted_initial_secondary_virtual_energy_history": "unweighted initial secondary penalty",
        "unweighted_final_secondary_virtual_energy_history": "unweighted final secondary penalty",
        "unweighted_initial_quality_energy_history": "unweighted initial quality energy",
        "unweighted_final_quality_energy_history": "unweighted final quality energy",
        "worst_quality_history": "worst quality",
        "min_relative_volume_history": "min. volume (relative to average)",
        "virtual_weight_scaling_history": "penalty weight",
    }

    # Plot 2: middle (log scale)
    mid_fields = {
        "unweighted_final_quality_gradient_norm_history": "final quality gradient norm",
        "unweighted_final_primary_virtual_gradient_norm_history": "final primary virtual gradient norm",
        "unweighted_final_secondary_virtual_gradient_norm_history": "final secondary virtual gradient norm",
        "final_total_gradient_max_coeff_history": "final total gradient max coeff",
        "final_total_Hessian_max_coeff_history": "final total hessian max coeff",
    }

    # Plot 3: bottom left axis (linear scale)
    bottom_left_fields = {
        "max_dihedral_angle_diff_deg_history": "max dihedral angle diff (deg)",
        "avg_dihedral_angle_diff_deg_history": "avg dihedral angle diff (deg)",
        "wrapped_ratio_history": "snapped ratio (%)",
    }

    # Plot 3: bottom right axis (dotted, x marker)
    bottom_right_fields = {
        "newton_iters_history": "newton iterations",
        "opti_time_s_history": "optimization time (s)",
    }

    # Timing fields for the stacked bar panel
    timing_fields = [
        ("pre_processing_time_s",  "pre-processing"),
        ("total_opti_time_s",      "optimization"),
        ("total_remeshing_time_s", "remeshing"),
        ("total_upkeep_time_s",    "upkeep"),
    ]

    info_fields = [
        ("mesh_name",       "mesh"),
        ("n_init_vertices", "#vertices"),
        ("n_init_cells",    "#cells"),
        ("runtime_s",       "runtime (s)"),
        ("opti_result",     "result"),
    ]

    # --- Build info string ---
    info_parts = []
    for key, label in info_fields:
        val = data.get(key, "N/A")
        if isinstance(val, float):
            val = f"{val:.4g}"
        info_parts.append(f"{label}: {val}")
    info_str = "    |    ".join(info_parts)

    # --- Colors (tab20, assigned in order across all fields) ---
    palette = plt.cm.tab20.colors
    all_keys = list(top_fields) + list(mid_fields) + list(bottom_left_fields) + list(bottom_right_fields)
    field_colors = {key: palette[i % len(palette)] for i, key in enumerate(all_keys)}

    # Histogram colors matched to their corresponding curve colors
    energy_hist_color = field_colors["unweighted_final_primary_virtual_energy_history"]
    quality_hist_color = field_colors["worst_quality_history"]

    # --- Histogram data ---
    quality_history = data.get("full_quality_history", [])
    energy_history = data.get("full_final_virtual_energies_history", [])
    n_iters = len(quality_history)
    has_quality = n_iters > 0
    has_energy = len(energy_history) > 0
    n_hist_rows = int(has_energy) + int(has_quality)

    # --- Layout ---
    height_px = 1080       # total figure height in pixels; width = height * 16/9
    dpi = 150
    height_in = height_px / dpi
    width_in = height_in * 16 / 9

    n_cols = max(n_iters, 1)

    fig = plt.figure(figsize=(width_in, height_in), dpi=dpi)
    fig.subplots_adjust(top=0.95, bottom=0.05, left=0.01, right=0.95)

    # Top-level GridSpec: 2 columns — left timing panel and right content area
    timing_col_weight = 1   # narrow left panel
    content_col_weight = 9  # wide right content
    top_gs = gridspec.GridSpec(
        1, 2,
        figure=fig,
        width_ratios=[timing_col_weight, content_col_weight],
        wspace=0.08,
    )

    # =========================================================
    # LEFT PANEL: timing stacked bar
    # =========================================================
    ax_timing = fig.add_subplot(top_gs[0, 0])

    runtime_s = data.get("runtime_s", None)
    timing_colors = ["#7bafd4", "#f4a460", "#8fbc8f", "#d98080", "#c0c0c0"]
    timing_values = [data.get(key, 0.0) or 0.0 for key, _ in timing_fields]
    timing_labels = [label for _, label in timing_fields]

    # Reverse order: pre-processing at the top, upkeep at the bottom
    timing_values_rev = list(reversed(timing_values))
    timing_labels_rev = list(reversed(timing_labels))
    timing_colors_rev = list(reversed(timing_colors[:len(timing_values)]))

    bar_x = 0.5
    bar_width = 0.5
    bottom_val = 0.0
    segment_bottoms = []
    for val, color in zip(timing_values_rev, timing_colors_rev):
        ax_timing.bar(bar_x, val, width=bar_width, bottom=bottom_val,
                      color=color, edgecolor="white", linewidth=0.8)
        segment_bottoms.append(bottom_val)
        bottom_val += val

    total_height = bottom_val

    # Label on the left side of each segment, rotated 90°, and value on top of segment
    for seg_bottom, val, label in zip(segment_bottoms, timing_values_rev, timing_labels_rev):
        seg_mid = seg_bottom + val / 2
        # Left-side label, rotated 90°
        ax_timing.text(bar_x - bar_width / 2 - 0.04, seg_mid, label,
                       ha="right", va="center", fontsize=6, rotation=90,
                       transform=ax_timing.transData)
        # Value written on top of the segment (inside, near the top edge)
        if val > total_height * 0.04:   # only if segment is tall enough
            ax_timing.text(bar_x, seg_bottom + val - total_height * 0.005,
                           f"{val:.3g}s", ha="center", va="top", fontsize=5.5)

    # Total runtime label above the bar
    if isinstance(runtime_s, (int, float)):
        ax_timing.text(bar_x, total_height + total_height * 0.02,
                       f"Total: {runtime_s:.4g}s",
                       ha="center", va="bottom", fontsize=6.5, fontweight="bold")

    ax_timing.set_xlim(0, 1)
    ax_timing.set_ylim(0, total_height * 1.15)
    ax_timing.axis("off")

    # =========================================================
    # RIGHT CONTENT AREA: curves + histograms
    # =========================================================
    content_spec = top_gs[0, 1]

    if n_hist_rows > 0:
        hist_height_ratio = 1.0
        curves_weight = 6
        hists_weight = hist_height_ratio * n_hist_rows
        inner_content_gs = gridspec.GridSpecFromSubplotSpec(
            2, 1,
            subplot_spec=content_spec,
            height_ratios=[curves_weight, hists_weight],
            hspace=0.1,   # hspace2: space between bottom curve plot and histograms
        )
        inner_curves_gs = gridspec.GridSpecFromSubplotSpec(
            3, 1,
            subplot_spec=inner_content_gs[0],
            hspace=0.05,   # hspace1: space between the curve plots
        )
        ax_top = fig.add_subplot(inner_curves_gs[0])
        ax_mid = fig.add_subplot(inner_curves_gs[1])
        ax_bot = fig.add_subplot(inner_curves_gs[2])
        # Histogram axes will be added via fig.add_axes after layout is known (see below)
        energy_axes = []
        quality_axes = []

    else:
        # No histograms
        inner_curves_gs = gridspec.GridSpecFromSubplotSpec(
            3, 1,
            subplot_spec=content_spec,
            hspace=0.2,
        )
        ax_top = fig.add_subplot(inner_curves_gs[0])
        ax_mid = fig.add_subplot(inner_curves_gs[1])
        ax_bot = fig.add_subplot(inner_curves_gs[2])
        energy_axes = []
        quality_axes = []

    # --- Info text ---
    fig.text(0.5, 0.97, info_str, ha="center", va="center", fontsize=8.5,
             bbox=dict(boxstyle="round,pad=0.4", facecolor="#f0f0f0", edgecolor="#aaaaaa"))

    # --- no-shift keys ---
    no_shift_keys = {
        "wrapped_ratio_history",
        "max_dihedral_angle_diff_deg_history",
        "avg_dihedral_angle_diff_deg_history",
        "worst_quality_history",
        "min_relative_volume_history",
    }

    # Compute x range
    all_series_lengths = [
        len(data.get(k, [])) for k in all_keys if data.get(k)
    ]
    x_right = (max(all_series_lengths) if all_series_lengths else 10) - 0.5

    # =========================================================
    # TOP PLOT
    # =========================================================
    ax_top.set_yscale("log")
    ax_top.yaxis.grid(True, linestyle="--", linewidth=0.5, color="lightgray", alpha=0.7)
    ax_top.xaxis.grid(True, linestyle="--", linewidth=0.5, color="lightgray", alpha=0.7)
    ax_top.set_axisbelow(True)

    top_lines = []
    for key, label in top_fields.items():
        y = data.get(key, [])
        if y:
            color = field_colors[key]
            shift = 0 if key in no_shift_keys else 1
            x = list(range(shift, len(y) + shift))
            line, = ax_top.plot(x, y, marker="o", markersize=4, label=label, color=color, linewidth=1.8)
            top_lines.append(line)

    ax_top.legend(handles=top_lines, loc="upper left", fontsize=6, framealpha=0.85)
    ax_top.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    ax_top.set_xlim(-0.5, x_right)
    ax_top.tick_params(axis="x", which="both", bottom=False, labelbottom=False)

    # =========================================================
    # MIDDLE PLOT
    # =========================================================
    ax_mid.set_yscale("log")
    ax_mid.yaxis.grid(True, linestyle="--", linewidth=0.5, color="lightgray", alpha=0.7)
    ax_mid.xaxis.grid(True, linestyle="--", linewidth=0.5, color="lightgray", alpha=0.7)
    ax_mid.set_axisbelow(True)

    mid_lines = []
    for key, label in mid_fields.items():
        y = data.get(key, [])
        if y:
            color = field_colors[key]
            shift = 0 if key in no_shift_keys else 1
            x = list(range(shift, len(y) + shift))
            line, = ax_mid.plot(x, y, marker="o", markersize=4, label=label, color=color, linewidth=1.8)
            mid_lines.append(line)

    ax_mid.legend(handles=mid_lines, loc="upper left", fontsize=6, framealpha=0.85)
    ax_mid.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    ax_mid.set_xlim(-0.5, x_right)
    ax_mid.tick_params(axis="x", which="both", bottom=False, labelbottom=False)

    # =========================================================
    # BOTTOM PLOT
    # =========================================================
    ax_bot.yaxis.grid(True, linestyle="--", linewidth=0.5, color="lightgray", alpha=0.7)
    ax_bot.xaxis.grid(True, linestyle="--", linewidth=0.5, color="lightgray", alpha=0.7)
    ax_bot.set_axisbelow(True)

    bot_left_lines = []
    for key, label in bottom_left_fields.items():
        y = data.get(key, [])
        if y:
            color = field_colors[key]
            if key == "wrapped_ratio_history":
                y = [v * 100 for v in y]
            shift = 0 if key in no_shift_keys else 1
            x = list(range(shift, len(y) + shift))
            line, = ax_bot.plot(x, y, marker="s", markersize=4, label=label, color=color, linewidth=1.8, linestyle="--")
            bot_left_lines.append(line)

    ax_bot_r = ax_bot.twinx()
    ax_bot_r.spines["right"].set_color("black")
    ax_bot_r.tick_params(axis="y", colors="black")

    bot_right_lines = []
    for key, label in bottom_right_fields.items():
        y = data.get(key, [])
        if y:
            color = field_colors[key]
            shift = 0 if key in no_shift_keys else 1
            x = list(range(shift, len(y) + shift))
            line, = ax_bot_r.plot(x, y, marker="x", markersize=5, label=label, color=color, linewidth=1.8, linestyle=":")
            bot_right_lines.append(line)
            # Add cumulative sum curve for opti_time_s_history on the left axis
            if key == "opti_time_s_history":
                y_cum = np.cumsum(y)
                y_acc = list(y_cum / y_cum[-1] * 100)
                acc_color = palette[len(all_keys) % len(palette)]
                line_acc, = ax_bot.plot(x, y_acc, marker="s", markersize=5, label="cumulative optimization time (%)",
                                        color=acc_color, linewidth=1.8, linestyle="--")
                bot_left_lines.append(line_acc)

    left_legend = ax_bot.legend(handles=bot_left_lines, loc="upper left", fontsize=6, framealpha=0.85)
    ax_bot.add_artist(left_legend)
    ax_bot_r.legend(handles=bot_right_lines, loc="upper right", fontsize=6, framealpha=0.85)
    ax_bot.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    ax_bot.set_xlim(-0.5, x_right)

    # =========================================================
    # DEEP REMESHING ITERATION MARKERS
    # =========================================================
    deep_remeshing_iters = data.get("deep_remeshing_iterations", [])
    print("deep remeshing iterations: ", deep_remeshing_iters)
    for idx in deep_remeshing_iters:
        for ax in (ax_top, ax_mid, ax_bot):
            ax.axvline(x=idx, color="#1f77b4", linewidth=1.2, linestyle="-", alpha=0.6, zorder=3)

    # =========================================================
    # HISTOGRAMS — placed via coordinate transform for exact x-axis alignment
    # =========================================================
    def make_hist_axes(row_start_iter, n_histograms, hist_row_top_fig, hist_row_bot_fig):
        """
        Create n_histograms axes for one histogram row.
        Each histogram i corresponds to data iteration (row_start_iter + i).
        The left edge of histogram i is aligned with x-tick (row_start_iter + i) on ax_bot,
        and the right edge with x-tick (row_start_iter + i + 1).
        Returns a list of axes in order.
        """
        fig.canvas.draw()  # needed to resolve transforms
        trans = ax_bot.transData + fig.transFigure.inverted()
        axes = []
        for i in range(n_histograms):
            tick = row_start_iter + i
            x_left_fig,  _ = trans.transform((tick - 0.5, 0))
            x_right_fig, _ = trans.transform((tick + 0.5, 0))
            h_pad = 0.003  # horizontal padding between histograms (figure fraction)
            w = x_right_fig - x_left_fig - 2 * h_pad
            h = hist_row_top_fig - hist_row_bot_fig
            ax_h = fig.add_axes([x_left_fig + h_pad, hist_row_bot_fig, w, h])
            axes.append(ax_h)
        return axes

    def draw_histograms(axes, history, row_label, x_min, x_max, start_iter=0, color="steelblue", show_sum=False):
        """Fill a list of histogram axes with data. axes has no None entries."""
        bins = np.arange(x_min, x_max + 1)
        tick_positions = np.arange(x_min, x_max + 1)
        selected_ticks = [tick_positions[1], tick_positions[len(tick_positions) // 2], tick_positions[-2]]

        for i, ax_h in enumerate(axes):
            arr = history[start_iter + i]
            vals = np.array(arr, dtype=float)
            vals = vals[vals > 0]
            log_vals = np.clip(np.log10(vals), x_min, x_max) if len(vals) else np.array([])
            total_sum = vals.sum() if len(vals) else 0

            if len(log_vals):
                counts, _ = np.histogram(log_vals, bins=bins)
                total = counts.sum()
                norm_counts = counts / total if total > 0 else counts
                bin_width = bins[1] - bins[0]
                ax_h.bar(bins[:-1] + bin_width * (1 - 0.85) / 2, norm_counts, width=bin_width * 0.85,
                         color=color, edgecolor="black", linewidth=0.4, align="edge")

            ax_h.set_xlim(x_min, x_max)
            ax_h.tick_params(axis="x", labelsize=6)
            ax_h.tick_params(axis="y", labelsize=6)
            ax_h.set_xticks(selected_ticks)
            def fmt_x(v, _, x_min=x_min, x_max=x_max):
                iv = int(v)
                if iv == x_min:
                    return f"$\\leq10^{{{iv}}}$"
                elif iv == x_max:
                    return f"$\\geq10^{{{iv}}}$"
                return f"$10^{{{iv}}}$"
            ax_h.xaxis.set_major_formatter(ticker.FuncFormatter(fmt_x))

            ax_h.set_yscale("log")
            ax_h.set_ylim(1e-4, 2e0)
            ax_h.set_yticks([1e-3, 1e-2, 1e-1, 1e0])
            ax_h.yaxis.set_minor_locator(ticker.NullLocator())
            ax_h.xaxis.set_minor_locator(ticker.NullLocator())

            if i > 0:
                ax_h.set_yticklabels([])
                ax_h.tick_params(axis="x", labelbottom=False)
            else:
                ax_h.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"$10^{{{int(np.log10(v))}}}$"))
                ax_h.set_ylabel(row_label, fontsize=7)

            if show_sum:
                ax_h.set_xlabel(f"sum={total_sum:.3g}", fontsize=6, labelpad=2)

    if n_hist_rows > 0:
        # Compute vertical band for the histogram rows using inner_content_gs[1]
        fig.canvas.draw()
        hist_spec = inner_content_gs[1]
        hist_ss = hist_spec.get_position(fig)   # a Bbox in figure coordinates
        hist_fig_top = hist_ss.y1
        hist_fig_bot = hist_ss.y0
        hist_row_height = (hist_fig_top - hist_fig_bot) / n_hist_rows
        # small vertical padding between rows (fraction of row height)
        v_pad = 0.015

        if has_energy:
            row_idx = 0
            row_top = hist_fig_top - row_idx * hist_row_height - v_pad
            row_bot = hist_fig_top - (row_idx + 1) * hist_row_height + v_pad
            energy_axes = make_hist_axes(
                row_start_iter=1, n_histograms=n_iters - 1,
                hist_row_top_fig=row_top, hist_row_bot_fig=row_bot,
            )
            draw_histograms(energy_axes, energy_history, "final penalty\ncount",
                            x_min=-4, x_max=7, start_iter=0,
                            color=energy_hist_color, show_sum=False)

        if has_quality:
            row_idx = int(has_energy)
            row_top = hist_fig_top - row_idx * hist_row_height - v_pad
            row_bot = hist_fig_top - (row_idx + 1) * hist_row_height + v_pad
            quality_axes = make_hist_axes(
                row_start_iter=0, n_histograms=n_iters,
                hist_row_top_fig=row_top, hist_row_bot_fig=row_bot,
            )
            draw_histograms(quality_axes, quality_history, "quality\ncount",
                            x_min=-6, x_max=0, start_iter=0,
                            color=quality_hist_color)

    plt.savefig(args.plot_output_path, bbox_inches="tight")
    print(f"Plot saved to {args.plot_output_path}")
    if args.show_plots:
        plt.show()

if __name__ == "__main__":
    main()