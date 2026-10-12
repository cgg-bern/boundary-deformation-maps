import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os
import numpy as np
import matplotlib.colors as colors
import math
from matplotlib.ticker import FuncFormatter, FixedLocator

def signed_log10(v):
    SIGN_MASK = np.sign(v)
    return np.log10(np.abs(v) + 1e-9) * SIGN_MASK

def signed_log_formatter(val, pos=None):
    """Formatter for signed log scale."""
    if val == 0:
        return '0'
    sign = '-' if val < 0 else ''
    return f'{sign}{10**abs(val):.0e}'

def generate_plots(csv_filepath, max_height_px=1000, marker_size=50, signed_log_grid=False):

    try:
        df = pd.read_csv(csv_filepath, sep=';')
    except FileNotFoundError:
        print(f"Error: CSV file not found at '{csv_filepath}'")
        return
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        return

    x_col = 'parameter1'
    y_col = 'parameter2'
    facet_col = 'parameter3'
    deform_col = 'deformation'
    plot_fields = ['opti_iters', 'converged', 'final_SDE', 'volume']
    plot_fields = ['opti_iters', 'converged', 'final_SDE']

    unique_param3_values = sorted(df[facet_col].unique())
    num_facets = len(unique_param3_values)
    num_cols = num_facets  # Ensure we have enough columns
    num_rows_total = len(plot_fields)

    # Determine aspect ratio and calculate figure height
    aspect_ratio = num_cols / num_rows_total
    fig_width_inches = 14  # Base width, adjust as needed
    fig_height_inches_estimated = fig_width_inches / aspect_ratio

    # Convert estimated height to pixels (assuming a DPI of 100)
    estimated_height_px = fig_height_inches_estimated * 100

    # Scale height if it exceeds the maximum
    if estimated_height_px > max_height_px:
        scale_factor = max_height_px / estimated_height_px
        fig_height_inches = fig_height_inches_estimated * scale_factor
        fig_width_inches = fig_width_inches * scale_factor # Scale width proportionally
    else:
        fig_height_inches = fig_height_inches_estimated

    fig, all_axes = plt.subplots(num_rows_total, num_cols, figsize=(fig_width_inches * 1.3, fig_height_inches)) # Increased width

    if num_rows_total == 1 and num_cols == 1:
        all_axes = np.array([[all_axes]])
    elif num_rows_total == 1:
        all_axes = all_axes[np.newaxis, :]
    elif num_cols == 1:
        all_axes = all_axes[:, np.newaxis]

    fig.suptitle("Left: domain deformation, right: codomain deformation", fontsize=16, y=0.92) # Adjusted uptitle position

    triangle_vertices_x = [-0.5, 0.5, 0, -0.5]
    triangle_vertices_y = [-0.288675, -0.288675, 0.57735, -0.288675]

    for row_idx, value_col in enumerate(plot_fields):
        for col_idx, param3_value in enumerate(unique_param3_values):
            ax = all_axes[row_idx, col_idx]
            subset_df = df[df[facet_col] == param3_value].copy()

            filtered_df = subset_df[~((subset_df[x_col] <= 0) & (subset_df[deform_col] == 0)) &
                                    ~((subset_df[x_col] > 0) & (subset_df[deform_col] == 1))]

            deform_1_df = filtered_df[filtered_df[deform_col] == 1].copy()
            deform_0_df = filtered_df[filtered_df[deform_col] == 0].copy()

            x_data_1 = deform_1_df[x_col].values
            y_data_1 = deform_1_df[y_col].values
            x_data_0 = deform_0_df[x_col].values
            y_data_0 = deform_0_df[y_col].values
            non_positive_x_1 = deform_1_df.loc[deform_1_df['final_SDE'] <= 0, x_col].values
            non_positive_y_1 = deform_1_df.loc[deform_1_df['final_SDE'] <= 0, y_col].values
            positive_x_1 = deform_1_df.loc[deform_1_df['final_SDE'] > 0, x_col].values
            positive_y_1 = deform_1_df.loc[deform_1_df['final_SDE'] > 0, y_col].values
            non_positive_x_0 = deform_0_df.loc[deform_0_df['final_SDE'] <= 0, x_col].values
            non_positive_y_0 = deform_0_df.loc[deform_0_df['final_SDE'] <= 0, y_col].values
            positive_x_0 = deform_0_df.loc[deform_0_df['final_SDE'] > 0, x_col].values
            positive_y_0 = deform_0_df.loc[deform_0_df['final_SDE'] > 0, y_col].values

            ax.set_yticks([])  # Always remove y-axis ticks and labels
            ax.tick_params(axis='x', which='major', labelsize='small') # Make x-axis tick labels smaller

            if not filtered_df.empty and value_col in filtered_df.columns:
                if signed_log_grid:
                    x_data_1 = signed_log10(x_data_1)
                    y_data_1 = signed_log10(y_data_1)
                    x_data_0 = signed_log10(x_data_0)
                    y_data_0 = signed_log10(y_data_0)
                    positive_x_1 = signed_log10(positive_x_1)
                    positive_y_1 = signed_log10(positive_y_1)
                    positive_x_0 = signed_log10(positive_x_0)
                    positive_y_0 = signed_log10(positive_y_0)
                    non_positive_x_1 = signed_log10(non_positive_x_1)
                    non_positive_y_1 = signed_log10(non_positive_y_1)
                    non_positive_x_0 = signed_log10(non_positive_x_0)
                    non_positive_y_0 = signed_log10(non_positive_y_0)

                    # Define representative tick locations in the transformed space for x-axis
                    tick_positions_x = [np.min(x_data_1), 0, np.max(x_data_0)]
                    #tick_positions_x = [np.min(x_data_1), 0.5 * np.min(x_data_1), 0, 0.5 * np.max(x_data_0), np.max(x_data_0)]
                    ax.xaxis.set_major_locator(FixedLocator(tick_positions_x))
                    ax.xaxis.set_major_formatter(FuncFormatter(signed_log_formatter))


                if value_col == 'converged':
                    cmap = plt.cm.get_cmap('viridis', 2)
                    scatter1 = ax.scatter(x_data_1, y_data_1, c=deform_1_df[value_col], cmap=cmap, marker='o', s=marker_size, alpha=0.7, vmin=0, vmax=1)
                    scatter0 = ax.scatter(x_data_0, y_data_0, c=deform_0_df[value_col], cmap=cmap, marker='o', s=marker_size, alpha=0.7, vmin=0, vmax=1)
                    cbar = fig.colorbar(scatter1, ax=ax, ticks=[0, 1], shrink=0.6, location='left', pad=0.05)
                    cbar.ax.yaxis.set_label_position('left')
                    cbar.ax.yaxis.set_ticks_position('left')
                elif value_col == 'final_SDE':
                    # Log color scale for positive final_SDE values
                    positive_sde_1 = deform_1_df.loc[deform_1_df['final_SDE'] > 0, 'final_SDE'].values
                    positive_sde_0 = deform_0_df.loc[deform_0_df['final_SDE'] > 0, 'final_SDE'].values

                    min_val_1 = positive_sde_1.min()
                    max_val_1 = positive_sde_1.max()
                    norm_1 = colors.LogNorm(vmin=min_val_1, vmax=max_val_1)
                    cmap_1 = plt.cm.viridis
                    scatter1 = ax.scatter(positive_x_1, positive_y_1, c=positive_sde_1, cmap=cmap_1, norm=norm_1, marker='o', s=marker_size, alpha=0.7, edgecolors='none')
                    cbar_1 = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap_1, norm=norm_1), ax=ax, shrink=0.6, location='left', pad=0.05)
                    cbar_1.ax.yaxis.set_label_position('left')
                    cbar_1.ax.yaxis.set_ticks_position('left')
                    cbar_1.ax.yaxis.set_offset_position('right')

                    min_val_0 = positive_sde_0.min()
                    max_val_0 = positive_sde_0.max()
                    norm_0 = colors.LogNorm(vmin=min_val_0, vmax=max_val_0)
                    cmap_0 = plt.cm.plasma # Using a different colormap to distinguish
                    scatter0 = ax.scatter(positive_x_0, positive_y_0, c=positive_sde_0, cmap=cmap_0, norm=norm_0, marker='o', s=marker_size, alpha=0.7, edgecolors='none')
                    cbar_0 = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap_0, norm=norm_0), ax=ax, shrink=0.6, location='right', pad=0.05)
                    cbar_0.ax.yaxis.set_label_position('right')
                    cbar_0.ax.yaxis.set_ticks_position('right')
                    cbar_0.ax.yaxis.set_offset_position('left')

                    # Plot empty black circles for non-positive values
                    ax.scatter(non_positive_x_1, non_positive_y_1, marker='o', s=marker_size, facecolors='none', edgecolors='black', alpha=0.7)
                    ax.scatter(non_positive_x_0, non_positive_y_0, marker='o', s=marker_size, facecolors='none', edgecolors='black', alpha=0.7)

                elif value_col == 'volume':
                    # Log color scale for volume
                    volume_1 = deform_1_df[value_col].values
                    volume_0 = deform_0_df[value_col].values

                    min_val_1 = volume_1.min() if volume_1.size > 0 and volume_1.min() > 0 else 1e-6
                    max_val_1 = volume_1.max() if volume_1.size > 0 else 1
                    norm_1 = colors.LogNorm(vmin=min_val_1, vmax=max_val_1)
                    cmap_1 = plt.cm.viridis
                    scatter1 = ax.scatter(x_data_1, y_data_1, c=volume_1, cmap=cmap_1, norm=norm_1, marker='o', s=marker_size, alpha=0.7)
                    cbar_1 = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap_1, norm=norm_1), ax=ax, shrink=0.6, location='left', pad=0.05)
                    cbar_1.ax.yaxis.set_label_position('left')
                    cbar_1.ax.yaxis.set_ticks_position('left')
                    cbar_1.ax.yaxis.set_offset_position('right')

                    min_val_0 = volume_0.min() if volume_0.size > 0 and volume_0.min() > 0 else 1e-6
                    max_val_0 = volume_0.max() if volume_0.size > 0 else 1
                    norm_0 = colors.LogNorm(vmin=min_val_0, vmax=max_val_0)
                    cmap_0 = plt.cm.plasma # Using a different colormap to distinguish
                    scatter0 = ax.scatter(x_data_0, y_data_0, c=volume_0, cmap=cmap_0, norm=norm_0, marker='o', s=marker_size, alpha=0.7)
                    cbar_0 = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap_0, norm=norm_0), ax=ax, shrink=0.6, location='right', pad=0.05)
                    cbar_0.ax.yaxis.set_label_position('right')
                    cbar_0.ax.yaxis.set_ticks_position('right')
                    cbar_0.ax.yaxis.set_offset_position('left')
                else:
                    # Independent linear color scales for other fields
                    min_val_1 = deform_1_df[value_col].min()
                    max_val_1 = deform_1_df[value_col].max()
                    if(signed_log_grid):
                        norm_1 = colors.LogNorm(vmin=min_val_1, vmax=max_val_1)
                    else:
                        norm_1 = colors.Normalize(vmin=min_val_1, vmax=max_val_1)
                    cmap_1 = plt.cm.viridis
                    scatter1 = ax.scatter(x_data_1, y_data_1, c=deform_1_df[value_col], cmap=cmap_1, norm=norm_1, marker='o', s=marker_size, alpha=0.7)
                    cbar_1 = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap_1, norm=norm_1), ax=ax, shrink=0.6, location='left', pad=0.05)
                    cbar_1.ax.yaxis.set_label_position('left')
                    cbar_1.ax.yaxis.set_ticks_position('left')
                    cbar_1.ax.yaxis.set_offset_position('right')

                    min_val_0 = deform_0_df[value_col].min()
                    max_val_0 = deform_0_df[value_col].max()
                    norm_0 = colors.Normalize(vmin=min_val_0, vmax=max_val_0)
                    cmap_0 = plt.cm.plasma # Using a different colormap to distinguish
                    scatter0 = ax.scatter(x_data_0, y_data_0, c=deform_0_df[value_col], cmap=cmap_0, norm=norm_0, marker='o', s=marker_size, alpha=0.7)
                    cbar_0 = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap_0, norm=norm_0), ax=ax, shrink=0.6, location='right', pad=0.05)
                    cbar_0.ax.yaxis.set_label_position('right')
                    cbar_0.ax.yaxis.set_ticks_position('right')
                    cbar_0.ax.yaxis.set_offset_position('left')

                if not signed_log_grid:
                    ax.plot(np.array(triangle_vertices_x), np.array(triangle_vertices_y), 'k-', linewidth=2)

                ax.set_title(f'height = {param3_value:.1e}')
                ax.set_aspect('equal', adjustable='box')
                ax.set_xlabel('')
                ax.set_ylabel('')

            else:
                ax.set_title(f'height = {param3_value:.1e}')
                ax.set_aspect('equal', adjustable='box')
                ax.axis('off')

            if row_idx == 0:
                ax.set_title(f'height = {param3_value:.1e}')
            if col_idx == 0:
                ax.set_ylabel(value_col, rotation=90, labelpad=60)

    plt.tight_layout(rect=[0.05, 0.02, 0.95, 0.92], h_pad=0.5, w_pad=0.7)
    plt.subplots_adjust(left=0.15) # Adjust this value

    base_name = os.path.splitext(os.path.basename(csv_filepath))[0]
    output_filename = f"{base_name}_plots.png"

    plt.savefig(output_filename, dpi=300)
    plt.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate multi-row faceted scatter plots with adjusted left margin.")
    parser.add_argument('--csv-file', type=str, required=True, help='Path to the input CSV file.')
    parser.add_argument('--max-height', type=int, default=1000, help='Maximum height of the figure in pixels.')
    parser.add_argument('--marker-size', type=float, default=50.0, help='Size of the scatter plot markers.')
    parser.add_argument('--log-grid', action='store_true', help='Use log scale for x and y axes of all subplots.')
    parser.add_argument('--signed-log-grid', action='store_true', help='Use signed log-like transformation for x and y coordinates.')
    args = parser.parse_args()
    generate_plots(args.csv_file, args.max_height, args.marker_size, args.signed_log_grid)