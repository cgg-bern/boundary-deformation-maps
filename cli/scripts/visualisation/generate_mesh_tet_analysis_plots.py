import pandas as pd
import matplotlib.pyplot as plt
import argparse
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
import matplotlib.ticker as mticker

def plot_data(csv_file_path, output_png_path=None):
    """
    Loads a CSV file, generates scatter subplots and histograms, with log scaling/transformation
    for specified columns, and optionally saves the figure as a high-quality PNG.
    Log-scale histograms are configured to display a maximum of 6 tick labels.

    Args:
        csv_file_path (str): The path to the CSV file.
        output_png_path (str, optional): The path to save the figure as a PNG file (optional).
    """
    try:
        df = pd.read_csv(csv_file_path, sep=';')
    except FileNotFoundError:
        print(f"Error: CSV file not found at '{csv_file_path}'")
        return
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        return

    boundary_vertices = df['n_boundary_vertices'].unique()
    colors = plt.cm.get_cmap('viridis', len(boundary_vertices))  # Choose a colormap

    # Define figure and axes with a fixed 4x3 layout
    fig = Figure(figsize=(15, 12), dpi=300) # Increased DPI for higher quality
    canvas = FigureCanvasAgg(fig)
    gs = fig.add_gridspec(4, 3, height_ratios=[0.35, 0.35, 0.15, 0.15]) # Fixed 4 rows, 3 columns, height ratios

    scatter_axes = [fig.add_subplot(gs[i, j]) for i in range(1) for j in range(3)]
    hist_axes = [fig.add_subplot(gs[1 + i, j]) for i in range(3) for j in range(3)]

    plot_pairs = [
        ('initial_SDE', 'final_SDE'),
        ('initial_SDE', 'opti_iters'),
        ('final_SDE', 'opti_iters'),
        #('boundary_size_variation', 'final_SDE'),
        #('final_SDE', 'bbox_diagonal_norm')
    ]

    log_cols = ['initial_SDE', 'final_SDE', 'ref_SDE', 'volume', 'bbox_diagonal_norm', 'boundary_size_variation', 'volume_variation', 'height_variation', 'min_dihedral_angle_deg', 'max_dihedral_angle_deg', 'volume', 'min_interior_angle_deg', 'max_interior_angle_deg']

    hist_cols = ['initial_SDE', 'final_SDE', 'opti_iters', 'min_dihedral_angle_deg', 'max_dihedral_angle_deg', 'volume', 'min_interior_angle_deg', 'max_interior_angle_deg']

    for i, (x_col, y_col) in enumerate(plot_pairs):
        if i < len(scatter_axes): # Ensure we don't go beyond the defined scatter axes
            ax = scatter_axes[i]
            for j, n_vertices in enumerate(boundary_vertices):
                subset = df[df['n_boundary_vertices'] == n_vertices].copy()  # Create a copy
                x_data = subset[x_col]
                y_data = subset[y_col]

                # Apply log scale and handle non-positive values for scatter plots
                if x_col in log_cols:
                    x_data = np.where(x_data > 0, x_data, np.nan)
                if y_col in log_cols:
                    y_data = np.where(y_data > 0, y_data, np.nan)

                if np.any(x_data > 0):
                    ax.scatter(x_data, y_data, label=f'{n_vertices} vertices', color=colors(j), s=20)

            ax.set_xlabel(x_col)
            ax.set_ylabel(y_col)
            ax.set_title(f'Scatter of {y_col} vs {x_col}')
            ax.grid(True)
            if i == 0:
                ax.legend(title='n_boundary_vertices')

            # Set log scale for relevant axes in scatter plots
            if x_col in log_cols:
                ax.set_xscale('log', nonpositive='clip')
            if y_col in log_cols:
                ax.set_yscale('log', nonpositive='clip')

    for i, col in enumerate(hist_cols):
        ax = hist_axes[i]
        data = df[col].dropna()

        print("plotting histogram for ", col)
        # Calculate mean and standard deviation of the raw data
        mean_raw = np.nanmean(data)
        std_raw = np.nanstd(data)

        if col in log_cols:
            if col == "max_dihedral_angle_deg" or col == "max_interior_angle_deg":
                col = col + "_diff_to_180"
                data = 180.0 - data

            textstr = f'Mean: {mean_raw:.2e}\nStd Dev: {std_raw:.2e}'  # Scientific notation
            data_log = np.log10(data[data > 0])  # Apply log10 to positive values
            ax.hist(data_log, bins=20, color='skyblue', edgecolor='black')

            min_log = np.floor(data_log.min())
            max_log = np.ceil(data_log.max())
            num_ticks = min(6, int(max_log - min_log + 1)) # Ensure no more than 6 ticks
            ticks_log = np.linspace(min_log, max_log, num_ticks)
            labels_log = [f'$10^{{{int(t)}}}$' for t in ticks_log]
            ax.set_xticks(ticks_log)
            ax.set_xticklabels(labels_log)
            ax.set_yscale('log')

        else:
            textstr = f'Mean: {mean_raw:.2f}\nStd Dev: {std_raw:.2f}'  # Decimal notation
            ax.hist(data, bins=20, color='lightcoral', edgecolor='black')

        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=10, verticalalignment='top', bbox=props)

        ax.set_ylabel('Frequency')
        ax.set_title(f'Histogram of {col}')
        ax.grid(axis='y', alpha=0.75)

    fig.tight_layout()

    if output_png_path:
        try:
            canvas.draw()
            fig.savefig(output_png_path, dpi=300)
            width, height = fig.get_size_inches() * fig.dpi
            print(f"Figure saved as PNG: {output_png_path} (width: {int(width)} pixels, height: {int(height)} pixels)")
        except Exception as e:
            print(f"Error saving figure: {e}")

    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate scatter plots and histograms from a CSV file with log scaling.")
    parser.add_argument("--csv-file", required=True, help="Path to the CSV file.")
    parser.add_argument("--output-png", help="Path to save the figure as a PNG file (optional).")
    args = parser.parse_args()
    plot_data(args.csv_file, args.output_png)