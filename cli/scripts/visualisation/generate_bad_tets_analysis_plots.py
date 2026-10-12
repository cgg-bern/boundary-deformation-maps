import pandas as pd
import matplotlib.pyplot as plt
import argparse
import numpy as np
import matplotlib.gridspec as gridspec
import os  # Import the os module

def generate_loglog_plots(csv_path, param1_label='parameter1', param2_label='parameter2', param3_label='parameter3'):
    """
    Args:
        csv_path (str): The path to the input CSV file.
        param1_label (str, optional): Label for the 'parameter1' axis. Defaults to 'parameter1'.
        param2_label (str, optional): Label for the 'parameter2' axis. Defaults to 'parameter2'.
        param3_label (str, optional): Label for the 'parameter3' axis (not currently used in plots). Defaults to 'parameter3'.
    """
    try:
        df = pd.read_csv(csv_path, sep=';')
    except FileNotFoundError:
        print(f"Error: CSV file not found at '{csv_path}'.")
        return
    except Exception as e:
        print(f"An error occurred while reading the CSV file: {e}")
        return

    # Create the output filename
    base_filename = os.path.splitext(os.path.basename(csv_path))[0]  # Get filename without extension
    output_filename = f"{base_filename}_plot.png"

    # Create a figure with a 2x3 grid of subplots
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    color_opti_0 = 'tab:blue'
    color_opti_1 = 'tab:orange'
    color_param_sde_0 = 'tab:blue'
    color_param_sde_1 = 'tab:orange'
    color_row2 = 'tab:green'

    # Get unique deformation values
    deformation_values = sorted(df['deformation'].unique())
    deformation_labels = {0: 'codomain deformation', 1: 'domain deformation'}

    # Small replacement value for zeros
    replacement_value = 1e-30

    # --- Nested Subplots in axes[0, 0] with adjusted sizes ---
    main_ax = axes[0, 0]
    gs0 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=main_ax.get_subplotspec(), height_ratios=[7, 3])

    # Top nested subplot for opti_iters (70% height)
    ax_opti = fig.add_subplot(gs0[0, 0], sharex=main_ax)
    for i, deformation in enumerate(deformation_values):
        color = color_opti_0 if deformation == 0 else color_opti_1
        subset = df[df['deformation'] == deformation].copy()
        subset['parameter1'] = subset['parameter1'].replace(0, replacement_value)
        subset['opti_iters'] = subset['opti_iters'].replace(0, replacement_value)
        label = deformation_labels[deformation]
        ax_opti.loglog(subset['parameter1'], subset['opti_iters'], linestyle='-', label=label, color=color)
    ax_opti.set_ylabel('Optimization Iterations')
    ax_opti.grid(True, which="both")
    ax_opti.legend(loc='upper left', fontsize='small')
    plt.setp(ax_opti.get_xticklabels(), visible=False) # Hide x-axis labels for the top subplot

    # Bottom nested subplot for converged (30% height)
    ax_conv = fig.add_subplot(gs0[1, 0], sharex=main_ax)
    for i, deformation in enumerate(deformation_values):
        color = color_opti_0 if deformation == 0 else color_opti_1 # Use opti colors
        subset = df[df['deformation'] == deformation].copy()
        subset['parameter1'] = subset['parameter1'].replace(0, replacement_value)
        ax_conv.plot(subset['parameter1'], subset['converged'], linestyle='--', label=deformation_labels[deformation], color=color)
    ax_conv.set_xlabel(param1_label)
    ax_conv.set_ylabel('Converged')
    ax_conv.grid(True)
    ax_conv.set_yticks([0, 1])
    ax_conv.set_ylim(-0.1, 1.1)

    # Disable the original top-left subplot's axes (to avoid overlapping labels)
    main_ax.axis('off')

    # --- Plot 2 (Top-Middle): parameter1 vs final_SDE ---
    ax = axes[0, 1]
    for i, deformation in enumerate(deformation_values):
        color = color_param_sde_0 if deformation == 0 else color_param_sde_1
        subset = df[df['deformation'] == deformation].copy()
        # Identify original zero values
        zero_mask = subset['final_SDE'] == 0
        subset['final_SDE_replaced'] = subset['final_SDE'].replace(0, replacement_value)
        subset['parameter1_replaced'] = subset['parameter1'].replace(0, replacement_value)
        label = deformation_labels[deformation]
        ax.loglog(subset['parameter1_replaced'], subset['final_SDE_replaced'], linestyle='-', label=label, color=color)
        # Plot larger, black, empty circles at original zero values
        ax.plot(subset.loc[zero_mask, 'parameter1'], [replacement_value] * zero_mask.sum(), marker='o', linestyle='none', color='black', markersize=6, fillstyle='none')
    ax.set_xlabel(param1_label)
    ax.set_ylabel('Final SDE')
    ax.grid(True, which="both")
    ax.legend(loc='best', fontsize='small')

    # --- Remove the upper right subplot (axes[0, 2]) ---
    fig.delaxes(axes[0, 2])

    # --- Plot 3 (Bottom-Left): parameter1 vs initial_SDE ---
    ax = axes[1, 0]
    subset = df.copy()
    subset['parameter1'] = subset['parameter1'].replace(0, replacement_value)
    subset['initial_SDE'] = subset['initial_SDE'].replace(0, replacement_value)
    ax.loglog(subset['parameter1'], subset['initial_SDE'], linestyle='-', color=color_row2)
    ax.set_xlabel(param1_label)
    ax.set_ylabel('Initial SDE')
    ax.grid(True, which="both")
    ax.legend(loc='best', fontsize='small')

    # --- Plot 4 (Bottom-Middle): parameter1 vs volume ---
    ax = axes[1, 1]
    subset = df.copy()
    subset['parameter1'] = subset['parameter1'].replace(0, replacement_value)
    subset['volume'] = subset['volume'].replace(0, replacement_value)
    ax.loglog(subset['parameter1'], subset['volume'], linestyle='-', color=color_row2)
    ax.set_xlabel(param1_label)
    ax.set_ylabel('Volume')
    ax.grid(True, which="both")
    ax.legend(loc='best', fontsize='small')

    # --- Plot 5 (Bottom-Right): volume vs initial_SDE ---
    ax = axes[1, 2]
    subset = df.copy()
    subset['volume'] = subset['volume'].replace(0, replacement_value)
    subset['initial_SDE'] = subset['initial_SDE'].replace(0, replacement_value)
    ax.loglog(subset['volume'], subset['initial_SDE'], linestyle='-', color=color_row2)
    ax.set_xlabel('Volume')
    ax.set_ylabel('Initial SDE')
    ax.grid(True, which="both")
    ax.legend(loc='best', fontsize='small')

    # Remove all subplot titles
    for ax in axes.flat:
        ax.set_title('')

    fig.tight_layout(w_pad=3.0)  # Increase horizontal padding

    plt.savefig(output_filename)
    plt.show()

    print(f"Final log-log plots (with two curves per plot) have been generated and saved as '{output_filename}'")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate final log-log plots from a CSV file with two curves per plot.")
    parser.add_argument("--csv-file", required=True, help="The path to the input CSV file.")
    parser.add_argument("--param1-label", default='parameter1', help="Custom label for the 'parameter1' axis.")
    parser.add_argument("--param2-label", default='parameter2', help="Custom label for the 'parameter2' axis.")
    parser.add_argument("--param3-label", default='parameter3', help="Custom label for the 'parameter3' axis (not currently used).")
    args = parser.parse_args()

    generate_loglog_plots(args.csv_file, args.param1_label, args.param2_label, args.param3_label)