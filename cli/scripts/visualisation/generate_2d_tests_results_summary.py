import json
import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm  # Import colormap module
import numpy as np
import argparse  # For command-line arguments

def merge_json_files(directory):
    """
    Reads multiple JSON files from a directory and merges them into a Pandas DataFrame.

    Args:
        directory: The path to the directory containing the JSON files.

    Returns:
        A Pandas DataFrame containing the merged data, or None if an error occurs
        or no JSON files are found.  Prints error messages to help diagnose
        issues.
    """

    all_data = []
    json_files_found = False  # Flag to track if any JSON files were found

    try:
        for filename in os.listdir(directory):
            if filename.endswith(".json"):
                json_files_found = True
                filepath = os.path.join(directory, filename)

                try:
                    with open(filepath, 'r', encoding='utf-8') as f:  # Handle potential encoding issues
                        try:
                            data = json.load(f)
                            if isinstance(data, list): # Check if the JSON is a list of objects
                                all_data.extend(data)
                            elif isinstance(data, dict): # Check if the JSON is a single object
                                all_data.append(data)
                            else:
                                print(f"Warning: JSON file '{filename}' does not contain a list or dictionary. Skipping.")

                        except json.JSONDecodeError as e:
                            print(f"Error decoding JSON in file '{filename}': {e}")
                            return None # Or continue to the next file if you don't want to stop on error
                except FileNotFoundError:
                    print(f"Error: File '{filename}' not found (should not happen).")
                    return None
                except Exception as e:  # Catch other potential file reading errors
                    print(f"Error reading file '{filename}': {e}")
                    return None

        if not json_files_found:
            print(f"No JSON files found in the directory: {directory}")
            return None

        if not all_data: #Check if all_data is still empty after the loop
            print(f"All JSON files in {directory} were empty or contained invalid data.")
            return None

        df = pd.DataFrame(all_data)
        return df

    except FileNotFoundError:
        print(f"Error: Directory '{directory}' not found.")
        return None
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return None

def create_and_save_histograms(directory_path, output_path, show_plots):
    merged_df = merge_json_files(directory_path)

    if merged_df is not None:
        hist_fields = ["opti_result", "runtime_s", "vertex_ratio", "face_ratio", "solver_iters", "split_rounds", "total_collapse_count", "total_split_count"]
        hist_titles = ["optimization result", "runtime (s)", "vertex growth ratio", "face growth ratio", "total solver iteration count", "split rounds count", "total collapse count", "total split count"]
        #num_plots = 2 * len(hist_fields)
        line_fields = ["n_verts_history", "n_elems_history"]  # New line plot fields
        line_titles = ["#vertices", "#elements"]
        num_hist_plots = 2 * len(hist_fields)
        num_line_plots = 2 * len(line_fields)
        compression_modes = [0, 1]  # The compression modes to plot
        categories = ["extension", "compression"]

        nrows = len(compression_modes) * 2
        ncols = 5
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(22, 12), layout="constrained")

        cmap = cm.get_cmap('viridis', num_hist_plots // 2)

        title_offset = 0  # To keep track of title row index

        for c_idx, compression_mode in enumerate(compression_modes):
            if("compression_mode" in merged_df.index):
                filtered_df = merged_df[merged_df["compression_mode"] == compression_mode]  # Filter data
            else:
                filtered_df = merged_df

            # Calculate ratios
            filtered_df['vertex_ratio'] = filtered_df['n_final_vertices'] / filtered_df['n_init_vertices']
            filtered_df['face_ratio']   = filtered_df['n_final_faces']    / filtered_df['n_init_faces']

            title_row = c_idx * 2  # Calculate title row index
            title_y = (nrows + title_offset - 0.5) / (nrows + len(compression_modes))
            fig.text(0.5, 0.23 + (1-compression_mode) * 0.53 , f"{len(filtered_df)} {categories[c_idx]} cases", ha='center', va='center', fontsize=14)  # Add title

            compression_shift = compression_mode * ncols // 2
            for i, field in enumerate(hist_fields):
                if num_hist_plots > 1:
                    row = (i // (ncols-1))
                    col = (i % (ncols-1))
                    ax = axes[row + compression_shift, col]
                else:
                    ax = axes

                color = cmap(i)

                if field == "opti_result":
                    zero_percentage = (filtered_df[field] == 0).sum() / len(filtered_df[field]) * 100
                    ax.hist(filtered_df[field], bins=[-1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5], edgecolor='black', color=color, label=f"Success Rate: {zero_percentage:.1f}%")
                    ax.set_xticks([-1, 0, 1, 2, 3, 4])
                    ax.set_xticklabels(["error", "success", "reached max iterations", "setup error", "infinite initial energy", "validation failure"], rotation=30, ha='right')
                else:
                    ax.hist(filtered_df[field], edgecolor='black', color=color, label=f"Average: {filtered_df[field].mean():.2f}")  # Added label with average

                #ax.set_xlabel(field)
                #ax.set_ylabel("Frequency")
                ax.legend(loc='upper right')  # Show legend
                ax.set_title(hist_titles[i])

            #if num_plots > 1:
            #    for j in range(i + 1, nrows * ncols):
            #        row = (j // ncols)
            #        col = (j % ncols)
                    #axes[row, col].set_axis_off()

            # Line plots
            line_cmap = cm.get_cmap('plasma', num_line_plots + 1) # Colormap for line plots
            for l_idx, line_field in enumerate(line_fields):
                row = l_idx + compression_mode * 2
                col = 4
                ax = axes[row, col]

                line_color = line_cmap(row)
                for index, history in filtered_df[line_field].items(): # Iterate over each row and each array
                    ax.plot(history, color=line_color, alpha=0.5) # Plot all lines on the same subplot

                ax.set_xlabel("Iterations")
                #ax.set_ylabel(line_field)
                ax.set_title(line_titles[l_idx])
                if(compression_mode):
                    ax.set_yscale('log')  # Set y-axis to logarithmic scale


        #plt.subplots_adjust(top=91, bottom=0.1) # Add vertical spacing between subplots
        plt.tight_layout(h_pad=2)
        plt.savefig(output_path, dpi=300)  # Save as PNG with specified path and DPI
        if(show_plots):
            plt.show()
        plt.close(fig) # Close the figure to free memory
        print(f"Histograms saved to {output_path}")

    else:
        print("Could not create DataFrame, so cannot create histograms.")


if __name__ == "__main__":  # Ensures the code is only executed when the script is run directly
    parser = argparse.ArgumentParser(description="Merge JSON files and create histograms.")
    parser.add_argument("--directory_path", help="Path to the directory containing JSON files.")
    parser.add_argument("--output_path", help="Path to save the output PNG image.")
    parser.add_argument("--show_plots", action="store_true", help="states whether the plots should be shown directly, on top of being exported as an image")

    args = parser.parse_args()
    create_and_save_histograms(args.directory_path, args.output_path, args.show_plots)
