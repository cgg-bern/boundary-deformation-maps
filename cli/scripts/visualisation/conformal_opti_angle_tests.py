import argparse
import pandas as pd
import matplotlib.pyplot as plt
import os

def plot_csv_data(input_path, log_xscale=False, log_yscale=False, output_path=None):
    """
    Reads a CSV file without a header, and plots the data.

    Args:
        input_path (str): The path to the CSV file.
        log_scale (bool): If True, the x-axis will be on a log scale.
        output_path (str): The path to save the plot as a PNG file.
    """
    if not os.path.exists(input_path):
        print(f"Error: The file '{input_path}' was not found.")
        return

    try:
        # Read the CSV file without a header
        df = pd.read_csv(input_path, header=None)

        # The first column is the x-axis data
        x_data = df.iloc[:, 0]

        # The rest of the columns are the y-axis curves
        y_data = df.iloc[:, 1:]

        # Create the plot
        plt.figure(figsize=(10, 6))

        # Plot each curve
        for col in y_data.columns:
            plt.plot(x_data, y_data[col], label=f'edge {col}')


        plt.plot(x_data, x_data, color="black", label="theta line")
        plt.xlabel('theta')

        # Set the x-axis to log scale if the flag is present
        if log_xscale:
            plt.xscale('log')
            #plt.yscale('log')
        if log_yscale:
            plt.yscale('log')

        plt.ylabel('dihedral angles')
        plt.title('CSV Data Plot')
        plt.legend()
        plt.grid(True)
        plt.tight_layout()

        # Save the plot if an output path is provided
        if output_path:
            plt.savefig(output_path)
            print(f"Plot saved to '{output_path}'")

        plt.show()

    except Exception as e:
        print(f"An error occurred: {e}")


def main():
    """
    Parses command-line arguments and calls the plotting function.
    """
    parser = argparse.ArgumentParser(description='Plot data from a CSV file.')
    parser.add_argument('--input-path', type=str, required=True,
                        help='Path to the CSV file to plot.')
    parser.add_argument('--log-x', action='store_true',
                        help='Use a logarithmic scale for the x-axis.')
    parser.add_argument('--log-y', action='store_true',
                        help='Use a logarithmic scale for the y-axis.')
    parser.add_argument('--output-path', type=str,
                        help='Path to a .png file to save the plot to.')

    args = parser.parse_args()
    plot_csv_data(args.input_path, args.log_x, args.log_y, args.output_path)


if __name__ == '__main__':
    main()