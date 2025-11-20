import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def plot_benchmark_results(file_path):
    """
    Reads benchmark results from a TSV file and plots boxplots for
    Elapsed Time and Max RSS.
    """

    # 1. Load the dataset
    # We use sep='\t' because the file extension is .tsv (Tab Separated Values)
    try:
        df = pd.read_csv(file_path, sep='\t')
    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
        return

    # Optional: Convert KB to MB for better readability
    df['Max_RSS_MB'] = df['Max_RSS_KB'] / 1024

    # 2. Set up the visual style
    sns.set_theme(style="whitegrid")

    # Create a figure with two subplots side-by-side
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # 3. Plot Elapsed Time (Subplot 1)
    sns.boxplot(
        data=df,
        x='Tool_Name',
        y='Elapsed_Time_s',
        ax=axes[0],
        palette="viridis",
        hue='Tool_Name',
        legend=False
    )
    axes[0].set_title('Execution Time Distribution', fontsize=14)
    axes[0].set_ylabel('Elapsed Time (seconds)', fontsize=12)
    axes[0].set_xlabel('Tool', fontsize=12)

    # 4. Plot Memory Usage (Subplot 2)
    sns.boxplot(
        data=df,
        x='Tool_Name',
        y='Max_RSS_MB',
        ax=axes[1],
        palette="magma",
        hue='Tool_Name',
        legend=False
    )
    axes[1].set_title('Peak Memory Usage Distribution', fontsize=14)
    axes[1].set_ylabel('Max RSS (MB)', fontsize=12)
    axes[1].set_xlabel('Tool', fontsize=12)

    # 5. Final Adjustments
    plt.tight_layout()

    # Save the plot
    output_filename = 'benchmark_boxplots.png'
    plt.savefig(output_filename, dpi=300)
    print(f"Plot saved successfully as '{output_filename}'")

    # Show the plot (optional, requires a windowing environment)
    # plt.show()

if __name__ == "__main__":
    # Ensure the file name matches your uploaded file
    TSV_FILE = "benchmark_results_raw.tsv"
    plot_benchmark_results(TSV_FILE)