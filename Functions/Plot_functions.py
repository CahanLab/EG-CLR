import matplotlib.pyplot as plt
import random
random.seed(10)

### plot

def plot_peaks_per_gene(CLR_peaks):
    
    # Calculate the Number that each gene has non-zero values across all peaks (rows)
    non_zero_Number = (CLR_peaks != 0).sum(axis=0)

    # Create a plot
    plt.figure(figsize=(10, 6))
    plt.bar(non_zero_Number.index, non_zero_Number.values)
    
    # Set the labels and title
    plt.xlabel('Genes')
    plt.ylabel('Number of associated peaks')
    plt.title('Number of associated peaks for each gene')
    plt.xticks(rotation=90)
    
    # Show the plot
    plt.tight_layout()
    plt.show()
    
def plot_genes_per_peak(CLR_peaks):
    
     # Calculate the frequency that each gene has non-zero values across all peaks (rows)
    non_zero_Number = (CLR_peaks != 0).sum(axis=1)

    # Create a plot
    plt.figure(figsize=(10, 6))
    plt.bar(CLR_peaks.index, non_zero_Number.values)

    # Set the labels and title
    plt.xlabel('Peaks')
    plt.ylabel('Number of associated genes')
    plt.title('Number of associated genes for each peaks')

    # Show every 100th x-axis label
    plt.xticks(ticks=range(0, len(CLR_peaks.index), 100), labels=CLR_peaks.index[::100], rotation=90)

    # Show the plot
    plt.tight_layout()
    plt.show()