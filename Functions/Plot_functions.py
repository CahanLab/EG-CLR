import matplotlib.pyplot as plt
import random
import seaborn as sns
import numpy as np
random.seed(10)

### plot

# plot peaks per gene 
def plot_peaks_per_gene(Peak_Gene_df, peak_col = "Peak", gene_col = "Gene"):
    
    # (a) count unique genes per element
    genes_per_element = (
        Peak_Gene_df.groupby(peak_col)[gene_col].nunique()                    
    )

    # (b) turn that into a distribution:  x = k genes,  y = #elements with k genes
    dist_elems = (
        genes_per_element.value_counts()   # how many elements have k genes
                    .sort_index()        # sort so bars appear left‑to‑right
    )

    # (c) plot
    plt.figure()
    dist_elems.plot(kind='bar')
    plt.xlabel('Number of genes linked to a peak')
    plt.ylabel('Number of peaks')
    plt.title('Peaks grouped by their gene count')
    plt.tight_layout()
    plt.show()

    print("Average peak for each gene", genes_per_element.mean())

# plot genes per peak
def plot_genes_per_peak(Peak_Gene_df, peak_col = "Peak", gene_col = "Gene"):
    
    # (a) count unique genes per element
    elems_per_gene = (
        Peak_Gene_df.groupby(gene_col)[peak_col].nunique()                  
    )

    # (b) distribution:  x = k elements,  y = #genes with k elements
    dist_genes = elems_per_gene.value_counts().sort_index()

    # (c) plot
    plt.figure()
    dist_genes.plot(kind='bar')
    plt.xlabel('Number of peaks linked to a gene')
    plt.ylabel('Number of genes')
    plt.title('Genes grouped by their peak count')
    plt.tight_layout()
    plt.show()

    print("Average gene for each peak", elems_per_gene.mean())
    
# Plot two values in scatterplot with hue and optional line of best fit
def plot_values(prediction_df, x_col = "distance_to_tss", y_col="Values", hue_col="Significant", palette_discrete = True,  line_fit=True):
    
    # decide palette
    if palette_discrete:
        palette = {False: "grey", True: "crimson"}
    else:
        palette = "viridis"
    
    # graph 
    sns.scatterplot(
        data=prediction_df,
        x=x_col,
        y=y_col,
        hue=hue_col,
        palette=palette,
        alpha=0.75
    )

    # overall line of best‑fit
    if line_fit:
        x = prediction_df[x_col].to_numpy()
        y = prediction_df[y_col].to_numpy()
        slope, intercept = np.polyfit(x, y, 1)
        x_fit = np.linspace(x.min(), x.max(), 200)
        y_fit = slope * x_fit + intercept
        plt.plot(x_fit, y_fit, "--", color="black", label=f"Best fit: y = {slope:.3f}x + {intercept:.3f}")

    # ucosmetics
    plt.xlabel(x_col)
    plt.ylabel(y_col)
    plt.title(f"{y_col} vs. {x_col}")
    plt.legend(title=hue_col)
    plt.tight_layout()
    plt.show()
    
# plot stacked bar chart 
def plot_stack_bar(prediction_df, bin = 10, value_col = "distance_to_tss", hue_col = "Significant"):

    sns.histplot(
        data=prediction_df,
        x=value_col,
        hue=hue_col,        #  
        bins= bin,  
        palette={False: "grey",True: "crimson"},
        #stat="probability",      # turn counts into fractions
        #multiple="fill",         # stack & normalise to 1
        #element="step"           # or "bars" for solid bars
    )

    plt.xlabel(value_col)
    plt.ylabel("Percentage")
    plt.title(f"Stacked bar chart of {value_col} by {hue_col}")
    plt.show()
    
    
    