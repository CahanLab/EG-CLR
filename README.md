# Enhancer-Gene Inference via CLR (EG-CLR)

# Workflow

<img width="512" alt="image" src="https://github.com/user-attachments/assets/7275480f-33ed-4827-ab8f-b4d208c1ba7c" />

# Environment:

```python

# Python Version
Python >= 3.10

# Essential Packages
anndata >= 0.11.4
numpy >= 2.1.3
pandas >= 2.2.3
scanpy >= 1.11.1
muon >= 0.1.7
matplotlib >= 3.10.1
scipy >= 1.15.2
joblib >= 1.4.2
pybedtools >= 0.12.0

```

Full package list is in environment.yml

# Data Input Requirement:

* Paired single-cell multiomics data (read count matrix in adata):
   * Single-cell RNA-seq
   * Single-cell ATAC-seq

* Peak annotation file (tsv format) for scATAC-seq (indicating promoter, enhancer, gene regions):
   * e.g. 10X ARC provided peak_annotation.tsv in reference 
      1. If a peak overlaps with promoter region (-1000bp, +100bp) of any TSS, it is annotated as a promoter peak of the gene.
      2. If a peak is within 200kb of the closest TSS, and if it is not a promoter peak of the gene of the closest TSS, it will be annotated as a distal peak of that gene.
      3. If a peak overlaps the body of a transcript, and it is not a promoter nor a distal peak of the gene, it will be annotated as a distal peak of that gene with distance set as zero.
      4. If a peak has not been mapped to any gene at the step, it will be annotated as an intergenic peak without a gene symbol assigned.
   * Any customized peak annotation in tsv format will be accepted.

# Output
