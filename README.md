# Enhancer-Gene Inference via CLR (EG-CLR)

# Workflow
<img width="655" alt="Screenshot 2025-03-18 at 09 45 40" src="https://github.com/user-attachments/assets/64f18a31-96ba-4682-bb14-92d96d9b8a3b" />

# Environment:
* Python >= 3.10


# Data input requirement: 
* Paired single-cell multiomics data (in count matrix format):
  * Single-cell RNA-seq 
  * Single-cell ATAC-seq
* Peak annotation file (tsv format) for scATAC-seq (indicating promoter, enhancer, gene regions): 
  * 10X ARC provided peak_annotation.tsv with the following rule:
  	1. If a peak overlaps with promoter region (-1000bp, +100bp) of any TSS, it is annotated as a promoter peak of the gene.
	  2. If a peak is within 200kb of the closest TSS, and if it is not a promoter peak of the gene of the closest TSS, it will be annotated as a distal peak of that gene.
	  3. If a peak overlaps the body of a transcript, and it is not a promoter nor a distal peak of the gene, it will be annotated as a distal peak of that gene with distance set as zero.
	  4. If a peak has not been mapped to any gene at the step, it will be annotated as an intergenic peak without a gene symbol assigned.
  * Any customized peak annotation in tsv format will be accepted. 
