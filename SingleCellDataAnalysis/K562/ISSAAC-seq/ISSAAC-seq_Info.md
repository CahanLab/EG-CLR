# ISSAAC-seq Info

Source: https://github.com/dbrg77/ISSAAC-seq/tree/main

Cell line: k562

Data type: paired single-cell RNA and single-cell ATAC seq

Year: 2022

Format: FASTQ 

Processed pipeline: snakemake file provided by ISSAAC-seq github (download and convert to barcode, feature, matrix files for RNA and ATAC)

Software needed: STAR, Chromap (because it is not 10x, no cellranger)

Reference Genome: hg38 (in orginal paper), hg19 (used in here fore benchmarking purpose)