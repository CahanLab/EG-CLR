# ISSAAC-seq Info

Source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140203

Cell line: k562

Data type: paired single-cell RNA and single-cell ATAC seq (one is drop-seq based, another is FACS based)

Year: 2022

Format: FASTQ

Processed pipeline: snakemake file provided by ISSAAC-seq github

Software needed: STAR, Chromap (orginal paper), or 10X cellranger in local or cloud

Reference Genome: hg38 (orignally used), hg19 (used for ABC benchmarking)