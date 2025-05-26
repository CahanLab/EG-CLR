# Share-seq Info

Source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140203

Cell line: k562

Data type: paired single-cell RNA and single-cell ATAC seq

Year: 2020

Format: Raw counts (the file is already a raw counts and fastq not provided because it is processed by Share-seq paper)

Processed pipeline: snakemake file provided by ISSAAC-seq github (download and convert to barcode, feature, matrix files for RNA and ATAC)

Software needed: STAR (only to do some QC on raw counts)

Reference Genome: hg19 (used in here fore benchmarking purpose)