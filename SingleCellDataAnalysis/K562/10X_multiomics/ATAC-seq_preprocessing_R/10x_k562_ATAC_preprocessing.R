library(Signac)
library(Seurat)
library(GenomicRanges)
library(GenomeInfoDb) 
library(ggplot2)
library(patchwork)

# This folder should contain barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz
matrix_dir <- "/Volumes/G-DRIVE mobile USB-C/Single-cell_data/K562/10x/ISSAACC-seq_generated/hg19_10xCloud_aligned_data/filtered_feature_bc_matrix"
fragments_path <- "/Volumes/G-DRIVE mobile USB-C/Single-cell_data/K562/10x/ISSAACC-seq_generated/hg19_10xCloud_aligned_data/atac_fragments.tsv.gz"


# Read 10x data
multiome_data <- Read10X(data.dir = matrix_dir)
dim(multiome_data[['Gene Expression']])
dim(multiome_data[['Peaks']])

# Create Seurat object by RNA
seurat_obj <- CreateSeuratObject(counts = multiome_data$`Gene Expression`, assay = "RNA")
head(seurat_obj$nFeature_RNA)

# Add fragments
chrom_assay <- CreateChromatinAssay(
  counts = multiome_data$Peaks,
  genome = 'hg19',
  sep = c(":", "-"),
  fragments = fragments_path,
  min.cells = 1
)

# Add ATAC to Seurat 
seurat_obj[["ATAC"]] <- chrom_assay

# Change the mode to ATAC
DefaultAssay(seurat_obj) <- "ATAC"


seurat_obj[["ATAC"]] 
granges(seurat_obj)


# keep standard chromosome
peaks.keep <- seqnames(granges(seurat_obj)) %in% standardChromosomes(granges(seurat_obj))
seurat_obj <- seurat_obj[as.vector(peaks.keep), ]

'''
# Basically no possible because hg19 is too old to have annotation file

# Add annotation
library(AnnotationHub)
ah <- AnnotationHub()

# Search for the Ensembl 98 EnsDb for Homo sapiens on AnnotationHub
query(ah, "EnsDb.Hsapiens.v75")
query(ah, c("EnsDb", "Hsapiens", "GRCh37")) 
'''


# compute nucleosome signal score per cell
seurat_obj <- NucleosomeSignal(object = seurat_obj)

# compute TSS enrichment score per cell
seurat_obj <- TSSEnrichment(object = seurat_obj)

# add fraction of reads in peaks
seurat_obj$pct_reads_in_peaks <- seurat_obj$peak_region_fragments / seurat_obj$passed_filters * 100

# add blacklist ratio
seurat_obj$blacklist_ratio <- FractionCountsInRegion(
  object = seurat_obj, 
  assay = 'peaks',
  regions = blacklist_hg19
)

