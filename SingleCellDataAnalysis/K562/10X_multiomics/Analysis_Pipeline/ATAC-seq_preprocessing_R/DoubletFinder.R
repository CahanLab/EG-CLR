library(Signac)
library(Seurat)
library(SeuratDisk)
library(GenomicRanges)
library(GenomeInfoDb) 
library(ggplot2)
library(patchwork)
library(sceasy)
library(reticulate)
library(SingleCellExperiment)
library(DoubletFinder)    # exposes paramSweep()


# load adata, convert to rds, load rds to seurat
h5ad_file = "/Volumes/G-DRIVE mobile USB-C/Single-cell_data/K562/10x/ISSAACC-seq_generated/hg19_10xCloud_aligned_data/processed_data/rna-seq_data_4.23.25.h5ad"
sceasy::convertFormat(h5ad_file, from="anndata", to="seurat", outFile='/Volumes/G-DRIVE mobile USB-C/Single-cell_data/K562/10x/ISSAACC-seq_generated/hg19_10xCloud_aligned_data/processed_data/rna-seq_data_4.23.25.rds')
my_object <- readRDS("/Volumes/G-DRIVE mobile USB-C/Single-cell_data/K562/10x/ISSAACC-seq_generated/hg19_10xCloud_aligned_data/processed_data/rna-seq_data_4.23.25.rds")


# inspect the first few rows for cell obs
head(my_object@meta.data)          # base R
head(my_object[["leiden"]], 10)    # Seurat helper


# inspect X
first_cell_id <- colnames(my_object)[1]                 # cell 1’s name
expr_vec      <- GetAssayData(my_object, slot = "counts")[, first_cell_id]
expr_vec

# inspect the first few rows for peaks var
assay <- my_object[["RNA"]]      # set assay 
head(assay@meta.features)


# doublet finder
my_object <- NormalizeData(my_object)
my_object <- FindVariableFeatures(my_object, selection.method = "vst", nfeatures = 2000)
my_object <- ScaleData(my_object)
my_object <- RunPCA(my_object)
my_object <- RunUMAP(my_object, dims = 1:10)

head(Embeddings(my_object, reduction = "umap"))
head(Embeddings(my_object, reduction = "pca"))


sweep.res.my_object <- paramSweep(my_object, PCs = 1:10, sct = FALSE)
sweep.my_object <- summarizeSweep(sweep.res.my_object, GT = FALSE)
bcmvn_my_object <- find.pK(sweep.my_object)

# Homotypic Doublet Proportion Estimate
annotations <- my_object@meta.data$leiden
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(my_object@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


# Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
my_object <- doubletFinder(my_object, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = NULL, sct = FALSE)
head(my_object@meta.data) 


group_col <- "log1p_total_counts"

my_object[[group_col]] <- factor(my_object[[group_col]][,1],
                                 levels = c("Singlet", "Doublet"))

DimPlot(
  my_object,
  reduction = "umap",
  group.by  = group_col,
  cols      = c("steelblue", "firebrick"),   # Singlet first, Doublet second
  pt.size   = 0.6                            # tweak point size to taste
) + ggtitle("DoubletFinder calls on UMAP")


FeaturePlot(
  my_object,
  features  = "total_counts",      # the metadata column we just created / confirmed
  reduction = "umap",
  pt.size   = 0.1,           # tweak point size
  cols     = c("yellow", "purple")          # keep original cell order
) + ggtitle("UMAP coloured by log1p(total counts)")
