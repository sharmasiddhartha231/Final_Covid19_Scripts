library(Seurat)
library(cowplot)
library(harmony)
library(ggplot2)
library(sctransform)
library(dplyr)
library(tidyverse)
library(SeuratDisk)
library(patchwork)

setwd("~/path/to/output")
getwd()

pbmc <- readRDS("./PooledSeuratRNA.rds")
pbmc

DefaultAssay(pbmc) <- "RNA"
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)


# store mitochondrial percentage in object meta data
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
# run sctransform
pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)
pbmc <- RunPCA(pbmc, verbose = FALSE)
pbmc <- RunUMAP(pbmc, dims = 1:30, verbose = FALSE)
pbmc <- FindNeighbors(pbmc, dims = 1:30, verbose = FALSE)
pbmc <- FindClusters(pbmc, verbose = FALSE)
p1 <- DimPlot(object = pbmc, reduction = "pca", pt.size = .1, group.by = "sample")
p2 <- DimPlot(object = pbmc, reduction = "umap", pt.size = .1, group.by = "sample")
pbmc <- RunHarmony(pbmc, "sample", plot_convergence = TRUE, assay.use = "SCT")
p3 <- DimPlot(object = pbmc, reduction = "harmony", pt.size = .1, group.by = "sample")
pbmc <- RunUMAP(pbmc, reduction = "harmony", dims = 1:20, reduction.name = "harmony.umap")
p4 <- DimPlot(object = pbmc, reduction = "harmony.umap", pt.size = .1, group.by = "sample")
pbmc <- FindNeighbors(pbmc, reduction = "harmony", dims = 1:20)
pbmc <- FindClusters(pbmc, verbose = FALSE)

## Initial UMAP plots
pdf("./UMAP_Plots.pdf", height = 15, width = 15)
p1 + p3
p2 + p4
dev.off()

saveRDS(pbmc, "./SeuratObject_SCT_NoAnnotations.rds")


## Load reference dataset for first set of annotations.
# Load the pre-processed scRNA-seq data for PBMCs
# Path to file: https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat
# The above object was used to annotate the snATAC samples. Please download to directory to use for annotation.

reference <- LoadH5Seurat("./pbmc_multimodal.h5seurat")
anchors <- FindTransferAnchors(
  reference = reference,
  query = pbmc,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50
)

pbmc <- MapQuery(
  anchorset = anchors,
  query = pbmc,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca",
  reduction.model = "wnn.umap"
)

p5 = DimPlot(pbmc, reduction = "harmony.umap", group.by = "predicted.celltype.l1", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
p6 = DimPlot(pbmc, reduction = "harmony.umap", group.by = "predicted.celltype.l2", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend()
p7 = DimPlot(pbmc, reduction = "umap", group.by = "predicted.celltype.l1", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
p8 = DimPlot(pbmc, reduction = "umap", group.by = "predicted.celltype.l2", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend()

## Save Final SCT annotated object
saveRDS(pbmc, "./SeuratObject_FinalQC_Processed.rds")

## Final UMAP Plots
pdf("./UMAP_Plots_WithAnnotations.pdf", height = 15, width = 15)
p5
p6
p7
p8
dev.off()

## For Differential Analysis
DefaultAssay(pbmc) <- "RNA"
## Change Idents to whichever Celltype annotations are to be used. Marker Annotations are the broad annotations we selected from marker genes.
Idents(pbmc) <- pbmc$MarkerAnnotations
table(Idents(pbmc))

##Change ident.1 and ident.2 to different condition to test between. group.by Status does the analysis between these health conditions. Change the subset.ident to celltype. 
## Following example tests for Differential Genes between Severe EC and Healthy CD14 monocyte cells at 0.25 logFC threshold.
## Filter genes based on p-value cutoffs after running differential analysis.
Mono_Sev_HD <- FindMarkers(pbmc, ident.1 = "Severe EC", ident.2 = "Healthy", group.by = "Status", subset.ident = "CD14 Monocytes", logfc.threshold = 0.25, min.pct = 0.1, assay = "RNA")

