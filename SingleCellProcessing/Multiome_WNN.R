library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)


## Please run these steps for every individual sample obtained from the Multiome Data for initial filtering and annotation.

inputdata.10x <- Read10X_h5("~/path to filtered_feature_bc_matrix.h5 file")

# extract RNA and ATAC data
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

# Create Seurat object
pbmc <- CreateSeuratObject(counts = rna_counts)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

frag.file <- "~/path to atac_fragments.tsv.gz file"
chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)
pbmc[["ATAC"]] <- chrom_assay

## To save Initial QC plots for Multiome Data
pdf("QC_Plots.pdf", height = 12, width = 12)
VlnPlot(pbmc, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
        log = TRUE, pt.size = 0) + NoLegend()
dev.off()

## Save Signac object with all cells. 
saveRDS(pbmc, file = "SeuratObject_No_QC.rds")

## Save all barcodes for each sample.
Barcodes_NoQC <- pbmc@active.ident
write.table(Barcodes_NoQC, file = "./Barcodes_NoQC.txt", sep = "\t", row.names = T)

##IMPORTANT. PLEASE SEE
## These are parameters provided from Seurat. Please refer to Supplementary tables for individual cutoffs. The cutoffs are different for each sample.

pbmc <- subset(
  x = pbmc,
  subset = nCount_ATAC < 7e4 &
    nCount_ATAC > 5e3 &
    nCount_RNA < 25000 &
    nCount_RNA > 1000 &
    percent.mt < 20
)

pbmc

## Save Signac object with initially filtered cells. 
saveRDS(pbmc, file = "SeuratObject_QC.rds")

## Save initially filtered barcodes for each sample.
Barcodes_QC <- pbmc@active.ident
write.table(Barcodes_QC, file = "./Barcodes_QC.txt", sep = "\t", row.names = T)

## IMPORTANT
## BEFORE RUNNING THE NEXT STEPS, PLEASE REMOVE THE CELLS IDENTIFIED BY THE snATAC DOUBLET DETECTOR TOOL (AMULET) AND SCRUBLET (scRNA DOUBLET DETECTION TOOL)
## WE RAN AMULET ON THE ATAC CELLS POST INITIAL FILTERING FROM SIGNAC/SEURAT PARAMETERS. WE RAN SCRUBLET ON INITIAL FILES OBTAINED FROM CELLRANGER PROCESSING.
## PLEASE REFER TO THE AMULET AND SCRUBLET SCRIPTS FOR FURTHER DETAILS.

## Add Amulet Doublets
Doublets <- read.table(file = "./DoubletBarcodes_01.txt", sep = "\t", header = F)
Doublets$V2 <- 1
colnames(Doublets) <- c("Barcode", "Doublet")

## Add Scrublet Doublets
Scrublets <- read.table(file = "./Outliers_Scrublet.txt", sep = "\t", header = F)
Scrublets$V2 <- 1
colnames(Scrublets) <- c("Barcode", "Doublet")

#Combine Doublets
AllDoublets <- rbind(Scrublets, Doublets)
AllDoublets <- unique(AllDoublets)

pbmc$Barcode <- rownames(pbmc@meta.data)
AllBarcodes <- pbmc$Barcode
AllBarcodes <- as.data.frame(AllBarcodes)
colnames(AllBarcodes) <- "Barcode"

Trial <- merge(AllBarcodes, AllDoublets, by = "Barcode", all.x = TRUE)
Trial[is.na(Trial)] <- 0
pbmc$Doublet <- Trial$Doublet

## Removes All Doublets
pbmc <- subset(pbmc, subset = Doublet == 0)
pbmc

# RNA analysis
DefaultAssay(pbmc) <- "RNA"
pbmc <- SCTransform(pbmc, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(pbmc) <- "ATAC"
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)
pbmc <- RunUMAP(pbmc, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

pbmc <- FindMultiModalNeighbors(pbmc, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
pbmc <- RunUMAP(pbmc, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
pbmc <- FindClusters(pbmc, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

# Load the pre-processed scRNA-seq data for PBMCs
# Path to file: https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat
# The above object was used to annotate the snATAC samples. Please download to directory to use for annotation.
reference <- LoadH5Seurat("~/path to pbmc_multimodal.h5seurat object")

##Save UMAP plots with annotations based on scRNA object.
pdf("UMAP_Plots.pdf", width = 15, height = 12)
DimPlot(object = reference, reduction = "wnn.umap", group.by = "celltype.l1", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
DimPlot(object = reference, reduction = "wnn.umap", group.by = "celltype.l2", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
dev.off()

pbmc <- SCTransform(pbmc, verbose = FALSE)

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
## Save Seurat object. Keep for further use.
saveRDS(pbmc, file = "./SeuratObject_QC_Annotated.rds")

## Save final set of barcodes for each sample.
Barcodes <- pbmc@active.ident
write.table(Barcodes, file = "./Barcodes.txt", sep = "\t", row.names = T)

