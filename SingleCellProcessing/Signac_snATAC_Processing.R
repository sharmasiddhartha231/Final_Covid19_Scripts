library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)
set.seed(1234)


## Please run the following steps for each snATAC sample file to individually process and annotate them using the reference dataset provided by Seurat.
counts <- Read10X_h5(filename = "~/path to filtered_peak_bc_matrix.h5 file")
metadata <- read.csv(
  file = "~/path to singlecell.csv file",
  header = TRUE,
  row.names = 1
)

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = '~ path to fragments.tsv.gz file',
  min.cells = 10,
  min.features = 200
)

pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)

pbmc
pbmc[['peaks']]
granges(pbmc)

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to UCSC style since the data was mapped to hg38
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

# add the gene information to the object
Annotation(pbmc) <- annotations

# compute nucleosome signal score per cell
pbmc <- NucleosomeSignal(object = pbmc)

# compute TSS enrichment score per cell
pbmc <- TSSEnrichment(object = pbmc, fast = FALSE)

# add blacklist ratio and fraction of reads in peaks
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100
pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments

pbmc$high.tss <- ifelse(pbmc$TSS.enrichment > 2, 'High', 'Low')
p1 <- TSSPlot(pbmc, group.by = 'high.tss') + NoLegend()

pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
p2 <- FragmentHistogram(object = pbmc, group.by = 'nucleosome_group')

p3 <- VlnPlot(
  object = pbmc,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)

## To save Initial QC plots for snATAC
pdf(file = "./QC_Plots.pdf", width = 15, height = 12)
p1
p2
p3
dev.off()

## Save Signac object with all cells. 
saveRDS(pbmc, file = "SeuratObject_No_QC.rds")

## Save all barcodes for each sample.
Barcodes_NoQC <- pbmc@active.ident
write.table(Barcodes_NoQC, file = "./Barcodes_NoQC.txt", sep = "\t", row.names = T)


##IMPORTANT. PLEASE SEE
## These are parameters provided from Signac. Please refer to Supplementary tables for individual cutoffs. The cutoffs are different for each sample.
pbmc <- subset(
  x = pbmc,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 20000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)
pbmc

## Save Signac object with initially filtered cells. 
saveRDS(pbmc, file = "SeuratObject_QC.rds")

## Save initially filtered barcodes for each sample.
Barcodes_QC <- pbmc@active.ident
write.table(Barcodes_QC, file = "./Barcodes_QC.txt", sep = "\t", row.names = T)

## IMPORTANT
## BEFORE RUNNING THE NEXT STEPS, PLEASE REMOVE THE CELLS IDENTIFIED BY THE snATAC DOUBLET DETECTOR TOOL (AMULET).
## WE RAN AMULET ON THE ATAC CELLS POST INITIAL FILTERING FROM SIGNAC/SEURAT PARAMETERS. 
## PLEASE REFER TO THE AMULET SCRIPT FOR FURTHER DETAILS.

## Add Amulet Doublets
Doublets <- read.table(file = "./DoubletBarcodes_01.txt", sep = "\t", header = F)
Doublets$V2 <- 1
colnames(Doublets) <- c("Barcode", "Doublet")

pbmc$Barcode <- rownames(pbmc@meta.data)
AllBarcodes <- pbmc$Barcode
AllBarcodes <- as.data.frame(AllBarcodes)
colnames(AllBarcodes) <- "Barcode"

Trial <- merge(AllBarcodes, Doublets, by = "Barcode", all.x = TRUE)
Trial[is.na(Trial)] <- 0
pbmc$Doublet <- Trial$Doublet

## Removes All Doublets
pbmc <- subset(pbmc, subset = Doublet == 0)
pbmc

pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)

pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)

gene.activities <- GeneActivity(pbmc)

# add the gene activity matrix to the Seurat object as a new assay and normalize it
pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)
pbmc <- NormalizeData(
  object = pbmc,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(pbmc$nCount_RNA)
)

DefaultAssay(pbmc) <- 'RNA'

# Load the pre-processed scRNA-seq data for PBMCs
# Path to file: https://signac-objects.s3.amazonaws.com/pbmc_10k_v3.rds
# The above object was used to annotate the snATAC samples. Please download to directory to use for annotation.
pbmc_rna <- readRDS("~/path to scRNA Object")

transfer.anchors <- FindTransferAnchors(
  reference = pbmc_rna,
  query = pbmc,
  reduction = 'cca'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = pbmc_rna$celltype,
  weight.reduction = pbmc[['lsi']],
  dims = 2:30
)

pbmc <- AddMetaData(object = pbmc, metadata = predicted.labels)

p2 <- DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("ATAC")
p5 <- DimPlot(pbmc, reduction = "umap", group.by = "predicted.id", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("ATAC") + NoLegend()

##Save UMAP plots with annotations based on scRNA object.
pdf("UMAP_Plots.pdf", width = 15, height = 12)
p2
p5
dev.off()

## Save Signac object. Keep for further use.
saveRDS(pbmc, file = "SeuratObject_QC_Annot.rds")


