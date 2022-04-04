library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(harmony)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)
set.seed(1234)

plan("multiprocess", workers = 4)
options(future.globals.maxSize = 100000 * 1024^2) # for 100 Gb RAM

setwd("~/path/to/output")

## From Run_IntegrateSamples_ATAC.sh
pbmc <- readRDS("PooledATAC.rds")

granges(pbmc)

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

# add the gene information to the object
Annotation(pbmc) <- annotations

# compute nucleosome signal score per cell
pbmc <- NucleosomeSignal(object = pbmc)

# compute TSS enrichment score per cell
pbmc <- TSSEnrichment(object = pbmc, fast = FALSE)

# add blacklist ratio and fraction of reads in peaks
#pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100
#pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments

pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)

pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)

saveRDS(pbmc, "./PooledSeurat_Final_NoHarmony.rds")

pbmc <- RunHarmony(
  object = pbmc,
  group.by.vars = 'sample',
  reduction = 'lsi',
  assay.use = 'ATAC',
  project.dim = FALSE
)

# re-compute the UMAP using corrected LSI embeddings
pbmc <- RunUMAP(pbmc, dims = 2:30, reduction = 'harmony', reduction.name = "harmony.umap")

saveRDS(pbmc, "./PooledSeurat_Final_WithHarmony.rds")

gene.activities <- GeneActivity(pbmc)

# add the gene activity matrix to the Seurat object as a new assay and normalize it
pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)
pbmc <- NormalizeData(
  object = pbmc,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(pbmc$nCount_RNA)
)

## Add Motif Data
# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
pbmc <- AddMotifs(
  object = pbmc,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

## Add Chromvar Data
pbmc <- RunChromVAR(
  object = pbmc,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

## Example list of motifs to run footprints on. You can add or remove based on requirement and preference.
Motifs <-c("FOSL1::JUNB","FOS::JUNB","FOSL2","FOSL1::JUND","FOSL2::JUN","FOSL2::JUND","FOSL1::JUN","FOS::JUND","FOSB::JUNB","BATF","FOS::JUN","FOSL2::JUNB","BATF3","JUN(var.2)","BATF::JUN","FOSL1","JUN::JUNB","JUND","JDP2","FOS","BACH1","JUNB","NFE2L1","JUN","JUND(var.2)","MAF::NFE2","NRF1","ATF3","ATF7","NFE2","BACH2","FOS::JUN(var.2)","FOSB::JUN","FOSL1::JUN(var.2)","FOSL2::JUN(var.2)","JUN::JUNB(var.2)","FOSB::JUNB(var.2)","FOSL2::JUNB(var.2)","JUNB(var.2)","FOSL1::JUND(var.2)","FOSL2::JUND(var.2)","ATF6","BACH2(var.2)","ATF2","CEBPA","CEBPB","CEBPD","KLF13","KLF14","KLF16","KLF10","KLF11","KLF15","KLF17","KLF2","IRF1","IRF2","IRF3","IRF4","STAT1","STAT3","STAT1::STAT2","MAFF","MAF","MAFA","SPI1","SPIB","SPIC","IRF5","IRF6","IRF7","IRF8","IRF9","NFKB1","NFKB2","CTCF","CTCFL","FOSL2", "FOS", "JUN", "BATF","MAF","CEBPA","KLF4","SPI1","GATA2","CEBPB","NR4A1","FLI1","HOXA10","ZBTB14","EGR3","KLF17","TFDP1","NRF1","NR5A1","KLF3","INSM1","KLF9","POU4F1","HINFP","TCF3","ZKSCAN5")

footprints = Motifs
footprints <- unique(footprints)
pbmc <- Footprint(
  object = pbmc,
  motif.name = footprints,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

## Save Final Object
saveRDS(pbmc, "./PooledSeurat_Final_WithGeneActivity.rds")


## For Differential Analysis
DefaultAssay(pbmc) <- "ATAC"
## Change Idents to whichever Celltype annotations are to be used. Marker Annotations are the broad annotations we selected from marker genes.
Idents(pbmc) <- pbmc$MarkerAnnotations

##Change ident.1 and ident.2 to different condition to test between. group.by Status does the analysis between these health conditions. Change the subset.ident to celltype. 
## Following example tests for Differential Genes between Severe EC and Healthy CD14 monocyte cells at 0.25 logFC threshold.
## Filter genes based on p-value cutoffs after running differential analysis. 
## We only reported differential peaks between a few celltypes and conditions within ATACseq data.
## For differential motif analysis, change DefaultAssay to chromvar and re run based on given conditions.
Mono_Sev_HD <- FindMarkers(pbmc, ident.1 = "Severe EC", ident.2 = "Healthy", group.by = "Status", subset.ident = "CD14 Monocytes", only.pos = TRUE, test.use = 'LR')


