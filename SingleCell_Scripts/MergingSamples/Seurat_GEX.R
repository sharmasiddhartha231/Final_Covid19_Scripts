library(Seurat)
library(Signac)
library(dplyr)
library(tidyverse)

## Run this script to create the GEX objects for each Multiome Sample. Then run Merge_Seurat_RNA.R to merge all of them together. 
## Run this in individual directory of each sample. 
## You'll need in the working directory the barcodes.tsv.gz, features.tsv.gz and the matrix.mtx.gz to create the initial Seurat object.
## From the Barcodes.txt file saved while individual processing, extract the barcodes and remove it from individual samples to remove low QC cells.

pbmc.data <- Read10X(data.dir = ".")
pbmc.data <- pbmc.data$`Gene Expression`
pbmc <- CreateSeuratObject(counts = pbmc.data)
pbmc

## Save Initial object.
saveRDS(pbmc, file = "SeuratObject.rds")

## Read in barcodes from Barcodes.txt file and use that list to filter out barcodes that were removed from individual QC
barcodes <- read.table(file = "./Barcodes.txt", sep = "\t", header = F)
barcodes$V2 <- 1
colnames(barcodes) <- c("Barcode","Keep")
pbmc$barcodes <- rownames(pbmc@meta.data)
AllBarcodes <- pbmc$barcodes
AllBarcodes <- as.data.frame(AllBarcodes)
colnames(AllBarcodes) <- "Barcode"
Trial <- merge(AllBarcodes, barcodes, by = "Barcode", all.x = TRUE)
Trial[is.na(Trial)] <- 0
pbmc$Final <- Trial$Keep
pbmc <- subset(pbmc, subset = Final == 1)
pbmc

## Save Final GEX object.
saveRDS(pbmc, file = "Sample_Name.rds")

## Run this script for each Multiome sample used to merge data together. Then run Merge_Seurat_RNA.R to create final merged object.