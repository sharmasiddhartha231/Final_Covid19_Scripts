#!/bin/sh

## Updated version of Amulet can be accessed here: https://github.com/UcarLab/AMULET
## For the analysis for snATAC and Multiome ATAC samples, the following steps were followed.

## For this step, extract the cell barcodes that were filtered out in the snATAC processing from the Signac (for snATAC) and Seurat WNN (for Multiome) and extract those cells from the singlecell csv file.
## Human autosomes can be found on the Github website.
## STEP 1
java -jar ~/path/to/ATAC-DoubletDetector-main/snATACOverlapCounter.jar ~/path/to/atac_possorted_bam.bam ~/path/to/singlecell.csv ~/path/to/human_autosomes.txt ~/path/to/output_directory

## From STEP 1, you will obtain a few output files. Use the overlaps.txt and Overlap_Summary.txt file in the next step to obtain the doublets.
## STEP 2
python3 ~/path/to/ATAC-DoubletDetector-main/ATACDoubletDetector.py --rfilter ~/path/to/ATAC-DoubletDetector-main/blacklist_repeats_segdups_rmsk_hg38.txt ~/path/to/overlaps.txt ~/path/to/Overlap_Summary.txt ~/path/to/output_directory

## STEP 2 provides a list of Doublet Barcodes file (DoubletBarcodes_01.txt/MultipletBarcodes_01.txt). 
## Use this text file to filter out cells from the Signac/Seurat object before proceeding with the remaining processing steps.
