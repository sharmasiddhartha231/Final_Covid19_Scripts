#!/bin/sh

## THIS SCRIPT HELPS REMOVE FILTERED OUT CELLS AND RUN BOTH AMULET AND SCRUBLET
## THIS IS DESIGNED FOR EACH MULTIOME SAMPLE. THIS SCRIPT CAN BE MODIFIED FOR THE snATAC SAMPLES AS WELL.
cd ~/path/to/individual_sample_directory

## FILES YOU NEED IN THIS DIRECTORY
## Barcodes_NoQC.txt, Barcodes_QC.txt - obtained from the Multiome_WNN.R/Signac_snATAC_Processing.R
## per_barcode_metrics.csv from Cellranger-arc processing. This file needs to be edited so that column 1,9,10 are the following: barcode, is_cell, is_cell_barcode respectively for the version of Amulet used in the filtering to Run. 
# Multiome per_barcode_metrics.csv are not provided in this format. They need to reformatted as such. the singlecell.csv file from snATAC samples are formatted accordingly so they do not need to be remodified.

## filtered_feature_bc_matrix directory from Cellranger-arc processing. 

## Steps to remove all outliers.
## Extract cells failing QC initial QC filters.
comm -23 Barcodes_NoQC.txt Barcodes_QC.txt > Differences

# Reformat per_barcode_metrics.csv for Amulet analysis. 
awk -F, '{ t = $4; $4 = $10; $10 = t; print; } ' OFS=, per_barcode_metrics.csv > AllBarcodes.csv
cp AllBarcodes.csv Temp.csv
tail -n+2 Temp.csv > Body.csv
head -1 Temp.csv > Header.csv
awk -F"," 'BEGIN { OFS = "," } {$9="cell_id"; print}' Header.csv > Headerout.csv
awk -F"," 'BEGIN { OFS = "," } {$9="None"; print}' Body.csv > Bodyout.csv
cat Headerout.csv Bodyout.csv > OutTemp.csv
rm Headerout.csv Header.csv Body.csv Bodyout.csv Temp.csv
tail -n+2 OutTemp.csv > A.csv
head -1 OutTemp.csv > Header
sort -t, -k10,10 -n -r A.csv > B.csv
awk -F, '{ if ($10 == "1") $9="_cell_"(NR-1)} 1' OFS=, B.csv > C.csv
sort -t, -r C.csv > D.csv
cat Header D.csv > E.csv
rm OutTemp.csv A.csv B.csv C.csv D.csv AllBarcodes.csv Header
mv E.csv AllBarcodes.csv
cut -f1 Differences > Outliers.txt
sed -i -e 's/"//g' Outliers.txt
awk -F "," 'FNR==NR{a[$1];next};!($1 in a)' Outliers.txt AllBarcodes.csv > AllBarcodes_OutliersRemoved.csv

## AllBarcodes_OutliersRemoved.csv is the csv file used in Amulet (snATAC Doublet Detection step). This file has the initial cells filtered out by first QC removed from the dataset. Thus Amulet is run on the filtered dataset.
## Functioning of Amulet explained in Amulet.sh provided in the directory.
java -jar ~/path/to/ATAC-DoubletDetector-main/snATACOverlapCounter.jar ~/path/to/atac_possorted_bam.bam ~/path/to/singlecell.csv ~/path/to/human_autosomes.txt ~/path/to/output_directory
python3 ~/path/to/ATAC-DoubletDetector-main/ATACDoubletDetector.py --rfilter ~/path/to/ATAC-DoubletDetector-main/blacklist_repeats_segdups_rmsk_hg38.txt ~/path/to/overlaps.txt ~/path/to/Overlap_Summary.txt ~/path/to/output_directory

## Extract Doublet Barcodes detected by Amulet and remove them from singlecell.csv file to create new csv file.
awk -F "," 'FNR==NR{a[$1];next};!($1 in a)' DoubletBarcodes_01.txt AllBarcodes_OutliersRemoved.csv > AllBarcodes_OutliersandDoubletRemoved.csv

## For scrublet, extract the files from filtered_feature_bc_matrix and get the Gene Expression data in Genes.tsv for scrublet to run.
gunzip ./filtered_feature_bc_matrix/*
cp ./filtered_feature_bc_matrix/matrix.mtx .
cp ./filtered_feature_bc_matrix/features.tsv .
cp ./filtered_feature_bc_matrix/barcodes.tsv .
awk -F'\t' '{print>$3}' features.tsv
mv Gene\ Expression Genes.tsv
mv Peaks Peaks.tsv  

## Run Scrublet.
python3 /projects/ucar-lab/Siddhartha/Covid_19/snATACseq/MultiOme/2nd_multiome/MultiomeRuns/scrublet_script.py 

## Remove Doublets detected from Scrublet from barcodes.tsv from filtered_feature_bc_matrix directory. 
paste barcodes.tsv predicted_doublet_mask.txt > Temp
awk -F'\t' '$2 == "True"' Temp > Temp1
cut -f1 Temp1 > Outliers_Scrublet.txt
rm Temp Temp1

## Extract the 2nd set of doublets from singlecell.csv file to obtain final csv file. 
awk -F "," 'FNR==NR{a[$1];next};!($1 in a)' Outliers_Scrublet.txt AllBarcodes_OutliersandDoubletRemoved.csv > AllBarcodes_OutliersandDoubletRemovedTwice.csv	

## Keep DoubletBarcodes_01.txt and Outliers_Scrublet.txt for further analysis.
## You can skip the reformatting and scrublet steps for snATAC samples. Please check singlecell files in case they are formatted differently in which case you will need to reformat them.
