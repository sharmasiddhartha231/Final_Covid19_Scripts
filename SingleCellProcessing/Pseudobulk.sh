#!/bin/sh

## To create pseudobulk profiles run the following command. The pseudobulk tool has been provided in the Software packet.
## Run this for each sample.
## The bam file is obtained from the cellranger analysis.
## The Singlecell csv file (AllBarcodes_OutliersandDoubletRemovedTwice.csv) can be obtained from running the Edit.sh script.
## The Markers.txt file contains the cluster numbers for each cell. This splits the file into different bam files labelling them as clusterbam_0.bam, clusterbam_1.bam and so on (Given number of clusters)
java -jar snATACClusteringTools.jar splitbams ~/path/to/atac_possorted_bam.bam ~/path/to/AllBarcodes_OutliersandDoubletRemovedTwice.csv ~/path/to/Markers.txt ~/path/to/output

## For each pseudobulked output bam file, run the following steps.
samtools sort -o Sorted_clusterbam_{cluster_number}.bam clusterbam_{cluster_number}.bam
samtools index Sorted_clusterbam_{cluster_number}.bam
macs2 callpeak -t Sorted_clusterbam_{cluster_number}.bam -f BAMPE --outdir ./ -n Sorted_clusterbam_{cluster_number}_bampe_peaks -g 'hs'
