# Covid19_Scripts

The following scripts were used for the processing of samples used in the analysis for the following publication:
**Epigenetic Memory of COVID-19 in Innate Immune Cells and Their Progenitors**


The scripts have been split into two major directories: Bulk Data (**BulkData_Scripts**) and Single Cell Data (**SingleCell_Scripts**)

Within **BulkData_Scripts**, we have the following:

**BulkDataProcessing**

Above subdirectory has **BulkATAC_Processing.sh** which is used to process all bulk ATAC samples. The script is run using **sbatch_bulk.sh** for multiple samples together. **Subread_Frip.sh** is used to calculate FRiP scores which is used as a Quality check. 

**DifferentialAnalysis_Bulk**

The above subdirectory contains scripts for runnning differential analysis. **dbind.R** is used to calculate the consensus matrix using all bulk samples and is run using **dbind.sh**. **dbind.R** utilizes a csv file to access all samples for which the consensus needs to be calculated. An example file (**CD14_hg38.csv**) has been provided as an example. **DifferentialAnalysis_cinaR.R** is used to run the Differential Analysis step using the consensus matrix from the previous step. The **PCA.R** and **TimeSeries.R** scripts build a PCA plot and conduct trend analysis and visualizations from the differential data.

**Enrichment**

The scripts in the above subdirectory conduct Hypergeometric enrichment analyses from the differential data. **Enrichment_Analyses.R** runs the enrichment by calling the **HPEA.R** script. The remaining files in the directory contain the enrichment genesets for different datasets.


Within **SingleCell_Scripts**, we have the following:

**CellrangerProcessing**

The scripts in the above subdirectory process the single cell fastq files using Cellranger-ATAC for snATAC samples and Cellranger-ARC for Multiomic Samples. **snATAC_Cellranger.sh** runs cellranger-ATAC for snATAC samples and **Multiome_Cellranger_Example.sh** runs Cellranger-ARC for Multiomic samples. Cellranger-ARC requires a library file to run and an examples file (**Example_library.csv**) has been provided for the same. 

Run these scripts for each sample to complete initial processing.

**IndividualSampleProcessing**

The scripts in the above subdirectory contains script for individual processing of the snATAC and Multiomic samples. The **Signac_snATAC_Processing.R** runs the Signac pipeline for the snATAC samples and the **Multiome_WNN.R** runs the Seurat WNN pipeline for the Multiome samples. In between, for Quality checks, we used scrublet (**scrublet_script.py**) to remove RNA doublets and AMULET/snATAC_DoubletDetector (**Amulet.sh**) to remove ATAC doublets. The **Edit.sh** script incorporates the Amulet and scrublet script and output a **AllBarcodes_OutliersandDoubletRemovedTwice.csv** file which is required for pseudobulk processes. This script was written for the Multiome samples but can be modified to exclude certain steps to make it run for snATAC samples.

For snATAC_DoubletDetector version used in Doublet Removal for Single cell data, download the software from given link https://drive.google.com/drive/folders/1LTvF61PjECwkk4K-zGxH8EgWLNO6xW2L?usp=sharing

**MergingSamples**

The scripts in the above subdirectory merges RNA and ATAC profiles from the Multiome samples only post conducting QC for each individual sample. **Seurat_GEX.R** is run for each sample to create RNA profiles for each Multiome sample and these profiles are merged using **Merge_Seurat_RNA.R**. For ATAC profiles, we run **IntegrateSamples_ATAC.R** using **Run_IntegrateSamples_ATAC.sh** to recall peaks for each sample and recount fragments and do the merging. IntegrateSamples_ATAC.R requires a text file containing the location of bed files and fragment files for the merging step. An example text file has been provided for the same (**FinalSampleMergeSignac.txt**).

**MergedDataProcessing**

The scripts in the above subdirectory run the analysis for the Merged RNA and ATAC objects. The **QC_DORC_Sai.R** script took the pooled Seurat objects and ran an additional QC step based on the ATAC quality and removed a few more cells based on that. It runs **ATACRNAcorr.R** to calculate correlation between the RNA and ATAC object which is used to calculate the DORC matrix. The **Merged_ATAC_Processing.R** and **Merged_RNA_Processing.R** then run the Signac and Seurat pipelines respectively for the final set of cells. Those scripts also include the parameters for Single Cell Differential Analysis, Motif Analysis, UMAP visualization steps, footprinting and chromvar matrix calculation.

**Pseudobulk**

The scripts in the above subdirectory creates pseudobulk ATAC profiles from the Multiome and snATAC samples. We created these profiles for HSPC and CD14 Monocytes to analyze these samples alongside the available bulk ATAC data. **Pseudobulk.sh** creates the pseudobulk profiles and does peak calling for the same. It utilizes two other files to run (**Markers.txt** and **AllBarcodes_OutliersandDoubletRemovedTwice.csv**). **HINT_SPARK.R** conducts motif analyses for each clinical condition by pooling these pseudobulk profiles per condition and also creates bedgraph files which are used for visualization using SPARK.
