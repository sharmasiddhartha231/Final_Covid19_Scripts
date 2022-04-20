#!/bin/sh

#SBATCH --job-name=Multiome_ScATAC_Batch_2_jcov49-2
#SBATCH --ntasks=1
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=50000

## Showing example run of Multiome processing using cellranger-arc
## Please change reference and libraries directory according to where they are stored on your system.
## Please change path of cellranger to where it is installed on your system.
## Example library file provided for said sample
## Example library csv file requires full path in 2 steps. For example, if the RNA libraries are stored in the following directory: 
# /projects/ucar-lab/Siddhartha/Covid_19/snATACseq/MultiOme/2nd_multiome/2nd_Multiome_Package_ATAC/Josefowicz-JC-9745_2020_12_10/jcov42, 
#then the fastqs column will be: /projects/ucar-lab/Siddhartha/Covid_19/snATACseq/MultiOme/2nd_multiome/2nd_Multiome_Package_ATAC/Josefowicz-JC-9745_2020_12_10
#and the sample column will be: jcov42
#For library type, specify Chromatin Accessibility or Gene Expression for ATAC and RNA libraries.

## We merged certain samples coming from the same patients. To do that step, specify each RNA library to be combined together and ATAC libraries to be combined after the RNA libraries in the libraries csv file.

cellranger-arc-1.0.0/cellranger-arc count --id=jcov42 \
				                                  --reference=/projects/ucar-lab/Siddhartha/Covid_19/snATACseq/MultiOme/refdata-cellranger-arc-GRCh38-2020-A \
					                                --libraries=/projects/ucar-lab/Siddhartha/Covid_19/snATACseq/MultiOme/2nd_multiome/MultiomeRuns/Example_library.csv \
					                                --localcores=16 \
                                          --localmem=50
