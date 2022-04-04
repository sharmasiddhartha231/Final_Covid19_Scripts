#!/bin/sh

## The 2 snATAC samples were processed using Cellranger-ATAC 1.2.0 using hg38 as reference sample.
## Only 2 snATAC samples are used in the cohort. 

/home/sharms/cellranger-atac-1.2.0/cellranger-atac count --id=HDf \
                   					                             --reference=/home/sharms/refdata-cellranger-atac-GRCh38-1.2.0/ \
                  					                             --fastqs=/projects/ucar-lab/Siddhartha/Covid_19/snATACseq/snATAC/Covid19/outs/fastq_path/H5YFFBGXF \
   					                                             --sample=HDf

/home/sharms/cellranger-atac-1.2.0/cellranger-atac count --id=jcov31 \
                                                         --reference=/home/sharms/refdata-cellranger-atac-GRCh38-1.2.0/ \
                                                         --fastqs=/projects/ucar-lab/Siddhartha/Covid_19/snATACseq/snATAC/Covid19/outs/fastq_path/H5YFFBGXF \
                                                         --sample=jcov31

