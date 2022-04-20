#!/bin/sh

## The following code creates pseudobulk bam files for each clinical condition for Motif Analysis using HINT and creates bedgraph files for Spark visualization.

## Pseudobulk bam files for each clinical condition. Merge each pseudobulk profile from the Multiomic profiles.

cd /path/to/output

## Repeat for each clinical condition samples.The Individual samples can be created using 
samtools merge {Condition}.bam {Individual_1}.bam {Individual_2}.bam ... {Individual_n}.bam
samtools sort -o Sorted_{Condition}.bam {Condition}.bam
samtools index Sorted_{Condition}.bam
macs2 callpeak -t Sorted_{Condition}.bam -f BAMPE --outdir ./ -n Sorted_{Condition}_bampe_peaks -g 'hs'


## This step creates the bam files
bamCoverage -b Sorted_{Condition}.bam -o Sorted_{Condition}.bdg -bs 1 -of bedgraph

## HINT Steps
rgt-hint footprinting --atac-seq --paired-end --organism=hg38 --output-location=./ --output-prefix={Condition} Sorted_{Condition}.bam Sorted_{Condition}_bampe_peaks.narrowPeak 

## Run this step after running the footprinting step for each condition sample. It needs two files to run a comparative analysis.
rgt-motifanalysis matching --organism=hg38 --input-files {Condition_1}.bed {Condition_2}.bed
rgt-hint differential --organism=hg38 --bc --nc 30 --mpbs-files=./match/{Condition_1}_mpbs.bed,./match/{Condition_2}_mpbs.bed --reads-files=Sorted_{Condition_1}.bam,Sorted_{Condition_2}.bam --conditions={Condition_1},{Condition_2} --output-location={Condition_1}_{Condition_2}
