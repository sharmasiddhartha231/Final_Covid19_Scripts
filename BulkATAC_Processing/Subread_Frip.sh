#!/bin/sh

## Need to install subread for checking Frip Scores
## PeakList.txt contains full path of peak files for samples.
## Bamlist.txt contains full path of final sorted indexed bam files for all samples.
## The following script runs a batch job for all samples specified in the BamList.txt and PeakList.txt

while IFS= read -r p && IFS= read -r q <&3; do
	echo ${p}
	echo ${q}
	F=$(basename ${q})
	F1="${F%.*}"
	echo ${F1}
	awk 'OFS="\t" {print $1"-"$2+1"-"$3, $1, $2+1, $3, "+"}' ${p} > foo.saf
	featureCounts -p -a foo.saf -F SAF -o ${F1}.txt ${q} 
	rm foo.saf	
	mv ${F1}.txt Subreads
done <PeakList.txt 3<BamList.txt


