#!/bin/sh

usage()
{
  echo "Usage: BulkATAC_Processing.sh [-a Input Fastq File 1 - Full Path] [-b Input Fastq File 2 - Full Path] [-o Full path of the Output Directory]";
}
if [[ $1 == "" ]]; then
    echo "No options provided.";
    usage;
    exit;
else
  while getopts "a:b:o:" OPTIONS; do
    case $OPTIONS in
      a) inputdir1=$OPTARG;;
      b) inputdir2=$OPTARG;;
      o) outputdir=$OPTARG;;
      *) echo "Incorrect Options Provided"
         usage;
         exit 1;;
    esac
  done

  ## Preferably requires the file to be in .fastq format. Works with .fastq.gz format files as well.
  ## FileNameWithExtension gives the entire filename without the entire path.
  ## FileNameWithoutExtension removes the extension. Keeps the .fastq in case a zipped file is provided.
  FileNameWithExtension1=$(basename $inputdir1)
  FileNameWithExtension2=$(basename $inputdir2)
  FileNameWithoutExtension1="${FileNameWithExtension1%.*}"
  FileNameWithoutExtension2="${FileNameWithExtension2%.*}"
  FNA1="${FileNameWithoutExtension1%.*}"
  FNA2="${FileNameWithoutExtension2%.*}"

  echo $FileNameWithExtension1
  echo $FileNameWithoutExtension1
  echo $FileNameWithExtension2
  echo $FileNameWithoutExtension2
  echo $FNA1
  echo $FNA2

  ## CommonName keeps the Filename uptil the substring that is common for both the read files.
  ## Makes a directory by the common name to store the output files.
  CommonName=$(printf "%s\n%s\n" "$FileNameWithoutExtension1" "$FileNameWithoutExtension2" | sed -e 'N;s/^\(.*\).*\n\1.*$/\1/')
  echo $CommonName
  echo $outputdir

  ## Make all the Directories for a sample within the specified output directory name.
  mkdir -p ${outputdir}
  mkdir -p ${outputdir}/${CommonName}
  mkdir -p ${outputdir}/${CommonName}/FastQC
  mkdir -p ${outputdir}/${CommonName}/Trimmomatic
  #mkdir -p ${outputdir}/${CommonName}/Trimmomatic/adapterTrimmed
  mkdir -p ${outputdir}/${CommonName}/bwa
  #mkdir -p ${outputdir}/${CommonName}/Broad_Peaks
  mkdir -p ${outputdir}/${CommonName}/Bampe_Peaks
  mkdir -p ${outputdir}/${CommonName}/IGV
  #mkdir -p ${outputdir}/${CommonName}/BedGraphs
  #cd ${outputdir}
  ## Generates QC stats for both paired end files. Just to see how your data looks.
  ## Example of good quality QC from FastQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/good_sequence_short_fastqc.html
  ## Stores QC output in FastQC directory
  /path/to/FastQC/0.11.3/fastqc -t 2 --noextract $inputdir1 -o ${outputdir}/${CommonName}/FastQC
  /path/to/FastQC/0.11.3/fastqc -t 2 --noextract $inputdir2 -o ${outputdir}/${CommonName}/FastQC

  ## Read Trimming Step before alignment
  java -jar /path/to/Trimmomatic/0.33/trimmomatic-0.33.jar PE -threads 2 $inputdir1 $inputdir2 ${outputdir}/${CommonName}/Trimmomatic/Paired_${FileNameWithExtension1} ${outputdir}/${CommonName}/Trimmomatic/UnPaired_${FileNameWithExtension1} ${outputdir}/${CommonName}/Trimmomatic/Paired_${FileNameWithExtension2} ${outputdir}/${CommonName}/Trimmomatic/Unpaired_${FileNameWithExtension2} TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 ILLUMINACLIP:/path/to/Trimmomatic/0.33/adapters/TruSeq3-PE.fa:2:30:10
  #trim_adapters ${outputdir}/${CommonName}/Trimmomatic/Paired_${FileNameWithExtension1} ${outputdir}/${CommonName}/Trimmomatic/Paired_${FileNameWithExtension2}
  
  ## Alignment step. Converts to Sam file. Make sure to specify the right genome.
  /path/to/bwa/0.7.12/bin/bwa mem -M /path/to/BWAIndex/hg38/hg38.fa ${outputdir}/${CommonName}/Trimmomatic/Paired_${FileNameWithExtension1} ${outputdir}/${CommonName}/Trimmomatic/Paired_${FileNameWithExtension2} > ${outputdir}/${CommonName}/bwa/${CommonName}.sam

  ## Quality checks. SortSam sorts input file by coordinate.
  ## MarkDuplicates locates and tags duplicate reads in input file.
  ## CollectInsertSizeMetrics provides useful metrics for validating library construction including the insert size distribution and read orientation of paired-end libraries.
  ## ATAC_BAM_shifter_gappedAlign.pl is for adjusting input files for the 9bp binding site on Tn5.
  java -Xms1g -Xmx4g -jar /path/to/picard/1.95/SortSam.jar INPUT=${outputdir}/${CommonName}/bwa/${CommonName}.sam OUTPUT=${outputdir}/${CommonName}/bwa/${CommonName}.sorted.sam SO=coordinate
  java -Xms1g -Xmx4g -jar /path/to/picard/1.95/MarkDuplicates.jar INPUT=${outputdir}/${CommonName}/bwa/${CommonName}.sorted.sam OUTPUT=${outputdir}/${CommonName}/bwa/${CommonName}.sorted_rmdup.sam METRICS_FILE=${outputdir}/${CommonName}/bwa/${CommonName}.metrics.txt
  java -Xms1g -Xmx4g -jar /path/to/picard/1.95/CollectInsertSizeMetrics.jar METRIC_ACCUMULATION_LEVEL=ALL_READS OUTPUT=${outputdir}/${CommonName}/bwa/${CommonName}_insertSize.txt HISTOGRAM_FILE=${outputdir}/${CommonName}/bwa/${CommonName}_insertSize.pdf INPUT=${outputdir}/${CommonName}/bwa/${CommonName}.sorted_rmdup.sam
  perl /path/to/HumanCardiomyocytes/ATACSeq/auyar/ATAC_BAM_shifter_gappedAlign.pl ${outputdir}/${CommonName}/bwa/${CommonName}.sorted_rmdup.sam ${outputdir}/${CommonName}/bwa/${CommonName}.sorted_rmdup_shifted
  
  ## Undo this to remove all sam files created uptil this point
  rm ${outputdir}/${CommonName}/bwa/${CommonName}*sam
  
  ## Sort and Index the bam file before calling peaks.
  ## bedtools bamtobed is a conversion utility that converts sequence alignments in BAM format into BED, BED12, and/or BEDPE records.
  #samtools view -Sb -o ${CommonName}/bwa/${CommonName}.sorted_rmdup.bam ${CommonName}/bwa/${CommonName}.sorted_rmdup.sam
  samtools sort -o ${outputdir}/${CommonName}/bwa/${CommonName}.sorted_rmdup_shifted_sorted.bam ${outputdir}/${CommonName}/bwa/${CommonName}.sorted_rmdup_shifted.bam
  samtools index ${outputdir}/${CommonName}/bwa/${CommonName}.sorted_rmdup_shifted_sorted.bam
  bedtools bamtobed -i ${outputdir}/${CommonName}/bwa/${CommonName}.sorted_rmdup_shifted_sorted.bam > ${outputdir}/${CommonName}/bwa/${CommonName}.sorted_rmdup_shifted_sorted.bed

  ## Undo this to remove the unsorted bam files.
  rm ${outputdir}/${CommonName}/bwa/${CommonName}.sorted_rmdup_shifted.bam

  ## IGV browser creation
  igvtools count ${outputdir}/${CommonName}/bwa/${CommonName}.sorted_rmdup_shifted_sorted.bam ${outputdir}/${CommonName}/IGV/${CommonName}.sorted_rmdup_shifted_sorted.tdf hg38

  ## Peak calling using macs2. BAMPE is preferred for paired end reads.
  # macs2 callpeak -t ${outputdir}/${CommonName}/bwa/${CommonName}.sorted_rmdup_shifted_sorted.bam -f BAM --outdir ${outputdir}/${CommonName}/Broad_Peaks -n ${CommonName}.sorted_rmdup_shifted_sorted_broad_peaks -g 'hs' --nomodel --shift -100 --extsize 200 --broad -B
  macs2 callpeak -t ${outputdir}/${CommonName}/bwa/${CommonName}.sorted_rmdup_shifted_sorted.bam -f BAMPE --outdir ${outputdir}/${CommonName}/Bampe_Peaks -n ${CommonName}.sorted_rmdup_shifted_sorted_bampe_peaks -g 'hs'
  
 fi

## Tools you will need to install: samtools, igvtools, macs2, bedtools, atactk
## Conda installation for most of these should be straigtforward.
   #conda create -n ATAC_Prep openjdk=11 python=3.9
   #conda activate ATAC_Prep
   #java -version
   #conda install -c bioconda igvtools bedtools samtools macs2
   ##conda create -n Trimming python=3.5.2
   ##conda activate Trimming
   ##conda install -c bioconda atactk
## You can install the remaining packages within this environment. Just make sure to activate it before running the script. 
## atactk does not work for python3.7 and above. It works with python2.7 and python3.5. This is the adapter trimming step.
## macs2 does not have support for python3.5.2 and atactk only works with python3.5.2 and below. 
## One way to get around to it is to create a conda environment with python3.5.2 and use it till the trimming step. Then switch over to the normal base setup and run the rest of the script. 

## You can copy over the following folders into your directory and make changes accordingly. Otherwise, it should run from my directory as well.
  # /path/to/picard/1.95/ (For Picard Steps.)
  # /path/to/HumanCardiomyocytes/ATACSeq/auyar/ (for the ATAC_BAM_shifter_gappedAlign.pl script to account for the Tn5 shifting)
  # /path/to/bwa/0.7.12/bin/bwa (Alignment step).
  # /path/to/Trimmomatic/0.33/ (Trimming step).
  # /path/to/Trimmomatic/0.33/adapters/TruSeq3-PE.fa:2:30:10 (adapter trimmers)
  # /path/to/FastQC/0.11.3/ (Initial QC).
  # /path/to/BWAIndex/hg38 (Alignment Genome for hg38)
  # I downloaded the reference genome from: http://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/
  # To index the genome, I used: bwa index -a bwtsw hg38.fa, samtools faidx hg38.fa,java -jar /path/to/picard/1.95/CreateSequenceDictionary REFERENCE=hg38.fa OUTPUT=hg38.dict