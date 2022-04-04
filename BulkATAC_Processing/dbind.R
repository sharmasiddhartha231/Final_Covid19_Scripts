#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#The following script calculates consensus peaks between given set of samples. IT requires a csv file as input containing paths of final sorted indexed bam files and peak files (Example provided alongside code.)

library(DiffBind)

data_path <- paste0 (args[1], '/CD14_hg38.csv')
DBdata <- dba(sampleSheet=data_path,config=data.frame(RunParallel=TRUE))
DBdata <- dba.count(DBdata, minOverlap=2, bParallel=TRUE, score = DBA_SCORE_READS)

##Saves the consensus count matrix to as R object
save(DBdata, file = paste0 (args[1],"/CD14.RData"))

Activation <- dba.peakset(DBdata, bRetrieve=TRUE, DataType=DBA_DATA_FRAME)
##Saves the consensus count matrix to text file
write.table(Activation, file = "path/to/output/CD14.txt", sep = '\t', row.names = F)


