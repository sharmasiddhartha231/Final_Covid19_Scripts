library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)
library(hdf5r)
library(GenomicRanges)
library(future)

args = commandArgs(trailingOnly = TRUE)

infile = args[1] #rows = samples, cols = id, peak, singlecell, fragmentfile
infile
projectname = args[2] #What the project name should be
projectname
outfile = args[3] #the output rds file.
outfile

plan("multiprocess", workers = 4)
options(future.globals.maxSize = 50000 * 1024^2) 

readInput <- function(infile){
  inputdata <- read.table(infile)
  rv <- list(inputdata[,1], inputdata[,2], inputdata[,3], inputdata[,4])
  return(rv)
}

getPeakUnion <- function(ids, peakfiles){
  peakobjs <- c()
  for(i in 1:length(ids)) {
      peakobjs <- append(peakobjs, makeGRangesFromDataFrame(read.table(
      file = peakfiles[i],
      col.names = c("chr", "start", "end")
    )))
  }
  peakobjs <- peakobjs[grepl("chr", peakobjs@seqnames),]
  print(peakobjs)
  
  combinedpeaks <- reduce(x=peakobjs)
  
  peakwidths <- width(combinedpeaks)
  combinedpeaks <- combinedpeaks[peakwidths  < 10000 & peakwidths > 20]
  
  return(combinedpeaks)
}

getFilteredSeuratObject <- function(id, singlecellcsv, fragmentfile, combinedpeaks) {
  
  metadata <- read.table(
    file = singlecellcsv,
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )[-1, ]
  metadata = metadata[metadata$is_cell == 1,]
  metadata

  fragobj <- CreateFragmentObject(
    path = fragmentfile,
    cells = rownames(metadata)
  )
  
  countmatrix <- FeatureMatrix(
    fragments = fragobj,
    features = combinedpeaks,
    cells = rownames(metadata)
  )
  
  chromassay <- CreateChromatinAssay(countmatrix, fragments = fragobj)
  rv <- CreateSeuratObject(chromassay, assay = "ATAC")
  
  return(rv)
}

getMergedSeuratObject <- function(ids, peakfiles, singlecellfiles, fragmentfiles, projectname) {
  combinedpeaks <- getPeakUnion(ids, peakfiles)
  
  print(combinedpeaks)
  
  firstobject <- getFilteredSeuratObject(ids[1], singlecellfiles[1], fragmentfiles[1], combinedpeaks)
  
  restobjects <- c()
  for(i in 2:length(ids)) {
    restobjects <- c(restobjects, getFilteredSeuratObject(ids[i], singlecellfiles[i], fragmentfiles[i], combinedpeaks))
  }
  
  rv <- merge(firstobject, y = restobjects, add.cell.ids = ids, project = projectname)
  
  return(rv)
}

inputdata = readInput(infile)

ids <- inputdata[[1]]
peakfiles <- inputdata[[2]]
singlecellfiles <- inputdata[[3]]
fragmentfiles <- inputdata[[4]]


mergedobject = getMergedSeuratObject(ids, peakfiles, singlecellfiles, fragmentfiles, projectname)
saveRDS(mergedobject, outfile)
