setwd("/Volumes/GoogleDrive/My Drive/000 Broad Data/20210405_covid_v3/Rcode/")
# Load packages -----------------------------------------------------------
options(timeout = 300)
pacman::p_load(
  BiocParallel, chromVAR, motifmatchr, BSgenome.Hsapiens.UCSC.hg19, BSgenome.Mmusculus.UCSC.mm10, 
  gplots, RColorBrewer, chromVARmotifs, BuenColors, mclust, scater,Seurat, 
  dplyr, tidyr, data.table, ggplot2,SummarizedExperiment, preprocessCore, 
  reshape, pheatmap,Rtsne, spatstat, Matrix, proxy, densityClust, 
  SummarizedExperiment, ggrepel,BuenRTools, FNN, ggpubr, NetworkToolbox, # chromVARxx
  BuenColors, MASS, viridis, tidyverse, umap, igraph, cisTopic, mgcv, 
  Rsamtools, irlba, DoubletFinder, ComplexHeatmap, seriation, rhdf5)
mytheme2 <- theme(text = element_text(size=24*0.2),
                  axis.text = element_text(size=20*0.2), 
                  axis.line = element_line(colour = 'black', size = 1*0.2),
                  axis.ticks = element_line(colour = 'black', size = 1*0.2),
                  axis.ticks.length=unit(.2*0.2, "cm"),
                  legend.text=element_text(size=16*0.2),
                  legend.title=element_text(size=16*0.2),
                  legend.key.height = unit(.15, "cm"),
                  legend.key.width = unit(.15, "cm"))
data("human_pwms_v1"); data("mouse_pwms_v1"); 
sout <- sapply(strsplit(names(human_pwms_v1), split = "_"), function(s) c(s[3]))
human_pwms_v2 <- human_pwms_v1[match(unique(sout), sout)]
sout <- sapply(strsplit(names(mouse_pwms_v1), split = "_"), function(s) c(s[3]))
mouse_pwms_v2 <- mouse_pwms_v1[match(unique(sout), sout)]
save(human_pwms_v2, file = 'human_pwms_v2.rda'); save(mouse_pwms_v2,  file = 'mouse_pwms_v2.rda')
register(MulticoreParam(1))
tsnecols = c("#E31A1C","#FFD700","#771122","#777711","#1F78B4","#68228B","#AAAA44",
             "#60CC52","pink","#DDDD77","#774411","#AA7744","#AA4455","#117744",
             "#000080","#44AA77","#AA4488","#DDAA77")
tSNEpath <- "/Users/saima/Documents/tSNE"
umap.defaults$min_dist <- 0.3
register(MulticoreParam(2))
options(stringsAsFactors = FALSE)
dir.create("../ATAC")


setOldClass(Classes = 'package_version')
setClassUnion(name = 'AnyMatrix', members = c("matrix", "dgCMatrix"))
setClassUnion(name = 'OptionalCharacter', members = c('NULL', 'character'))

# load atac and RNA data ----------------------------------------------------------

seuset <- readRDS("../raw/PooledSeuratRNA.rds")
# filter RNA
seuset <- subset(x = seuset, subset = nFeature_RNA > 300 & nFeature_RNA < 7000 & percent.mt < 25)


# Load ATAC data (covid peakset) ----------------------------------------------------------
covid.peak <- getPeaks(filename = "covid.peak.20210521.bed", sort_peaks = F)
# saveRDS(covid.peak, "covid.peak.rds")

read10x.multiome.atac <- function(dir, sample, cellid="", peaks=covid.peak){
  if (cellid == ''){
    counts <- getCountsFromFrags(paste(dir, "/atac_fragments.tsv.gz", sep=""), peaks)
  } else {
    counts <- getCountsFromFrags(paste(dir, "/atac_fragments.tsv.gz", sep=""), peaks, barcodeList = cellid)
  }
  # counts$original <- colnames(counts)
  colnames(counts) <- paste(sample, colnames(counts), sep=".")
  counts$sample <- sample
  return(counts)
}

## only exact barcode exits in RNA data
table(seuset$sample)
sample <- "HD_AD"; cells <- seuset@meta.data$original[seuset$sample == sample]
bulk1 <- read10x.multiome.atac(dir = paste("../raw/", sample, sep=""), sample = sample, cellid = cells, peaks=covid.peak)
sample <- "HDa"; cells <- seuset@meta.data$original[seuset$sample == sample]
bulk2 <- read10x.multiome.atac(dir = paste("../raw/", sample, sep=""), sample = sample, cellid = cells, peaks=covid.peak)
sample <- "HDj"; cells <- seuset@meta.data$original[seuset$sample == sample]
bulk3 <- read10x.multiome.atac(dir = paste("../raw/", sample, sep=""), sample = sample, cellid = cells, peaks=covid.peak)
sample <- "HDu"; cells <- seuset@meta.data$original[seuset$sample == sample]
bulk4 <- read10x.multiome.atac(dir = paste("../raw/", sample, sep=""), sample = sample, cellid = cells, peaks=covid.peak)
sample <- "HDv"; cells <- seuset@meta.data$original[seuset$sample == sample]
bulk5 <- read10x.multiome.atac(dir = paste("../raw/", sample, sep=""), sample = sample, cellid = cells, peaks=covid.peak)
sample <- "jcov21"; cells <- seuset@meta.data$original[seuset$sample == sample]
bulk6 <- read10x.multiome.atac(dir = paste("../raw/", sample, sep=""), sample = sample, cellid = cells, peaks=covid.peak)
sample <- "jcov42_20"; cells <- seuset@meta.data$original[seuset$sample == sample]
bulk7 <- read10x.multiome.atac(dir = paste("../raw/", sample, sep=""), sample = sample, cellid = cells, peaks=covid.peak)
sample <- "jcov42"; cells <- seuset@meta.data$original[seuset$sample == sample]
bulk8 <- read10x.multiome.atac(dir = paste("../raw/", sample, sep=""), sample = sample, cellid = cells, peaks=covid.peak)
sample <- "jcov49_2"; cells <- seuset@meta.data$original[seuset$sample == sample]
bulk9 <- read10x.multiome.atac(dir = paste("../raw/", sample, sep=""), sample = sample, cellid = cells, peaks=covid.peak)
sample <- "jcov52"; cells <- seuset@meta.data$original[seuset$sample == sample]
bulk10 <- read10x.multiome.atac(dir = paste("../raw/", sample, sep=""), sample = sample, cellid = cells, peaks=covid.peak)
sample <- "jcov56_9"; cells <- seuset@meta.data$original[seuset$sample == sample]
bulk11 <- read10x.multiome.atac(dir = paste("../raw/", sample, sep=""), sample = sample, cellid = cells, peaks=covid.peak)
sample <- "jcov60_6"; cells <- seuset@meta.data$original[seuset$sample == sample]
bulk12 <- read10x.multiome.atac(dir = paste("../raw/", sample, sep=""), sample = sample, cellid = cells, peaks=covid.peak)
sample <- "jcov69_8"; cells <- seuset@meta.data$original[seuset$sample == sample]
bulk13 <- read10x.multiome.atac(dir = paste("../raw/", sample, sep=""), sample = sample, cellid = cells, peaks=covid.peak)
sample <- "jcov72_5"; cells <- seuset@meta.data$original[seuset$sample == sample]
bulk14 <- read10x.multiome.atac(dir = paste("../raw/", sample, sep=""), sample = sample, cellid = cells, peaks=covid.peak)
sample <- "jcov92"; cells <- seuset@meta.data$original[seuset$sample == sample]
bulk15 <- read10x.multiome.atac(dir = paste("../raw/", sample, sep=""), sample = sample, cellid = cells, peaks=covid.peak)
sample <- "jcov93"; cells <- seuset@meta.data$original[seuset$sample == sample]
bulk16 <- read10x.multiome.atac(dir = paste("../raw/", sample, sep=""), sample = sample, cellid = cells, peaks=covid.peak)
sample <- "jcov101_2"; cells <- seuset@meta.data$original[seuset$sample == sample]
bulk17 <- read10x.multiome.atac(dir = paste("../raw/", sample, sep=""), sample = sample, cellid = cells, peaks=covid.peak)
sample <- "jcov101"; cells <- seuset@meta.data$original[seuset$sample == sample]
bulk18 <- read10x.multiome.atac(dir = paste("../raw/", sample, sep=""), sample = sample, cellid = cells, peaks=covid.peak)
sample <- "jcov110"; cells <- seuset@meta.data$original[seuset$sample == sample]
bulk19 <- read10x.multiome.atac(dir = paste("../raw/", sample, sep=""), sample = sample, cellid = cells, peaks=covid.peak)
sample <- "jcov111"; cells <- seuset@meta.data$original[seuset$sample == sample]
bulk20 <- read10x.multiome.atac(dir = paste("../raw/", sample, sep=""), sample = sample, cellid = cells, peaks=covid.peak)
sample <- "jcov114"; cells <- seuset@meta.data$original[seuset$sample == sample]
bulk21 <- read10x.multiome.atac(dir = paste("../raw/", sample, sep=""), sample = sample, cellid = cells, peaks=covid.peak)
sample <- "jcov115"; cells <- seuset@meta.data$original[seuset$sample == sample]
bulk22 <- read10x.multiome.atac(dir = paste("../raw/", sample, sep=""), sample = sample, cellid = cells, peaks=covid.peak)
sample <- "jcov117"; cells <- seuset@meta.data$original[seuset$sample == sample]
bulk23 <- read10x.multiome.atac(dir = paste("../raw/", sample, sep=""), sample = sample, cellid = cells, peaks=covid.peak)
sample <- "jcov119"; cells <- seuset@meta.data$original[seuset$sample == sample]
bulk24 <- read10x.multiome.atac(dir = paste("../raw/", sample, sep=""), sample = sample, cellid = cells, peaks=covid.peak)
sample <- "jcov122"; cells <- seuset@meta.data$original[seuset$sample == sample]
bulk25 <- read10x.multiome.atac(dir = paste("../raw/", sample, sep=""), sample = sample, cellid = cells, peaks=covid.peak)
sample <- "jcov124"; cells <- seuset@meta.data$original[seuset$sample == sample]
bulk26 <- read10x.multiome.atac(dir = paste("../raw/", sample, sep=""), sample = sample, cellid = cells, peaks=covid.peak)
sample <- "jcov127"; cells <- seuset@meta.data$original[seuset$sample == sample]
bulk27 <- read10x.multiome.atac(dir = paste("../raw/", sample, sep=""), sample = sample, cellid = cells, peaks=covid.peak)
sample <- "lgtd8m"; cells <- seuset@meta.data$original[seuset$sample == sample]
bulk28 <- read10x.multiome.atac(dir = paste("../raw/", sample, sep=""), sample = sample, cellid = cells, peaks=covid.peak)
sample <- "lgtd17"; cells <- seuset@meta.data$original[seuset$sample == sample]
bulk29 <- read10x.multiome.atac(dir = paste("../raw/", sample, sep=""), sample = sample, cellid = cells, peaks=covid.peak)
sample <- "lgtd18"; cells <- seuset@meta.data$original[seuset$sample == sample]
bulk30 <- read10x.multiome.atac(dir = paste("../raw/", sample, sep=""), sample = sample, cellid = cells, peaks=covid.peak)


counts.raw.covid <- cbind(bulk1, bulk2, bulk3,  bulk4,
                         bulk5, bulk6,  bulk7,  bulk8,
                         bulk9, bulk10, bulk11, bulk12,
                         bulk13,bulk14, bulk15, bulk16, 
                         bulk17,bulk18, bulk19, bulk20,
                         bulk21,bulk22, bulk23, bulk24,
                         bulk25,bulk26, bulk27, bulk28,
                         bulk29,bulk30)
table(counts.raw.covid$sample)
# saveRDS(counts.raw.covid, "counts.raw.covid.rds")

# double check intersect of ATAC and RNA object
atac.se <- counts.raw.covid
Intersect <- intersect(colnames(seuset), colnames(atac.se)); length(Intersect)

# filtering sample on ATAC
minpeak=0.3; mindepth=1000; readCut=500
filtering_plot  <- filterSamplesPlot(atac.se, min_depth = mindepth, min_in_peaks = minpeak, use_plotly = F); filtering_plot
index1 <- filterSamples(atac.se, min_depth = mindepth, min_in_peaks = minpeak, shiny=F, ix_return = T)
atac.se <- atac.se[, index1]; atac.se
seuset <- subset(seuset, cells =  colnames(atac.se))


# Cis peak-gene association ---------------------------------------------------------------------

RNA <- seuset@assays$SCT@counts
RNA <- RNA[rowSums(RNA) > 0, ]
RNA <- RNA[!(rownames(RNA) %like% "Gm" | rownames(RNA) %like% "Rik"), ]
dim(RNA)

library(BSgenome.Hsapiens.UCSC.hg38)

source("/Users/saima/Google Drive/000 Broad Data/ATACRNAcorr.R")

SE.filt <- atac.se; rm(atac.se)
RNAmat <- RNA; rm(seuset); rm(RNA)

RNAmat <- RNAmat[rowSums(RNAmat) > 0, ]
SE.filt <- SE.filt[rowSums(counts(SE.filt))> 0, ]

SE.filt <- addGCBias(SE.filt, genome = BSgenome.Hsapiens.UCSC.hg38)
ATACmat <- counts(SE.filt)
rownames(ATACmat) <- paste("Peak",1:nrow(ATACmat), sep="")
dir.create("../ATACRNA/")
dir.create("../ATACRNA/corFiles/")
dim(ATACmat); dim(RNAmat)


## pre-run with 100 interations, and 50k window
runTest <- TRUE
if(runTest){
  # Keep genes that have annotation and are in RNA matrix
  hg38TSS <- BuenRTools::hg38TSSRanges
  names(hg38TSS) <- as.character(hg38TSS$gene_name)
  
  # Checking in case some genes in RNA don't overlap our TSS annotations
  genesToKeep <- intersect(names(hg38TSS),rownames(RNAmat))
  cat("Num genes overlapping hg38 TSS annotation and RNA matrix being considered:\n")
  print(length(genesToKeep))
  
  # Match gene order in RNA matrix and TSS ranges
  RNAmat <- RNAmat[genesToKeep,]
  hg38TSS <- hg38TSS[genesToKeep]
  
  # Pad TSS by this much *either side*
  windowSize <- 50000
  
  cat("Taking +/-",windowSize,"bp around each TSS to determine peak-gene overlaps ..\n")
  TSSflank <- GenomicRanges::flank(hg38TSS, width = windowSize, both = TRUE)
  saveRDS(TSSflank, "../ATACRNA/TSSflank.rds")
  peakRanges <- granges(SE.filt)
  
  # Get peak summit
  cat("Taking peak summits from peak windows ..\n")
  peakSummits <- resize(peakRanges,width = 1,fix = "center")
  
  # Find overlap of all peaks to all genes given window
  # Subject is Peaks, query is Gene
  cat("Finding overlapping peak-gene pairs ..\n")
  genePeakOv <- findOverlaps(query = TSSflank,subject = peakSummits)
  numPairs <- length(genePeakOv)
  
  cat("Found ",numPairs,"total gene-peak pairs for", windowSize, "TSS window ..\n")
  
  cat("Number of peak summits that overlap any gene TSS window:\n")
  print(length(unique(subjectHits(genePeakOv))))
  cat("Number of gene TSS windows that overlap any peak summit:\n")
  print(length(unique(queryHits(genePeakOv))))
  
  saveRDS(genePeakOv,paste0("../ATACRNA/genePeakOv_",as.character(windowSize*2/1000),"_kb.rds"))
  
  # Also save the hg38 TSS and peakRanges for these hits (useful for filtering/indexing/subsetting)
  saveRDS(hg38TSS,"../ATACRNA/hg38TSS_filt_GRanges.rds")
  saveRDS(peakRanges,"../ATACRNA/peak_filt_GRanges.rds")
  
  # For each gene, determine observed correlation of each overlapping peak to its associated gene (gene expression)
  
  # For each of those genes, also determine correlation based on background peaks and save
  # Do this in parallel, and get data frame of gene-peak-pearson values
  # Fetch background peaks for each peak tested (i.e. that has overlap in window with gene)
  set.seed(42)
  cat("Fetching background peaks ..\n")
  
  bg <- getBackgroundPeaks(SE.filt,niterations=100)
  # bg <- readRDS("bg.rds")
  saveRDS(bg,"bg.rds")
  
  cat("Computing gene-peak correlations ..\n")
  
  numCores <- 1
  pairsPerChunk <- 500
  
  # Folder where chunked files will be saved
  outDir <- "../ATACRNA/corFiles/"
  
  # This defines the outer (larger chunks)
  largeChunkSize <- 10000
  
  # If for any reason the workers fail, resume from where it failed by specifying the starting point here
  # You can keep track of this by looking at the files that are successfully written after each chunk is run (see below)
  
  #startingPoint <- 120001
  startingPoint <- 1
  chunkStarts <- seq(startingPoint, numPairs, largeChunkSize)
  chunkEnds <- chunkStarts + largeChunkSize -1
  chunkEnds[length(chunkEnds)] <- numPairs
  
  
  for(i in 1:length(chunkStarts)){
    cat("Running pairs: ",chunkStarts[i], "to",chunkEnds[i],"\n")
    # This fill further chunk this up and run in parallel, saving the merged output ObsCor
    ObsCor <- PeakGeneCor(ATAC = ATACmat,
                          RNA = RNAmat,
                          OV = genePeakOv[chunkStarts[i]:chunkEnds[i]],
                          chunkSize = pairsPerChunk,
                          ncores = numCores,
                          bg = bg)
    gc()
    # Add time stamp and save to file
    st=format(Sys.time(), "%Y-%m-%d_%H:%M")
    
    # Write chunk to txt file, so it can be loaded and merged later for analysis
    chunkFile <- paste0(outDir,"hg38TSS_",as.character(windowSize*2/1000),"kb_GenePeakCorr_spearman_",chunkStarts[i],"_",chunkEnds[i],"_",st,".txt")
    cat("Saving output to file: ",chunkFile,"\n")
    write.table(ObsCor,chunkFile,sep="\t",quote=FALSE)
  }
  
  cat("Finished!\n")
}

getZscore <- function(x){
  pop_sd <- sd(x[4:length(x)])
  pop_mean <- mean(x[4:length(x)])
  zscore <- data.frame(
    Gene=x[1],
    Peak=x[2],
    Zscore=(x[3] - pop_mean) / pop_sd,
    Pvalue=1 - pnorm(x[3], pop_mean, pop_sd),
    Corr=x[3])
  return(zscore)
}

files <- list.files("../ATACRNA/corFiles/"); files
count=1
for (file in files){
  print(file)
  temp <- read.table(paste("../ATACRNA/corFiles/", files[count], sep=""))
  temp2 <- apply(temp, 1, getZscore)
  df <- data.frame(matrix(unlist(temp2), nrow=length(temp2), byrow=T))
  colnames(df) <- colnames(temp2[1]$`1`)
  if (count == 1){
    Cis.pre <- df
  } else {
    Cis.pre <- rbind(Cis.pre, df)
  }
  count=count+1
}

head(Cis.pre); dim(Cis.pre)
# TSSflank <- readRDS("../ATACRNA/TSSflank.rds")
Cis.pre$genename <- TSSflank$gene_name[Cis.pre$Gene]
# saveRDS(Cis.pre, "Cis.pre.rds")

Cis.pre.sig <- Cis.pre[Cis.pre$Pvalue < 0.05, ]
dim(Cis.pre.sig)


## filter duplicatd Cis
temp <- Cis.pre.sig %>% group_by(Peak)  %>% filter(Pvalue == min(Pvalue)) %>% filter(Corr == max(Corr)); temp
temp <- temp[!duplicated(temp$Peak), ]
dim(temp)
Cis.pre.sig <- as.data.frame(temp); dim(Cis.pre.sig)

peaks <- data.frame(chr=seqnames(SE.filt), start=start(SE.filt), end=end(SE.filt))
Cis.pre.sig$chr <- peaks$chr[Cis.pre.sig$Peak]
Cis.pre.sig$start <- peaks$start[Cis.pre.sig$Peak]
Cis.pre.sig$end <- peaks$end[Cis.pre.sig$Peak]
head(Cis.pre.sig)

Cis.pre.sig$genename <- factor(Cis.pre.sig$genename, levels = TSSflank$gene_name)

# saveRDS(Cis.pre.sig, "Cis.pre.sig.rds")
# write.csv(Cis.pre.sig, "Cis.pre.sig.csv", quote = F)

plot(sort(table(Cis.pre.sig$genename)), ylab="No. peak-gene association")
tail(sort(table(Cis.pre.sig$genename)), 20)
temp <- table(Cis.pre.sig$genename)
temp["LEF1"]

# write.table(as.data.frame((tail(sort(table(Cis.pre.sig$genename)), 20))), "topDORCs.txt", quote=F)

temp <- as.data.frame(sort(table(Cis.pre.sig$genename)))
head(temp)
interst.gene <- tail(temp, 60)$Var1
temp$index <- 1:nrow(temp)
temp$interst.gene <- ""
temp$interst.gene[temp$Var1 %in% interst.gene]  <- as.character(temp$Var1[temp$Var1 %in% interst.gene])
unique(temp$interst.gene)

ggplot(temp, aes(x=index, y=Freq)) + geom_point(color="darkred", size=1, stroke=0) +
  geom_line(color="darkred", size=0.1) + 
  geom_text_repel(aes(label = interst.gene), nudge_x = -200, size= 3, segment.size = 0.2, point.padding = 0.5, box.padding = 0.1) +
  ylab("Number of correlated peaks") +
  xlab("Rank") +
  mytheme2 +
  pretty_plot(fontsize = 10) +
  L_border()
ggsave("../ATAC/Jplot.png", dpi = 1200, width = 12, height = 9, units = "in")
ggsave("../ATAC/Jplot.pdf", width = 12, height = 9, units = "in", useDingbats=FALSE)


sum(table(Cis.pre.sig$genename) >= 8)
DORCgene <- names(table(Cis.pre.sig$genename)[table(Cis.pre.sig$genename) >= 8])


Cis.pre.sig$DORC.gene <- FALSE
Cis.pre.sig$DORC.gene[Cis.pre.sig$genename %in% DORCgene] <- TRUE
head(Cis.pre.sig)
write.table(Cis.pre.sig, "Cis.pre.sig.csv", quote = F, sep="\t")

# saveRDS(DORCgene, "DORCgene.rds")

# calculate DORC score -------------------------------------------------------
## filter peaks with multiple genes, choose minimium p value
Cis.target.sig.filter <- Cis.pre.sig[Cis.pre.sig$genename %in% DORCgene, ]
length(unique(Cis.target.sig.filter$Peak))

# saveRDS(Cis.target.sig.filter, "Cis.target.sig.filter.rds")

# sum of super correlators 50k window
SEsum <- data.frame(matrix(0,nrow = length(unique(Cis.target.sig.filter$genename)), ncol = ncol(counts(atac.se))))
colnames(SEsum) <- colnames(atac.se)
rownames(SEsum) <- unique(Cis.target.sig.filter$genename)
SEsum[1:3, 1:3]


## patallel
registerDoParallel(4)
library(foreach)
library(parallel)

ATAC <- counts(atac.se) #[ , 1:100]
count=0; temp <- foreach (i=unique(Cis.target.sig.filter$genename), .combine=cbind) %do% {
  count=count+1
  i <- as.character(i)
  print(paste(count,i, sep=":"))
  colSums(ATAC[Cis.target.sig.filter$Peak[Cis.target.sig.filter$genename == i], ])
}

SEsum <- t(temp)
rownames(SEsum) <- unique(Cis.target.sig.filter$genename)
SEsum.norm = sweep(SEsum, 2, colMeans(counts(atac.se)), `/`)

# saveRDS(SEsum, "SEsum.rds")
# saveRDS(SEsum.norm, "SEsum.norm.rds")


# smooth DORC over ATAC 50 KNNs
lsi <- readRDS("lsi.rds")
lsi <- lsi[colnames(atac.se), ]
k = 50 # smaller k, get more clusters
temp <- lsi
knn.lsi = get.knn(lsi[, ], k = k)
dim(knn.lsi$nn.index)
# saveRDS(knn.lsi, "knn.lsi.rds")


library("parallel")
registerDoParallel(4)
SEsum_smooth <- foreach (j=1:nrow(SEsum.norm), .combine=rbind) %do% {
  print(j)
  SE.smooth <- as.numeric(SEsum.norm[j, ])
  SE <- SE.smooth; for (i in 1:ncol(SEsum.norm)){ SE.smooth[i] <- mean(SE[knn.lsi$nn.index[i,1:50]])}
  SE.smooth
}

colnames(SEsum_smooth) <- colnames(SEsum.norm)
rownames(SEsum_smooth) <- rownames(SEsum.norm)
# saveRDS(SEsum_smooth, "SEsum_smooth.rds")


