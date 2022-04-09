library(cinaR)
library(TCseq)
library(pheatmap)
library(ggplot2)
library(tidyr)
library(dplyr)
library(cinaRgenesets)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(RColorBrewer)
library(reshape2)
library(ggforce)
library(ggplot2)
library(paletteer)
library(pheatmap)

setwd("~/path/to.output")

#Specify number of clusters you want to break it down into
Clusters <- 12
#Create a folder by 10_Clusters in the working directory to store output.
#Create subdirectories in the output folder.
dir.create(paste(Clusters,"_Clusters", sep = ""))
dir.create(paste(Clusters,"_Clusters/Motif_Enrichment", sep = ""))
dir.create(paste(Clusters,"_Clusters/Enrichment", sep = ""))

#Read Raw Counts file
Counts <- read.table(file = "~/path/to/ConsensusPeakSet.txt", sep = '\t', header = T)

## Next steps arrange the peak data according to our clinical conditions used as contrast and remove certain peaks. Can change accordingly.
#CD14 <- subset(CD14, nchar(as.character(CHR)) <= 5)
#CD14 <- data.frame(CD14[,c(1:40,42:59,41)])
#CD14 <- CD14[,-59]

## Can be found in Supplementary Tables.
Metadata <- read_excel("~/path/to/Metadata.xlsx")

contrast <- Metadata$Contrast
age <- Metadata$Age
sex <- Metadata$Sex

## Runs Differential Analysis and Annotations using Chipseeker.
## Saves the newly created consensus peaks and differential peaks in a R object.
results <- cinaR(Counts, contrast, reference.genome = "hg38", batch.correction = T, additional.covariates = age, DA.fdr.threshold = 0.1, enrichment.FDR.cutoff = 0.1)
cp <- results$cp

## Extract peaks for Mild vs Healthy, Severe vs Healthy and ICU vs Healthy for the remaining steps. You can add or remove select peaks according to your requirement.
DA.res <- results$DA.peaks[c("Mild_HD", "Sev_HD", "Cont_HD")]

## We ran the analysis for the Severe vs Healthy peaks with FDR less than 0.1 (Significance cutoff). Can be modified according to need.
peak.locs <- results$DA.peaks$Sev_HD$FDR < 0.1

logFCs <- sapply(DA.res, function(x){
  res <- x[peak.locs,"logFC"]
  names(res) <- x[peak.locs,]$Row.names
  return(res)
})

# Added zero columns for healthy samples
logFCs<- cbind(HD = 0, logFCs)
colnames(logFCs) <- strsplit(colnames(logFCs), split = "_", fixed = T) %>% sapply(function(x){x[1]})
#logFCs <- as.data.frame(logFCs)
tca <- timeclust(logFCs, algo = "cm", k = Clusters, standardize = TRUE)

p <- timeclustplot(tca, value = "Normalized LFC", cols = 3)
p
print(p)
save(cp, DA.res, tca, p, file = (paste0(Clusters, "_Clusters/CD14_Severe_Healthy_peaks_results.rda")))
pdf(paste0(Clusters, "_Clusters/CD14_Severe_Healthy_peaks_clustering.pdf"))
print(p)
dev.off()

all.peaks <- DA.res$Cont_HD
DA.Sev_Hea <- DA.res$Sev_HD
k = Clusters # number of clusters
CD14.cluster.genes <- NULL
for (i in c(1:k)){
  CD14.cluster.genes <- rbind(CD14.cluster.genes, 
                              cbind(Cluster = i,
                                    all.peaks[all.peaks$Row.names %in% names(tca@cluster[tca@cluster == i]),
                                              c("Row.names","gene_name")])
  )
  file.to.write <- all.peaks[all.peaks$Row.names %in% names(tca@cluster[tca@cluster == i]),c(1:4,6)]
  file.to.write$strand <- "+"
  colnames(file.to.write) <- c("UniqueID", "Chr", "Start", "End", "Strand")
  write.table(file.to.write, file = paste0(Clusters, "_Clusters/Motif_Enrichment/Cluster", i, "_DA_Peaks.txt"), 
              quote = F, row.names = F, sep = "\t")
  annot.file <- merge(file.to.write, DA.Sev_Hea, by.x = "UniqueID", by.y = "Row.names")
  write.table(annot.file, file = paste0(Clusters, "_Clusters/Enrichment/Cluster", i, "_DA_Annotations.txt"), 
              quote = F, row.names = F, sep = "\t")
}

background.peaks <- DA.Sev_Hea[abs(DA.Sev_Hea$logFC) < 0.2,c(1:4,6)]
background.peaks$Strand <- "+"
colnames(background.peaks) <-  c("UniqueID", "Chr", "Start", "End", "Strand")
write.table(background.peaks, file = paste0(Clusters, "_Clusters/Motif_Enrichment/CD14_Background_peaks.txt"), 
            quote = F, row.names = F, sep = "\t")

pm.data <- tca@data
pm.meta <- data.frame(Cluster = tca@cluster)

pm.order <- order(pm.meta$Cluster)
breaksList = seq(-2, 2, by = .001)

pdf(paste0(Clusters, "_Clusters/Pheatmap_CD14_Severe_Healthy_Clustering.pdf"), useDingbats = FALSE, width = 5, height = 5)
pheatmap(pm.data[pm.order,], cluster_rows = FALSE, cluster_cols = FALSE, annotation_row = pm.meta,
         show_rownames = FALSE, breaks = breaksList, 
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(length(breaksList)), border_color = NA)
dev.off()

# temporal pattern peaks chromatin accessibility heatmap

# reverse the order for visualization
cp <- cp[,ncol(cp):1]
age_cd14 <- rev(age_cd14)
contrast_cd14 <- rev(contrast_cd14)

## S is SingleCell. B is Bulk. The vector counts number of bulk and pseudobulk samples for each condition in the analysis for the heatmap.
exp.type <- c(rep("S",6), rep("B", 11), rep("S", 4), rep("B", 4), 
              rep("S",3), rep("B",4), rep("S", 4), rep("B", 12))

exp.type[exp.type == "S"] <- "snATAC-seq"
exp.type[exp.type == "B"] <- "bulk ATAC-seq"

ann_cols = data.frame (Age = age_cd14, Status = contrast_cd14, Experiment = exp.type)
rownames(ann_cols) <- colnames(cp)

cluster.peaks <- cp[names(tca@cluster[order(tca@cluster)]),]

all.peaks <- DA.res$Sev_HD

all.peaks[all.peaks$Row.names %in% rownames(cluster.peaks),c("Row.names", "gene_name")]

ann_rows = data.frame(Cluster = tca@cluster)


palette1 <- paletteer_d("ggsci::nrc_npg")
palette2 <- paletteer_d("awtools::a_palette")


ann_colors = list(
  Experiment = c("bulk ATAC-seq" = palette2[3], "snATAC-seq" = palette2[7]),
  Status = c("HD" = palette2[2], "Mild" = palette2[4], "Sev" = palette2[6], "Cont" = palette2[8]),
  Cluster = c("1" = palette1[1], "2" = palette1[2], "3" = palette1[3], "4" ="#f7941d")
)

label.row <- rep("", nrow(cluster.peaks))
names(label.row) <- rownames(cluster.peaks)

breaksList = seq(-2, 2, by = .001)

pmap.cd14 <- pheatmap(cluster.peaks, cluster_cols = FALSE, scale = "row", annotation_col = ann_cols, cluster_rows = FALSE,
                      color = colorRampPalette(rev(brewer.pal(n = 9, name = "PuOr")))(length(breaksList)),
                      breaks = breaksList, gaps_col = c(17, 25, 32), show_rownames = F, show_colnames = F, annotation_row = ann_rows,
                      annotation_colors = ann_colors)

pdf(paste0(Clusters, "_Clusters/CD14_accessibility_levels_DA_Severe_vs_Healthy.pdf"), height = 5, width = 7)
print(pmap.cd14)
dev.off()

## BOXPLOT PANEL FOR CLUSTER 
for (i in 1:Clusters){
  
  cluster.peaks <- tca@data[tca@cluster == i,]
  cluster.melt <- melt(cluster.peaks)
  
  ggplot(cluster.melt, aes(x = Var2, y= value, color = Var2)) +
    geom_line(aes(group = Var1), alpha = 0.1, color = "darkgray") + 
    geom_violin() + geom_sina(size = 1) + scale_color_manual(values = ann_colors$Status) +
    theme_minimal(base_size = 22) + labs(title = paste0("Cluster", i, "(CD14+)"), color = "Status") +
    xlab("") + ylab("") + theme(legend.position = "none")
  ggsave(paste0(Clusters, "_Clusters/CD14_Severe_Healhty_peaks_cluster",i,"_boxplot.pdf"), width = 6, height = 5)
}

