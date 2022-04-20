setwd("~/path/to/output")
# Hypergeometric P-value Enrichment Analyses

# Required libraries
library(dplyr)
library(purrr)
library(cinaR)
library(cinaRgenesets)
library(ggplot2)
library(plyr)
#Load HPEA
source("~/path/to/HPEA.R")

## The following steps are for running the Wikipathways/ImmuneModules dataset. Please comment out everything until Fold change and p-value cutoffs and change geneset name accordingly 
## Change name to vp2008 to run for for Immune Modules.
## Change name to wp to run for for Wikipathways.
load("~/path/to/enrichment_analysis.Rdata")
geneset.name <- "wp"

# Selected modules 
# There are multiple genesets consisting of multiple modules.
# We select one geneset per test and correct the p-values for it.
geneset <- selected_genesets[[geneset.name]]

# Module annotations
geneset.label <- selected_genesets_labels[[geneset.name]]

# Merge annotations via modules
geneset.merged <- merge(geneset, geneset.label, by = "Module.ID")

## The following steps are for running the Activated Immune dataset. Please comment out after loading the HPEA.R script until Fold change and p-value cutoffs  and change geneset name accordingly 
##ACTIVATED IMMUNE
# Load genesets
load("~/DiffBind/Peaks/DiffBind_Results/activated.immune.rda")

# Inflammation module from Virginia Pasqual (2008)
geneset.name <- "activated_immune"

A1 <- as.data.frame(activated.immune$Monocyte_LPS)
colnames(A1) <- "GeneName"
A1$Module.ID <- "Module_1"
A1$Module.Name <- "Monocyte_LPS"

A2 <- as.data.frame(activated.immune$Bcell_AntiCD3)
colnames(A2) <- "GeneName"
A2$Module.ID <- "Module_2"
A2$Module.Name <- "Bcell_AntiCD3"

A3 <- as.data.frame(activated.immune$NK_AntiCD3)
colnames(A3) <- "GeneName"
A3$Module.ID <- "Module_3"
A3$Module.Name <- "NK_AntiCD3"

A4 <- as.data.frame(activated.immune$CD8T_Naive_AntiCD3)
colnames(A4) <- "GeneName"
A4$Module.ID <- "Module_4"
A4$Module.Name <- "CD8T_Naive_AntiCD3"

A5 <- as.data.frame(activated.immune$CD8T_Mem_AntiCD3)
colnames(A5) <- "GeneName"
A5$Module.ID <- "Module_5"
A5$Module.Name <- "CD8T_Mem_AntiCD3"

A6 <- as.data.frame(activated.immune$CD4T_Naive_AntiCD3)
colnames(A6) <- "GeneName"
A6$Module.ID <- "Module_6"
A6$Module.Name <- "CD4T_Naive_AntiCD3"

A7 <- as.data.frame(activated.immune$CD4T_Mem_AntiCD3)
colnames(A7) <- "GeneName"
A7$Module.ID <- "Module_7"
A7$Module.Name <- "CD4T_Mem_AntiCD3"

geneset.merged <- rbind(A1, A2, A3, A4, A5, A6, A7)
rm(A1)
rm(A2)
rm(A3)
rm(A4)
rm(A5)
rm(A6)
rm(A7)

## The following steps are for running the Immune Wikipathways dataset. Please comment out after loading the HPEA.R script until Fold change and p-value cutoffs  and change geneset name accordingly 
#WIKI IMMUNE
results <- wiki.immune
R <- as.data.frame(as.matrix(results))
class(R)
R$Module.Name <- rownames(R)
R$Module.ID <- "Module_1"
for (i in 1:nrow(R)){
  R$Module.ID[i] = paste("Module",i, sep = "_")  
}
R <- unnest(R, V1)
geneset.merged <- as.data.frame(R)
colnames(geneset.merged) <- c("GeneName", "Module.Name", "Module.ID")
geneset.name <- "wiki.immune"
rm(R)
rm(results)

## The following steps are for running the Interferon dataset. Please comment out after loading the HPEA.R script until Fold change and p-value cutoffs  and change geneset name accordingly 
##INTERFERON
R <- read.table(file = "~/DiffBind/Peaks/DiffBind_Results/Interferon.txt", sep = "\t", header = F)
colnames(R) <- c("GeneName", "Module.Name")

R$Module.ID <- "Module_1"
for (i in 1:nrow(R)){
  if (R$Module.Name[i] == "New_IFNa-induced") {
    R$Module.ID[i] = "Module_2"
  }
  if (R$Module.Name[i] == "New_IFNg-induced") {
    R$Module.ID[i] = "Module_3"
  }
  if (R$Module.Name[i] == "Type_I") {
    R$Module.ID[i] = "Module_4"
  }
  if (R$Module.Name[i] == "Type_I_GAMMA") {
    R$Module.ID[i] = "Module_5"
  }
}
geneset.merged <- as.data.frame(R)
geneset.name <- "interferon"


# Fold change and p-value cutoffs
cutoff.fc <- 0
cutoff.p  <- 0.1
cutoff.tss<- 50e3

# Contrast file names
filelist = list.files(pattern = "*txt") 
filelist 

contrast.list <- (filelist) %>% strsplit("_", fixed = T) %>% 
  sapply(function(x){paste(x[1:4], collapse = "_")}) %>% 
  strsplit(".", fixed = T) %>% sapply(function(x){x[1]})
contrast.list


## This function reads in all the text files present in the given directory. Please save the text files accordingly in one place and name them by the conditions the differential analysis was run in.
## Our files were saved as CellType_Condition1_Condition2.txt
## Update the contrast.list function above if the filenames are saved differently.

HPEA.list <- lapply(filelist, function(file){
  
  contrast <- file %>% 
    strsplit("_", fixed = T) %>% 
    sapply(function(x){paste(x[1:2], collapse = "_")})
  
  # read the DE results
  peaks <- read.csv(file, sep = "\t")
  peaks$geneSymbol <- peaks$gene_name
  # filter out the genes bigger than tss cutoff
  peaks.filtered <- peaks[abs(peaks$distanceToTSS) < cutoff.tss,]
  
  # order according to tss
  peaks.filtered <- peaks.filtered %>% arrange(abs(distanceToTSS))
  
  # remove duplicated gene names
  peaks.filtered <- peaks.filtered[!duplicated(peaks.filtered$geneSymbol),]
  
  #peaks.filtered$logFC <- peaks.filtered$avg_logFC
  # For RNA-seq analyses you need up and down regulated genes (and opening/closing peaks for ATACseq).
  # The sign of this regulation is determined via fold-change (FC) or log(FC) of the genes, 
  # which are calculated with differential analyses. Then, these genes are filtered further
  # with respect to their p-values (or adjusted p-values). 
  
  # DE genes up-regulated
  DE.genes.ur <- peaks.filtered[peaks.filtered$logFC > cutoff.fc & peaks.filtered$FDR < cutoff.p,]
  # DE genes down-regulated
  DE.genes.dr <- peaks.filtered[peaks.filtered$logFC < cutoff.fc & peaks.filtered$FDR < cutoff.p,]
  
  # If both up and down regulated genes do not exist skip the contrast
  if(nrow(DE.genes.ur) > 0 | nrow(DE.genes.dr) > 0){
    
    gene.names.ur <- DE.genes.ur$geneSymbol
    gene.names.dr <- DE.genes.dr$geneSymbol
    
    result.ur <- HPEA(geneset.merged, gene.names.ur, geneset.name = geneset.name, contrast = contrast)
    result.dr <- HPEA(geneset.merged, gene.names.dr, geneset.name = geneset.name, contrast = contrast)
    
    return(list(UR = result.ur, DR = result.dr))
  }
  return(NULL)
  
})

names(HPEA.list) <- contrast.list

# Merge UR/DR lists
HPEA.list <- lapply(HPEA.list, function(x){
  if (!is.null(x)){
    rbind(
      cbind(x[["UR"]], Status = "Up"), 
      cbind(x[["DR"]], Status = "Down")
    )
  }
})

# Make it a table
HPEA.table <- map_df(HPEA.list, ~as.data.frame(.x), .id="contrast")

# Remove empty rows
HPEA.table.filtered <- HPEA.table[HPEA.table$p < 1,]
write.csv(HPEA.table.filtered, file = paste0("./", geneset.name, "_hpea-list.csv"))
df.plot <- HPEA.table.filtered


filter.pathways <- TRUE
fdr.cutoff <- 0.1

if(filter.pathways){
  if (sum(df.plot$adj.p < fdr.cutoff) == 0){
    stop("You can't filter because there are no pathways to be displayed!")
  }
  df.plot <- subset(df.plot, adj.p < fdr.cutoff)
}

df.plot$grouping = "Healthy"

## Make sure the contrast of the loop matches the labels of the contrasts in the text files used in the differential analysis.
## If the names of the contrasts are different, please rename the loop conditions accordingly. Refer to df.plot dataframe when naming the conditions.

for (i in 1:nrow(df.plot)) {
  if (df.plot$contrast[i] == "Sev_Mild" || df.plot$contrast[i] == "Sev_F_Mild" || df.plot$contrast[i] == "Mild_F_Mild") {
    df.plot$grouping[i] <- "Mild EC"
  }
  if (df.plot$contrast[i] == "Mild_ICU_Cont" || df.plot$contrast[i] == "Mild_F_ICU_Cont" || df.plot$contrast[i] == "Sev_ICU_Cont" || df.plot$contrast[i] == "Sev_F_ICU_Cont") {
    df.plot$grouping[i] <- "ICU Control"
  }
  if (df.plot$contrast[i] == "Sev_Mild_F" || df.plot$contrast[i] == "Sev_F_Mild_F") {
    df.plot$grouping[i] <- "Mild LC"
  }
  if (df.plot$contrast[i] == "Sev_F_Sev") {
    df.plot$grouping[i] <- "Severe EC"
  }
}
## If the names of the contrasts are different, please rename the factor conditions accordingly. Refer to df.plot dataframe when naming the conditions.
df.plot$contrast <- factor(df.plot$contrast, levels=c("Mild_HD", "Mild_F_HD", "Sev_HD", "Sev_F_HD", "ICU_Cont_HD", "Mild_ICU_Cont", "Mild_F_ICU_Cont", "Sev_ICU_Cont", "Sev_F_ICU_Cont", "Sev_F_Sev", "Mild_F_Mild", "Sev_Mild", "Sev_Mild_F", "Sev_F_Mild", "Sev_F_Mild_F"))

df.plot$grouping <- factor(df.plot$grouping, levels=c("Healthy", "ICU Control", "Mild EC", "Mild LC", "Severe EC", "Severe LC"))
df.plot$Status <- factor(df.plot$Status, levels=c("Down", "Up"))

# create ggplot
plot.dot <- ggplot2::ggplot(df.plot,
                            ggplot2::aes(x = contrast,
                                         y = module.name,
                                         size = ifelse(adj.p < fdr.cutoff, -log(adj.p), NA),
                                         color = Status))

plot.dot <- plot.dot + ggplot2::geom_point()

plot.dot <- plot.dot + facet_grid(~grouping,scales='free', space = "free")
plot.dot <- plot.dot + ggplot2::labs(x = "Contrast",
                                     y = "Pathways",
                                     color = "Sign",
                                     size = "-log10(adj.p)",
                                     caption = paste0("FDR < ", fdr.cutoff))
color_values <- color_values[c(4,7)]
plot.dot <- plot.dot + ggplot2::scale_color_manual(values = color_values)

pdf("~/path/to/Dotplot.pdf", height = 40, width = 15)
plot.dot
dev.off()