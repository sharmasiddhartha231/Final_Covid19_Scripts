library(cinaR)
library(cinaRgenesets)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(readxl)
setwd("~/path/to/output")

Metadata <- read_excel("~/path/to/Metadata.xlsx")

contrast <- Metadata$Contrast

## Type is whether data is bulkATAC or pseudobulk ATAC.
type <- Metadata$Type

## Use normalized consensus peaks obtained after running cinaR on the dataset.
## To run PCA on select differential peaks, load in differential peak data and use those peak counts instead of all the consensus peaks for PCA.
cp <- read.table("~/path/to/ConsensusPeaks.txt", sep = '\t', header = T)

##PCA
#silence CRAN build NOTES
PC1 <- PC2 <- NULL

cp <- cp 

sample.names = NULL
if(is.null(sample.names)){
  sample.names <- colnames(cp)
} else{
  if(length(sample.names) != ncol(cp)){
    stop("The length of `sample.names` should be equal to number of samples.")
  }
}
  
# eliminate NaN values before-hand if there is any.
pca <- stats::prcomp(t(stats::na.omit(cp)), center = TRUE)
  
d  <- round(pca$sdev^2/sum(pca$sdev^2)*100, digits=1)
xl <- sprintf("PC 1: %.1f %%", d[1])
yl <- sprintf("PC 2: %.1f %%", d[2])
  
overlaid.info = contrast
plot.df <- data.frame(PC1 = as.numeric(pca$x[,1]),
                      PC2 = as.numeric(pca$x[,2]),
                      overlaid.info = overlaid.info,
                      names = sample.names
)
  
plot.pca <- ggplot2::ggplot(plot.df, ggplot2::aes(PC1, PC2, color = overlaid.info, shape = factor(type))) +
  ggplot2::geom_point(size = 4, alpha = 1) +
  ggplot2::labs(x=xl,y=yl) +
  ggplot2::theme_minimal() +
  ggplot2::labs(color = "Status", shape = "Data Type") +
  ggplot2::coord_fixed(ratio = 1) +
  ggplot2::theme_light() + 
  ggplot2::scale_shape(solid = T)

#Color Scales based on clinical conditions.  
colr <- c("#27AAE1", "#26428B", "#008D48", "#61BE79", "#CE1D1E", "#F6BD31")

if (typeof(overlaid.info) %in% c("character", "factor")){
  plot.pca <- plot.pca +
    ggplot2::scale_color_manual(values = colr)
}

show.names = F  
if(show.names){
  plot.pca <- plot.pca + ggrepel::geom_text_repel(ggplot2::aes(label = names))
}

pdf("~/Path/to/PCA.pdf", height = 10, width = 10)
plot.pca  
dev.off()
