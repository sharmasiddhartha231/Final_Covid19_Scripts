library(cinaR)
library(cinaRgenesets)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(tidyverse)
library(readxl)
setwd("~/path/to/output")

## Use consensus peaks obtained from Diffbind Script.
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
cinaR_results <- cinaR(Counts, contrast, reference.genome = "hg38", batch.correction = T, additional.covariates = age, DA.fdr.threshold = 0.1, enrichment.FDR.cutoff = 0.1)

## Save R object. 
saveRDS(cinaR_results, "DifferentialResults.rds")
