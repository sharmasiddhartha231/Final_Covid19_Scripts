## COpy each individual object from Seurat_GEX.R in one directory.

setwd("~/path/to/GEX_Objects")
library(Seurat)

HDu <- readRDS("./HDu.rds")
HDa <- readRDS("./HDa.rds")
HDv <- readRDS("./HDv.rds")
HDj <- readRDS("./HDj.rds")
HD_AD <- readRDS("./HD_AD.rds")

jcov101 <- readRDS("./jcov101.rds")
jcov101_2 <- readRDS("./jcov101-2.rds")
jcov110 <- readRDS("./jcov110.rds")
jcov111 <- readRDS("./jcov111.rds")
jcov114 <- readRDS("./jcov114.rds")
jcov115 <- readRDS("./jcov115.rds")
jcov117 <- readRDS("./jcov117.rds")
jcov119 <- readRDS("./jcov119.rds")
jcov122 <- readRDS("./jcov122.rds")
jcov124 <- readRDS("./jcov124.rds")
jcov127 <- readRDS("./jcov127.rds")

jcov21 <- readRDS("./jcov21.rds")
jcov42 <- readRDS("./jcov42.rds")
jcov42_20 <- readRDS("./jcov42-20.rds")
jcov49_2 <- readRDS("./jcov49-2.rds")
jcov52 <- readRDS("./jcov52.rds")
jcov56_9 <- readRDS("./jcov56-9.rds")
jcov60_6 <- readRDS("./jcov60-6.rds")
jcov69_8 <- readRDS("./jcov69-8.rds")
jcov72_5 <- readRDS("./jcov72-5.rds")
jcov92 <- readRDS("./jcov92.rds")
jcov93 <- readRDS("./jcov93.rds")

lgtd17 <- readRDS("./lgtd17.rds")
lgtd18 <- readRDS("./lgtd18.rds")
lgtd8m <- readRDS("./lgtd8m.rds")

SeuratRNA <- merge(HDu, y = c(HDa, HDv, HDj, HD_AD, jcov101, jcov101_2, jcov110, jcov111, jcov114, jcov115, jcov117, jcov119, jcov122, jcov124, jcov127, jcov21, jcov42, jcov42_20, jcov49_2, jcov52, jcov56_9, jcov60_6, jcov69_8, jcov72_5, jcov92, jcov93, lgtd17, lgtd18, lgtd8m), add.cell.ids = c("HDu", "HDa", "HDv", "HDj", "HD_AD", "jcov101", "jcov101_2", "jcov110", "jcov111", "jcov114", "jcov115", "jcov11
7", "jcov119", "jcov122", "jcov124", "jcov127", "jcov21", "jcov42", "jcov42_20", "jcov49_2", "jcov52", "jcov56_9", "jcov60_6", "jcov69_8", "jcov72_5", "jcov92", "jcov93", "lgtd17", "lgtd18", "lgtd8m"), project = "Covid_19")

## Save Final Unprocessed Object for RNA
saveRDS(SeuratRNA, file = "./PooledSeuratRNA.rds")
