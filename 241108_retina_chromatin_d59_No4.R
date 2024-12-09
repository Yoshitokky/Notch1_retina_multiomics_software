# GSE183684

library(Signac)
library(Seurat)
library(GenomicRanges)
library(ggplot2)
library(patchwork)

# RNA_d59 <- readRDS("241102_retina_chromatin_d59_RNAseq_No3.rds")
multiome_d59 <- readRDS("241108_retina_chromatin_d59_ATACseq_No3.rds")

# change back to working with peaks instead of gene activities
DefaultAssay(multiome_d59) <- 'peaks'

Idents(object = multiome_d59) <- "predicted.id"

multiome_d59 <- SortIdents(multiome_d59)

CoveragePlot(
  object = multiome_d59,
  region = "NOTCH1",
  extend.upstream = 1000,
  extend.downstream = 1000
)

CoveragePlot(
  object = multiome_d59,
  region = "NOTCH2",
  extend.upstream = 1000,
  extend.downstream = 1000
)

CoveragePlot(
  object = multiome_d59,
  region = "NOTCH3",
  extend.upstream = 1000,
  extend.downstream = 1000
)

CoveragePlot(
  object = multiome_d59,
  region = "NOTCH4",
  extend.upstream = 1000,
  extend.downstream = 1000
)