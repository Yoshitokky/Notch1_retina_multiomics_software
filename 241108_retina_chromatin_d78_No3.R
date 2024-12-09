# GSE183684

library(Signac)
library(Seurat)
library(GenomicRanges)
library(ggplot2)
library(patchwork)

set.seed(1)

multiome_d78 <- readRDS("241108_retina_chromatin_d78_No2.rds")

multiome_d78 <- RunTFIDF(multiome_d78)
multiome_d78 <- FindTopFeatures(multiome_d78, min.cutoff = 'q0')
multiome_d78 <- RunSVD(multiome_d78)

DepthCor(multiome_d78)

multiome_d78 <- RunUMAP(object = multiome_d78, reduction = 'lsi', dims = 2:30)
multiome_d78 <- FindNeighbors(object = multiome_d78, reduction = 'lsi', dims = 2:30)
multiome_d78 <- FindClusters(object = multiome_d78, verbose = FALSE, algorithm = 3)
DimPlot(object = multiome_d78, label = TRUE) + NoLegend()

gene.activities <- GeneActivity(multiome_d78)

# add the gene activity matrix to the Seurat object as a new assay and normalize it
multiome_d78[['RNA']] <- CreateAssayObject(counts = gene.activities)
multiome_d78 <- NormalizeData(
  object = multiome_d78,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(multiome_d78$nCount_RNA)
)

DefaultAssay(multiome_d78) <- 'RNA'

FeaturePlot(
  object = multiome_d78,
  features = c('NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)

# RNAseq

RNA_d78 <- readRDS("241108_GSE183684_D78.rds")

new.cluster.ids <- c("RPC1", "RPC2", "RPC3", 
                     "ELAVL2+/4+ RPC1", "ELAVL2+/4+ RPC2",
                     "ONECUT1+/2+/LHX1+ RPC", "ONECUT1+/2+ RPC1",
                     "ELAVL2+/4+ RPC3", "ONECUT1+/MEIS2+ RPC1",
                     "ONECUT1+/2+/PTF1A+ RPC", "ONECUT1+/MEIS2+ RPC1", 
                     "ONECUT1+/MEIS2+ RPC2", "ONECUT1+/MEIS2+ RPC3",
                     "ONECUT1+/2+ RPC2", "Proliferating RPC1", "Proliferating RPC2",
                     "VSX2++ RPC", "EEF1A1+ RPC")
names(new.cluster.ids) <- levels(RNA_d78)
RNA_d78 <- RenameIdents(RNA_d78, new.cluster.ids)
RNA_d78$celltype <- Idents(RNA_d78)

DimPlot(RNA_d78, reduction = "umap", label = TRUE, repel = TRUE)

RNA_d78 <- UpdateSeuratObject(RNA_d78)

transfer.anchors <- FindTransferAnchors(
  reference = RNA_d78,
  query = multiome_d78,
  reduction = 'cca'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = RNA_d78$celltype,
  weight.reduction = multiome_d78[['lsi']],
  dims = 2:30
)

multiome_d78 <- AddMetaData(object = multiome_d78, metadata = predicted.labels)

cluster_colors <- c(
  "RPC1" = "#1f77b4", 
  "RPC2" = "#ff7f0e", 
  "RPC3" = "#2ca02c", 
  "ELAVL2+/4+ RPC1" = "#d62728", 
  "ELAVL2+/4+ RPC2" = "#9467bd", 
  "ONECUT1+/2+/LHX1+ RPC" = "#8c564b", 
  "ONECUT1+/2+ RPC1" = "#e377c2", 
  "ELAVL2+/4+ RPC3" = "#7f7f7f", 
  "ONECUT1+/MEIS2+ RPC1" = "#bcbd22", 
  "ONECUT1+/2+/PTF1A+ RPC" = "#17becf", 
  "ONECUT1+/MEIS2+ RPC1" = "#f44336", 
  "ONECUT1+/MEIS2+ RPC2" = "#2196f3", 
  "ONECUT1+/MEIS2+ RPC3" = "#009688", 
  "ONECUT1+/2+ RPC2" = "#FF1493", 
  "Proliferating RPC1" = "#ff69b4", 
  "Proliferating RPC2" = "#1f77b4", 
  "VSX2++ RPC" = "#ff9800", 
  "EEF1A1+ RPC" = "#8a2be2"
)

plot1 <- DimPlot(
  object = RNA_d78,
  group.by = 'celltype',
  label = TRUE,
  repel = TRUE,
  cols = cluster_colors) + NoLegend() + ggtitle('scRNA-seq') + 
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20))

plot2 <- DimPlot(
  object = multiome_d78,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE,
  cols = cluster_colors) + NoLegend() + ggtitle('scATAC-seq') + 
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20))

plot1
plot2

saveRDS(object = RNA_d78, file = "241108_retina_chromatin_d78_RNAseq_No3.rds")
saveRDS(object = multiome_d78, file = "241108_retina_chromatin_d78_ATACseq_No3.rds")