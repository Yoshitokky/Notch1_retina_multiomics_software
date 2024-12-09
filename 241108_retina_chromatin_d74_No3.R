# GSE183684

library(Signac)
library(Seurat)
library(GenomicRanges)
library(ggplot2)
library(patchwork)

set.seed(1)

multiome_d74 <- readRDS("241108_retina_chromatin_d74_No2.rds")

multiome_d74 <- RunTFIDF(multiome_d74)
multiome_d74 <- FindTopFeatures(multiome_d74, min.cutoff = 'q0')
multiome_d74 <- RunSVD(multiome_d74)

DepthCor(multiome_d74)

multiome_d74 <- RunUMAP(object = multiome_d74, reduction = 'lsi', dims = 2:30)
multiome_d74 <- FindNeighbors(object = multiome_d74, reduction = 'lsi', dims = 2:30)
multiome_d74 <- FindClusters(object = multiome_d74, verbose = FALSE, algorithm = 3)
DimPlot(object = multiome_d74, label = TRUE) + NoLegend()

gene.activities <- GeneActivity(multiome_d74)

# add the gene activity matrix to the Seurat object as a new assay and normalize it
multiome_d74[['RNA']] <- CreateAssayObject(counts = gene.activities)
multiome_d74 <- NormalizeData(
  object = multiome_d74,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(multiome_d74$nCount_RNA)
)

DefaultAssay(multiome_d74) <- 'RNA'

FeaturePlot(
  object = multiome_d74,
  features = c('NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)

# RNAseq

RNA_d74 <- readRDS("241106_GSE183684_D74.rds")

new.cluster.ids <- c("RPC1", "ELAVL2+/4+ RPC1", "RPC2", 
                     "RPC3", "ELAVL2+/4+ RPC2",
                     "ONECUT1+/2+/LHX1+/PTF1A+ RPC1", 
                     "ELAVL2+/4+ RPC3", "ONECUT1+/MEIS2+ RPC1",
                     "RPC4", "ONECUT2+ RPC", "ONECUT1+/2+/LHX1+/PTF1A+ RPC2", 
                     "RPC5", "ONECUT1+/MEIS2+ RPC2", "VSX2++ RPC", "RPC6")
names(new.cluster.ids) <- levels(RNA_d74)
RNA_d74 <- RenameIdents(RNA_d74, new.cluster.ids)
RNA_d74$celltype <- Idents(RNA_d74)

DimPlot(RNA_d74, reduction = "umap", label = TRUE, repel = TRUE)

RNA_d74 <- UpdateSeuratObject(RNA_d74)

transfer.anchors <- FindTransferAnchors(
  reference = RNA_d74,
  query = multiome_d74,
  reduction = 'cca'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = RNA_d74$celltype,
  weight.reduction = multiome_d74[['lsi']],
  dims = 2:30
)

multiome_d74 <- AddMetaData(object = multiome_d74, metadata = predicted.labels)

# RNA_d74 <- readRDS("241108_retina_chromatin_d74_RNAseq_No3.rds")
# multiome_d74 <- readRDS("241108_retina_chromatin_d74_ATACseq_No3.rds")

cluster_colors <- c(
  "RPC1" = "#1f77b4", 
  "ELAVL2+/4+ RPC1" = "#ff7f0e", 
  "RPC2" = "#2ca02c", 
  "RPC3" = "#d62728", 
  "ELAVL2+/4+ RPC2" = "#9467bd", 
  "ONECUT1+/2+/LHX1+/PTF1A+ RPC1" = "#8c564b", 
  "ELAVL2+/4+ RPC3" = "#e377c2", 
  "ONECUT1+/MEIS2+ RPC1" = "#7f7f7f", 
  "RPC4" = "#bcbd22", 
  "ONECUT2+ RPC" = "#17becf", 
  "ONECUT1+/2+/LHX1+/PTF1A+ RPC2" = "#f44336", 
  "RPC5" = "#2196f3", 
  "ONECUT1+/MEIS2+ RPC2" = "#009688", 
  "VSX2++ RPC" = "#FF1493", 
  "RPC6" = "#ff69b4"
)

plot1 <- DimPlot(
  object = RNA_d74,
  group.by = 'celltype',
  label = TRUE,
  repel = TRUE,
  cols = cluster_colors) + NoLegend() + ggtitle('scRNA-seq') + 
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20))


plot2 <- DimPlot(
  object = multiome_d74,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE,
  cols = cluster_colors) + NoLegend() + ggtitle('scATAC-seq') + 
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20))

plot1
plot2

saveRDS(object = RNA_d74, file = "241108_retina_chromatin_d74_RNAseq_No3.rds")
saveRDS(object = multiome_d74, file = "241108_retina_chromatin_d74_ATACseq_No3.rds")