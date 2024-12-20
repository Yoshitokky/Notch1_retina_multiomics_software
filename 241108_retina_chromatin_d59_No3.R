# GSE183684

library(Signac)
library(Seurat)
library(GenomicRanges)
library(ggplot2)
library(patchwork)

set.seed(1)

multiome_d59 <- readRDS("241108_retina_chromatin_d59_No2.rds")

multiome_d59 <- RunTFIDF(multiome_d59)
multiome_d59 <- FindTopFeatures(multiome_d59, min.cutoff = 'q0')
multiome_d59 <- RunSVD(multiome_d59)

DepthCor(multiome_d59)

multiome_d59 <- RunUMAP(object = multiome_d59, reduction = 'lsi', dims = 2:30)
multiome_d59 <- FindNeighbors(object = multiome_d59, reduction = 'lsi', dims = 2:30)
multiome_d59 <- FindClusters(object = multiome_d59, verbose = FALSE, algorithm = 3)
DimPlot(object = multiome_d59, label = TRUE) + NoLegend()

gene.activities <- GeneActivity(multiome_d59)

# add the gene activity matrix to the Seurat object as a new assay and normalize it
multiome_d59[['RNA']] <- CreateAssayObject(counts = gene.activities)
multiome_d59 <- NormalizeData(
  object = multiome_d59,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(multiome_d59$nCount_RNA)
)

DefaultAssay(multiome_d59) <- 'RNA'

FeaturePlot(
  object = multiome_d59,
  features = c('NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)

# RNAseq

RNA_d59 <- readRDS("241108_GSE183684_D59_No0.rds")

new.cluster.ids <- c("EEF1A1+ RPC", "RPC1", "ONECUT1+/2+ RPC1", 
                     "ELAVL2+/4+ RPC1", "Proliferating RPC", "ELAVL2+/4+ RPC2", 
                     "ELAVL2+/4+ RPC3", "ONECUT1+/2+ RPC1",
                     "Transient RPC", "ELAVL2+/4+ RPC4", "RPC2", "RPC3", 
                     "ONECUT1+/2+ RPC2", "Glial cell")
names(new.cluster.ids) <- levels(RNA_d59)
RNA_d59 <- RenameIdents(RNA_d59, new.cluster.ids)
RNA_d59$celltype <- Idents(RNA_d59)

DimPlot(RNA_d59, reduction = "umap", label = TRUE, repel = TRUE)

RNA_d59 <- UpdateSeuratObject(RNA_d59)

transfer.anchors <- FindTransferAnchors(
  reference = RNA_d59,
  query = multiome_d59,
  reduction = 'cca'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = RNA_d59$celltype,
  weight.reduction = multiome_d59[['lsi']],
  dims = 2:30
)

multiome_d59 <- AddMetaData(object = multiome_d59, metadata = predicted.labels)

cluster_colors <- c(
  "EEF1A1+ RPC" = "#1f77b4", 
  "RPC1" = "#ff7f0e", 
  "ONECUT1+/2+ RPC1" = "#2ca02c", 
  "ELAVL2+/4+ RPC1" = "#d62728", 
  "Proliferating RPC" = "#9467bd", 
  "ELAVL2+/4+ RPC2" = "#8c564b", 
  "ELAVL2+/4+ RPC3" = "#e377c2", 
  "ONECUT1+/2+ RPC1" = "#7f7f7f", 
  "Transient RPC" = "#bcbd22", 
  "ELAVL2+/4+ RPC4" = "#17becf", 
  "RPC2" = "#f44336", 
  "RPC3" = "#2196f3", 
  "ONECUT1+/2+ RPC2" = "#009688", 
  "Glial cell" = "#FF1493"
)

plot1 <- DimPlot(
  object = RNA_d59,
  group.by = 'celltype',
  label = TRUE,
  repel = TRUE,
  cols = cluster_colors) + NoLegend() + ggtitle('scRNA-seq') + 
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20))

plot2 <- DimPlot(
  object = multiome_d59,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE,
  cols = cluster_colors) + NoLegend() + ggtitle('scATAC-seq') + 
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20))

plot1
plot2

saveRDS(object = RNA_d59, file = "241108_retina_chromatin_d59_RNAseq_No3.rds")
saveRDS(object = multiome_d59, file = "241108_retina_chromatin_d59_ATACseq_No3.rds")