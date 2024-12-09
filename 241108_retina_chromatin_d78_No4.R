# GSE183684

library(Signac)
library(Seurat)
library(GenomicRanges)
library(ggplot2)
library(patchwork)

# RNA_d78 <- readRDS("241108_retina_chromatin_d78_RNAseq_No3.rds")
multiome_d78 <- readRDS("241108_retina_chromatin_d78_ATACseq_No3.rds")

# change back to working with peaks instead of gene activities
DefaultAssay(multiome_d78) <- 'peaks'

Idents(object = multiome_d78) <- "predicted.id"

multiome_d78 <- SortIdents(multiome_d78)

CoveragePlot(
  object = multiome_d78,
  region = "NOTCH1",
  extend.upstream = 1000,
  extend.downstream = 1000
)

CoveragePlot(
  object = multiome_d78,
  region = "NOTCH2",
  extend.upstream = 1000,
  extend.downstream = 1000
)

CoveragePlot(
  object = multiome_d78,
  region = "NOTCH3",
  extend.upstream = 1000,
  extend.downstream = 1000
)

CoveragePlot(
  object = multiome_d78,
  region = "NOTCH4",
  extend.upstream = 1000,
  extend.downstream = 1000
)

# sessionInfo()
# R version 4.4.1 (2024-06-14)
# Platform: x86_64-pc-linux-gnu
# Running under: Ubuntu 22.04.4 LTS
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.10.0 
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
# [3] LC_TIME=ja_JP.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=ja_JP.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=ja_JP.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=ja_JP.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: Asia/Tokyo
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] patchwork_1.3.0      ggplot2_3.5.1        GenomicRanges_1.56.2
# [4] GenomeInfoDb_1.40.1  IRanges_2.38.1       S4Vectors_0.42.1    
# [7] BiocGenerics_0.50.0  Seurat_5.1.0         SeuratObject_5.0.2  
# [10] sp_2.1-4             Signac_1.14.0       
# 
# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3      rstudioapi_0.17.1       jsonlite_1.8.9         
# [4] magrittr_2.0.3          spatstat.utils_3.1-0    farver_2.1.2           
# [7] zlibbioc_1.50.0         vctrs_0.6.5             ROCR_1.0-11            
# [10] spatstat.explore_3.3-3  Rsamtools_2.20.0        RcppRoll_0.3.1         
# [13] htmltools_0.5.8.1       sctransform_0.4.1       parallelly_1.38.0      
# [16] KernSmooth_2.23-24      htmlwidgets_1.6.4       ica_1.0-3              
# [19] plyr_1.8.9              plotly_4.10.4           zoo_1.8-12             
# [22] igraph_2.1.1            mime_0.12               lifecycle_1.0.4        
# [25] pkgconfig_2.0.3         Matrix_1.6-5            R6_2.5.1               
# [28] fastmap_1.2.0           GenomeInfoDbData_1.2.12 fitdistrplus_1.2-1     
# [31] future_1.34.0           shiny_1.9.1             digest_0.6.37          
# [34] colorspace_2.1-1        tensor_1.5              RSpectra_0.16-2        
# [37] irlba_2.3.5.1           labeling_0.4.3          progressr_0.15.0       
# [40] fansi_1.0.6             spatstat.sparse_3.1-0   httr_1.4.7             
# [43] polyclip_1.10-7         abind_1.4-8             compiler_4.4.1         
# [46] withr_3.0.2             BiocParallel_1.38.0     fastDummies_1.7.4      
# [49] MASS_7.3-61             tools_4.4.1             lmtest_0.9-40          
# [52] httpuv_1.6.15           future.apply_1.11.3     goftest_1.2-3          
# [55] glue_1.8.0              nlme_3.1-165            promises_1.3.0         
# [58] grid_4.4.1              Rtsne_0.17              cluster_2.1.6          
# [61] reshape2_1.4.4          generics_0.1.3          gtable_0.3.6           
# [64] spatstat.data_3.1-2     tidyr_1.3.1             data.table_1.16.2      
# [67] utf8_1.2.4              XVector_0.44.0          spatstat.geom_3.3-3    
# [70] RcppAnnoy_0.0.22        ggrepel_0.9.6           RANN_2.6.2             
# [73] pillar_1.9.0            stringr_1.5.1           spam_2.11-0            
# [76] RcppHNSW_0.6.0          later_1.3.2             splines_4.4.1          
# [79] dplyr_1.1.4             lattice_0.22-5          survival_3.7-0         
# [82] deldir_2.0-4            tidyselect_1.2.1        Biostrings_2.72.1      
# [85] miniUI_0.1.1.1          pbapply_1.7-2           gridExtra_2.3          
# [88] scattermore_1.2         matrixStats_1.4.1       stringi_1.8.4          
# [91] UCSC.utils_1.0.0        lazyeval_0.2.2          codetools_0.2-19       
# [94] tibble_3.2.1            cli_3.6.3               uwot_0.2.2             
# [97] xtable_1.8-4            reticulate_1.39.0       munsell_0.5.1          
# [100] Rcpp_1.0.13             globals_0.16.3          spatstat.random_3.3-2  
# [103] png_0.1-8               spatstat.univar_3.0-1   parallel_4.4.1         
# [106] dotCall64_1.2           bitops_1.0-9            listenv_0.9.1          
# [109] viridisLite_0.4.2       scales_1.3.0            ggridges_0.5.6         
# [112] leiden_0.4.3.1          purrr_1.0.2             crayon_1.5.3           
# [115] rlang_1.1.4             cowplot_1.1.3           fastmatch_1.1-4