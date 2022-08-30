# 05_manual_annot

## Objectives

1.  Manual Annotation via Seurat

### load data and make seurat object

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(Seurat)
```

    ## Attaching SeuratObject

    ## Attaching sp

``` r
library(patchwork)
library(ggplot2)


source("../tools/spongy_panda/export_gdcmatrix.R")

data.id = "gw10"

data <- readRDS("../../data/gse165388_processed/gw10_seuratobject.rds")
```

## Data scaling (for PCA)

``` r
all.genes <- rownames(data)
data <- ScaleData(data, features = all.genes)
```

    ## Centering and scaling data matrix

``` r
data <- FindVariableFeatures(object = data)

feat <- VariableFeatures(object = data)
```

## Export

``` r
save_gdcmatrix(
  GetAssayData(data)[feat, ],
  file = paste0("./", data.id, "feat_matrix")
)

saveRDS(
  data,
  file = paste0("./", data.id, "feat_seurat_pbj.rds")
)
```