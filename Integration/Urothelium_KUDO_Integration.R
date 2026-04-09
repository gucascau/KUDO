#!/usr/bin/env Rscript
library(Seurat)
library(gridExtra)
library(ggplot2)
library(tidyverse)
library(data.table)
#install.packages('R.utils')
library(R.utils)
library(reshape2)
library(stringr)
library(harmony)
library(cowplot)

library(patchwork)
#library(decoupleR)

# Extra libraries
library(dplyr)
library(pheatmap)

library(patchwork)
# we use v3 instead of v5
#options(Seurat.object.assay.version = "v5")

options(future.globals.maxSize = 1e12)



sessionInfo()
setwd("/home/gdjacksonlab/lab/xxw004/UUO/Results/KUDOs/Integration/IntegratingUUO/")

UUODir<- c("/home/gdjacksonlab/lab/xxw004/UUO/Datasets/")
KUDODir <- c("/home/gdjacksonlab/lab/xxw004/UUO/Results/KUDOs/Integration/DoubletFinder/")

scrna <-readRDS(file = "KudoUUOComparision_Scrna_PCA.rds")
# then use the Harmony to remove batches
scrna <- RunHarmony(scrna, group.by.vars = 'sample', dims.use = 1:20, plot_convergence= FALSE, reduction = "pca", assay.use = "SCT", reduction.save = "harmony")
#saveRDS(scrna, file = "KudoUUOComparision_Scrna_Harmony.rds")
scrna <- RunUMAP(scrna, dims = 1:20,reduction = "harmony")
#saveRDS(scrna, file = "KudoUUOComparision_Scrna_RunUMAP.rds")
scrna <- FindNeighbors(scrna, dims = 1:20, reduction = "harmony")
#saveRDS(scrna, file = "KudoUUOComparision_Scrna_FindNeighbor.rds")
scrna <- FindClusters(scrna, resolution = 0.2)

DimPlot(scrna, group.by = "sample")

saveRDS(scrna, file = "KudoUUOComparision_Scrna.rds")
