#!/usr/bin/env Rscript
# ---
# title: "Urothelium_KUDO_UrotheliumOrganoid_UUOandIRI.qmd"
# author: "Xin Wang"
# date: "2026-04-06"
# email: xin.wang@nationwidechildrens.org
# output: html_document

# Description:
#   Method Part1: Select the Urothelium cells in the object
#   Method reference: https://satijalab.org/seurat/articles/integration_introduction.html

# Dataset Description: 

#   uni-IRI and UUO on wild type adult male mice: which enthanized at multiple timepoints
#   as the disease progress (0.6 house and 2, 7. 14 28 days post uni-IRI or 0,2,4, 6, 10,14 days post UUO), n=2 for each timepoint)
#   GEO: GSE190887
#   GEO reference: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi
#   Reference: Li H, Dixon EE, Wu H, Humphreys BD. Comprehensive single-cell transcriptional profiling defines shared and unique epithelial injury responses during kidney fibrosis. Cell Metab 2022 Dec 6;34(12):1977-1998.e9.
  
#   Platform: Illumina NovaSeq 6000
  
#   snRNA-seq provied another strategy for performing single cell transcriptomics where individual nuclei instead of cells are captured and sequenced.
#    Advantages: snRNA-seq over scRNA-seq is that they don't need dissociation to preseve cell integrity during sample preparation.
#    Disadvantages: loss of transcripts are primarily located in the cytoplasm potentially limited the availiablity of biological signals for gene with little nuclear localization.


#    Urothelium organoid dataset: GSE131909_RAW
#  Reference: Santos, C.P., Lapi, E., Martínez de Villarreal, J. et al. Urothelial organoids originating from Cd49fhigh mouse stem cells display Notch-dependent differentiation capacity. Nat Commun 10, 4407 (2019). https://doi.org/10.1038/s41467-019-12307-1

# ---

# ```{r}
# module purge
# ml GCC/9.3.0
# ml OpenMPI/4.0.3
# ml R/4.4.0
# R
# capabilities()
# plot(1:20)
# ```

## Step 1: Load the packages
## Please install all the packages before loading them.
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
#library(decoupleR)

# Extra libraries
library(dplyr)
library(pheatmap)

library(patchwork)
# Here we use V5
#options(Seurat.object.assay.version = "v5")

options(future.globals.maxSize = 1e12)


# remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)
# remotes::install_github("satijalab/seurat-data", "seurat5", quiet = TRUE)
# remotes::install_github("satijalab/azimuth", "seurat5", quiet = TRUE)
# remotes::install_github("satijalab/seurat-wrappers", "seurat5", quiet = TRUE)
# remotes::install_github("stuart-lab/signac", "seurat5", quiet = TRUE)
# #remotes::install_github("bnprks/BPCells", quiet = TRUE)
# remotes::install_github("bnprks/BPCells", quiet = TRUE)
# # Enter commands in R (or R studio, if installed)




## Step 2: Set up the workspace environment:


# set up the working directory
Outdir<- c("/home/gdjacksonlab/lab/xxw004/UUO/Results/KUDOs/Integration/IntegratingUUO_Organoid/")

setwd("/home/gdjacksonlab/lab/xxw004/UUO/Results/KUDOs/Integration/IntegratingUUO_Organoid/")

# create the output directory folder for biomarkers, we will save the biomarkers in this folder
# 

BMdir<- paste0(Outdir,"Biomarkers/")
CAbdir<- paste0(Outdir,"CorrelationAnalysis/")
FEdir<- paste0(Outdir,"FunctionalEnrichment/")


dir.create(BMdir, showWarnings = F)
dir.create(CAbdir, showWarnings = F)
dir.create(FEdir, showWarnings = F)



# # set up the input datasets
# # 1.Urothelium organoid dataset: GSE131909_RAW
# # Ref Santos, C.P., Lapi, E., Martínez de Villarreal, J. et al. Urothelial organoids originating from Cd49fhigh mouse stem cells display Notch-dependent differentiation capacity. Nat Commun 10, 4407 (2019). https://doi.org/10.1038/s41467-019-12307-1
# # We have preprocessed the dataset and we will use the preprocessed dataset for the integration: remove doublets and perform the batch correction using harmony, and we will use the batch corrected dataset for the integration
# UroOrganoidInputDir <- c("/home/gdbecknelllab/xxw004/gdjacksonlab/UUO/Results/UroBladderOrganoid/IntegratedUroBladder/")

# # 2. UUO and IRI dataset: GSE190887_RAW
# # ref: Li H, Dixon EE, Wu H, Humphreys BD. Comprehensive single-cell transcriptional profiling defines shared and unique epithelial injury responses during kidney fibrosis. Cell Metab 2022 Dec 6;34(12):1977-1998.e9.
# UUODir<- c("/home/gdbecknelllab/xxw004/gdjacksonlab/UUO/Datasets/Mouse/UUO_CellMeta2022/")

# # 3. We also use the mouse kidney organoid dataset from KUDO as the reference for the comparison 
# # set up the correlation folders : scrna_merged_withdoublets.rds, there are four different condition, Yoda, LOF, GOF and Vehicle
# KUDODir <- c("/home/gdjacksonlab/lab/xxw004/UUO/Results/KUDOs/Integration/DoubletFinder/")


#####################################################################################
# Sample Integration
#####################################################################################

# 1. read the UUO and IRI dataset
## read the preprocessing
# UUOProjectObject_F<- readRDS(paste0(UUODir,"UUP_Urothelium_ObjectMatrixTest.rds"))

# MetaAnnotation<- read.csv(paste0(UUODir,"GSE190887_meta_cell_type_sample.csv"), header = T, row.names = 1)
# rownames(MetaAnnotation)<-str_replace_all(rownames(MetaAnnotation), "\\.","_")
# MetaAnnotation %>% head()
# UUOProjectObject_F@meta.data
# UUOProjectObject_F<-AddMetaData(object = UUOProjectObject_F, metadata = MetaAnnotation)
# ## we also add the meta data inside
# UUOProjectObject_F@meta.data %>% head()

# ## add the mitchondiron percentage
# #Test<- as.Seurat(UUOProject0,counts = "counts",project = "SingleCellExperiment",)
# UUOProjectObject_F[['percent.mt']] <- PercentageFeatureSet(UUOProjectObject_F, pattern = "^mt[-\\.]")

# ## add the ribosome genes
# UUOProjectObject_F[['percent.rRNA']] <- PercentageFeatureSet(UUOProjectObject_F, pattern = "^Rp[sl][[:digit:]]")

# Idents(UUOProjectObject_F) <- UUOProjectObject_F@meta.data$celltype0421

# # We only check the expression of Upk genes and Krt genes
# Idents(UUOProjectObject_F)

# UUOProjectObject_F@meta.data <- UUOProjectObject_F@meta.data %>% mutate(Group= case_when(grepl("UUO",sample) ~"UUO", 
#                                                  grepl("Health",sample) ~ "Health",
#                                                 grepl("IRI",sample) ~"IRI",
#                                                 .default = as.character(sample)
#                                                  ))
# # we only select the Urothelium that in Health and UUO model
# UUOProjectObject_F_UUOHealth <- subset (UUOProjectObject_F, subset =(Group == "Health" | Group== "UUO") )

# # We only select the Urothelium that in Health and UUO model, and we will use this dataset for the integration with KUDO dataset
# UUOProjectObject_F_UUOHealth <- subset (UUOProjectObject_F, subset =(Group == "Health" | Group== "UUO") & celltype0421 == "Uro" )


# saveRDS(UUOProjectObject_F_UUOHealth, file = "UUOProjectObject_F_UUOHealth_UrotheliumOnly.rds")
# #Urothelium <-subset(UUOProjectObject_F_UUOHealth, subset = Upk3a > 2 | Upk1a >2 |Upk1b >2 |Upk2 > 2 | Krt5 > 2 | Krt14 > 2 | Krt15 > 2 | Krt17 > 2 | Krt18 > 2 | Krt19 > 2)

# UUOProjectObject_F_UUOHealth<- readRDS(file = "UUOProjectObject_F_UUOHealth_UrotheliumOnly.rds")
# # 2. read the KUDO datasets
# KudoObj<- readRDS(paste0(KUDODir,"scrna_merged_withdoublets.rds"))

# KudoObj@meta.data %>% head()
# # we filtered the doublets in the KUDO dataset, and we will use the filtered dataset for the integration
# KudoObj_filtered <- subset(KudoObj, doublet_id == "Singlet")

# rm(KudoObj)

# # 3. read the Urothelium organoid dataset
# UroOrganoidObj <- readRDS(paste0(UroOrganoidInputDir,"UroBladderOrganoid_singlecell_doublet_harmony_v0618.RDS"))

# # we compared the KUDO with Mouse Urothelium
# # merge the Kudo Single cell and Mouse Urothelium
# # delete the objects:

# scrna <-
#   merge(
#     x = KudoObj_filtered,
#     y = c(UUOProjectObject_F_UUOHealth,UroOrganoidObj),
#     project = "KudoUUOUroBladderComparision"
#   )

# saveRDS(scrna, file = "KudoUUOComparision_Scrna_Merged.rds")

#####################################################################################
### Part 2: Preprocessing the merged dataset
#####################################################################################

scrna <- readRDS(file = "KudoUUOComparision_Scrna_Merged.rds")

# add the meta data for the UUO and KUDO
scrna@meta.data %>% head()
table(scrna@meta.data$sample)
table(scrna@meta.data$orig.ident)
# change the sample name NA as the KUDOVehicle
# add the sample and groups
scrna@meta.data <- scrna@meta.data %>% mutate(Sample = case_when(grepl("Urotherlium",orig.ident) ~ sample, 
                                 TRUE ~ as.character(orig.ident)), 
                                Group = case_when(grepl("GOF",orig.ident) ~ "KUDO", 
                                 grepl("LOF",orig.ident) ~ "KUDO",
                                 grepl("Yoda",orig.ident) ~ "KUDO",
                                 grepl("Vehicle",orig.ident) ~ "KUDO",
                                 grepl("Urotherlium",orig.ident) ~ "UUO",
                                 grepl("GSM3827175_NMU_O_P",orig.ident) ~ "UroBladderOrganoid",
                                 grepl("GSM3827176_NMU_O_D",orig.ident) ~ "UroBladderOrganoid",
                                TRUE ~ as.character(orig.ident)))

# check the number of cells in each sample(x, size, replace = FALSE, prob = NULL)
table((scrna$Group))
table((scrna$Sample))
length(scrna$Group)
nrow (scrna@meta.data)
scrna <- JoinLayers(scrna)  # ensure layers are joined before splitting
# we run out of the memories, therefore we split by sam[ple , SCTtransform each sample, and then merge them together, and then run the harmony for the batch correction
#scrna <- SplitObject(scrna, split.by = "Sample")
scrna[["RNA"]] <- split(scrna[["RNA"]], f = scrna$Sample)
scrna <- NormalizeData(scrna)
scrna <- FindVariableFeatures(scrna)
scrna <- ScaleData(scrna)
scrna <- RunPCA(scrna)

scrna <- FindNeighbors(scrna, dims = 1:30, reduction = "pca")
scrna <- FindClusters(scrna, resolution = 2, cluster.name = "unintegrated_clusters")
scrna <- RunUMAP(scrna, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
saveRDS(scrna, file = "KudoUUOComparision_Scrna_Merged_uumap.unintegrated.rds")


# integrating with harmony -- method 1: we use the harmony wrapper in Seurat, which is more memory efficient than the original harmony method, and we will use the PCA reduction for the integration, and then we will use the harmony reduction for the downstream analysis
scrna <- IntegrateLayers(
  object = scrna, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)
scrna <- IntegrateLayers(
  object = scrna, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE
)
scrna <- IntegrateLayers(
  object = scrna, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = FALSE
)

# scrna <- IntegrateLayers(
#   object = scrna, method = scVIIntegration,
#   new.reduction = "integrated.scvi",
#   conda_env = "/home/gdbecknelllab/xxw004/miniconda3/envs/scvi-env", verbose = FALSE
# )
scrna <- FindNeighbors(scrna, reduction = "harmony", dims = 1:30)
scrna <- FindClusters(scrna, resolution = 2, cluster.name = "harmony_clusters")
scrna <- RunUMAP(scrna, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
#scrna <- JoinLayers(scrna)
saveRDS(scrna, file = "KudoUUOComparision_Scrna_Merged_harmonyIntegration.rds")

# CCA integration method, which is the original method used in the Seurat v3, and we will use the PCA reduction for the integration, and then we will use the CCA reduction for the downstream analysis

scrna <- FindNeighbors(scrna, reduction = "integrated.cca", dims = 1:30)
scrna <- FindClusters(scrna, resolution = 2, cluster.name = "cca_clusters")
scrna <- RunUMAP(scrna, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")
#scrna <- JoinLayers(scrna)
saveRDS(scrna, file = "KudoUUOComparision_Scrna_Merged_ccaIntegration.rds")

# we used RPCAIntegration method which is more memory efficient than the original RPCA method, and we will use the PCA reduction for the integration, and then we will use the RPCA reduction for the downstream analysis

scrna <- FindNeighbors(scrna, reduction = "integrated.rpca", dims = 1:30)
scrna <- FindClusters(scrna, resolution = 2, cluster.name = "rpca_clusters")
scrna <- RunUMAP(scrna, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")
#scrna <- JoinLayers(scrna)
saveRDS(scrna, file = "KudoUUOComparision_Scrna_Merged_rpcaIntegration.rds")

# we also use mmn
scrna <- IntegrateLayers(
  object = scrna, method = FastMNNIntegration,
  new.reduction = "integrated.mnn",
  verbose = FALSE
)
scrna <- FindNeighbors(scrna, reduction = "integrated.mnn", dims = 1:30)
scrna <- FindClusters(scrna, resolution = 2, cluster.name = "mnn_clusters")
scrna <- RunUMAP(scrna, reduction = "integrated.mnn", dims = 1:30, reduction.name = "umap.mnn")
#scrna <- JoinLayers(scrna)
saveRDS(scrna, file = "KudoUUOComparision_Scrna_Merged_mnnIntegration.rds")



# scrna <- FindNeighbors(scrna, reduction = "integrated.scvi", dims = 1:30)
# scrna <- FindClusters(scrna, resolution = 2, cluster.name = "scvi_clusters")
# scrna <- RunUMAP(scrna, reduction = "integrated.scvi", dims = 1:30, reduction.name = "umap.scvi")
# scrna <- JoinLayers(scrna)
# saveRDS(scrna, file = "KudoUUOComparision_Scrna_Merged_scviIntegration.rds")

# # Step 2: SCTransform per sample (memory-safe)
# scrna <- lapply(scrna, function(x) {
#   SCTransform(
#     x,
#     vars.to.regress = "percent.mt",
#     method = "glmGamPoi",
#     verbose = FALSE
#   )
# })
# # Step 3: Integration (critical for your design)
# features <- SelectIntegrationFeatures(objs, nfeatures = 3000)

# objs <- PrepSCTIntegration(objs, anchor.features = features)

# anchors <- FindIntegrationAnchors(
#   objs,
#   normalization.method = "SCT",
#   anchor.features = features
# )

# scrna <- IntegrateData(
#   anchors,
#   normalization.method = "SCT"
# )


# # Step 4: PCA + Harmony (optional)
# scrna <- RunPCA(scrna, features = VariableFeatures(object = scrna),assay = "SCT",npcs=20)

# scrna <- RunHarmony(scrna, group.by.vars = 'Sample', dims.use = 1:20, plot_convergence= FALSE, reduction = "pca", assay.use = "SCT", reduction.save = "harmony")

# # step 5 downstream
# scrna <- RunUMAP(scrna, dims = 1:20, reduction = "harmony")
# scrna <- FindNeighbors(scrna, dims = 1:20, reduction = "harmony")
# scrna <- FindClusters(scrna, resolution = 0.2)

# #scrna <- scrna %>%  
# #  NormalizeData() %>% 
# #  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
# #  ScaleData() %>% 
# #  SCTransform(vars.to.regress = c("percent.mt"))

# ## find most vaible features across sample to integrate


# # scrna <-  RunPCA(scrna, features = VariableFeatures(object = scrna),assay = "SCT",npcs=20)
# # saveRDS(scrna, file = "KudoUUOComparision_Scrna_PCA.rds")
# # #scrna <-readRDS(file = "KudoUUOComparision_Scrna_PCA.rds")
# # # then use the Harmony to remove batches
# # scrna <- RunHarmony(scrna, group.by.vars = 'Sample', dims.use = 1:20, plot_convergence= FALSE, reduction = "pca", assay.use = "SCT", reduction.save = "harmony")
# # #saveRDS(scrna, file = "KudoUUOComparision_Scrna_Harmony.rds")
# # scrna <- RunUMAP(scrna, dims = 1:20,reduction = "harmony")
# # #saveRDS(scrna, file = "KudoUUOComparision_Scrna_RunUMAP.rds")
# # scrna <- FindNeighbors(scrna, dims = 1:20, reduction = "harmony")
# # #saveRDS(scrna, file = "KudoUUOComparision_Scrna_FindNeighbor.rds")
# # scrna <- FindClusters(scrna, resolution = 0.2)


