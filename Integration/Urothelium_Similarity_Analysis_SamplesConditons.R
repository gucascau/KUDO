#!/usr/bin/env Rscript
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

library(patchwork)
library(decoupleR)

# Extra libraries
library(dplyr)
library(pheatmap)

library(patchwork)

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

sessionInfo()

## Step 2: Set up the workspace environment:

setwd("/home/gdbecknelllab/xxw004/gdjacksonlab/UUO/Results/KUDOs/Integration/IntegratingUUO_Organoid/")
KUDODir <- c("/home/gdbecknelllab/xxw004/gdjacksonlab/UUO/Results/KUDOs/Integration/IntegratingUUO_Organoid/")

# Make directories for the results
# annotation directory
dir.create(paste0(KUDODir,"Annotation/"))
# DEG directory
dir.create(paste0(KUDODir,"DEG/"))
# Feature directory
dir.create(paste0(KUDODir,"Feature/"))
# Cell cell communication directory
dir.create(paste0(KUDODir,"CellCommunication/"))

# Our Final Integrated Object

# HarmonyIntegration
AllUrothelium<- readRDS("KudoUUOComparision_Scrna_Merged_RunHarmony.rds")

#AllUrothelium <- JoinLayers(AllUrothelium)


# ─── Pairwise correlation matrix for all 12 samples ───────────────────────────

# Step 1: Split AllUrothelium by Sample into a named list of Seurat objects
sample_levels <- c("Vehicle", "Yoda", "LOF", "GOF", "Health",
                   "UUO_Day2", "UUO_Day4", "UUO_Day6", "UUO_Day10", "UUO_Day14",
                   "GSM3827175_NMU_O_P", "GSM3827176_NMU_O_D")

DefaultAssay(AllUrothelium) <- "RNA"
AllUrothelium <- JoinLayers(AllUrothelium)

sample_list <- lapply(sample_levels, function(s) {
  subset(AllUrothelium, subset = Sample == s)
})
names(sample_list) <- sample_levels

# Step 2: Compute per-sample average expression vectors
avg_expr_list <- lapply(sample_list, function(obj) {
  Matrix::rowMeans(GetAssayData(obj, assay = "RNA", layer = "data"))
})

# Step 3: Find genes expressed (mean > 0) in ALL 12 samples
shared_genes <- Reduce(intersect, lapply(avg_expr_list, function(avg) {
  names(avg[avg > 0])
}))
message("Shared expressed genes across all 12 samples: ", length(shared_genes))

# Step 4: Build a gene × sample matrix restricted to shared genes
expr_matrix <- do.call(cbind, lapply(avg_expr_list, function(avg) avg[shared_genes]))
colnames(expr_matrix) <- sample_levels

# Step 5: Compute pairwise Spearman correlation matrix
cor_matrix <- cor(expr_matrix, method = "spearman")

# Step 6: Save the numeric matrix
write.csv(cor_matrix, "AllUrothelium_12Sample_SpearmanCorrelation.csv")

# Step 7: Heatmap
library(pheatmap)
pdf("AllUrothelium_12Sample_SpearmanCorrelation_Heatmap.pdf", width = 8, height = 7)
pheatmap(cor_matrix,
         color        = colorRampPalette(c("white", "#2166ac"))(100),
         breaks       = seq(0.8, 1, length.out = 101),
         display_numbers = TRUE,
         number_format   = "%.3f",
         fontsize_number = 7,
         cluster_rows = TRUE, cluster_cols = TRUE,
         main   = "Spearman Correlation — 12 Samples (shared expressed genes)")

dev.off()


### this section is for the marker gene expression
markers_mouse_expanded <- c(
  # Basal
  "Krt5","Krt14","Trp63","Itga6","Krt17","Col17a1","Sox2","Ngfr","Cd44","Lgr5","Pdpn","Epcam",
  # Intermediate
  "Krt8","Krt18","Foxa1","Gata3","Pparg","Grhl3","Upk1a","Upk1b","Upk2","Upk3a","S100a1","Cdh1","Tjp1",
  # Umbrella
  "Krt20","Shh","Plk2","Umod","Lypd3","Cldn8","Cldn4","Ocln"
)

Idents(AllUrothelium) <- "Sample"
UrothelialMarkersDotPlot<- DotPlot(AllUrothelium, features = markers_mouse_expanded) + RotatedAxis()
ggsave("UrothelialMarkersDotPlot_Urothelim_KUDO_MouseH_BySamples.pdf", plot=UrothelialMarkersDotPlot, height=3, width=18)



Idents(AllUrothelium) <- "Group"
UrothelialMarkersDotPlot<- DotPlot(AllUrothelium, features = markers_mouse_expanded) + RotatedAxis()
ggsave("UrothelialMarkersDotPlot_Urothelim_KUDO_MouseH_ByGroup.pdf", plot=UrothelialMarkersDotPlot, height=3, width=19)

markers_mouse_expanded_withothers<- c("Krt5","Krt14","Trp63",
  "Krt8","Krt18","Foxa1",
  "Upk1a","Upk1b","Upk2","Upk3b","Upk3a","Krt20",
  "Lrp2","Slc34a1","Slc13a3", # Proximal tubule cells (PT)
  "Slc12a3","Pvalb","Wnk1",# Distal connecting tublue (DCT)
"Umod","Slc12a1","Cldn10", # Thick ascending limb cells (TAL)
                      "Klk1","Slc8a1","Calb1", # Connecting tubule (CNT)
                     "Slc14a2","Fst","Bst1", # descending thin  limb cells (DTL)
                     #"Akr1b1","Sh3gl3","Prox1", # ascending thin limb cells (ATL)
                     #"Enox1","Thsd4", #Macula densa
                      "Aqp2","Hsd11b2","Scnn1g", # principal cells (PC)
                    "Atp6v1g3","Aqp6","Slc26a7", # Intercalated cells
                    "Cdh5","Igfbp3","Pecam1", # Endothelium
                    "Col3a1","Fbn1","Lum","Col1a1", # Fibroblast
                    "Ireb2","Alas2", # Plasmacytoid dentritic cells
                    "S100a8","S100a9","Il1b", # Neutrophil
                    "C1qa","C1qb","Aif1", # Macrophage
                    "Cxcr6","Cd247","Cd3e", # T cells
                    "Igkc","Cd79a","Cd79b", # B cells
                    "Ccl5","Nkg7","Cd7" # NK/Cd8 cells
  
  )
Idents(AllUrothelium) <- "Sample"
KidneyMarkersDotPlot<- DotPlot(AllUrothelium, features = markers_mouse_expanded_withothers) + RotatedAxis()
ggsave("KidneyMarkersDotPlot_Urothelim_KUDO_MouseH_Reclustered.pdf", plot=KidneyMarkersDotPlot, height=4, width=15)

# check the average expression of the markers in the mouse urothelium and organoid urothelium

avg_exp <- AverageExpression(
  object = AllUrothelium,
  features = markers_mouse_expanded_withothers,
  group.by = "Sample"   # or "seurat_clusters", "condition"
)
avg_exp%>% head()
# set up the order of the columns in the same order as the markers
aveg_exp_df <- as.data.frame(avg_exp$RNA)


# identify missing genes
missing_genes <- setdiff(markers_mouse_expanded_withothers, rownames(aveg_exp_df))
if (length(missing_genes) > 0) {
  zero_mat <- matrix(
    0,
    nrow = length(missing_genes),
    ncol = ncol(aveg_exp_df),
    dimnames = list(missing_genes, colnames(aveg_exp_df))
  )
  
  aveg_exp_df <- rbind(aveg_exp_df, zero_mat)
}

aveg_exp_df$Gene <- rownames(aveg_exp_df)

# change to the long type
aveg_exp_df_long <- aveg_exp_df %>%
  pivot_longer(cols = -Gene, names_to = "Sample", values_to = "AverageExpression")

library(ggplot2)
# set the order of the genes in the same order as the markers
aveg_exp_df_long$Gene <- factor(aveg_exp_df_long$Gene, levels = markers_mouse_expanded_withothers)

AvereageExpressionBarPlot <-
ggplot(aveg_exp_df_long, aes(x = Gene, y = AverageExpression, fill= Sample)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  theme_classic(base_size = 14) +
  labs(y = "Average Expression", x = NULL) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  )
ggsave("AvereageExpressionBarPlot_Urothelim_KUDO_MouseH_Reclustered.pdf", plot=AvereageExpressionBarPlot, height=4, width=15)
