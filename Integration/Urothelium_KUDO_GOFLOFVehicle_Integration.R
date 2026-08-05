#!/usr/bin/env Rscript
## Please install all the packages before loading them.
sessionInfo()
library(sf)
library(dplyr)
library(Seurat)
library(patchwork)

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
library(devtools)

library("sctransform")
library("ggthemes")

# BiocManager::install("limma")
# BiocManager::install("fgsea")
# BiocManager::install("msigdbr")
library(limma)
library(nichenetr)
library(tidyverse)
library(dplyr)
library("data.table")

library(umap)
library(patchwork)
library(cowplot)

library("ggsci")
library("ggplot2")
library("gridExtra")

library(harmony)
library(tidyverse)
library(ggplot2)
library("AnnotationDbi")
library("org.Hs.eg.db")
library(here)
library(future)
library(future.callr)

#remotes::install_github(repo ='chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)

# GSEA / gene set tools
library(msigdbr)
library(fgsea)
library(pheatmap)
library(ggpubr)

options(future.globals.maxSize = 1e12)

Outdir <- c("/home/gdbecknelllab/xxw004/gdjacksonlab/UUO/Results/KUDOs/Integration/GOFLOFVehicle/")
dir.create(Outdir, showWarnings = FALSE, recursive = TRUE)
setwd(Outdir)

Indir <- c("/home/gdbecknelllab/xxw004/gdjacksonlab/UUO/Results/KUDOs/Preprocess/")

# we only keep the three conditions of interest -- Yoda is intentionally excluded
Conditions <- c("Vehicle", "GOF", "LOF")

QCdir <- paste0(Outdir, "QC/")
DoubletDir <- paste0(Outdir, "DoubletFinder/")
IntegratedDir <- paste0(Outdir, "IntegratedH/")
CandidateMarkersDir <- paste0(Outdir, "CandidateMarkers/")
LineageDir <- paste0(Outdir, "LineageProportions/")
NotchDir <- paste0(Outdir, "NotchSignaling/")
GSEADir <- paste0(Outdir, "GSEA_ProlifDiff/")

for (d in c(QCdir, DoubletDir, IntegratedDir, CandidateMarkersDir, LineageDir, NotchDir, GSEADir)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

set.seed(10000)

scrna.list <- list()

for (i in Conditions) {
  scrna.list[[i]] <- ReadMtx(
    mtx = paste0(Indir, i, "/raw_matrix/matrix.mtx.gz"),
    features = paste0(Indir, i, "/raw_matrix/features.tsv.gz"),
    cells = paste0(Indir, i, "/raw_matrix/barcodes.tsv.gz")
  )

  scrna.list[[i]] <- as(scrna.list[[i]], "dgCMatrix")
  scrna.list[[i]] <- CreateSeuratObject(
    counts = scrna.list[[i]],
    project = i,
    min.cells = 3,
    min.features = 200
  )

  scrna.list[[i]] <- RenameCells(scrna.list[[i]], add.cell.id = i)
  assign(i, scrna.list[[i]])
  scrna.list[[i]]$DataSet <- rep(i, length(scrna.list[[i]]$orig.ident))

  scrna.list[[i]]$percent.mt <- PercentageFeatureSet(scrna.list[[i]], pattern = "^mt-")
  scrna.list[[i]]$rDNA <- PercentageFeatureSet(scrna.list[[i]], pattern = "^Rp[sl][[:digit:]]")

  scrna.list[[i]] <- subset(
    scrna.list[[i]],
    subset = nFeature_RNA > 200 &
      nFeature_RNA < 8000 &
      nCount_RNA < 16000 & percent.mt < 20 & rDNA < 40
  )

  scrna.list[[i]] <- NormalizeData(scrna.list[[i]])
  scrna.list[[i]] <- ScaleData(scrna.list[[i]])
  scrna.list[[i]] <- FindVariableFeatures(scrna.list[[i]], selection.method = "vst", nfeatures = 2000)
  scrna.list[[i]] <- RunPCA(scrna.list[[i]])
  scrna.list[[i]] <- FindNeighbors(scrna.list[[i]], dims = 1:20)
  scrna.list[[i]] <- RunUMAP(scrna.list[[i]], dims = 1:20)
}

setwd(QCdir)

scrna <- merge(
  x = scrna.list[[1]],
  y = c(scrna.list[[2]], scrna.list[[3]]),
  project = "KudoGOFLOFVehicle"
)

levels(as.factor(scrna@meta.data$DataSet))

scrna@meta.data$DataSet <- ordered(
  factor(scrna@meta.data$DataSet),
  levels = c("Vehicle", "GOF", "LOF")
)

scrna[["percent.mt"]] <- PercentageFeatureSet(scrna, pattern = "^mt-")
scrna[["rDNA"]] <- PercentageFeatureSet(scrna, pattern = "^Rp[sl][[:digit:]]")

# 1. QC plots restricted to GOF/LOF/Vehicle
QCReport <- VlnPlot(
  scrna,
  features = c("rDNA", "percent.mt", "nCount_RNA", "nFeature_RNA"),
  split.by = "DataSet",
  pt.size = 0,
  ncol = 1
)
ggsave(filename = "GOFLOFVehicle_SampleQC_reports.pdf", plot = QCReport, height = 12, width = 12)

pdf("GOFLOFVehicle_Feature.pdf", height = 4, width = 6)
VlnPlot(scrna, features = c("nFeature_RNA"), group.by = "DataSet", pt.size = 0)
dev.off()

pdf("GOFLOFVehicle_nCounts.pdf", height = 4, width = 6)
VlnPlot(scrna, features = c("nCount_RNA"), group.by = "DataSet", pt.size = 0)
dev.off()

pdf("GOFLOFVehicle_mtPercent.pdf", height = 4, width = 6)
VlnPlot(scrna, features = c("percent.mt"), group.by = "DataSet", pt.size = 0, y.max = 80)
dev.off()

pdf("GOFLOFVehicle_rDNA.pdf", height = 4, width = 6)
VlnPlot(scrna, features = c("rDNA"), group.by = "DataSet", pt.size = 0)
dev.off()

scrna <- subset(
  scrna,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 8000 &
    nCount_RNA < 16000 & percent.mt < 20 & rDNA < 40
)

setwd(DoubletDir)

FindDoublets <- function(library_id, seurat_aggregate) {
  seurat_obj <- subset(seurat_aggregate, idents = library_id)
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  seurat_obj <- RunPCA(seurat_obj)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)

  sweep.res.list_kidney <- paramSweep(seurat_obj, PCs = 1:20, sct = F)
  sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
  bcmvn_kidney <- find.pK(sweep.stats_kidney)
  pK <- bcmvn_kidney %>%
    dplyr::filter(BCmetric == max(BCmetric)) %>%
    dplyr::select(pK)
  pK <- as.numeric(as.character(pK[[1]]))
  seurat_doublets <- doubletFinder(
    seurat_obj, PCs = 1:20, pN = 0.25, pK = pK,
    nExp = round(0.05 * length(seurat_obj@active.ident)),
    reuse.pANN = FALSE, sct = F
  )

  DF.class <- names(seurat_doublets@meta.data) %>% str_subset("DF.classifications")
  pANN <- names(seurat_doublets@meta.data) %>% str_subset("pANN")

  p1 <- ggplot(bcmvn_kidney, aes(x = pK, y = BCmetric)) +
    geom_bar(stat = "identity") +
    ggtitle(paste0("pKmax=", pK)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p2 <- DimPlot(seurat_doublets, group.by = DF.class)
  p3 <- FeaturePlot(seurat_doublets, features = pANN)

  ggsave(filename = paste0(library_id, "_pKmax_distribution.pdf"), p1, height = 5, width = 6)
  ggsave(filename = paste0(library_id, "_Dimplotdoublet_distribution.pdf"), p2, height = 5, width = 5)
  ggsave(filename = paste0(library_id, "_Featureplotdoublet_distribution.pdf"), p3, height = 5, width = 5)

  df_doublet_barcodes <- as.data.frame(cbind(rownames(seurat_doublets@meta.data), seurat_doublets@meta.data[[DF.class]]))
  return(df_doublet_barcodes)
}

Idents(scrna) <- scrna$DataSet
orig.ident <- levels(Idents(scrna))

list.doublet.bc <- lapply(orig.ident, function(x) { FindDoublets(x, seurat_aggregate = scrna) })

doublet_id <- list.doublet.bc %>%
  bind_rows() %>%
  dplyr::rename("doublet_id" = "V2") %>%
  tibble::column_to_rownames(var = "V1")

table(doublet_id)
scrna <- AddMetaData(scrna, doublet_id)

saveRDS(scrna, file = "GOFLOFVehicle_merged_withdoublets.rds")

scrna <- subset(scrna, subset = (doublet_id == "Singlet"))

setwd(IntegratedDir)

scrna <- scrna %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst") %>%
  ScaleData()

scrna <- RunPCA(scrna, verbose = TRUE, npcs = 30)
scrna <- RunHarmony(scrna, group.by.vars = c("DataSet"), reduction.save = "harmony")
scrna <- RunUMAP(scrna, dims = 1:30, reduction = "harmony")
scrna <- FindNeighbors(scrna, dims = 1:30, reduction = "harmony")
scrna <- FindClusters(scrna, resolution = 0.4, algorithm = 1)

DimplotScrna <- DimPlot(scrna, reduction = "umap", label = T)
DimplotScrnaSplit <- DimPlot(scrna, reduction = "umap", split.by = "DataSet", label = T)

ggsave(filename = "GOFLOFVehicle_Dimplot.pdf", plot = DimplotScrna, width = 6, height = 5)
ggsave(filename = "GOFLOFVehicle_DimplotSplit.pdf", plot = DimplotScrnaSplit, width = 15, height = 5)

saveRDS(scrna, file = "Kudo_GOFLOFVehicle_harmony.RDS")

scrna <- FindClusters(scrna, resolution = 1.5, algorithm = 1)
saveRDS(scrna, file = "Kudo_GOFLOFVehicle_harmony_res1.5.RDS")

Idents(scrna) <- scrna$DataSet

setwd(CandidateMarkersDir)

# basal / intermediate / umbrella urothelial lineage markers used throughout this repo
urothelium.markers <- c(
  "Krt5", "Krt14", "Trp63", "Itga6", "Krt17",
  "Krt8", "Krt18", "Foxa1", "Gata3", "Pparg",
  "Upk1a", "Upk1b", "Upk2", "Upk3b", "Upk3a", "Krt20"
)

for (cond in Conditions) {
  obj <- subset(scrna, subset = DataSet == cond)
  p <- DotPlot(obj, features = urothelium.markers) + RotatedAxis() + ggtitle(cond)
  ggsave(paste0(cond, "_UrotheliumMarkers_DotPlot.pdf"), plot = p, height = 4, width = 8)
}

# combined view, split by condition, for direct comparison
CombinedUroDotPlot <- DotPlot(scrna, features = urothelium.markers, group.by = "DataSet") +
  RotatedAxis()
ggsave("GOFLOFVehicle_UrotheliumMarkers_DotPlot_Combined.pdf", plot = CombinedUroDotPlot, height = 4, width = 8)

# nephron-segment and immune markers used elsewhere in this repo to check for
# non-urothelial (kidney) cell contamination in the organoid
markers_mouse_expanded_withothers <- c(
  "Krt5", "Krt14", "Trp63",
  "Krt8", "Krt18", "Foxa1",
  "Upk1a", "Upk1b", "Upk2", "Upk3b", "Upk3a", "Krt20",
  "Lrp2", "Slc34a1", "Slc13a3", # Proximal tubule cells (PT)
  "Slc12a3", "Pvalb", "Wnk1", # Distal connecting tubule (DCT)
  "Umod", "Slc12a1", "Cldn10", # Thick ascending limb cells (TAL)
  "Klk1", "Slc8a1", "Calb1", # Connecting tubule (CNT)
  "Slc14a2", "Fst", "Bst1", # descending thin limb cells (DTL)
  "Aqp2", "Hsd11b2", "Scnn1g", # principal cells (PC)
  "Atp6v1g3", "Aqp6", "Slc26a7", # Intercalated cells
  "Cdh5", "Igfbp3", "Pecam1", # Endothelium
  "Col3a1", "Fbn1", "Lum", "Col1a1", # Fibroblast
  "Ireb2", "Alas2", # Plasmacytoid dendritic cells
  "S100a8", "S100a9", "Il1b", # Neutrophil
  "C1qa", "C1qb", "Aif1", # Macrophage
  "Cxcr6", "Cd247", "Cd3e", # T cells
  "Igkc", "Cd79a", "Cd79b", # B cells
  "Ccl5", "Nkg7", "Cd7" # NK/Cd8 cells
)

avg_exp <- AverageExpression(
  object = scrna,
  features = markers_mouse_expanded_withothers,
  group.by = "DataSet"
)

aveg_exp_df <- as.data.frame(avg_exp$RNA)

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
aveg_exp_df$Gene <- factor(aveg_exp_df$Gene, levels = markers_mouse_expanded_withothers)

aveg_exp_df_long <- aveg_exp_df %>%
  pivot_longer(cols = -Gene, names_to = "DataSet", values_to = "AverageExpression")
aveg_exp_df_long$DataSet <- factor(aveg_exp_df_long$DataSet, levels = Conditions)

KidneyMarkersBarPlot <- ggplot(aveg_exp_df_long, aes(x = Gene, y = AverageExpression, fill = DataSet)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_npg() +
  theme_classic(base_size = 14) +
  labs(y = "Average Expression", x = NULL) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  )
ggsave("GOFLOFVehicle_KidneyMarkersBarPlot.pdf", plot = KidneyMarkersBarPlot, height = 5, width = 16)
write.csv(aveg_exp_df_long, "GOFLOFVehicle_KidneyMarkers_AverageExpression.csv", row.names = FALSE)

setwd(LineageDir)

lineage.expr <- FetchData(
  scrna,
  vars = c("Upk1a", "Upk1b", "Upk2", "Upk3a", "Krt5", "Krt14", "DataSet")
)

# Upk+ takes priority over Krt5+ when a cell expresses both lineages of markers,
# matching the >2 (log-normalized) threshold convention already used in this repo
lineage.expr <- lineage.expr %>%
  mutate(Lineage = case_when(
    Upk1a > 2 | Upk1b > 2 | Upk2 > 2 | Upk3a > 2 ~ "Upk+",
    Krt5 > 2 | Krt14 > 2 ~ "Krt5+",
    TRUE ~ "Other"
  ))

scrna$Lineage <- lineage.expr$Lineage

lineage.counts <- lineage.expr %>%
  count(DataSet, Lineage) %>%
  group_by(DataSet) %>%
  mutate(Proportion = n / sum(n)) %>%
  ungroup()

write.csv(lineage.counts, "GOFLOFVehicle_Lineage_Proportions.csv", row.names = FALSE)

# sort conditions by their Upk+ proportion (descending) and fix the lineage stacking order
condition_order <- lineage.counts %>%
  filter(Lineage == "Upk+") %>%
  arrange(desc(Proportion)) %>%
  pull(DataSet)

lineage.counts$DataSet <- factor(lineage.counts$DataSet, levels = condition_order)
lineage.counts$Lineage <- factor(lineage.counts$Lineage, levels = c("Upk+", "Krt5+", "Other"))

LineageProportionPlot <- ggplot(lineage.counts, aes(x = DataSet, y = Proportion, fill = Lineage)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_npg() +
  theme_classic(base_size = 14) +
  labs(y = "Proportion of cells", x = NULL)
ggsave("GOFLOFVehicle_LineageProportion_Stacked.pdf", plot = LineageProportionPlot, height = 5, width = 5)

LineageProportionDodge <- ggplot(lineage.counts, aes(x = reorder(DataSet, -Proportion), y = Proportion, fill = Lineage)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_npg() +
  theme_classic(base_size = 14) +
  labs(y = "Proportion of cells", x = NULL)
ggsave("GOFLOFVehicle_LineageProportion_Dodge.pdf", plot = LineageProportionDodge, height = 5, width = 6)

# pairwise proportion tests for Upk+ and Krt5+ fractions between conditions
lineage.wide <- lineage.expr %>% count(DataSet, Lineage) %>% pivot_wider(names_from = Lineage, values_from = n, values_fill = 0)
lineage.wide$Total <- rowSums(lineage.wide[, c("Upk+", "Krt5+", "Other")])

pairwise_prop_test <- function(lineage_col) {
  pairs <- combn(as.character(lineage.wide$DataSet), 2, simplify = FALSE)
  res <- lapply(pairs, function(p) {
    x <- lineage.wide[lineage.wide$DataSet %in% p, ]
    test <- prop.test(x[[lineage_col]], x$Total)
    data.frame(Comparison = paste(p, collapse = " vs "), Lineage = lineage_col, p.value = test$p.value)
  })
  bind_rows(res)
}

lineage_stats <- bind_rows(pairwise_prop_test("Upk+"), pairwise_prop_test("Krt5+"))
lineage_stats$p.adj <- p.adjust(lineage_stats$p.value, method = "BH")
write.csv(lineage_stats, "GOFLOFVehicle_Lineage_PairwiseProportionTests.csv", row.names = FALSE)

setwd(NotchDir)

# MSigDB C2:CP:WIKIPATHWAYS human gene set, converted to mouse orthologs with nichenetr
# to match this Seurat object's mouse gene symbols
notch.human.geneset <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:WIKIPATHWAYS") %>%
  dplyr::filter(gs_name == "WP_NOTCH_SIGNALING_PATHWAY") %>%
  dplyr::pull(gene_symbol) %>%
  unique()

notch.mouse.genes <- convert_human_to_mouse_symbols(notch.human.geneset)
notch.mouse.genes <- notch.mouse.genes[!is.na(notch.mouse.genes)]
notch.mouse.genes <- intersect(notch.mouse.genes, rownames(scrna))

length(notch.mouse.genes)
write.csv(data.frame(Gene = notch.mouse.genes), "WP_NOTCH_SIGNALING_PATHWAY_MouseGenes.csv", row.names = FALSE)

# average expression heatmap of Notch pathway genes across the three conditions
notch.avg <- AverageExpression(scrna, features = notch.mouse.genes, group.by = "DataSet")$RNA
notch.avg <- notch.avg[, Conditions]
notch.avg.mat <- as.matrix(notch.avg)
# drop genes with zero variance across the three conditions (can't be z-scored, only 3 points)
notch.avg.mat <- notch.avg.mat[apply(notch.avg.mat, 1, sd) > 0, , drop = FALSE]
notch.z <- t(scale(t(notch.avg.mat)))

pdf("GOFLOFVehicle_NotchSignaling_Heatmap.pdf", height = 12, width = 5)
pheatmap(
  notch.z,
  cluster_cols = FALSE,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  main = "WP_NOTCH_SIGNALING_PATHWAY (z-score of average expression)"
)
dev.off()

# per-cell Notch pathway module score, compared pairwise across conditions
scrna <- AddModuleScore(scrna, features = list(notch.mouse.genes), name = "NotchScore")

NotchScoreViolin <- VlnPlot(scrna, features = "NotchScore1", group.by = "DataSet", pt.size = 0) +
  stat_compare_means(comparisons = combn(Conditions, 2, simplify = FALSE), method = "wilcox.test") +
  ggtitle("Notch signaling module score")
ggsave("GOFLOFVehicle_NotchModuleScore_Violin.pdf", plot = NotchScoreViolin, height = 5, width = 5)

notch.score.stats <- pairwise.wilcox.test(scrna$NotchScore1, scrna$DataSet, p.adjust.method = "BH")
capture.output(notch.score.stats, file = "GOFLOFVehicle_NotchModuleScore_PairwiseWilcox.txt")

setwd(GSEADir)

# GO:BP gene sets from MSigDB (human), converted to mouse orthologs, restricted to
# proliferation- and differentiation-related terms
gobp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")

proliferation.sets <- gobp %>% dplyr::filter(grepl("PROLIFERATION", gs_name))
differentiation.sets <- gobp %>% dplyr::filter(grepl("DIFFERENTIATION", gs_name))
gsea.sets <- bind_rows(proliferation.sets, differentiation.sets)

pathway.list <- split(gsea.sets$gene_symbol, gsea.sets$gs_name) %>%
  lapply(function(g) {
    mouse.g <- convert_human_to_mouse_symbols(unique(g))
    mouse.g[!is.na(mouse.g)]
  })

Idents(scrna) <- scrna$DataSet

run_pairwise_gsea <- function(ident1, ident2) {
  markers <- FindMarkers(
    scrna,
    ident.1 = ident1,
    ident.2 = ident2,
    group.by = "DataSet",
    logfc.threshold = 0,
    min.pct = 0.05
  )
  markers$gene <- rownames(markers)
  markers <- markers[!is.na(markers$avg_log2FC), ]

  ranks <- markers$avg_log2FC
  names(ranks) <- markers$gene
  ranks <- sort(ranks, decreasing = TRUE)

  fgseaRes <- fgsea(pathways = pathway.list, stats = ranks, minSize = 10, maxSize = 500, eps = 0)
  fgseaRes$Comparison <- paste0(ident1, "_vs_", ident2)
  fgseaRes[order(fgseaRes$pval), ]
}

comparisons <- list(
  c("GOF", "Vehicle"),
  c("LOF", "Vehicle"),
  c("GOF", "LOF")
)

gsea.results <- lapply(comparisons, function(p) run_pairwise_gsea(p[1], p[2]))
gsea.results.df <- bind_rows(gsea.results)

# leadingEdge is a list-column; flatten it for a clean CSV export
gsea.results.export <- gsea.results.df
gsea.results.export$leadingEdge <- vapply(gsea.results.export$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
write.csv(gsea.results.export, "GOFLOFVehicle_GSEA_ProlifDiff_AllResults.csv", row.names = FALSE)

# summary dot plot of NES across comparisons for the top significant terms in either direction
top.terms <- gsea.results.df %>%
  filter(padj < 0.05) %>%
  group_by(pathway) %>%
  filter(any(padj < 0.05)) %>%
  ungroup() %>%
  distinct(pathway) %>%
  pull(pathway)

plot.data <- gsea.results.df %>% filter(pathway %in% top.terms)
plot.data$Comparison <- factor(plot.data$Comparison, levels = c("GOF_vs_Vehicle", "LOF_vs_Vehicle", "GOF_vs_LOF"))

GSEASummaryDotPlot <- ggplot(plot.data, aes(x = Comparison, y = reorder(pathway, NES), size = -log10(padj), color = NES)) +
  geom_point() +
  scale_color_gradient2(low = "navy", mid = "white", high = "firebrick3", midpoint = 0) +
  theme_bw(base_size = 10) +
  labs(x = NULL, y = NULL, color = "NES", size = "-log10(FDR)")
ggsave("GOFLOFVehicle_GSEA_ProlifDiff_SummaryDotPlot.pdf", plot = GSEASummaryDotPlot,
       height = max(4, 0.25 * length(top.terms)), width = 9)

# individual enrichment plots for the top proliferation and differentiation term per comparison
for (comp in unique(gsea.results.df$Comparison)) {
  res.comp <- gsea.results.df %>% filter(Comparison == comp)
  top.prolif <- res.comp %>% filter(grepl("PROLIFERATION", pathway)) %>% arrange(pval) %>% slice(1)
  top.diff <- res.comp %>% filter(grepl("DIFFERENTIATION", pathway)) %>% arrange(pval) %>% slice(1)

  ident1 <- strsplit(comp, "_vs_")[[1]][1]
  ident2 <- strsplit(comp, "_vs_")[[1]][2]
  markers <- FindMarkers(scrna, ident.1 = ident1, ident.2 = ident2, group.by = "DataSet", logfc.threshold = 0, min.pct = 0.05)
  ranks <- sort(setNames(markers$avg_log2FC, rownames(markers)), decreasing = TRUE)

  if (nrow(top.prolif) == 1) {
    pdf(paste0("GOFLOFVehicle_GSEA_", comp, "_", top.prolif$pathway, "_EnrichmentPlot.pdf"), height = 4, width = 6)
    print(plotEnrichment(pathway.list[[top.prolif$pathway]], ranks) + ggtitle(paste(comp, top.prolif$pathway)))
    dev.off()
  }
  if (nrow(top.diff) == 1) {
    pdf(paste0("GOFLOFVehicle_GSEA_", comp, "_", top.diff$pathway, "_EnrichmentPlot.pdf"), height = 4, width = 6)
    print(plotEnrichment(pathway.list[[top.diff$pathway]], ranks) + ggtitle(paste(comp, top.diff$pathway)))
    dev.off()
  }
}

setwd(Outdir)
saveRDS(scrna, file = "Urothelium_KUDO_GOFLOFVehicle_Final.rds")
