#!/usr/bin/env Rscript
# =============================================================================
# Download_GEO_GEOquery.R
# Author: Xin Wang — xin.wang@nationwidechildrens.org
# Purpose: Download and inspect GEO supplementary files for the renal pelvis
#          urothelium atlas publication using GEOquery.
#          Run AFTER Download_PublicDatasets.sh or as a standalone alternative.
#
# Datasets:
#   GSE234788 — Li et al. Cell Metab 2024  (SHARE-seq human kidney, Priority 1)
#   GSE157079 — GUDMAP Human Bladder       (Priority 3)
#   GSE129845 — Bhatt et al. JCI 2019      (Priority 4)
#   GSE119531 — GUDMAP Mouse Urinary Tract (Priority 4)
# =============================================================================

# ── Install / load packages ───────────────────────────────────────────────────
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("GEOquery",    quietly = TRUE)) BiocManager::install("GEOquery")

library(GEOquery)
library(tidyverse)

options(timeout = 3600)   # 1-hour timeout for large files

# ── Paths ─────────────────────────────────────────────────────────────────────
DATASETS_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets"

dataset_dirs <- list(
  GSE234788 = file.path(DATASETS_DIR, "Human/Li_CellMetab2024_SHAREseq_GSE234788"),
  GSE157079 = file.path(DATASETS_DIR, "Human/GUDMAP_HumanBladder_GSE157079"),
  GSE129845 = file.path(DATASETS_DIR, "Human/Bhatt_JCI2019_HumanMouseBladder_GSE129845"),
  GSE119531 = file.path(DATASETS_DIR, "Mouse/GUDMAP_MouseUrinaryTract_GSE119531")
)

# Ensure directories exist
invisible(lapply(dataset_dirs, dir.create, showWarnings = FALSE, recursive = TRUE))

# ── Helper: fetch metadata + supplementary files for one accession ─────────────
fetch_geo_dataset <- function(accession, outdir) {
  message("\n", strrep("=", 60))
  message("Fetching: ", accession)
  message("Destination: ", outdir)
  message(strrep("=", 60))

  # --- 1. Series metadata (GSE-level) ----------------------------------------
  tryCatch({
    gse <- getGEO(accession, destdir = outdir, GSEMatrix = TRUE, AnnotGPL = FALSE)
    message("[", accession, "] Metadata fetched. Samples: ",
            length(gse[[1]]@phenoData@data[[1]]))

    # Save phenotype table
    pheno_df <- pData(gse[[1]])
    write.csv(pheno_df,
              file = file.path(outdir, paste0(accession, "_phenoData.csv")),
              row.names = TRUE)
    message("[", accession, "] phenoData saved.")
  }, error = function(e) {
    message("[", accession, "] WARNING — metadata fetch failed: ", e$message)
  })

  # --- 2. Supplementary files (count matrices, metadata, etc.) ----------------
  tryCatch({
    suppl_files <- getGEOSuppFiles(accession, baseDir = outdir, fetch_files = TRUE)
    message("[", accession, "] Supplementary files downloaded:")
    print(suppl_files)
  }, error = function(e) {
    message("[", accession, "] WARNING — supplementary file download failed: ", e$message)
    message("  Try running Download_PublicDatasets.sh (wget-based) as fallback.")
  })

  message("[", accession, "] Done.\n")
}

# ── Download each dataset ──────────────────────────────────────────────────────

# Priority 1: Li et al. Cell Metab 2024 — SHARE-seq human kidney
# Contains Uro1 (ureter) and Uro2 (renal pelvis) urothelial clusters
# NOTE: Verify accession GSE234788 in the paper's Data Availability section
#       https://www.cell.com/cell-metabolism/fulltext/S1550-4131(24)00061-5
fetch_geo_dataset("GSE234788", dataset_dirs$GSE234788)

# Priority 3: GUDMAP Human Bladder Cell Atlas
# 29,834 cells; scRNA-seq + snRNA-seq; 4 bladder regions
fetch_geo_dataset("GSE157079", dataset_dirs$GSE157079)

# Priority 4: Bhatt et al. JCI 2019 — Human + Mouse Bladder
# Cross-species urothelial layer reference
fetch_geo_dataset("GSE129845", dataset_dirs$GSE129845)

# Priority 4: GUDMAP Mouse Urinary Tract Atlas
# Within-species comparison: renal pelvis vs. ureter vs. bladder
fetch_geo_dataset("GSE119531", dataset_dirs$GSE119531)

# ── Inspect downloaded files ───────────────────────────────────────────────────
message("\n", strrep("=", 60))
message("Downloaded file summary:")
message(strrep("=", 60))

for (acc in names(dataset_dirs)) {
  d <- dataset_dirs[[acc]]
  files <- list.files(d, recursive = TRUE, full.names = FALSE)
  message("\n[", acc, "] ", length(files), " file(s) in ", d)
  if (length(files) > 0) {
    sizes <- file.size(file.path(d, files))
    df <- data.frame(
      File = files,
      Size_MB = round(sizes / 1e6, 2)
    )
    print(df)
  }
}

# ── Manual download instructions for non-GEO datasets ─────────────────────────
message("\n", strrep("=", 60))
message("MANUAL DOWNLOAD REQUIRED:")
message(strrep("=", 60))
message("
[Priority 2] KPMP snRNA-seq v1.0
  DOI  : 10.48698/6k3k-q779
  URL  : https://www.kpmp.org/doi-collection/10-48698-6k3k-q779
  Dest : ", DATASETS_DIR, "/Human/KPMP_snRNAseq_v1
  Steps: Visit URL, click Download, move files to Dest above.

[Priority 2] KPMP Multiome v1.3 (snRNA-seq + snATAC-seq)
  DOI  : 10.48698/92nk-e805
  URL  : https://www.kpmp.org/doi-collection/10-48698-92nk-e805
  Dest : ", DATASETS_DIR, "/Human/KPMP_Multiome_v1.3

[Priority 3] Lake et al. Nature 2023 — Multimodal Human Kidney Atlas
  URL  : https://cellxgene.cziscience.com  (search 'Lake kidney 2023')
  Dest : ", DATASETS_DIR, "/Human/Lake_Nature2023_KidneyAtlas
  Steps: Download .h5ad file from CellxGene Discover.
         Convert to Seurat with SeuratDisk::Convert() if needed:
           library(SeuratDisk)
           Convert('file.h5ad', dest = 'h5seurat', overwrite = TRUE)
           obj <- LoadH5Seurat('file.h5seurat')
")

message("Script complete: ", Sys.time())
