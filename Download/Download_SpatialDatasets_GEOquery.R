#!/usr/bin/env Rscript
# =============================================================================
# Download_SpatialDatasets_GEOquery.R
# Author: Xin Wang — xin.wang@nationwidechildrens.org
# Purpose: Download and inspect spatial transcriptomics datasets for
#          validation of renal pelvis urothelial identity programs.
#
# Already downloaded (skipped here):
#   SP-3: GSE282059 — CosMx 6000-plex Human Kidney UUO + Healthy
#   SP-4: GSE253439 — CosMx 980-plex Human Kidney IgA/MCD/Pyelonephritis
#
# This script handles:
#   SP-1: GSE269884 -- Xuanyuan et al. Nat Commun 2025 — Mouse Xenium + Visium
#         (accession TBD — verify in paper)
#   SP-2: GSE211785 -- Muto et al. Nat Genet 2024     — Human Visium + CosMx
#         (accession TBD — verify in paper)
# =============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("GEOquery",    quietly = TRUE)) BiocManager::install("GEOquery")

library(GEOquery)
library(tidyverse)

options(timeout = 7200)

# ── Paths ─────────────────────────────────────────────────────────────────────
DATASETS_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets"

spatial_dirs <- list(
  # Replace accession keys below once confirmed from paper Data Availability
  SP3_GSE282059 = file.path(DATASETS_DIR, "Human/GSE282059"),
  SP4_GSE253439 = file.path(DATASETS_DIR, "Human/GSE253439"),
  SP1_Mouse_Xenium_Visium = file.path(DATASETS_DIR,
      "Mouse/Xuanyuan_NatCommun2025_MouseKidney_Xenium_Visium"),
  SP2_Human_Visium_CosMx  = file.path(DATASETS_DIR,
      "Human/Muto_NatGenet2024_HumanKidney_Visium_CosMx")
)

invisible(lapply(spatial_dirs, dir.create, showWarnings = FALSE, recursive = TRUE))

# ── Inspect already-downloaded datasets ───────────────────────────────────────

# SP-3: CosMx 6000-plex Human Kidney UUO + Healthy
message("\n", strrep("=", 60))
message("SP-3: GSE282059 — CosMx 6000-plex Human Kidney UUO + Healthy")
message("Status: ALREADY DOWNLOADED")
message(strrep("=", 60))

sp3_files <- list.files(spatial_dirs$SP3_GSE282059, full.names = FALSE)
message("Files present (", length(sp3_files), "):")
print(sp3_files)
message("
CosMx file format guide:
  *_exprMat_file.csv   — gene expression matrix (cells × genes)
  *_metadata_file.csv  — cell metadata (cell ID, FOV, x/y coordinates)
  *_tx_file.csv        — per-transcript table (gene, x, y, cell)
  *_fov_positions_file.csv — field-of-view coordinates
  *_polygons.csv       — cell boundary polygons

Load in R with Seurat:
  library(Seurat)
  nanostring_dir <- '", spatial_dirs$SP3_GSE282059, "'
  tx   <- read.csv(file.path(nanostring_dir, 'GSE282059_2B_FFPE_tx_file.csv'))
  meta <- read.csv(file.path(nanostring_dir, 'GSE282059_2B_FFPE_metadata_file.csv'))
  expr <- read.csv(file.path(nanostring_dir, 'GSE282059_2B_FFPE_exprMat_file.csv'), row.names=1)
  obj  <- CreateSeuratObject(counts = t(expr), meta.data = meta)
")

# SP-4: CosMx 980-plex Human Kidney Biopsies
message("\n", strrep("=", 60))
message("SP-4: GSE253439 — CosMx 980-plex Human Kidney IgA/MCD/Pyelonephritis")
message("Status: ALREADY DOWNLOADED — Seurat object available")
message(strrep("=", 60))

sp4_files <- list.files(spatial_dirs$SP4_GSE253439, full.names = FALSE)
message("Files present (", length(sp4_files), "):")
print(sp4_files)
message("
Load pre-built Seurat object directly:
  library(Seurat)
  cosmx_obj <- readRDS('",
  file.path(spatial_dirs$SP4_GSE253439, "GSE253439_seurat_object_cosmx.rds.gz"),
  "')
Note: This is a pre-processed Seurat object — ready for immediate analysis.
")

# ── SP-1: Xuanyuan et al. Mouse Xenium + Visium ───────────────────────────────
message("\n", strrep("=", 60))
message("SP-1: Xuanyuan et al. Nat Commun 2025 — Mouse Xenium + Visium")
message(strrep("=", 60))
message("
ACTION REQUIRED: Verify GEO accession number.
  1. Go to: https://www.nature.com/articles/s41467-025-62599-9
  2. Find 'Data Availability' section
  3. Note the GEO accession (e.g., GSExxxxxx)
  4. Uncomment and update the code below, then re-run this script.
")

# After confirming accession, uncomment:
SP1_ACCESSION <- "GSE269884"   # <-- replace with confirmed accession
tryCatch({
  gse <- getGEO(SP1_ACCESSION,
                destdir = spatial_dirs$SP1_Mouse_Xenium_Visium,
                GSEMatrix = TRUE)
  pheno <- pData(gse[[1]])
  write.csv(pheno,
    file.path(spatial_dirs$SP1_Mouse_Xenium_Visium,
              paste0(SP1_ACCESSION, "_phenoData.csv")))
  getGEOSuppFiles(SP1_ACCESSION,
                  baseDir = spatial_dirs$SP1_Mouse_Xenium_Visium)
  message("[SP-1] Download complete.")
}, error = function(e) message("[SP-1] Error: ", e$message))

# ── SP-2: Muto et al. Human Visium + CosMx ────────────────────────────────────
message("\n", strrep("=", 60))
message("SP-2: Muto et al. Nat Genet 2024 — Human Visium + CosMx")
message(strrep("=", 60))
message("
ACTION REQUIRED: Verify GEO accession number.
  1. Go to: https://www.nature.com/articles/s41588-024-01802-x
  2. Find 'Data Availability' section
  3. Note the GEO accession (e.g., GSExxxxxx)
  4. Uncomment and update the code below, then re-run this script.
")

# After confirming accession, uncomment:
SP2_ACCESSION <- "GSE211785"   # <-- replace with confirmed accession
tryCatch({
  gse <- getGEO(SP2_ACCESSION,
                destdir = spatial_dirs$SP2_Human_Visium_CosMx,
                GSEMatrix = TRUE)
  pheno <- pData(gse[[1]])
  write.csv(pheno,
    file.path(spatial_dirs$SP2_Human_Visium_CosMx,
              paste0(SP2_ACCESSION, "_phenoData.csv")))
  getGEOSuppFiles(SP2_ACCESSION,
                  baseDir = spatial_dirs$SP2_Human_Visium_CosMx)
  message("[SP-2] Download complete.")
}, error = function(e) message("[SP-2] Error: ", e$message))

# ── Analysis roadmap for spatial datasets ─────────────────────────────────────
message("\n", strrep("=", 60))
message("SPATIAL ANALYSIS ROADMAP:")
message(strrep("=", 60))
message("
SP-3 (GSE282059 — Human CosMx 6000-plex UUO):
  Manuscript role: Figure 7 spatial validation
  Key analyses:
    1. Load CosMx data into Seurat
    2. Annotate urothelial cells using renal pelvis markers
       (UPK1A, UPK2, UPK3A, KRT20, PPARG, FOXA1)
    3. Compare renal pelvis identity score: UUO vs. Healthy
    4. Visualize spatial distribution of renal pelvis TF regulon activity
    5. Show spatial loss of PPARG/FOXA1 signal in UUO tissue sections

SP-4 (GSE253439 — Human CosMx 980-plex Pyelonephritis):
  Manuscript role: Figure 7 extended / Supplementary
  Key analyses:
    1. Load pre-built Seurat RDS object
    2. Focus on pyelonephritis samples (upper urinary tract infection)
    3. Score renal pelvis identity signature
    4. Compare to IgA / MCD controls

SP-1 (Mouse Xenium + Visium — pending accession):
  Manuscript role: Figure 2 / Figure 5 spatial validation
  Key analyses:
    1. Spatially map renal pelvis-specific marker genes in mouse kidney
    2. Confirm that site-specific genes are expressed in papilla/pelvis region
    3. Use Xenium single-cell data to validate cell type assignments
    4. Overlay TF regulon scores onto Visium spatial coordinates

SP-2 (Human Visium + CosMx — pending accession):
  Manuscript role: Figure 4 cross-species spatial validation
  Key analyses:
    1. Identify Uro2-equivalent cells in spatial context
    2. Map PPARG/FOXA1 regulon activity in renal pelvis region
    3. Compare spatial TF programs: healthy vs. CKD

SP-5 (KPMP Visium — manual download):
  Manuscript role: Figure 7 clinical validation
  Key analyses:
    1. Download from atlas.kpmp.org
    2. Score renal pelvis signature on Visium spots
    3. Correlate with histological fibrosis score
")

message("Script complete: ", Sys.time())
