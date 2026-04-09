#!/usr/bin/env Rscript
# =============================================================================
# Download_DevelopmentalDatasets_GEOquery.R
# Author: Xin Wang — xin.wang@nationwidechildrens.org
# Purpose: Download and inspect developmental kidney / urinary tract datasets
#          using GEOquery. Provides metadata inspection + file inventory.
#
# Scientific rationale:
#   The overlapping basal cells between KUDO (renal pelvis) and bladder
#   organoid represent a shared embryonic progenitor (ureteric bud lineage).
#   These datasets establish when renal pelvis identity diverges from bladder
#   (cloacal) identity during development.
#
# Datasets:
#   DEV-1: Lindström et al. Cell Rep 2018  (GSE108291) — Mouse kidney dev
#   DEV-2: Cao et al. Science 2019         (GSE119945) — Mouse organogenesis
#   DEV-3: Menon et al. Dev Cell 2022      (GSE185477) — Human fetal kidney
# =============================================================================

# ── Packages ──────────────────────────────────────────────────────────────────
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("GEOquery",    quietly = TRUE)) BiocManager::install("GEOquery")

library(GEOquery)
library(tidyverse)

options(timeout = 7200)   # 2-hour timeout — MOCA is very large

# ── Paths ─────────────────────────────────────────────────────────────────────
DATASETS_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets"

dev_dirs <- list(
  GSE108291 = file.path(DATASETS_DIR, "Mouse/Lindstrom_CellRep2018_KidneyDev_GSE108291"),
  GSE119945 = file.path(DATASETS_DIR, "Mouse/Cao_Science2019_MouseOrganogenesis_GSE119945"),
  GSE185477 = file.path(DATASETS_DIR, "Human/Menon_DevCell2022_FetalKidney_GSE185477")
)

invisible(lapply(dev_dirs, dir.create, showWarnings = FALSE, recursive = TRUE))

# ── Dataset metadata ───────────────────────────────────────────────────────────
dev_meta <- list(
  GSE108291 = list(
    name     = "Lindström et al. Mouse Kidney Development",
    journal  = "Cell Reports 2018",
    priority = "DEV-1 — PRIMARY developmental dataset",
    use      = "Ureteric bud lineage → renal pelvis progenitors",
    note     = "Extract ureteric tip/trunk cells; compare TF programs to KUDO"
  ),
  GSE119945 = list(
    name     = "Cao et al. Mouse Organogenesis Cell Atlas (MOCA)",
    journal  = "Science 2019",
    priority = "DEV-2 — Supporting broad context",
    use      = "Whole-embryo: trace ureteric bud emergence at E9.5–E13.5",
    note     = "VERY LARGE (~2M cells). Subset to kidney/urinary tract after download."
  ),
  GSE185477 = list(
    name     = "Menon et al. Human Fetal Kidney",
    journal  = "Developmental Cell 2022",
    priority = "DEV-3 — Human cross-species developmental validation",
    use      = "Human ureteric bud / collecting duct / early urothelial cells (wk 8-23)",
    note     = "Cross-species closure: mouse dev (DEV-1) → human dev (DEV-3) → adult human (Li 2024)"
  )
)

# ── Helper: fetch metadata + supplementary files ───────────────────────────────
fetch_geo_dataset <- function(accession, outdir, meta) {
  message("\n", strrep("=", 60))
  message("Dataset  : ", meta$name)
  message("Accession: ", accession, " (", meta$journal, ")")
  message("Priority : ", meta$priority)
  message("Use      : ", meta$use)
  message("Note     : ", meta$note)
  message("Dest     : ", outdir)
  message(strrep("=", 60))

  # --- 1. Metadata ---
  tryCatch({
    gse <- getGEO(accession, destdir = outdir, GSEMatrix = TRUE, AnnotGPL = FALSE)
    pheno_df <- pData(gse[[1]])

    write.csv(pheno_df,
              file.path(outdir, paste0(accession, "_phenoData.csv")),
              row.names = TRUE)

    message("[", accession, "] Samples: ", nrow(pheno_df))
    message("[", accession, "] phenoData saved.")

    # Print sample summary
    if ("title" %in% colnames(pheno_df)) {
      message("[", accession, "] Sample titles:")
      print(head(pheno_df$title, 10))
    }

  }, error = function(e) {
    message("[", accession, "] WARNING — metadata failed: ", e$message)
  })

  # --- 2. Supplementary files ---
  tryCatch({
    message("[", accession, "] Downloading supplementary files...")
    suppl <- getGEOSuppFiles(accession, baseDir = outdir, fetch_files = TRUE)
    message("[", accession, "] Files downloaded:")
    print(suppl)
  }, error = function(e) {
    message("[", accession, "] WARNING — suppl download failed: ", e$message)
    message("  Fallback: run Download_DevelopmentalDatasets.sh (wget-based)")
  })

  message("[", accession, "] Done.\n")
}

# ── DEV-1: Lindström et al. — Mouse Kidney Development ────────────────────────
fetch_geo_dataset("GSE108291", dev_dirs$GSE108291, dev_meta$GSE108291)

# ── DEV-2: Cao et al. — Mouse Organogenesis Cell Atlas ────────────────────────
# WARNING: This is a very large dataset. If storage is limited, comment this
# out and download only the processed count matrix manually from:
# https://oncoscape.v3.sttrcancer.org/atlas.gs.washington.edu.mouse.rna/downloads
message("\n>>> DEV-2 WARNING: MOCA is ~10–20 GB. Proceeding...")
fetch_geo_dataset("GSE119945", dev_dirs$GSE119945, dev_meta$GSE119945)

# ── DEV-3: Menon et al. — Human Fetal Kidney ──────────────────────────────────
fetch_geo_dataset("GSE185477", dev_dirs$GSE185477, dev_meta$GSE185477)

# ── File inventory ─────────────────────────────────────────────────────────────
message("\n", strrep("=", 60))
message("Downloaded file inventory:")
message(strrep("=", 60))

for (acc in names(dev_dirs)) {
  d <- dev_dirs[[acc]]
  files <- list.files(d, recursive = TRUE, full.names = FALSE)
  total_mb <- sum(file.size(file.path(d, files)), na.rm = TRUE) / 1e6
  message("\n[", acc, "] ", length(files), " file(s) | Total: ",
          round(total_mb, 1), " MB | Dir: ", d)
  if (length(files) > 0 && length(files) <= 20) {
    df <- data.frame(
      File    = files,
      Size_MB = round(file.size(file.path(d, files)) / 1e6, 2)
    )
    print(df)
  } else if (length(files) > 20) {
    message("  (", length(files), " files — too many to list. Check directory.)")
  }
}

# ── Next steps ─────────────────────────────────────────────────────────────────
message("\n", strrep("=", 60))
message("NEXT STEPS after download:")
message(strrep("=", 60))
message("
DEV-1 (GSE108291 — Lindström):
  → Load count matrix into Seurat
  → Subset to ureteric bud / ureteric tip / trunk cells
  → Run decoupleR TF activity
  → Compare TF regulons to KUDO Vehicle (adult renal pelvis)
  → Key question: are adult renal pelvis TFs already active in ureteric bud?

DEV-2 (GSE119945 — Cao MOCA):
  → Subset to kidney lineage cells (annotated as 'Kidney' or 'Urinary tract')
  → Use as broad context for ureteric bud emergence timing
  → Focus on E11.5–E13.5 window (ureteric bud branching stage)

DEV-3 (GSE185477 — Menon human fetal):
  → Subset to ureteric bud / collecting duct / urothelial cells
  → Cross-species comparison: mouse ureteric bud TFs vs. human fetal TFs
  → Bridge: human fetal (DEV-3) → human adult Uro2 (Li et al. 2024)
  → Closes the developmental → adult → disease loop in human

DEV-4 (GUDMAP portal — manual):
  → Browse: gudmap.org → Gene Expression
  → Key genes: PPARG, FOXA1, GATA3, TP63, KRT5, UPK1A, KRT20
  → Filter: renal pelvis + ureter + bladder across E11.5 to adult
  → Download ISH images for manuscript supplementary figures
")

message("Script complete: ", Sys.time())
