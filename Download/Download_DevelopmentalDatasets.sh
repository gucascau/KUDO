#!/bin/bash
#SBATCH --job-name=Download_DevDatasets
#SBATCH --account=gdjacksonlab
#SBATCH --partition=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=slurm_download_dev_%j.out
#SBATCH --error=slurm_download_dev_%j.err

# =============================================================================
# Download_DevelopmentalDatasets.sh
# Author: Xin Wang — xin.wang@nationwidechildrens.org
# Purpose: Download developmental kidney / urinary tract single-cell datasets
#          to anchor the embryonic origin of renal pelvis urothelial identity.
#
# Scientific rationale:
#   The overlapping basal cells between KUDO (renal pelvis) and bladder organoid
#   represent a shared embryonic progenitor derived from the ureteric bud.
#   These datasets establish when and how renal pelvis identity diverges from
#   bladder (cloacal/urogenital sinus) identity during development.
#
# Datasets:
#   DEV-1: Lindström et al. Cell Rep 2018    (GSE108291) — Mouse kidney dev
#   DEV-2: Cao et al. Science 2019           (GSE119945) — Mouse organogenesis
#   DEV-3: Menon et al. Dev Cell 2022        (GSE185477) — Human fetal kidney
#   DEV-4: GUDMAP Development Portal         (manual)    — ISH/IHC validation
# =============================================================================

set -euo pipefail

# ── Paths ─────────────────────────────────────────────────────────────────────
DATASETS_DIR="/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets"
MOUSE_DIR="${DATASETS_DIR}/Mouse"
HUMAN_DIR="${DATASETS_DIR}/Human"

# ── Create developmental dataset folders ──────────────────────────────────────
echo "[$(date)] Creating developmental dataset directories..."

mkdir -p "${MOUSE_DIR}/Lindstrom_CellRep2018_KidneyDev_GSE108291"
mkdir -p "${MOUSE_DIR}/Cao_Science2019_MouseOrganogenesis_GSE119945"
mkdir -p "${HUMAN_DIR}/Menon_DevCell2022_FetalKidney_GSE185477"
mkdir -p "${DATASETS_DIR}/GUDMAP_UrinaryTractDev_Portal"

echo "[$(date)] Directories created."

# ── Helper: download GEO supplementary files via FTP ──────────────────────────
download_geo_suppl() {
    local ACCESSION=$1
    local OUTDIR=$2
    local PREFIX="${ACCESSION:0:${#ACCESSION}-3}nnn"

    local FTP_URL="https://ftp.ncbi.nlm.nih.gov/geo/series/${PREFIX}/${ACCESSION}/suppl/"

    echo ""
    echo "============================================================"
    echo "[$(date)] Downloading ${ACCESSION}"
    echo "  FTP : ${FTP_URL}"
    echo "  Dest: ${OUTDIR}"
    echo "============================================================"

    wget \
        --recursive \
        --no-parent \
        --no-host-directories \
        --cut-dirs=6 \
        --directory-prefix="${OUTDIR}" \
        --reject "index.html*" \
        --tries=5 \
        --wait=2 \
        "${FTP_URL}" \
        && echo "[$(date)] ${ACCESSION} complete." \
        || echo "[$(date)] WARNING: ${ACCESSION} may be incomplete — check ${OUTDIR}"
}

# ── DEV-1: Lindström et al. Cell Rep 2018 — Mouse Kidney Development ──────────
# Most directly relevant: ureteric bud / ureteric tip / trunk lineages
# Ureteric tip → collecting duct → renal pelvis progenitors
# Use to show shared basal cells = remnant of ureteric bud progenitor
echo ""
echo ">>> DEV-1: Lindström et al. Mouse Kidney Development (GSE108291)"
download_geo_suppl \
    "GSE108291" \
    "${MOUSE_DIR}/Lindstrom_CellRep2018_KidneyDev_GSE108291"

# ── DEV-2: Cao et al. Science 2019 — Mouse Organogenesis Cell Atlas (MOCA) ────
# Whole-embryo atlas E9.5–E13.5 (~2M cells)
# Subset to kidney/urinary tract lineage after download
# WARNING: Very large dataset — supplementary files may be multi-GB
echo ""
echo ">>> DEV-2: Cao et al. Mouse Organogenesis Cell Atlas (GSE119945)"
echo "    NOTE: This is a very large dataset. Consider downloading only"
echo "    the processed count matrix and metadata, not raw FASTQ."
download_geo_suppl \
    "GSE119945" \
    "${MOUSE_DIR}/Cao_Science2019_MouseOrganogenesis_GSE119945"

# ── DEV-3: Menon et al. Dev Cell 2022 — Human Fetal Kidney ───────────────────
# Human counterpart to Lindström
# Captures ureteric bud + collecting duct + early urothelial cells (wk 8–23)
# Cross-species developmental validation
echo ""
echo ">>> DEV-3: Menon et al. Human Fetal Kidney (GSE185477)"
download_geo_suppl \
    "GSE185477" \
    "${HUMAN_DIR}/Menon_DevCell2022_FetalKidney_GSE185477"

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo "[$(date)] GEO developmental downloads complete."
echo ""
echo "MANUAL DOWNLOAD REQUIRED:"
echo ""
echo "  [DEV-4] GUDMAP Mouse Urinary Tract Development Portal"
echo "    URL  : https://www.gudmap.org"
echo "    Dest : ${DATASETS_DIR}/GUDMAP_UrinaryTractDev_Portal"
echo "    Steps:"
echo "      1. Go to gudmap.org → Data → Gene Expression"
echo "      2. Search genes of interest: PPARG, FOXA1, GATA3, TP63, KRT5"
echo "      3. Filter by anatomy: renal pelvis, ureter, bladder"
echo "      4. Filter by stage: E11.5, E13.5, E15.5, P0, Adult"
echo "      5. Download ISH/IHC images and expression tables"
echo "      6. Save to: ${DATASETS_DIR}/GUDMAP_UrinaryTractDev_Portal"
echo ""
echo "    Key gene panels to download for manuscript figures:"
echo "      Basal TFs  : TP63, KRT5, KRT14, SOX2"
echo "      Renal pelvis TFs: PPARG, FOXA1, GATA3, GRHL3"
echo "      Differentiation: UPK1A, UPK2, UPK3A, KRT20"
echo "============================================================"

echo "[$(date)] Script finished."
