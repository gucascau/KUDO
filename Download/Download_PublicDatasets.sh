#!/bin/bash
#SBATCH --job-name=Download_PublicDatasets
#SBATCH --account=gdjacksonlab
#SBATCH --partition=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=slurm_download_%j.out
#SBATCH --error=slurm_download_%j.err

# =============================================================================
# Download_PublicDatasets.sh
# Author: Xin Wang — xin.wang@nationwidechildrens.org
# Purpose: Download all public datasets for the renal pelvis urothelium
#          single-cell atlas publication.
#
# Datasets:
#   Priority 1: Li et al. Cell Metab 2024 (GSE234788) — SHARE-seq human kidney
#   Priority 2: KPMP snRNA-seq v1.0 + Multiome v1.3  — see manual instructions
#   Priority 3: Lake et al. Nature 2023               — via CellxGene portal
#   Priority 3: GUDMAP Human Bladder (GSE157079)
#   Priority 4: Bhatt et al. JCI 2019 (GSE129845)    — human + mouse bladder
#   Priority 4: GUDMAP Mouse Urinary Tract (GSE119531)
# =============================================================================

set -euo pipefail

# ── Paths ─────────────────────────────────────────────────────────────────────
DATASETS_DIR="/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets"
HUMAN_DIR="${DATASETS_DIR}/Human"
MOUSE_DIR="${DATASETS_DIR}/Mouse"

# ── Create dataset-specific folders ───────────────────────────────────────────
echo "[$(date)] Creating dataset directories..."

mkdir -p "${HUMAN_DIR}/Li_CellMetab2024_SHAREseq_GSE234788"
mkdir -p "${HUMAN_DIR}/KPMP_snRNAseq_v1"
mkdir -p "${HUMAN_DIR}/KPMP_Multiome_v1.3"
mkdir -p "${HUMAN_DIR}/Lake_Nature2023_KidneyAtlas"
mkdir -p "${HUMAN_DIR}/GUDMAP_HumanBladder_GSE157079"
mkdir -p "${HUMAN_DIR}/Bhatt_JCI2019_HumanMouseBladder_GSE129845"
mkdir -p "${MOUSE_DIR}/GUDMAP_MouseUrinaryTract_GSE119531"

echo "[$(date)] Directories created."

# ── Helper function: download GEO supplementary files via FTP ─────────────────
download_geo_suppl() {
    local ACCESSION=$1
    local OUTDIR=$2
    local PREFIX="${ACCESSION:0:${#ACCESSION}-3}nnn"   # e.g. GSE234788 → GSE234nnn

    local FTP_URL="https://ftp.ncbi.nlm.nih.gov/geo/series/${PREFIX}/${ACCESSION}/suppl/"

    echo ""
    echo "============================================================"
    echo "[$(date)] Downloading ${ACCESSION} supplementary files"
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
        && echo "[$(date)] ${ACCESSION} download complete." \
        || echo "[$(date)] WARNING: ${ACCESSION} download may be incomplete — check ${OUTDIR}"
}

# ── Priority 1: Li et al. Cell Metab 2024 — SHARE-seq human kidney ────────────
# Explicitly identifies Uro1 (ureter) and Uro2 (renal pelvis) urothelial cells
# Contains simultaneous snRNA-seq + snATAC-seq from 54 human kidney samples
# NOTE: Verify accession GSE234788 in paper's Data Availability section:
#       https://www.cell.com/cell-metabolism/fulltext/S1550-4131(24)00061-5
download_geo_suppl \
    "GSE234788" \
    "${HUMAN_DIR}/Li_CellMetab2024_SHAREseq_GSE234788"

# ── Priority 3: GUDMAP Human Bladder Cell Atlas ────────────────────────────────
# 29,834 cells; scRNA-seq + snRNA-seq from 4 bladder regions
# Use for bladder vs. renal pelvis site-specificity comparison
download_geo_suppl \
    "GSE157079" \
    "${HUMAN_DIR}/GUDMAP_HumanBladder_GSE157079"

# ── Priority 4: Bhatt et al. JCI 2019 — Human + Mouse Bladder Map ─────────────
# Cross-species ortholog framework; well-cited urothelial layer reference
download_geo_suppl \
    "GSE129845" \
    "${HUMAN_DIR}/Bhatt_JCI2019_HumanMouseBladder_GSE129845"

# ── Priority 4: GUDMAP Mouse Urinary Tract Atlas ──────────────────────────────
# Mouse urinary tract scRNA-seq across anatomical segments
# Allows within-species comparison: renal pelvis vs. ureter vs. bladder
download_geo_suppl \
    "GSE119531" \
    "${MOUSE_DIR}/GUDMAP_MouseUrinaryTract_GSE119531"

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo "[$(date)] GEO downloads complete."
echo ""
echo "MANUAL DOWNLOAD REQUIRED for the following datasets:"
echo ""
echo "  [Priority 2] KPMP snRNA-seq v1.0"
echo "    DOI  : 10.48698/6k3k-q779"
echo "    URL  : https://www.kpmp.org/doi-collection/10-48698-6k3k-q779"
echo "    Dest : ${HUMAN_DIR}/KPMP_snRNAseq_v1"
echo "    Steps: 1. Visit URL above"
echo "           2. Click 'Download' on the dataset files listed"
echo "           3. Move downloaded files to Dest folder above"
echo ""
echo "  [Priority 2] KPMP Multiome v1.3 (snRNA-seq + snATAC-seq)"
echo "    DOI  : 10.48698/92nk-e805"
echo "    URL  : https://www.kpmp.org/doi-collection/10-48698-92nk-e805"
echo "    Dest : ${HUMAN_DIR}/KPMP_Multiome_v1.3"
echo "    Steps: Same as above"
echo ""
echo "  [Priority 3] Lake et al. Nature 2023 — Multimodal Human Kidney Atlas"
echo "    URL  : https://cellxgene.cziscience.com  (search 'Lake kidney 2023')"
echo "    URL  : https://www.nature.com/articles/s41586-023-05769-3 (Data Avail.)"
echo "    Dest : ${HUMAN_DIR}/Lake_Nature2023_KidneyAtlas"
echo "    Steps: Download .h5ad from CellxGene or follow paper's data link"
echo ""
echo "  [Supplementary] Human Protein Atlas"
echo "    URL  : https://www.proteinatlas.org"
echo "    Use  : Download IHC images for UPK1A / UPK2 / KRT5 / KRT20 / KRT14"
echo "           Filter: Tissue = kidney, urinary bladder, ureter"
echo ""
echo "  [Supplementary] GUDMAP Expression Portal"
echo "    URL  : https://www.gudmap.org"
echo "    Use  : Browse ISH/IHC for PPARG / FOXA1 / GATA3 / TP63 in urinary tract"
echo "============================================================"
