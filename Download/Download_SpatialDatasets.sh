#!/bin/bash
#SBATCH --job-name=Download_SpatialDatasets
#SBATCH --account=gdjacksonlab
#SBATCH --partition=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=slurm_spatial_%j.out
#SBATCH --error=slurm_spatial_%j.err

# =============================================================================
# Download_SpatialDatasets.sh
# Author: Xin Wang — xin.wang@nationwidechildrens.org
# Purpose: Download spatial transcriptomics datasets for validation of
#          renal pelvis urothelial identity in mouse and human tissue.
#
# NOTE: SP-3 (GSE282059) and SP-4 (GSE253439) are ALREADY DOWNLOADED.
#       This script handles the remaining spatial datasets.
#
# Spatial datasets:
#   SP-1: Xuanyuan et al. Nat Commun 2025   — Mouse Xenium + Visium (IRI)
#   SP-2: Muto et al. Nat Genet 2024        — Human Visium + CosMx (CKD)
#   SP-5: KPMP Visium                        — Human Visium (manual)
#   SP-6: 10x Visium HD Mouse Kidney FFPE   — Mouse reference (direct)
# =============================================================================

set -euo pipefail

DATASETS_DIR="/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets"
HUMAN_DIR="${DATASETS_DIR}/Human"
MOUSE_DIR="${DATASETS_DIR}/Mouse"
SPATIAL_DIR="${DATASETS_DIR}/Spatial"

# ── Create spatial dataset folders ────────────────────────────────────────────
echo "[$(date)] Creating spatial dataset directories..."

mkdir -p "${MOUSE_DIR}/Xuanyuan_NatCommun2025_MouseKidney_Xenium_Visium"
mkdir -p "${HUMAN_DIR}/Muto_NatGenet2024_HumanKidney_Visium_CosMx"
mkdir -p "${HUMAN_DIR}/KPMP_Visium_Spatial"
mkdir -p "${MOUSE_DIR}/10xGenomics_VisiumHD_MouseKidney_FFPE"
mkdir -p "${SPATIAL_DIR}/Already_Downloaded"

# Symlink already-downloaded spatial datasets for easy access
ln -sfn "${HUMAN_DIR}/GSE282059" \
    "${SPATIAL_DIR}/Already_Downloaded/SP3_CosMx6000_HumanKidney_UUO_GSE282059" 2>/dev/null || true
ln -sfn "${HUMAN_DIR}/GSE253439" \
    "${SPATIAL_DIR}/Already_Downloaded/SP4_CosMx980_HumanKidney_IgA_MCD_Pyelonephritis_GSE253439" 2>/dev/null || true

echo "[$(date)] Directories and symlinks created."

# ── Helper: GEO FTP download ──────────────────────────────────────────────────
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

    wget --recursive --no-parent --no-host-directories \
         --cut-dirs=6 --directory-prefix="${OUTDIR}" \
         --reject "index.html*" --tries=5 --wait=2 \
         "${FTP_URL}" \
        && echo "[$(date)] ${ACCESSION} complete." \
        || echo "[$(date)] WARNING: ${ACCESSION} may be incomplete."
}

# ── SP-1: Xuanyuan et al. Nat Commun 2025 — Mouse Xenium + Visium ─────────────
# Xenium (300-gene panel) + Visium CytAssist (whole transcriptome)
# 12 mouse kidneys × 6 IRI timepoints; >1 million cells
# Use to spatially validate mouse renal pelvis marker genes and TF programs
# NOTE: Verify GEO accession at:
#       https://www.nature.com/articles/s41467-025-62599-9 (Data Availability)
SP1_ACCESSION="VERIFY_IN_PAPER"   # <-- replace with confirmed GEO accession
SP1_DIR="${MOUSE_DIR}/Xuanyuan_NatCommun2025_MouseKidney_Xenium_Visium"

echo ""
echo ">>> SP-1: Xuanyuan et al. Nat Commun 2025 — Mouse Xenium + Visium"
echo "    IMPORTANT: Verify the GEO accession number before running."
echo "    Check Data Availability at:"
echo "    https://www.nature.com/articles/s41467-025-62599-9"
echo "    Then replace SP1_ACCESSION in this script and re-run."
echo "    Dest: ${SP1_DIR}"

# Uncomment after verifying accession:
# download_geo_suppl "${SP1_ACCESSION}" "${SP1_DIR}"

# ── SP-2: Muto et al. Nat Genet 2024 — Human Visium + CosMx ──────────────────
# 81 samples: snRNA + snATAC + Visium + CosMx
# Healthy + diabetic + hypertensive CKD
# Use to spatially validate renal pelvis TF programs in human CKD tissue
# NOTE: Verify GEO accession at:
#       https://www.nature.com/articles/s41588-024-01802-x (Data Availability)
SP2_ACCESSION="VERIFY_IN_PAPER"   # <-- replace with confirmed GEO accession
SP2_DIR="${HUMAN_DIR}/Muto_NatGenet2024_HumanKidney_Visium_CosMx"

echo ""
echo ">>> SP-2: Muto et al. Nat Genet 2024 — Human Visium + CosMx"
echo "    IMPORTANT: Verify the GEO accession number before running."
echo "    Check Data Availability at:"
echo "    https://www.nature.com/articles/s41588-024-01802-x"
echo "    Then replace SP2_ACCESSION in this script and re-run."
echo "    Dest: ${SP2_DIR}"

# Uncomment after verifying accession:
# download_geo_suppl "${SP2_ACCESSION}" "${SP2_DIR}"

# ── SP-6: 10x Genomics Visium HD Mouse Kidney FFPE ────────────────────────────
# High-resolution reference dataset (2µm — near single-cell)
# Free download from 10x Genomics website
# Use as mouse kidney spatial reference for renal pelvis region boundaries
echo ""
echo ">>> SP-6: 10x Genomics Visium HD Mouse Kidney FFPE"
SP6_DIR="${MOUSE_DIR}/10xGenomics_VisiumHD_MouseKidney_FFPE"

VISIUMHD_BASE="https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Mouse_Kidney"

echo "  Downloading Visium HD Mouse Kidney processed files..."
wget -P "${SP6_DIR}" --tries=3 \
    "${VISIUMHD_BASE}/Visium_HD_Mouse_Kidney_filtered_feature_bc_matrix.h5" \
    && echo "  filtered_feature_bc_matrix.h5 downloaded." \
    || echo "  WARNING: Check URL at 10xgenomics.com/datasets"

wget -P "${SP6_DIR}" --tries=3 \
    "${VISIUMHD_BASE}/Visium_HD_Mouse_Kidney_spatial.tar.gz" \
    && echo "  spatial.tar.gz downloaded." \
    || echo "  WARNING: Check URL at 10xgenomics.com/datasets"

# Extract spatial files
if [ -f "${SP6_DIR}/Visium_HD_Mouse_Kidney_spatial.tar.gz" ]; then
    tar -xzf "${SP6_DIR}/Visium_HD_Mouse_Kidney_spatial.tar.gz" -C "${SP6_DIR}"
    echo "  Spatial files extracted."
fi

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo "[$(date)] Spatial download script complete."
echo ""
echo "STATUS SUMMARY:"
echo ""
echo "  ALREADY DOWNLOADED (ready to analyze):"
echo "    SP-3: GSE282059 — CosMx 6000-plex Human Kidney UUO + Healthy"
echo "          ${HUMAN_DIR}/GSE282059"
echo "    SP-4: GSE253439 — CosMx 980-plex Human Kidney IgA/MCD/Pyelonephritis"
echo "          ${HUMAN_DIR}/GSE253439"
echo ""
echo "  DOWNLOADED NOW:"
echo "    SP-6: 10x Visium HD Mouse Kidney FFPE"
echo "          ${SP6_DIR}"
echo ""
echo "  MANUAL STEPS REQUIRED:"
echo ""
echo "  [SP-1] Xuanyuan et al. Nat Commun 2025 — Mouse Xenium + Visium"
echo "    Step 1: Visit https://www.nature.com/articles/s41467-025-62599-9"
echo "    Step 2: Find GEO accession in Data Availability section"
echo "    Step 3: Update SP1_ACCESSION in this script and re-run"
echo "    Dest  : ${SP1_DIR}"
echo ""
echo "  [SP-2] Muto et al. Nat Genet 2024 — Human Visium + CosMx"
echo "    Step 1: Visit https://www.nature.com/articles/s41588-024-01802-x"
echo "    Step 2: Find GEO accession in Data Availability section"
echo "    Step 3: Update SP2_ACCESSION in this script and re-run"
echo "    Dest  : ${SP2_DIR}"
echo ""
echo "  [SP-5] KPMP Visium Spatial Transcriptomics"
echo "    Step 1: Visit https://atlas.kpmp.org/repository/"
echo "    Step 2: Filter Data Type = Visium"
echo "    Step 3: Download relevant tissue files"
echo "    Dest  : ${HUMAN_DIR}/KPMP_Visium_Spatial"
echo "============================================================"
