#!/bin/bash
#SBATCH --job-name=Download_SpatialGEOquery
#SBATCH --account=gdjacksonlab
#SBATCH --partition=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=slurm_spatial_geoquery_%j.out
#SBATCH --error=slurm_spatial_geoquery_%j.err

module purge
ml GCC/9.3.0
ml OpenMPI/4.0.3
ml R/4.4.0

SCRIPT_DIR="/vast0/home/gdjacksonlab/lab/xxw004/UUO/Scripts/KUDO/Download"
Rscript "${SCRIPT_DIR}/Download_SpatialDatasets_GEOquery.R"
