#!/bin/bash
#SBATCH --time=30:00:00
#SBATCH --partition=himem
#SBATCH --cpus-per-task=24
#SBATCH --mail-type=FAIL,REQUEUE
#SBATCH --account=gdjacksonlab
#SBATCH --mail-user=xin.wang@nationwidechildrens.org
#SBATCH --job-name=LargeD
#SBATCH --output=slurm_%j.out
set -e

# loading the module

module purge
ml GCC/9.3.0
ml OpenMPI/4.0.3
ml R/4.4.0

CWD=/home/gdbecknelllab/xxw004/gdjacksonlab/UUO/Scripts/KUDO/Integration/

cd $CWD
Rscript  $CWD/Urothelium_KUDO_GOFLOFVehicle_Integration.R
