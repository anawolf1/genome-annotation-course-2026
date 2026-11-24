#!/bin/bash
#SBATCH --job-name=riparian
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Usage: sbatch scripts/18-run_riparian.sh                      (first pass)
#        sbatch scripts/18-run_riparian.sh "scaffold_3,scaffold_7"   (flip reversed scaffolds)
set -euo pipefail
WORKDIR=/data/users/$USER/genome_annotation_course
COURSEDIR=/data/courses/genome_annotation_resources/2026
GENESPACE_DIR="$WORKDIR/output/comparisons/GENESPACE"
FOCAL=Rab_R1
FLIPS="${1:-}"
cd "$WORKDIR"

apptainer exec --cleanenv --bind /data --bind "$WORKDIR" --pwd "$WORKDIR" \
  "$COURSEDIR/containers/genespace_latest.sif" \
  Rscript "$WORKDIR/scripts/18-riparian.R" "$GENESPACE_DIR" "$FOCAL" "$FLIPS"
