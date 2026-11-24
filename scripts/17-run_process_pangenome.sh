#!/bin/bash
#SBATCH --job-name=process_pangenome
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# BEFORE submitting (once):
#   cp $COURSEDIR/scripts/*process_pangenome.R scripts/17-process_pangenome.R
#   then edit in it:  wd <- ".../output/comparisons/GENESPACE" ;  focal_genome <- "Rab_R1"
# Covers Step 4 (core/accessory/specific) and Step 5 (copy-variable core entries).
set -euo pipefail
WORKDIR=/data/users/$USER/genome_annotation_course
COURSEDIR=/data/courses/genome_annotation_resources/2026
ANALYSIS_SIF="$COURSEDIR/containers/Annotation_analysis_R4.4.sif"
TMPDIR="$WORKDIR/output/comparisons/Annotation_R_tmp"
mkdir -p "$TMPDIR"
cd "$WORKDIR"

apptainer exec --cleanenv \
  --bind /data --bind "$WORKDIR" \
  --bind "$TMPDIR":/temp --env TMPDIR=/temp --pwd "$WORKDIR" \
  "$ANALYSIS_SIF" Rscript "$WORKDIR/scripts/17-process_pangenome.R"