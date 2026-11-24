#!/bin/bash
#SBATCH --job-name=LTR_identity_plot
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Plot LTR identity (proxy for insertion age) of intact LTR-RTs per TEsorter clade
WORKDIR=/data/users/$USER/genome_annotation_course
COURSEDIR=/data/courses/genome_annotation_resources/2026
cd "$WORKDIR/output/LTR_full_length"
mkdir -p plots
apptainer exec -C --bind "$WORKDIR" --bind "$COURSEDIR" --pwd "$PWD" \
  "$COURSEDIR/containers/Annotation_analysis_R4.4.sif" \
  Rscript --vanilla "$WORKDIR/scripts/03-full_length_LTRs_identity_copy.R"
