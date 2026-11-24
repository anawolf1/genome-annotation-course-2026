#!/bin/bash
#SBATCH --job-name=TE_plots
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# TE density along contigs (circlize) and TE divergence landscape per superfamily
set -euo pipefail
WORKDIR=/data/users/$USER/genome_annotation_course
COURSEDIR=/data/courses/genome_annotation_resources/2026
SIF=$COURSEDIR/containers/Annotation_analysis_R4.4.sif
cd $WORKDIR/output/EDTA
mkdir -p plots

run() { apptainer exec -C --bind $WORKDIR --bind $COURSEDIR --pwd "$PWD" $SIF "$@"; }

run samtools faidx assembly.fasta          # -> assembly.fasta.fai (contig lengths)
run Rscript --vanilla $WORKDIR/scripts/04-annotation_circlize.R
mv -f 02-TE_density.pdf plots/TE_density_circlize.pdf
run Rscript --vanilla $WORKDIR/scripts/06-plot_div.R
