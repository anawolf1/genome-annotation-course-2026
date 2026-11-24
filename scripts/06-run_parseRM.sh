#!/bin/bash
#SBATCH --job-name=parseRM
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Bin RepeatMasker divergence (1% bins) per TE superfamily, excluding 0% self-matches
set -euo pipefail
WORKDIR=/data/users/$USER/genome_annotation_course
COURSEDIR=/data/courses/genome_annotation_resources/2026
genome=assembly.fasta
cd $WORKDIR/output/EDTA/$genome.mod.EDTA.anno

awk 'NR == 1 || $2 != 0' $genome.mod.EDTA.RM.out > $genome.mod.EDTA.RM.nonzero.out

apptainer exec -C --bind $WORKDIR --bind $COURSEDIR --pwd "$PWD" \
  $COURSEDIR/containers/Annotation_analysis_R4.4.sif \
  perl $COURSEDIR/scripts/05-parseRM.pl -i "$PWD/$genome.mod.EDTA.RM.nonzero.out" -l 50,1 -v
