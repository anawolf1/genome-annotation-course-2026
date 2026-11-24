#!/bin/bash
#SBATCH --job-name=TEsorter_TElib
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Classify all Copia and Gypsy sequences of the EDTA TE library into clades (REXdb plant)
set -euo pipefail
WORKDIR=/data/users/$USER/genome_annotation_course
COURSEDIR=/data/courses/genome_annotation_resources/2026
TELIB=$WORKDIR/output/EDTA/assembly.fasta.mod.EDTA.TElib.fa
OUT=$WORKDIR/output/TEsorter_TElib
mkdir -p $OUT && cd $OUT

R_SIF=$COURSEDIR/containers/Annotation_analysis_R4.4.sif
TS_SIF=$COURSEDIR/containers/TEsorter_1.3.0.sif

apptainer exec -C --bind $WORKDIR $R_SIF seqkit grep -r -p "Copia" $TELIB > Copia_sequences.fa
apptainer exec -C --bind $WORKDIR $R_SIF seqkit grep -r -p "Gypsy" $TELIB > Gypsy_sequences.fa

for sf in Copia Gypsy; do
  apptainer exec --bind $WORKDIR --pwd $OUT $TS_SIF \
    TEsorter ${sf}_sequences.fa -db rexdb-plant -p $SLURM_CPUS_PER_TASK
done
