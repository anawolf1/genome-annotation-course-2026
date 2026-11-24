#!/bin/bash
#SBATCH --job-name=busco_maker
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=03:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# BUSCO completeness (brassicales_odb10) of the longest protein and transcript per gene
set -euo pipefail
FINALDIR=/data/users/$USER/genome_annotation_course/output/MAKER/final
protein=assembly.all.maker.proteins.fasta
transcript=assembly.all.maker.transcripts.fasta
cd $FINALDIR
module load BUSCO/5.4.2-foss-2021a
export PYTHONNOUSERSITE=1
busco -i ${protein}.renamed.longest.fasta -l brassicales_odb10 -o busco_proteins -m proteins -c 4 -f
busco -i ${transcript}.renamed.longest.fasta -l brassicales_odb10 -o busco_transcripts -m transcriptome -c 4 -f
