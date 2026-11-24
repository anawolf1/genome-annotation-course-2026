#!/bin/bash
#SBATCH --job-name=maker_merge
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
# Merge per-contig MAKER results into genome-wide GFF3 and FASTA files
WORKDIR=/data/users/$USER/genome_annotation_course
SIF=/data/courses/genome_annotation_resources/2026/containers/MAKER_3.01.03.sif
cd $WORKDIR/output/MAKER
IDX=$(ls -d *.maker.output)/$(basename $(ls -d *.maker.output) .maker.output)_master_datastore_index.log
apptainer exec --bind /data $SIF gff3_merge -s -d $IDX > assembly.all.maker.gff
apptainer exec --bind /data $SIF gff3_merge -n -s -d $IDX > assembly.all.maker.noseq.gff
apptainer exec --bind /data $SIF fasta_merge -d $IDX -o assembly
