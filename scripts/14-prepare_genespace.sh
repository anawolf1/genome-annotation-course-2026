#!/bin/bash
#SBATCH --job-name=genespace_prep
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00
#SBATCH --output=logs/%x_%j.out 
#SBATCH --error=logs/%x_%j.err

# Build 5 BED/FASTA pairs for GENESPACE: own accession + TAIR10 + 3 accessions from Lian et al.
set -euo pipefail
WORKDIR=/data/users/$USER/genome_annotation_course
COURSEDIR=/data/courses/genome_annotation_resources/2026
ACCESSION=Rab_R1                 # GENESPACE does not allow "-" in names
OTHERS=(Etna_2 Altai_5 Est_0)      # three accessions != your own
FINALDIR="$WORKDIR/output/MAKER/final"
GENES="$FINALDIR/filtered.genes.renamed.gff3"
LONGEST="$FINALDIR/assembly.all.maker.proteins.fasta.renamed.longest.fasta"
GS="$WORKDIR/output/comparisons/GENESPACE"
mkdir -p "$GS/bed" "$GS/peptide" "$GS/tmp"

# Own accession: BED is 0-based (GFF start - 1); FASTA IDs = gene IDs (strip -RA and descriptions)
awk -F'\t' 'BEGIN{OFS="\t"} $3=="gene"{split($9,a,";"); split(a[1],b,"="); print $1,$4-1,$5,b[2]}' \
  "$GENES" > "$GS/bed/${ACCESSION}.bed"
sed '/^>/s/ .*//' "$LONGEST" | sed '/^>/s/-R[A-Z]*$//' > "$GS/peptide/${ACCESSION}.fa"

# TAIR10 reference
cp "$COURSEDIR/data/TAIR10.bed" "$GS/bed/TAIR10.bed"
cp "$COURSEDIR/data/TAIR10.fa"  "$GS/peptide/TAIR10.fa"

# Other accessions (already underscore-named, 0-based BED)
for o in "${OTHERS[@]}"; do
  for f in bed/${o}.bed peptide/${o}.fa; do
    src="$COURSEDIR/data/Lian_et_al/GENESPACE/$f"
    if [[ ! -f "$src" ]]; then
      echo "Missing $src. Available:"; ls "$COURSEDIR/data/Lian_et_al/GENESPACE/bed"; exit 1
    fi
    cp "$src" "$GS/$f"
  done
done

# Validate: no duplicates, BED IDs == FASTA IDs
fail=0
for g in "$ACCESSION" TAIR10 "${OTHERS[@]}"; do
  cut -f4 "$GS/bed/$g.bed" | sort > "$GS/tmp/$g.bed.ids"
  grep '>' "$GS/peptide/$g.fa" | sed 's/^>//; s/ .*//' | sort > "$GS/tmp/$g.fa.ids"
  dup_b=$(uniq -d "$GS/tmp/$g.bed.ids" | wc -l)
  dup_f=$(uniq -d "$GS/tmp/$g.fa.ids"  | wc -l)
  only_b=$(comm -23 "$GS/tmp/$g.bed.ids" "$GS/tmp/$g.fa.ids" | wc -l)
  only_f=$(comm -13 "$GS/tmp/$g.bed.ids" "$GS/tmp/$g.fa.ids" | wc -l)
  echo "$g: genes=$(wc -l < "$GS/tmp/$g.bed.ids") dup_bed=$dup_b dup_fa=$dup_f only_in_bed=$only_b only_in_fa=$only_f"
  (( dup_b + dup_f + only_b + only_f == 0 )) || fail=1
done
(( fail == 0 )) && echo "ALL IDs MATCH" || { echo "ID MISMATCH: fix before running GENESPACE"; exit 1; }