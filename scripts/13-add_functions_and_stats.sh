#!/bin/bash
#SBATCH --job-name=uniprot_functions
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Add UniProt functions to FASTA/GFF3, then summary stats and FLC (AT5G10140) lookup
set -euo pipefail
WORKDIR=/data/users/$USER/genome_annotation_course
COURSEDIR=/data/courses/genome_annotation_resources/2026
FINALDIR="$WORKDIR/output/MAKER/final"
OUTDIR="$WORKDIR/output/comparisons"
PROTEINS="$FINALDIR/assembly.all.maker.proteins.fasta.renamed.filtered.fasta"
UNIPROT="$COURSEDIR/data/uniprot_viridiplantae_reviewed.fa"
GENES="$FINALDIR/filtered.genes.renamed.gff3"
MAKER_SIF="$COURSEDIR/containers/MAKER_3.01.03.sif"
ANALYSIS_SIF="$COURSEDIR/containers/Annotation_analysis_R4.4.sif"
mkdir -p "$OUTDIR/tmp"
cd "$OUTDIR"

mk()  { apptainer exec --cleanenv --env TMPDIR=/tmp --bind /data --bind "$WORKDIR" --bind "$OUTDIR/tmp":/tmp --pwd "$OUTDIR" "$MAKER_SIF" "$@"; }
ana() { apptainer exec --cleanenv --bind "$COURSEDIR" --bind "$WORKDIR" --pwd "$OUTDIR" "$ANALYSIS_SIF" "$@"; }

# Make UniProt IDs in the best-hit table match the FASTA header tokens (sp|ACC|NAME)
grep '>' "$UNIPROT" | awk -F'\t' -v OFS='\t' 'NR==FNR{t=$0; sub(/ .*/,"",t); sub(/^>/,"",t); split(t,a,"|"); full[a[2]]=t; next}
  {if($2 in full) $2=full[$2]; print}' - blastp_uniprot.besthits.tsv > blastp_uniprot.besthits.fixed.tsv
echo "best hits remapped: $(awk -F'\t' '$2 ~ /^sp\|/' blastp_uniprot.besthits.fixed.tsv | wc -l) of $(wc -l < blastp_uniprot.besthits.tsv)"

mk maker_functional_fasta "$UNIPROT" blastp_uniprot.besthits.fixed.tsv "$PROTEINS" > maker_proteins.filtered.Uniprot.fasta
mk maker_functional_gff   "$UNIPROT" blastp_uniprot.besthits.fixed.tsv "$GENES"    > filtered.genes.Uniprot.gff3 || echo "WARNING: maker_functional_gff failed" >&2

# ---- Guiding questions: hit proportions and length bias ----
{
  total=$(grep -c '>' "$PROTEINS")
  uni=$(wc -l < blastp_uniprot.besthits.tsv)
  tair=$(wc -l < blastp_TAIR10.besthits.tsv)
  unk=$(grep '>' maker_proteins.filtered.Uniprot.fasta | grep -c "Protein of unknown function" || true)
  unc=$(grep '>' maker_proteins.filtered.Uniprot.fasta | grep -ci "Similar to.*uncharacterized" || true)
  echo "filtered proteins:            $total"
  echo "with UniProt hit:             $uni ($(awk -v a=$uni -v b=$total 'BEGIN{printf "%.1f", 100*a/b}')%)"
  echo "with TAIR10 hit:              $tair ($(awk -v a=$tair -v b=$total 'BEGIN{printf "%.1f", 100*a/b}')%)"
  echo "'Protein of unknown function': $unk"
  echo "hits to 'uncharacterized':    $unc"
  echo "--- protein length with vs without UniProt hit (aa) ---"
  ana seqkit fx2tab -n -i -l "$PROTEINS" | \
    awk -F'\t' 'NR==FNR{h[$1]=1; next}
      {g=($1 in h)?"with_hit":"no_hit"; n[g]++; s[g]+=$2; if($2<100) sh[g]++}
      END{for(g in n) printf "%s: n=%d mean_len=%.0f short(<100aa)=%d (%.1f%%)\n", g, n[g], s[g]/n[g], sh[g], 100*sh[g]/n[g]}' \
    blastp_uniprot.besthits.tsv -
} > uniprot_hit_summary.txt
cat uniprot_hit_summary.txt

# ---- FLC (AT5G10140) ----
{
  echo "## Top TAIR10 hits to AT5G10140 (query, subject, %id, alnlen, evalue, bits)"
  awk -F'\t' '$2 ~ /^AT5G10140/ {print $1,$2,$3,$4,$11,$12}' OFS='\t' blastp_TAIR10.tsv | sort -k6,6gr | head -5
  echo "## Queries whose BEST TAIR10 hit is AT5G10140"
  ids=$(awk -F'\t' '$2 ~ /^AT5G10140/ {print $1}' blastp_TAIR10.besthits.tsv)
  echo "${ids:-none}"
  for id in $ids; do
    echo "## $id: mRNA line"
    grep -P "\tmRNA\t.*ID=${id};" "$GENES" || true
    echo "## $id: InterProScan (Pfam/GO)"
    grep -P "^${id}\t" "$FINALDIR/output.iprscan" || true
  done
} > FLC_candidates.txt
cat FLC_candidates.txt