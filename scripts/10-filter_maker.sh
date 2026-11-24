#!/bin/bash
#SBATCH --job-name=maker_filter
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Week 7 (handout section 6): rename gene models, Pfam/GO with InterProScan, AED distribution,
# quality filter, gene features only, filtered + longest-isoform FASTAs, AGAT statistics
set -euo pipefail
WORKDIR=/data/users/$USER/genome_annotation_course
COURSEDIR=/data/courses/genome_annotation_resources/2026
OUTDIR=$WORKDIR/output/MAKER
FINALDIR=$OUTDIR/final
MAKER_SIF=$COURSEDIR/containers/MAKER_3.01.03.sif
IPR_SIF=$COURSEDIR/containers/InterProScan_5.70.sif
ANALYSIS_SIF=$COURSEDIR/containers/Annotation_analysis_R4.4.sif
AGAT_SIF=$COURSEDIR/containers/AGAT_v1.4.1.sif
protein=assembly.all.maker.proteins.fasta
transcript=assembly.all.maker.transcripts.fasta
gff=assembly.all.maker.noseq.gff
prefix=RAB

mkdir -p $FINALDIR/tmp $OUTDIR/tmp
cd $FINALDIR
mk()  { apptainer exec --cleanenv --env TMPDIR=/tmp --bind $OUTDIR/tmp:/tmp --bind $COURSEDIR --bind $WORKDIR --pwd $FINALDIR $MAKER_SIF "$@"; }
ana() { apptainer exec --cleanenv --bind $COURSEDIR --bind $WORKDIR --pwd $FINALDIR $ANALYSIS_SIF "$@"; }

# 6.1 fresh copies of the merged files, then consistent IDs (RAB0000001-RA ...)
cp $OUTDIR/$gff ${gff}.renamed.gff
cp $OUTDIR/$protein ${protein}.renamed.fasta
cp $OUTDIR/$transcript ${transcript}.renamed.fasta
mk maker_map_ids --prefix $prefix --justify 7 ${gff}.renamed.gff > id.map
mk map_gff_ids id.map ${gff}.renamed.gff
mk map_fasta_ids id.map ${protein}.renamed.fasta
mk map_fasta_ids id.map ${transcript}.renamed.fasta

# 6.2 Pfam domains + GO terms
apptainer exec --cleanenv --pwd $FINALDIR \
  --bind $COURSEDIR/data/interproscan-5.70-102.0/data:/opt/interproscan/data \
  --bind $COURSEDIR --bind $WORKDIR --bind $FINALDIR/tmp:/temp \
  $IPR_SIF /opt/interproscan/interproscan.sh \
  -appl pfam --disable-precalc -f TSV --goterms --iprlookup --seqtype p \
  -i ${protein}.renamed.fasta -o output.iprscan

# 6.3 add InterProScan results to the GFF
mk ipr_update_gff ${gff}.renamed.gff output.iprscan > ${gff}.renamed.iprscan.gff

# 6.4 cumulative AED distribution (bins of 0.025)
mk perl /usr/local/bin/AED_cdf_generator.pl -b 0.025 ${gff}.renamed.gff > assembly.all.maker.renamed.gff.AED.txt

# 6.5 keep transcripts with AED < 1 and/or a Pfam domain
mk perl /usr/local/bin/quality_filter.pl -s ${gff}.renamed.iprscan.gff > ${gff}_iprscan_quality_filtered.gff

# 6.6 keep only gene features and their subfeatures
grep -P "\tgene\t|\tCDS\t|\texon\t|\tfive_prime_UTR\t|\tthree_prime_UTR\t|\tmRNA\t" \
  ${gff}_iprscan_quality_filtered.gff > filtered.genes.renamed.gff3
sed -i '1i##gff-version 3' filtered.genes.renamed.gff3

# 6.7 transcript and protein FASTAs of the retained models
grep -P "\tmRNA\t" filtered.genes.renamed.gff3 | awk '{print $9}' | cut -d ';' -f1 | sed 's/ID=//g' > list.txt
ana seqkit grep -f list.txt ${transcript}.renamed.fasta -o ${transcript}.renamed.filtered.fasta
ana seqkit grep -f list.txt ${protein}.renamed.fasta -o ${protein}.renamed.filtered.fasta

# DIY: longest isoform per gene (gene ID = everything before -RA/-RB...), input for BUSCO
for f in $protein $transcript; do
  ana seqkit fx2tab -n -i -l ${f}.renamed.filtered.fasta | \
    awk -F'\t' '{g=$1; sub(/[-.]R[A-Z]+$/,"",g); if(!(g in L) || $2+0>L[g]){L[g]=$2+0; B[g]=$1}} END{for(g in B) print B[g]}' > ${f}.longest.ids
  ana seqkit grep -f ${f}.longest.ids ${f}.renamed.filtered.fasta -o ${f}.renamed.longest.fasta
done

# AGAT statistics of the final annotation
apptainer exec --cleanenv --bind $COURSEDIR --bind $WORKDIR --pwd $FINALDIR $AGAT_SIF \
  agat_sp_statistics.pl --gff filtered.genes.renamed.gff3 --output annotation.stat

# Summary numbers for the README / guiding questions
{
  echo "genes before filter: $(awk -F'\t' '$3=="gene"' ${gff}.renamed.gff | wc -l)"
  echo "genes after filter:  $(awk -F'\t' '$3=="gene"' filtered.genes.renamed.gff3 | wc -l)"
  echo "mRNAs after filter:  $(awk -F'\t' '$3=="mRNA"' filtered.genes.renamed.gff3 | wc -l)"
  echo "longest proteins:    $(grep -c '>' ${protein}.renamed.longest.fasta)"
  # share of transcripts with a Pfam domain, low vs high AED
  awk -F'\t' '$3=="mRNA"{a=1; if(match($9,/_AED=[0-9.]+/)) a=substr($9,RSTART+5,RLENGTH-5)+0;
      c=(a<0.5)?"AED<0.5 ":"AED>=0.5"; n[c]++; if($9 ~ /Pfam:/) k[c]++}
      END{for(c in n) printf "%s: %d transcripts, %d with Pfam (%.1f%%)\n", c, n[c], k[c], 100*k[c]/n[c]}' \
      ${gff}.renamed.iprscan.gff
} > summary.txt
cat summary.txt
