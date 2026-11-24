#!/bin/bash
# Genome-wide abundance (copies, bp) of Copia/Gypsy clades:
# TEsorter clade of each library sequence (cls.tsv) joined to annotated copies in EDTA TEanno.gff3
WORKDIR=/data/users/$USER/genome_annotation_course
CLS="$WORKDIR/output/TEsorter_TElib/Copia_sequences.fa.rexdb-plant.cls.tsv $WORKDIR/output/TEsorter_TElib/Gypsy_sequences.fa.rexdb-plant.cls.tsv"
GFF=$WORKDIR/output/EDTA/assembly.fasta.mod.EDTA.TEanno.gff3
OUT=$WORKDIR/output/TEsorter_TElib/clade_abundance.tsv

cat $CLS | grep -v "^#" | \
awk -F'\t' 'NR==FNR { split($1,a,"#"); clade[a[1]]=$3"\t"$4; next }
  $0 !~ /^#/ && $3 !~ /repeat_region|long_terminal_repeat|target_site_duplication/ {
    if (match($9,/Name=[^;]+/)) { n=substr($9,RSTART+5,RLENGTH-5)
      if (n in clade) { c[clade[n]]++; bp[clade[n]]+=$5-$4+1 } } }
  END { print "Superfamily\tClade\tCopies\tbp"; for (k in c) print k"\t"c[k]"\t"bp[k] }' \
  - $GFF | sort -t$'\t' -k1,1 -k4,4nr > $OUT
cat $OUT
