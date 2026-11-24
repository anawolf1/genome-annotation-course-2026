# Organization and annotation of eukaryote genomes — *Arabidopsis thaliana* Rab-R1

**Course:** Organization and annotation of Eukaryote genomes (UE-SBL.30004), University of Fribourg / Bern, autumn 2026
**Student:** Ana Wolf Ruiz
**Accession:** Rab-R1 (Madeira genetic group in Lian et al. 2024)
**Cluster:** IBU HPC (SLURM, partition `pibu_el8`), project directory `/data/users/awolfruiz/genome_annotation_course`

## Input data

| Input | Source | Path in project |
|---|---|---|
| Genome assembly | Flye assembly from the Assembly course (147.2 Mb, 170 contigs, BUSCO 100%, QV 61.6) | `data/assemblies/assembly.fasta` |
| Transcriptome | Trinity assembly (accession Sha, shared course RNA-seq) | `data/transcripts/trinity.Trinity.fasta` |
| CDS for EDTA | TAIR10 representative CDS (course directory) | `/data/courses/genome_annotation_resources/2026/data/` |
| Proteins for MAKER | TAIR10 representative peptides + UniProt Viridiplantae (reviewed) | course directory |

The Flye assembly was chosen over hifiasm and LJA because its size matched the k-mer genome size estimate (GenomeScope 144.6 Mb) and it had the best BUSCO completeness and Merqury QV.

## Project structure

    genome_annotation_course/
    ├── data/{assemblies,transcripts}
    ├── scripts/            # numbered in execution order
    ├── logs/               # SLURM .out/.err
    ├── output/
    │   ├── EDTA/                  # TE library and genome-wide TE annotation
    │   ├── LTR_full_length/       # TEsorter on intact LTR-RTs + identity plots
    │   ├── TEsorter_TElib/        # TEsorter on Copia/Gypsy sequences of the TE library
    │   └── MAKER/                 # gene annotation
    ├── maker_*.ctl         # MAKER control files (edited maker_opts.ctl)
    └── README.md

## Scripts (execution order)

| Script | Tool | Purpose |
|---|---|---|
| `01-run_EDTA.sh` | EDTA 2.3.0 (`--sensitive 1 --anno 1`, TAIR10 CDS) | De novo TE library + genome-wide TE annotation |
| `02-run_TEsorter_full_length_LTRs.sh` | TEsorter 1.3.0 (REXdb-plant) | Clade classification of intact LTR retrotransposons |
| `03-run_LTR_identity_plot.sh` → `03-full_length_LTRs_identity_copy.R` | R | LTR identity (insertion age) per clade |
| `04-run_maker.sh` | MAKER 3.01.03 (MPI, 20 CPUs) | Gene annotation (ESTs, proteins, Augustus `arabidopsis`, EDTA TE library for masking) |
| `05-run_TEsorter_TElib.sh` | SeqKit + TEsorter | Clade classification of all Copia/Gypsy sequences in the TE library |
| `06-run_parseRM.sh` | parseRM.pl | RepeatMasker divergence binned in 1% classes (0% self-matches removed) |
| `07-run_TE_plots.sh` → `04-annotation_circlize.R`, `06-plot_div.R` | samtools, R (circlize, ggplot2) | TE density along contigs; TE divergence landscape |
| `08-clade_abundance.sh` | awk | Copies and bp per Copia/Gypsy clade (TEsorter clades of library families joined to `TEanno.gff3`) |
| `09-merge_maker.sh` | MAKER `gff3_merge`, `fasta_merge` | Merge the 170 per-contig MAKER results into genome-wide GFF3 and FASTA |

Scripts 01, 02, 04 and the R scripts `03-full_length_LTRs_identity_copy.R`, `04-annotation_circlize.R`, `06-plot_div.R` are course scripts adapted to this project (paths, resources, input names); `05-parseRM.pl` is called from the course directory. Scripts 03-run, 05, 06, 07, 08 and 09 were written for this project, following the commands in the course handout where available.

## Results

### 1. TE content (EDTA, `assembly.fasta.mod.EDTA.TEanno.sum`)

| Superfamily | bp masked | % of genome |
|---|---|---|
| LTR / unknown | 9,930,749 | 6.74 |
| LTR / Gypsy | 6,278,442 | 4.26 |
| Helitron | 4,236,031 | 2.88 |
| LTR / Copia | 2,033,594 | 1.38 |
| TIR / Mutator | 1,296,761 | 0.88 |
| TIR / CACTA | 1,201,281 | 0.82 |
| LINE / L1 | 763,902 | 0.52 |
| TIR / PIF-Harbinger | 435,068 | 0.30 |
| TIR / hAT | 111,929 | 0.08 |
| TIR / Tc1-Mariner | 102,998 | 0.07 |
| SINE / tRNA | 8,711 | 0.01 |
| Repeat fragments | 814,555 | 0.55 |
| **Total interspersed** | **27,214,021** | **18.48** |

**Interpretation:** LTR retrotransposons are the largest TE fraction (~12.4%), with Gypsy about 3× more abundant than Copia, followed by Helitrons. This matches Lian et al. (2024), where LTR-RTs and Helitrons were the largest TE fractions. The total (27 Mb) is higher than the 15–18 Mb reported by Lian et al., mainly because of the large "LTR/unknown" class (96,754 mostly short fragments not assigned to Copia or Gypsy). Lian et al. used a pan-TE library from 69 genomes, so absolute values are not directly comparable.

### 2. Intact LTR retrotransposons: clades and insertion age

Plots: `output/LTR_full_length/plots/01_LTR_Copia_Gypsy_cladelevel.png`, `02_LTR_Copia_Gypsy_global.png`

Identity between the two LTRs of an intact element estimates insertion age (identical at insertion). Ages: T = K / 2r with K = 1 − identity and r = 8.22 × 10⁻⁹ substitutions/site/year (Kagale et al. 2014).

| LTR identity | Approx. age |
|---|---|
| 0.995–1.00 | < 0.3 My |
| 0.99 | ~0.6 My |
| 0.975 | ~1.5 My |
| 0.95 | ~3.0 My |

Number of intact elements per clade (TEsorter, `output/LTR_full_length/intact_LTR.fa.rexdb-plant.cls.tsv`):

| Copia clade | n | Identity pattern | Gypsy clade | n | Identity pattern |
|---|---|---|---|---|---|
| Ale | 76 | peak at 1.00 (youngest) | Athila | 39 | broad, 0.93–1.00 |
| Tork | 28 | peak ~0.99 | Reina | 20 | 0.96–1.00 |
| Bianca | 23 | broad, 0.90–0.99 (older) | Retand | 20 | peak ~0.975 |
| Ivana | 22 | peak at 1.00 | Tekay | 17 | 0.96–0.99 |
| SIRE | 16 | 0.96–1.00 | CRM | 14 | 0.93–1.00 |
| Alesia / Ikeros / TAR / unknown | 3 / 3 / 1 / 4 | — | unknown / mixture | 4 / 2 | — |
| **Total Copia** | **176** | | **Total Gypsy** | **116** | |

**Interpretation:** Although Gypsy occupies ~3× more sequence than Copia, there are more intact Copia elements, and they are younger: Ale and Ivana have most copies at 99–100% identity (< ~0.6 My), indicating recent activity. Ale alone accounts for 43% of intact Copia, so the young Copia signal is driven mainly by one clade. Intact Gypsy elements peak at slightly lower identity (~0.975, ~1.5 My) and have almost no copies at 1.00. The Copia and Gypsy clades found match those reported as dominant in Brassicales (Copia: Ale, Bianca, Ivana, SIRE, TAR, Tork; Gypsy: Athila, CRM, Reina, Retand, Tekay). Caveat: intact elements are biased toward young insertions, since older copies are mutated, truncated or reduced to solo-LTRs and no longer detected as intact.

### 3. Genomic distribution of Gypsy and Copia

Plot: `output/EDTA/plots/TE_density_circlize.pdf` (20 largest contigs)

**Interpretation:** Gypsy elements are concentrated in a few dense blocks (e.g. contig_153, contig_64, contig_78, contig_224), often at contig ends, whereas Copia elements are dispersed along most contigs. This is consistent with Gypsy (particularly Athila) accumulating in (peri)centromeric heterochromatin, where the assembly also tends to break, and Copia inserting in chromosome arms. Because the assembly is contig-level, centromere positions are inferred, not observed.

### 4. TE divergence landscape

Plot: `output/EDTA/plots/TE_landscape_divergence.pdf`

| Divergence | Approx. age |
|---|---|
| 1% | ~0.6 My |
| 4% | ~2.4 My |
| 8% | ~4.9 My |
| 13% | ~7.9 My |
| 30% | ~18 My |

**Interpretation:** Gypsy dominates the landscape with a broad peak at ~4–13% divergence (maximum ~8%, ≈ 5 My), indicating an older amplification wave; most Gypsy sequence is therefore older than its intact copies. Copia is enriched at the lowest divergence (0–4%), consistent with its young intact copies. DNA transposons (mainly Mutator and CACTA) show low, continuous accumulation across divergence classes. Caveats: ages use uncorrected divergence (no Jukes–Cantor / CpG correction) and a constant neutral rate, so older ages are underestimated; Helitrons, LINEs and LTR/unknown are not shown by the course plotting script.

### 5. Genome-wide abundance of Copia and Gypsy clades

Script: `scripts/08-clade_abundance.sh` (TEsorter clades of the EDTA library joined to all annotated copies in `TEanno.gff3`); output `output/TEsorter_TElib/clade_abundance.tsv`.

| Superfamily | Clade | Annotated copies | bp | Intact copies | Copies per intact |
|---|---|---|---|---|---|
| Gypsy | Athila | 2,895 | 4,230,190 | 39 | 74 |
| Gypsy | Retand | 886 | 1,061,802 | 20 | 44 |
| Gypsy | Tekay | 349 | 643,382 | 17 | 21 |
| Gypsy | CRM | 162 | 239,130 | 14 | 12 |
| Gypsy | Reina | 181 | 179,362 | 20 | 9 |
| Copia | Ale | 602 | 835,323 | 76 | 8 |
| Copia | Bianca | 209 | 319,023 | 23 | 9 |
| Copia | Tork | 214 | 308,584 | 28 | 8 |
| Copia | Ivana | 165 | 293,101 | 22 | 8 |
| Copia | SIRE | 116 | 219,652 | 16 | 7 |
| Copia | Ikeros / TAR / Angela / Alesia | 32 / 21 / 49 / 5 | 39,523 / 19,669 / 18,680 / 16,153 | 3 / 1 / 0 / 3 | — |

Clade totals (Copia 2.07 Mb, Gypsy 6.35 Mb) agree with the EDTA superfamily summary.

**Interpretation:** Athila is the most abundant LTR clade (4.2 Mb, two-thirds of Gypsy sequence), but it is present mostly as fragments: ~74 annotated copies per intact element, compared with ~8 for Ale and the other Copia clades. This supports the conclusions from LTR identity and the divergence landscape: Gypsy (especially Athila) accumulated in an older amplification wave and has since decayed, while Copia (especially Ale) is dominated by recent, still-intact insertions.

## Gene annotation (in progress)

MAKER was run with the Flye assembly, the Trinity transcriptome (`est2genome=1`), TAIR10 + UniProt proteins (`protein2genome=1`), Augustus (`arabidopsis` model), the EDTA TE library (`rmlib`) and PTREP20 repeat proteins. It finished all 170 contigs in 5 h 05 min; per-contig results were merged with `scripts/09-merge_maker.sh`.

| Metric (before filtering) | Rab-R1 | Lian et al. 2024 |
|---|---|---|
| Gene models | 27,411 | 27,246–28,989 protein-coding genes per accession (Col-0: 27,445) |
| Transcripts / proteins | 30,165 | — |
| Transcripts with AED < 0.5 | 29,492 (97.8%) | — |

**Interpretation:** the raw gene count lies within the range reported for *A. thaliana* accessions, suggesting that repeat masking with the EDTA library prevented TE-derived genes from inflating the annotation. Nearly all models are well supported by evidence (AED < 0.5); this is expected because `keep_preds=0` discards ab initio predictions without transcript or protein support. Filtering, BUSCO and functional annotation follow.

## Problems encountered and how they were handled

- **Course account not available:** `sbatch` rejected `--account=class-sbl-30004-ws-2026` ("Invalid account or account/partition combination"); `sacctmgr` showed my user is not associated with it. The `--account` and `--reservation` lines were commented out (`##SBATCH`) and jobs were run in the normal `pibu_el8` queue under my personal account.
- **EDTA out of memory:** the EDTA job (96 GB) was OOM-killed after 26 h, during the final post-annotation statistics step (`TE.fa.stat*`). All outputs used here (`TElib.fa`, `TEanno.gff3`, `TEanno.sum`, `intact.gff3`, `RM.out`) had already been written, so EDTA was not rerun.
- **`apptainer` not available on the login node:** all container steps, including R plotting, were run as small SLURM jobs.
- **Slurm controller outage (27.10):** `sbatch` failed with "Unable to contact slurm controller"; the MAKER merge was submitted once the scheduler was back.
- **parseRM input not found:** with `apptainer exec -C` the container did not start in the working directory; fixed with `--pwd "$PWD"` and an absolute input path.

## To do

- MAKER: filter gene models, AED distribution plot, BUSCO, functional annotation (InterProScan)
- GENESPACE / OrthoFinder comparison with other accessions; FLC characterization

## References

- Lian Q. et al. (2024) A pan-genome of 69 *Arabidopsis thaliana* accessions reveals a conserved genome structure throughout the global species range. *Nature Genetics* 56:982–991.
- Ou S. et al. (2019) Benchmarking transposable element annotation methods for creation of a streamlined, comprehensive pipeline. *Genome Biology* 20:275 (EDTA).
- Zhang R.-G. et al. (2022) TEsorter: an accurate and fast method to classify LTR-retrotransposons in plant genomes. *Horticulture Research* 9:uhac017.
- Neumann P. et al. (2019) Systematic survey of plant LTR-retrotransposons elucidates phylogenetic relationships of their polyprotein domains and provides a reference for element classification. *Mobile DNA* 10:1.
- Wicker T. et al. (2007) A unified classification system for eukaryotic transposable elements. *Nature Reviews Genetics* 8:973–982.
- Kagale S. et al. (2014) Polyploid evolution of the Brassicaceae during the Cenozoic era. *Plant Cell* 26:2777–2791.
