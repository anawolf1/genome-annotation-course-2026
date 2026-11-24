library(GENESPACE)
args <- commandArgs(trailingOnly = TRUE)
wd <- args[1]
focal_genome <- args[2]
# optional 3rd arg: comma-separated scaffolds of the focal genome to flip for display
flip_chrs <- if (length(args) >= 3 && nzchar(args[3])) strsplit(args[3], ",")[[1]] else NULL

load(file.path(wd, "results", "gsParams.rda"))   # loads gsParam
outdir <- file.path(wd, "processed_riparian")
dir.create(outdir, showWarnings = FALSE)

# 1) focal vs TAIR10, scaffolds ordered by synteny
plot_riparian(gsParam, genomeIDs = c("TAIR10", focal_genome), refGenome = "TAIR10",
              reorderBySynteny = TRUE,
              pdfFile = file.path(outdir, "my_accession_vs_TAIR10.pdf"))

flips <- if (!is.null(flip_chrs)) data.frame(genome = focal_genome, chr = flip_chrs) else NULL

# 2) oriented version (flips optional)
plot_riparian(gsParam, genomeIDs = c("TAIR10", focal_genome), refGenome = "TAIR10",
              reorderBySynteny = TRUE, syntenyWeight = 1, invertTheseChrs = flips,
              pdfFile = file.path(outdir, "my_accession_vs_TAIR10_oriented.pdf"))

# 3) all five genomes
plot_riparian(gsParam, refGenome = "TAIR10", reorderBySynteny = TRUE, syntenyWeight = 1,
              invertTheseChrs = flips, scalePlotWidth = 1.5,
              pdfFile = file.path(outdir, "all_genomes_TAIR10.pdf"))
