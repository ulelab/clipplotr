# Script to plot gene structure given a GTF
# A. M. Chakrabarti
# Last updated: 28th June 2019

library(optparse)

# TODO: Nobby to convert this list into an optparse object for final script.

opt <- list(xlinks = "iclip1.bedgraph.gz,iclip2.bedgraph.gz,iclip3.bedgraph.gz",
	gtf = "~/Ule/ref/gencode.v27.annotation.gtf.gz",
	region = "1:123:123:+",
	gene = "ENSG",
	smoothing = "gaussian",
	normalisation = "libSize",
	output = "plot.pdf")

# ==========
# Part 0 - Get relevant regions in bedgraphs
# ==========

library(rtracklayer)

region.gr <- GRanges(seqnames = Rle(sapply(strsplit(opt$region), ":"), "[[", 1),
                     ranges = IRanges(start = Rle(sapply(strsplit(opt$region), ":"), "[[", 2),
                                      end = Rle(sapply(strsplit(opt$region), ":"), "[[", 3)),
                     strand = Rle(Rle(sapply(strsplit(opt$region), ":"), "[[", 4)))



# ==========
# Part 1 - top half: normalised and smoothed tracks
# ==========

library(ggplot2)
library(ggthemes)
library(cowplot)
library(smoother)
library(zoo)

# ==========
# Part 2 - bottom half: gene structures
# ==========

library(ggbio)
library(GenomicFeatures)

gtf <- rtracklayer::import.gff2(opt$gtf)
test <- gtf[gtf$gene_name == "CAMK2A"]
regions.gr <- test[test$type == "gene"]


TxDb <- makeTxDbFromGFF(opt$gtf)
# p.annot <- autoplot(TxDb, which = regions.gr)

annot.gr <- biovizBase::crunch(TxDb, which = regions.gr)
annot.grl <- split(annot.gr, annot.gr$tx_name)

# Need to trim exons so they don't overlap utrs to ensure they aren't overplotted
annot.grl <- GRangesList(lapply(annot.grl, function(x) {
  
  utr <- x[x$type == "utr"]
  exon <- x[x$type == "exon"]
  rest <- x[!x$type %in% c("utr", "exon")]

  exon <- setdiff(exon, utr, ignore.strand = TRUE)  
  mcols(exon)$type <- "exon"
  
  return(sort(c(rest, exon, utr)))
  
}))


p.annot <- ggplot(data = annot.grl) +
  geom_alignment(cds.rect.h = 0.25)


# ==========
# Part 3 - combine plots
# ==========

save_plot(plot_grid(p.top, p.annot, align = "hv", axis = "tlbr", nrow = 2, rel_heights, c(1, 2)), filename = opt$output)