# Script to plot gene structure given a GTF
# A. M. Chakrabarti
# Last updated: 28th June 2019

library(optparse)

# TODO: Nobby to convert this list into an optparse object for final script.

opt <- list(xlinks = "../tdp43_studentExercise/tardbp-esc-m-p2lox-gfp-tdp43-20151212_trimmed_single.bedgraph ../tdp43_studentExercise/tardbp-ngfp-esc-m-p2lox-gfp-tdp43-20151212_trimmed_single.bedgraph", #from https://imaps.genialis.com/iclip
           # track_names = "tdp43_1 tdp43_2", # or can be NULL 
	gtf = "../hg38_regions/gencode.v30.primary_assembly.annotation.gtf", #"~/Ule/ref/gencode.v27.annotation.gtf.gz"
	region = "chr3:35754106:35856276:+",
	gene = "ENSMUSG00000037400", #ID or name "Atp11b"
	smoothing = "gaussian", #gaussian or rollmean
	smoothing_window = 10, #both types of smoothing require a window
	normalisation = "libsize", #libsize or maxpeak or none
	output = "plot.pdf")

# ==========
# Part 0 - Get relevant regions in bedgraphs
# ==========

library(rtracklayer)

region.gr <- GRanges(seqnames = Rle(sapply(strsplit(opt$region, ":"), "[[", 1)),
                     ranges = IRanges(start = as.integer(sapply(strsplit(opt$region, ":"), "[[", 2)),
                                      end = as.integer(sapply(strsplit(opt$region, ":"), "[[", 3))),
                     strand = Rle(sapply(strsplit(opt$region, ":"), "[[", 4)))


# ==========
# Part 1 - top half: normalised and smoothed tracks
# ==========

library(ggplot2)
library(ggthemes)
library(cowplot)
library(smoother)
library(zoo)
library(dplyr)

# Read in xlinks
xlinks <- strsplit(opt$xlinks, " ")[[1]] %>% lapply(.,import, format="bedGRaph")
libSizes <- lapply(xlinks, function(x){sum(abs(x$score))})

# Subset for region
xlinks <- lapply(xlinks,subsetByOverlaps, region.gr)

# Names of bedgraphs : If name is supplied use that, if not then generate a name from the file name
if (!is.null(opt$track_names)) {
  track_names=strsplit(opt$track_names, " ")[[1]]
} else {
  track_names = lapply(strsplit(opt$xlinks, " ")[[1]], function(x) gsub(".bedgraph", "", basename(x)))
}
# Create dataframe for plotting
xlinks_df <- lapply(xlinks, as.data.frame)
for (i in 1:length(xlinks_df)){
  xlinks_df[[i]]$sample <- track_names[[i]]
  xlinks_df[[i]]$libSize <- libSizes[[i]]
}
xl_df <- dplyr::bind_rows(xlinks_df)
xl_df <- as.data.frame(xl_df)

# Do the normalisation
xl_df <- as.data.frame(switch(opt$normalisation, "libsize"=xl_df %>% dplyr::group_by(sample) %>% dplyr::mutate(norm=score/libSize),
       "maxpeak"=xl_df %>% dplyr::group_by(sample) %>% dplyr::mutate(norm=score/max(score)),
       "none"=xl_df %>% dplyr::group_by(sample) %>% dplyr::mutate(norm=score)))
# Do the smoothing
xl_df <- as.data.frame(switch(opt$smoothing, "gaussian"=xl_df %>% dplyr::group_by(sample) %>% dplyr::mutate(smoothed= smth.gaussian(norm, window=opt$smoothing_window)),
                              "rollmean"=xl_df %>% dplyr::group_by(sample) %>% dplyr::mutate(smoothed= rollmean(norm, opt$smoothing_window, fill=0))))
# Plot
ggplot(xl_df,aes(x=start,y=smoothed, group=sample, color=sample)) +
 geom_line() +
  theme_few() +
  ylab("smoothed cDNAs") +
  xlab("chr coordinate")


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
