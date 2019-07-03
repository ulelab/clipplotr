# Script to plot gene structure given a GTF
# A. M. Chakrabarti
# Last updated: 28th June 2019

library(optparse)

option_list <- list(make_option(c("-x", "--xlinks"), action = "store", type = "character", help = "Input iCLIP bedgraphs (space separated)"),
                    make_option(c("-g", "--gtf"), action = "store", type = "character", help = "Reference gtf (Gencode)"),
                    make_option(c("-r", "--region"), action = "store", type = "character", help = "Region of interest as chr3:35754106:35856276:+ or gene as ENSMUSG00000037400 or Atp11b"),
                    make_option(c("-n", "--normalisation"), action = "store", type = "character", help = "Normalisation options: none, maxpeak, libsize", default = "libsize"),
                    make_option(c("-s", "--smoothing"), action = "store", type = "character", help = "Normalisation options: none, gaussian, rollmean", default = "gaussian"),
                    make_option(c("-w", "--smoothing_window"), action = "store", type = "integer", help = "Smoothing window", default = 10),
                    make_option(c("-o", "--output", action = "store", type = "character", help = "Output plot filename")))
opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# opt <- list(xlinks = "../tdp43_studentExercise/tardbp-esc-m-p2lox-gfp-tdp43-20151212_trimmed_single.bedgraph ../tdp43_studentExercise/tardbp-ngfp-esc-m-p2lox-gfp-tdp43-20151212_trimmed_single.bedgraph", #from https://imaps.genialis.com/iclip
#            # track_names = "tdp43_1 tdp43_2", # or can be NULL 
# 	gtf = "../hg38_regions/gencode.v30.primary_assembly.annotation.gtf", #"~/Ule/ref/gencode.v27.annotation.gtf.gz"
# 	region = "chr3:35754106:35856276:+",
# 	gene = "ENSMUSG00000037400", #ID or name "Atp11b"
# 	smoothing = "gaussian", #gaussian or rollmean
# 	smoothing_window = 10, #both types of smoothing require a window
# 	normalisation = "libsize", #libsize or maxpeak or none
# 	output = "plot.pdf")

setwd("~/Ule/charlotte/CLIPplotR")

opt <- list(xlinks = "tardbp-esc-m-p2lox-gfp-tdp43-20151212_trimmed_single.bedgraph tardbp-ngfp-esc-m-p2lox-gfp-tdp43-20151212_trimmed_single.bedgraph", #from https://imaps.genialis.com/iclip
            # track_names = "tdp43_1 tdp43_2", # or can be NULL 
            gtf = "gencode.vM22.annotation.gtf.gz", #"~/Ule/ref/gencode.v27.annotation.gtf.gz"
            region = "chr3:35754106:35856276:+",
            # gene = "ENSMUSG00000037400", #ID or name "Atp11b"
            # gene = "Atp11b", #ID or name "Atp11b"
            smoothing = "none", #gaussian or rollmean or none
            smoothing_window = 10, #both types of smoothing require a window
            normalisation = "none", #libsize or maxpeak or none
            output = "plot.pdf")

# Actually just have it so that region can be a gene id or name

# ==========
# Part 0 - Get relevant regions in bedgraphs
# ==========

library(rtracklayer)
library(GenomicFeatures)

# Create rosetta and TxDb
message("Making annotation from GTF")
gtf <- rtracklayer::import.gff2(opt$gtf)
genes.gr <- gtf[gtf$type == "gene"]
TxDb <- makeTxDbFromGFF(opt$gtf)

# Define region and plot title
if(is.null(opt$region)) {
  
  stop("Need to supply a region or a gene id or name")
  
} else if(grepl("^ENS", opt$region)) {
  
  region.gr <- genes.gr[grepl(opt$gene, genes.gr$gene_id)]
  
} else if(grepl("^[A-Z]", opt$gene)) {
  
  region.gr <- genes.gr[grepl(opt$gene, genes.gr$gene_name)]
  
} else {
  
  region.gr <- GRanges(seqnames = Rle(sapply(strsplit(opt$region, ":"), "[[", 1)),
                       ranges = IRanges(start = as.integer(sapply(strsplit(opt$region, ":"), "[[", 2)),
                                        end = as.integer(sapply(strsplit(opt$region, ":"), "[[", 3))),
                       strand = Rle(sapply(strsplit(opt$region, ":"), "[[", 4)))
  
}


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
xlinks <- strsplit(opt$xlinks, " ")[[1]] %>% lapply(.,import, format="bedGraph")
libSizes <- lapply(xlinks, function(x){sum(abs(x$score))})

# Subset for region
xlinks <- lapply(xlinks, subsetByOverlaps, region.gr)

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
# 
# I get this warning here:
# Warning messages:
#   1: In bind_rows_(x, .id) : Unequal factor levels: coercing to character
# 2: In bind_rows_(x, .id) :
#   binding character and factor vector, coercing into character vector
# 3: In bind_rows_(x, .id) :
#   binding character and factor vector, coercing into character vector

xl_df <- as.data.frame(xl_df) # Isn't it already a data.frame?

# Do the normalisation
xl_df <- as.data.frame(switch(opt$normalisation, "libsize"=xl_df %>% dplyr::group_by(sample) %>% dplyr::mutate(norm=score/libSize),
       "maxpeak"=xl_df %>% dplyr::group_by(sample) %>% dplyr::mutate(norm=score/max(score)),
       "none"=xl_df %>% dplyr::group_by(sample) %>% dplyr::mutate(norm=score)))

# TODO: maybe convert to XPM for easier reading

# Do the smoothing


xl_df <- as.data.frame(switch(opt$smoothing, 
                              "gaussian"=xl_df %>% dplyr::group_by(sample) %>% dplyr::mutate(smoothed= smth.gaussian(norm, window=opt$smoothing_window)),
                              "rollmean"=xl_df %>% dplyr::group_by(sample) %>% dplyr::mutate(smoothed= rollmean(norm, opt$smoothing_window, fill=0)),
                              "none"=xl_df %>% dplyr::group_by(sample) %>% dplyr::mutate(smoothed= norm)))

# TODO: I think for the smoothing to work properly, the vector needs to be without breaks, 1 value per coordinate and include positions with counts of 0

# TODO: Also maybe smooth + and - strands separately... Or only keep xlinks on the same strand as the gene/region. Maybe latter. Actually, should filter the bedgraph properly using subsetByOverlaps

# Plot
p.iclip <- ggplot(xl_df,aes(x=start,y=smoothed, group=sample, color=sample)) +
  geom_line() +
  labs(title = opt$region,
       x = "Coordinate",
       y = "Normalised crosslinks")+
  theme_cowplot() + theme(legend.position = "top")

# TODO: remove x axis and x axis label once we are happy the two plots line up (just keep on bottom annotation plot)

# ==========
# Part 2 - bottom half: gene structures
# ==========

library(ggbio)

annot.gr <- biovizBase::crunch(TxDb, which = region.gr)
annot.grl <- split(annot.gr, annot.gr$tx_name)

# Need to trim exons so they don't overlap utrs to ensure they aren't overplotted
annot.grl <- GRangesList(lapply(annot.grl, function(x) {
  
  utr <- x[x$type == "utr"]
  exon <- x[x$type == "exon"]
  rest <- x[!x$type %in% c("utr", "exon")]

  exon <- GenomicRanges::setdiff(exon, utr, ignore.strand = TRUE)  
  mcols(exon)$type <- "exon"
  
  return(sort(c(rest, exon, utr)))
  
}))


p.annot <- ggplot(data = annot.grl) +
  geom_alignment(cds.rect.h = 0.25) +
  labs(x = "Coordinate")
  theme_cowplot()

p.annot <- p.annot@ggplot

# ==========
# Part 3 - combine plots
# ==========

# Interactive for now to test
plot_grid(p.iclip, p.annot, align = "hv", axis = "tlbr", nrow = 2, rel_heights = c(1, 2))

# save_plot(plot_grid(p.top, p.annot, align = "hv", axis = "tlbr", nrow = 2, rel_heights, c(1, 2)), filename = opt$output)
