#!/usr/bin/env Rscript

# Script to plot multiple CLIP tracks with gene structures
# A. M. Chakrabarti
# C. Capitanchik
# Last updated: 10th July 2019

# ==========
# Preamble
# ==========

CheckAndLoad <- function(package) {
  
  if(!suppressPackageStartupMessages(require(package, character.only = TRUE, quietly = TRUE))) {
    
    message("Installing ", package)
    install.packages(package, character.only = TRUE)
    suppressPackageStartupMessages(library(package, character.only = TRUE, quietly = TRUE))
    
  }
  
}

CheckAndLoad("optparse")

option_list <- list(make_option(c("-x", "--xlinks"), action = "store", type = "character", help = "Input iCLIP bedgraphs (space separated)"),
                    make_option(c("-l", "--label"), action = "store", type = "character", help = "iCLIP bedgraph labels (space separated)"),
                    make_option(c("-g", "--gtf"), action = "store", type = "character", help = "Reference gtf (Gencode)"),
                    make_option(c("-r", "--region"), action = "store", type = "character", help = "Region of interest as chr3:35754106:35856276:+ or gene as ENSMUSG00000037400 or Atp11b"),
                    make_option(c("-n", "--normalisation"), action = "store", type = "character", help = "Normalisation options: none, maxpeak, libsize [default %default]", default = "libsize"),
                    make_option(c("-s", "--smoothing"), action = "store", type = "character", help = "Smoothing options: none, rollmean, spline, gaussian [default %default]", default = "rollmean"),
                    make_option(c("-w", "--smoothing_window"), action = "store", type = "integer", help = "Smoothing window [default %default]", default = 100),
                    make_option(c("-o", "--output", action = "store", type = "character", help = "Output plot filename")))
opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Load CRAN packages
packages <- c("BiocManager", "ggplot2", "ggthemes", "cowplot", "smoother", "zoo", "dplyr")
for(package in packages) CheckAndLoad(package)

# Loas Bioconductor packages
biocpackages <- c("rtracklayer", "GenomicFeatures", "ggbio")
for(package in biocpackages) {
  
  if(!suppressPackageStartupMessages(require(package, character.only = TRUE, quietly = TRUE))) {
    
    message("Installing ", package)
    BiocManager::install(package)
    suppressPackageStartupMessages(library(package, character.only = TRUE, quietly = TRUE))
    
  }  
  
}

# print(opt)

# suppressPackageStartupMessages(library(rtracklayer))
# suppressPackageStartupMessages(library(GenomicFeatures))
# suppressPackageStartupMessages(library(ggplot2))
# suppressPackageStartupMessages(library(ggthemes))
# suppressPackageStartupMessages(library(cowplot))
# suppressPackageStartupMessages(library(smoother))
# suppressPackageStartupMessages(library(zoo))
# suppressPackageStartupMessages(library(dplyr))
# suppressPackageStartupMessages(library(ggbio))

# opt <- list(xlinks = "../tdp43_studentExercise/tardbp-esc-m-p2lox-gfp-tdp43-20151212_trimmed_single.bedgraph ../tdp43_studentExercise/tardbp-ngfp-esc-m-p2lox-gfp-tdp43-20151212_trimmed_single.bedgraph", #from https://imaps.genialis.com/iclip
#            # track_names = "tdp43_1 tdp43_2", # or can be NULL 
# 	gtf = "../hg38_regions/gencode.v30.primary_assembly.annotation.gtf", #"~/Ule/ref/gencode.v27.annotation.gtf.gz"
# 	region = "chr3:35754106:35856276:+",
# 	gene = "ENSMUSG00000037400", #ID or name "Atp11b"
# 	smoothing = "gaussian", #gaussian or rollmean
# 	smoothing_window = 10, #both types of smoothing require a window
# 	normalisation = "libsize", #libsize or maxpeak or none
# 	output = "plot.pdf")

# setwd("~/Ule/charlotte/CLIPplotR")
# 
# opt <- list(xlinks = "tardbp-esc-m-p2lox-gfp-tdp43-20151212_trimmed_single.bedgraph tardbp-ngfp-esc-m-p2lox-gfp-tdp43-20151212_trimmed_single.bedgraph", #from https://imaps.genialis.com/iclip
#             # track_names = "tdp43_1 tdp43_2", # or can be NULL
#             gtf = "gencode.vM22.annotation.gtf.gz", #"~/Ule/ref/gencode.v27.annotation.gtf.gz"
#             # region = "chr3:35754106:35856276:+",
#             # gene = "ENSMUSG00000037400", #ID or name "Atp11b"
#             region = "Atp11b", #ID or name "Atp11b"
#             smoothing = "rollmean", #gaussian or rollmean or none
#             smoothing_window = 1000, #both types of smoothing require a window
#             normalisation = "none", #libsize or maxpeak or none
#             output = "plot.pdf")

# ==========
# Part 0 - Create annotation and get relevant regions in bedgraphs
# ==========

# Create rosetta and TxDb
message("Reading in annotation from GTF")

gtf <- rtracklayer::import.gff2(opt$gtf)
genes.gr <- gtf[gtf$type == "gene"]

# Create TxDb if doesn't already exist
if(file.exists(gsub(".gtf.gz|.gf", ".sqlite", opt$gtf))) {

  message("Loading pre-existing annotation database")    
  TxDb <- loadDb(gsub(".gtf.gz|.gf", ".sqlite", opt$gtf))
  
} else {

  message("Creating annotation database for future runs")
  TxDb <- makeTxDbFromGFF(opt$gtf)
  saveDb(TxDb, gsub(".gtf.gz|.gf", ".sqlite", opt$gtf))

}

# Define region
if(is.null(opt$region)) {
  
  stop("Need to supply a region or a gene id or name")
  
} else if(grepl("^ENS", opt$region)) {
  
  region.gr <- genes.gr[grepl(opt$region, genes.gr$gene_id)]
  
} else if(grepl("^[A-Z]", opt$region)) {
  
  region.gr <- genes.gr[grepl(opt$region, genes.gr$gene_name)]
  
} else {
  
  region.gr <- GRanges(seqnames = Rle(sapply(strsplit(opt$region, ":"), "[[", 1)),
                       ranges = IRanges(start = as.integer(sapply(strsplit(opt$region, ":"), "[[", 2)),
                                        end = as.integer(sapply(strsplit(opt$region, ":"), "[[", 3))),
                       strand = Rle(sapply(strsplit(opt$region, ":"), "[[", 4)))
  
}

seqlevels(region.gr) <- as.character(unique(seqnames(region.gr))) # Cut down to one seqlevels for later comparision

# ==========
# Part 1 - top half: normalised and smoothed tracks
# ==========

ImportiMapsBedgraph <- function(bedgraph.file) {
  
  bg <- import.bedGraph(bedgraph.file)
  
  # Assign strand based on score
  strand(bg)[bg$score < 0] <- "-"
  strand(bg)[bg$score > 0] <- "+"
  
  # Convert scores to positives now that strands assigned
  bg$score <- abs(bg$score)
  
  return(bg)
  
}

# Read in xlinks
message("Loading bedgraphs")
xlinks <- strsplit(opt$xlinks, " ")[[1]] %>% lapply(., ImportiMapsBedgraph)
libSizes <- lapply(xlinks, function(x) { sum(abs(x$score)) })

# Subset for region and add in 0 count position
xlinks <- lapply(xlinks, function(x) {
  
  # Subset bedgraph to region
  xlinks.gr <- subsetByOverlaps(x, region.gr, ignore.strand = FALSE)
  
  # Get single position 0 counts
  zero.gr <- setdiff(region.gr, xlinks.gr)
  zero.gr <- unlist(tile(zero.gr, width = 1))
  zero.gr$score <- 0
  
  # Combine and sanity check
  xlinks.gr <- sort(c(xlinks.gr, zero.gr))
  
  stopifnot(all(width(xlinks.gr) == 1 & reduce(xlinks.gr) == region.gr))
  return(xlinks.gr)
  
})

# Names of bedgraphs : If name is supplied use that, if not then generate a name from the file name
if (!is.null(opt$track_names)) {
  
  track_names=strsplit(opt$label, " ")[[1]]
  
} else {
  
  track_names = lapply(strsplit(opt$xlinks, " ")[[1]], function(x) gsub(".bedgraph", "", basename(x)))
  
}

# Create dataframe for plotting
xlinks_df <- lapply(xlinks, as.data.frame)

for (i in 1:length(xlinks_df)) {
  
  xlinks_df[[i]]$sample <- track_names[[i]]
  xlinks_df[[i]]$libSize <- libSizes[[i]]
  
}

xl_df <- suppressWarnings(dplyr::bind_rows(xlinks_df))

# TODO: Fix this warning:
# Warning messages:
#   1: In bind_rows_(x, .id) : Unequal factor levels: coercing to character
# 2: In bind_rows_(x, .id) :
#   binding character and factor vector, coercing into character vector
# 3: In bind_rows_(x, .id) :
#   binding character and factor vector, coercing into character vector

xl_df <- as.data.frame(xl_df) # Isn't it already a data.frame?

# Do the normalisation
message("Normalising")
xl_df <- as.data.frame(switch(opt$normalisation, 
                              "libsize"=xl_df %>% dplyr::group_by(sample) %>% dplyr::mutate(norm=(score * 1e6)/libSize),
                              "maxpeak"=xl_df %>% dplyr::group_by(sample) %>% dplyr::mutate(norm=score/max(score)),
                              "none"=xl_df %>% dplyr::group_by(sample) %>% dplyr::mutate(norm=score)))


# Do the smoothing
message("Smoothing")
xl_df <- as.data.frame(switch(opt$smoothing, 
                              "gaussian"=xl_df %>% dplyr::group_by(sample) %>% dplyr::mutate(smoothed= smth.gaussian(norm, window=opt$smoothing_window)),
                              "rollmean"=xl_df %>% dplyr::group_by(sample) %>% dplyr::mutate(smoothed= rollmean(norm, opt$smoothing_window, fill=0)),
                              "none"=xl_df %>% dplyr::group_by(sample) %>% dplyr::mutate(smoothed= norm)))

# Plot top half
p.iclip <- ggplot(xl_df,aes(x=start,y=smoothed, group=sample, color=sample)) +
  geom_line() +
  labs(title = opt$region,
       x = "",
       y = "Crosslink signal",
       colour = "") +
  scale_colour_tableau(palette = "Tableau 10") +
  theme_cowplot() + theme(legend.position = "top")

# ==========
# Part 2 - bottom half: gene structures
# ==========

message("Creating annotation track")
annot.gr <- biovizBase::crunch(TxDb, which = region.gr)
annot.grl <- split(annot.gr, annot.gr$tx_name)

# Need to trim exons so they don't overlap utrs to ensure they aren't overplotted
annot.grl <- GRangesList(lapply(annot.grl, function(x) {

  exon <- x[x$type == "exon"]  
  utr <- x[x$type == "utr"]
  strand(utr) <- unique(strand(exon)) # Otherwise utr isn't stranded, but needed for setdiff later
  rest <- x[!x$type %in% c("utr", "exon")]

  exon <- GenomicRanges::setdiff(exon, utr, ignore.strand = FALSE)  
  mcols(exon)$type <- "exon"
  
  return(sort(c(rest, exon, utr)))
  
}))

# Add gene name to labels
gene_names <- gtf$gene_name[match(names(annot.grl), gtf$transcript_id)]
names(annot.grl) <- paste0(gene_names, " - ", names(annot.grl))

p.annot <- ggplot(data = annot.grl) +
  geom_alignment(cds.rect.h = 0.1, fill = "black", stat = "identity") +
  labs(x = "Coordinate")

# Select out ggplot object
p.annot <- p.annot@ggplot

# ==========
# Part 3 - combine plots
# ==========

# plot_grid(p.iclip, p.annot, align = "hv", axis = "tlbr", nrow = 2, rel_heights = c(1, 2))
ggsave(plot_grid(p.iclip, p.annot, align = "hv", axis = "tlbr", nrow = 2, rel_heights = c(1, 2)), width = 297, height = 210, units = "mm", filename = opt$output)
message("Completed")