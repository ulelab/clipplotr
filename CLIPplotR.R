#!/usr/bin/env Rscript

# Script to plot multiple CLIP tracks with gene structures
# A. M. Chakrabarti
# C. Capitanchik
# Last updated: 9th February 2020

# ==========
# Preamble
# ==========

CheckAndLoad <- function(package) {
  
  if(!suppressPackageStartupMessages(require(package, character.only = TRUE, quietly = TRUE))) {
    
    message("Installing ", package)
    install.packages(package, character.only = TRUE, repos = "https://cloud.r-project.org")
    suppressPackageStartupMessages(library(package, character.only = TRUE, quietly = TRUE))
    
  }
  
}

CheckAndLoad("optparse")

option_list <- list(make_option(c("-x", "--xlinks"), action = "store", type = "character", help = "Input iCLIP bedgraphs (space separated)"),
                    make_option(c("-l", "--label"), action = "store", type = "character", help = "iCLIP bedgraph labels (space separated)"),
                    make_option(c("-c", "--colours"), action = "store", type = "character", help = "iCLIP bedgraph colours (space separated)"),
                    make_option(c("-f", "--groups"), action = "store", type = "character", help = "Grouping of iCLIP bedgraphs for separate plots (space separated"),
                    make_option(c("-p", "--peaks"), action = "store", type = "character", help = "BED file of peaks (space separated"),                    
                    make_option(c("-g", "--gtf"), action = "store", type = "character", help = "Reference gtf (Gencode)"),
                    make_option(c("-r", "--region"), action = "store", type = "character", help = "Region of interest as chr3:35754106:35856276:+ or gene as ENSMUSG00000037400 or Atp11b"),
                    make_option(c("-n", "--normalisation"), action = "store", type = "character", help = "Normalisation options: none, maxpeak, libsize [default %default]", default = "libsize"),
                    make_option(c("-s", "--smoothing"), action = "store", type = "character", help = "Smoothing options: none, rollmean, spline, gaussian [default %default]", default = "rollmean"),
                    make_option(c("-w", "--smoothing_window"), action = "store", type = "integer", help = "Smoothing window [default %default]", default = 100),
                    make_option(c("-x", "--size_x"), action = "store", type = "integer", help = "Plot size in mm (x)", default = 300),
                    make_option(c("-y", "--size_y"), action = "store", type = "integer", help = "Plot size in mm (y)", default = 300),
                    make_option(c("-o", "--output", action = "store", type = "character", help = "Output plot filename")))
opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Load CRAN packages
packages <- c("BiocManager", "ggplot2", "ggthemes", "cowplot", "smoother", "zoo", "dplyr")
for(package in packages) CheckAndLoad(package)

# Load Bioconductor packages
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
#             region = "Ank3", #ID or name "Atp11b"
#             smoothing = "rollmean", #gaussian or rollmean or none
#             smoothing_window = 1000, #both types of smoothing require a window
#             normalisation = "none", #libsize or maxpeak or none
#             output = "plot.pdf")
# 
# opt <- list(region = "chr15:68206992:68210000:-", gtf = "gencode.v29.primary_assembly.annotation.gtf.gz")

# ==========
# Functions
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

SubsetBedgraph <- function(gr, selected.region.gr) {
  
  xlinks.gr <- unlist(tile(selected.region.gr, width = 1))
  
  ol <- findOverlaps(xlinks.gr, gr)
  xlinks.gr$score <- NA
  xlinks.gr[queryHits(ol)]$score <- gr[subjectHits(ol)]$score
  xlinks.gr$score[is.na(xlinks.gr$score)] <- 0
  
  return(xlinks.gr)
  
}

# ==========
# Part 0 - Create annotation and get relevant regions in bedgraphs
# ==========

# Create rosetta and TxDb
message("Reading in annotation from GTF")

gtf <- rtracklayer::import.gff2(opt$gtf)
genes.gr <- gtf[gtf$type == "gene"]

# Create TxDb if doesn't already exist
if(file.exists(gsub(".gtf.gz|.gtf", ".sqlite", opt$gtf))) {

  message("Loading pre-existing annotation database")    
  TxDb <- loadDb(gsub(".gtf.gz|.gtf", ".sqlite", opt$gtf))
  
} else {

  message("Creating annotation database for future runs")
  TxDb <- makeTxDbFromGFF(opt$gtf)
  saveDb(TxDb, gsub(".gtf.gz|.gtf", ".sqlite", opt$gtf))

}

# Define region
if(is.null(opt$region)) {
  
  stop("Need to supply a region or a gene id or name")
  
} else if(grepl("^ENS", opt$region)) {
  
  region.gr <- genes.gr[grepl(opt$region, genes.gr$gene_id)]
  
} else if(grepl(":", opt$region)) {
  
  region.gr <- GRanges(seqnames = Rle(sapply(strsplit(opt$region, ":"), "[[", 1)),
                       ranges = IRanges(start = as.integer(sapply(strsplit(opt$region, ":"), "[[", 2)),
                                        end = as.integer(sapply(strsplit(opt$region, ":"), "[[", 3))),
                       strand = Rle(sapply(strsplit(opt$region, ":"), "[[", 4)))
} else {
  # region.gr <- genes.gr[grepl(opt$region, genes.gr$gene_name)]
  region.gr <- genes.gr[opt$region == genes.gr$gene_name]  

}

seqlevels(region.gr) <- as.character(unique(seqnames(region.gr))) # Cut down to one seqlevels for later comparision

# ==========
# Part 1a - top half: normalised and smoothed tracks
# ==========

# Read in xlinks
message("Loading bedgraphs")
xlinks <- strsplit(opt$xlinks, " ")[[1]] %>% lapply(., ImportiMapsBedgraph)
libSizes <- lapply(xlinks, function(x) { sum(abs(x$score)) })

# Subset for region and add in 0 count position
xlinks <- lapply(xlinks, function(x) SubsetBedgraph(gr = x, selected.region.gr = region.gr))

# Names of bedgraphs : If name is supplied use that, if not then generate a name from the file name
if (!is.null(opt$label)) {
  
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

# Assign groups
if(!is.null(opt$groups)) {

  groups.df <- data.table(sample = track_names, grp = strsplit(opt$group, " ")[[1]])
  xl_df <- merge(xl_df, groups.df, by = sample)

}

# Plot top half
if(is.null(opt$colours)) {

p.iclip <- ggplot(xl_df,aes(x=start,y=smoothed, group=sample, color=sample)) +
  geom_line() +
  labs(title = opt$region,
       x = "",
       y = "Crosslink signal",
       colour = "") +
  scale_colour_tableau(palette = "Tableau 10") +
  theme_cowplot() + theme(legend.position = "top") +
      xlim(start(region.gr),end(region.gr))

  } else {

    cols <- strsplit(opt$colours, " ")[[1]]
    names(cols) <- track_names

    p.iclip <- ggplot(xl_df,aes(x=start,y=smoothed, group=sample, color=sample)) +
      geom_line() +
      labs(title = opt$region,
           x = "",
           y = "Crosslink signal",
           colour = "") +
      scale_colour_manual(values = cols) +
      theme_cowplot() + theme(legend.position = "top") +
      xlim(start(region.gr),end(region.gr))

}

# Facet if groups
if(!is.null(opt$groups)) {

  p.iclip <- p.iclip + facet_grid(grp ~ .)

}

# ==========
# Part 1b - top half: peak annotation
# ==========

# TODO: add in peak names as an option to supply

if(!is.null(opt$peaks)) {

  peak.files <- strsplit(opt$peaks, " ")[[1]]
  peaks.grl <- lapply(peaks.files, import.bed)
  peaks.grl <- lapply(peaks.grl, function(x) subsetByOverlaps(x, region.gr, ignore.strand = FALSE))

  peaks_df <- suppressWarnings(dplyr::bind_rows(lapply(peaks.grl, as.data.frame)))
  peaks_df$exp <- rep(gsub(".bed", "", basename(peaks.files)), elementNROWS(peaks.grl))
  peaks_df$centre <- with(peaks_df, start + width/2)

  if(nrow(peaks_df) == 0) {

    p.peaks <- ggplot() + theme_cowplot() + theme(axis.line = element_blank())

  } else {

  p.peaks <- ggplot(peaks_df, aes(x = centre, width = width, y = exp)) +
    geom_tile() +
    scale_y_discrete(breaks = gsub(".bed", "", basename(peaks.files)),
                     limits = gsub(".bed", "", basename(peaks.files))) +
    xlim(start(region.gr), end(region.gr)) +
    labs(y = "",
         x = "") +
    theme_cowplot()

  }

} else {

  p.peaks <- ggplot()

}

# ==========
# Part 2 - bottom half: gene structures
# ==========

message("Creating annotation track")

# # First deal with purely intronic regions, because biovisbase doesn't handle this
# if(is.null(exonsByOverlaps(TxDb, region.gr)$tx_id)){
#   overlapping_tscripts = transcriptsByOverlaps(TxDb, region.gr)$tx_id
#   overlapping_tscripts_name = transcriptsByOverlaps(TxDb, region.gr)$tx_name
#   introns_in_trscripts = intronsByTranscript(TxDb)[overlapping_tscripts,]
  
#   #Get overlapping introns
  
#   if (length(overlapping_tscripts)>1){
#     annot.gr=list()
#     for (i in 1:length(overlapping_tscripts)){
#       intron_reg = region.gr[subjectHits(findOverlaps(region.gr, introns_in_trscripts[[i]])),]
#       intron_reg$type="gap"
#       intron_reg$tx_id=overlapping_tscripts[[i]]
#       intron_reg$tx_name=overlapping_tscripts_name[[i]]}
#     annot.gr = c(final_int_reg,intron_reg)
#   } else {
#     intron_reg = region.gr[subjectHits(findOverlaps(region.gr, introns_in_trscripts)),]
#     intron_reg$type="gap"
#     intron_reg$tx_id=overlapping_tscripts
#     intron_reg$tx_name=overlapping_tscripts_name
#     annot.gr = intron_reg
#   }
# } else {
annot.gr <- biovizBase::crunch(TxDb, which = region.gr)
# }

if (is.na(annot.gr)){
  ggsave(p.iclip,height = 300, width = 300, units = "mm", filename = opt$output)
  message("Completed")
  quit(save="no")
}


annot.gr$tx_name <- as.character(annot.gr$tx_name)
annot.grl <- split(annot.gr, annot.gr$tx_name)


# Need to trim exons so they don't overlap utrs to ensure they aren't overplotted
# but only if there are exons 
annot.grl <- GRangesList(lapply(annot.grl, function(x) {
  if(any(grepl("exon",x$type))){
    exon <- x[x$type == "exon"]  
    utr <- x[x$type == "utr"]
    strand(utr) <- unique(strand(exon)) # Otherwise utr isn't stranded, but needed for setdiff later
    rest <- x[!x$type %in% c("utr", "exon")]
    
    exon <- GenomicRanges::setdiff(exon, utr, ignore.strand = FALSE)  
    mcols(exon)$type <- "exon"
    
    return(sort(c(rest, exon, utr)))} else {
      return(x)
  }
}))

# # For bug check
# annot.grl <- GRangesList(lapply(1:length(annot.grl), function(i) {
#   
#   message(i)
#   x <- annot.grl[[i]]
#   exon <- x[x$type == "exon"]  
#   utr <- x[x$type == "utr"]
#   strand(utr) <- unique(strand(exon)) # Otherwise utr isn't stranded, but needed for setdiff later
#   rest <- x[!x$type %in% c("utr", "exon")]
#   
#   exon <- GenomicRanges::setdiff(exon, utr, ignore.strand = FALSE)  
#   mcols(exon)$type <- "exon"
#   
#   return(sort(c(rest, exon, utr)))
#   
# }))

# Add gene name to labels
gene_names <- gtf$gene_name[match(names(annot.grl), gtf$transcript_id)]
names(annot.grl) <- paste0(gene_names, " - ", names(annot.grl))

p.annot <- ggplot(data = annot.grl) +
  geom_alignment(cds.rect.h = 0.1, length = unit(0.1, "cm"), fill = "black", stat = "identity") +
  labs(x = "Coordinate") +
      xlim(start(region.gr),end(region.gr)) +
  theme_cowplot() + theme(axis.line.y = element_blank())

# Select out ggplot object
p.annot <- p.annot@ggplot

# ==========
# Part 3 - combine plots
# ==========

# plot_grid(p.iclip, p.annot, align = "hv", axis = "tlbr", nrow = 2, rel_heights = c(1, 2))
ggsave(plot_grid(p.iclip, p.peaks, p.annot, align = "hv", axis = "tlbr", nrow = 3, rel_heights = c(1, 1, 2)), height = opt$size_y, width = opt$size_x, units = "mm", filename = opt$output)
message("Completed")
