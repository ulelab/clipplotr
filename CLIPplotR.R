#!/usr/bin/env Rscript

# Script to plot multiple CLIP tracks with gene structures
# A. M. Chakrabarti
# C. Capitanchik
# Last updated: 17th April 2020

# ==========
# Preamble
# ==========

suppressPackageStartupMessages(library(optparse))

option_list <- list(make_option(c("-x", "--xlinks"), action = "store", type = "character", help = "Input iCLIP bedgraphs (space separated)"),
                    make_option(c("-l", "--label"), action = "store", type = "character", help = "iCLIP bedgraph labels (space separated)"),
                    make_option(c("-c", "--colours"), action = "store", type = "character", help = "iCLIP bedgraph colours (space separated)"),
                    make_option(c("", "--groups"), action = "store", type = "character", help = "Grouping of iCLIP bedgraphs for separate plots (space separated"),
                    make_option(c("-p", "--peaks"), action = "store", type = "character", help = "BED file of peaks (space separated)"),                    
                    make_option(c("-g", "--gtf"), action = "store", type = "character", help = "Reference GTF file (Gencode)"),
                    make_option(c("-r", "--region"), action = "store", type = "character", help = "Region of interest as chr3:35754106:35856276:+ or gene as ENSMUSG00000037400 or Atp11b"),
                    make_option(c("-n", "--normalisation"), action = "store", type = "character", help = "Normalisation options: none, maxpeak, libsize [default %default]", default = "libsize"),
                    make_option(c("-s", "--smoothing"), action = "store", type = "character", help = "Smoothing options: none, rollmean, spline, gaussian [default %default]", default = "rollmean"),
                    make_option(c("-w", "--smoothing_window"), action = "store", type = "integer", help = "Smoothing window [default %default]", default = 100),
                    make_option(c("-a", "--annotation"), action = "store", type = "character", help = "Annotation options: original, gene, transcript, none [default %default]", default = "original"),
                    make_option(c("", "--size_x"), action = "store", type = "integer", help = "Plot size in mm (x) [default: %default]", default = 210),
                    make_option(c("", "--size_y"), action = "store", type = "integer", help = "Plot size in mm (y) [default: %default]", default = 297),
                    make_option(c("-o", "--output"), action = "store", type = "character", help = "Output plot filename"))

opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# print(opt)

# ==========
# Checks for minimum parameters
# ==========

if(is.null(opt$xlinks)) {

  message("ERROR: No iCLIP bedgraphs supplied")
  quit(save = "no")

}

if(is.null(opt$gtf)) {

  message("ERROR: No Reference GTF supplied")
  quit(save = "no")

}

if(is.null(opt$region)) {

  message("ERROR: No region defined")
  quit(save = "no")

}

if(is.null(opt$output)) {

  message("ERROR: No output defined")
  quit(save = "no")

}

# ==========
# Load libraries
# ==========

suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggthemes))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(smoother))
suppressPackageStartupMessages(library(zoo))
suppressPackageStartupMessages(library(data.table))

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

# opt <- list(xlinks = "tardbp-316del346-egfp-hek293-hd1-6-merged-20180829-ju.bedgraph.gz tardbp-316del346-egfp-hek293-hd2-5-merged-20180829-ju.bedgraph.gz tardbp-egfp-hek293-hd1-6-merged-20180829-ju.bedgraph.gz tardbp-egfp-hek293-hd2-5-merged-20180829-ju.bedgraph.gz",
#             peaks = "tdp43_binding_sites_genome_browser3.bed",
#             gtf = "~/Ule/charlotte/CLIPplotR/gencode.v29.primary_assembly.annotation.gtf.gz",
#             # region = "chr11:65503537:65505019:+",
#             region = "TARDBP",
#             output = "test_annot5.pdf",
#             normalisation = "libsize",
#             smoothing = "rollmean",
#             smoothing_window = 50,
#             annotation = "gene")
# 
# setwd("~/Downloads/martina_grouped_iclip_libs")

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
xlink.files <- strsplit(opt$xlinks, " ")[[1]]
xlinks <- lapply(xlink.files, ImportiMapsBedgraph)
libSizes <- lapply(xlinks, function(x) { sum(abs(x$score)) })

# Subset for region and add in 0 count position
xlinks <- lapply(xlinks, function(x) SubsetBedgraph(gr = x, selected.region.gr = region.gr))

# Names of bedgraphs : If name is supplied use that, if not then generate a name from the file name
if (!is.null(opt$label)) {
  
  track_names=strsplit(opt$label, " ")[[1]]
  
} else {
  
  track_names = lapply(strsplit(opt$xlinks, " ")[[1]], function(x) gsub(".bedgraph", "", basename(x)))
  
}

# Create datatable for plotting
xlinks.dt <- lapply(seq_along(xlinks), function(i) {
  
  dt <- as.data.table(xlinks[[i]])
  dt$sample <- track_names[[i]]
  dt$libSize <- libSizes[[i]]
  return(dt)
  
})
xlinks.dt <- rbindlist(xlinks.dt)  
setkey(xlinks.dt, sample)

# Do the normalisation
message("Normalising")
xlinks.dt[, norm := switch(opt$normalisation, 
                           "libsize" = (score * 1e6)/libSize,
                           "maxpeak" = score/max(score),
                           "none" = score), 
          by = sample]

# Do the smoothing
message("Smoothing")
xlinks.dt[, smoothed := switch(opt$smoothing, 
                               "gaussian" = smth.gaussian(norm, window = opt$smoothing_window),
                               "rollmean" = rollmean(norm, opt$smoothing_window, fill = 0),
                               "none" = norm), 
          by = sample]

# Assign groups
if(!is.null(opt$groups)) {

  groups.dt <- data.table(sample = track_names, grp = strsplit(opt$group, " ")[[1]])
  xlinks.dt <- merge(xlinks.dt, groups.dt, by = sample)

}

# Plot top half
if(is.null(opt$colours)) {

p.iclip <- ggplot(xlinks.dt, aes(x = start, y = smoothed, group = sample, color = sample)) +
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

    p.iclip <- ggplot(xlinks.dt, aes(x = start, y = smoothed, group = sample, color = sample)) +
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

  peaks.files <- strsplit(opt$peaks, " ")[[1]]
  peaks.grl <- lapply(peaks.files, import.bed)
  peaks.grl <- lapply(peaks.grl, function(x) subsetByOverlaps(x, region.gr, ignore.strand = FALSE))

  peaks.dt <- rbindlist(lapply(peaks.grl, as.data.table))
  peaks.dt$exp <- rep(gsub(".bed", "", basename(peaks.files)), elementNROWS(peaks.grl))
  peaks.dt$centre <- with(peaks.dt, start + width/2) # Replace this with data.table syntax, but need to check groups first

  if(nrow(peaks.dt) == 0) {

    p.peaks <- ggplot() + theme_cowplot() + theme(axis.line = element_blank())

  } else {

  p.peaks <- ggplot(peaks.dt, aes(x = centre, width = width, y = exp)) +
    geom_tile() +
    scale_y_discrete(breaks = gsub(".bed", "", basename(peaks.files)),
                     limits = gsub(".bed", "", basename(peaks.files))) +
    xlim(start(region.gr), end(region.gr)) +
    labs(y = "",
         x = "") +
    theme_cowplot()

  }

} else {

  p.peaks <- ggplot() + theme_cowplot() + theme(axis.line = element_blank())

}

# ==========
# Part 2 - bottom half: gene structures
# ==========

message("Creating annotation track")

rosetta.dt <- as.data.table(mcols(genes.gr))[, .(gene_id, gene_name)]
setkey(rosetta.dt, gene_id)

# ==========
# Plot gene kind of structure
# ==========

if(opt$annotation == "gene") {

  exons_g.grl <- exonsBy(TxDb, by = "gene")
  
  sel.genes <- subsetByOverlaps(genes.gr, region.gr)
  sel.exons <- exons(TxDb, filter = list(gene_id = sel.genes$gene_id), columns = "gene_id")
  
  # Region (for arrows)
  region.gr$gene_id <- NULL
  region.tiled <- rep(region.gr, length(sel.genes)) # Create one for each gene
  region.tiled <- tile(region.tiled, width = 500)
  
  sel.region.tiled <- GRangesList(lapply(seq_along(sel.genes), function(i) {
    gr <- subsetByOverlaps(region.tiled[[i]], sel.genes[i], type = "any")
    # Need to adjust for ends
    start(gr[1]) <- start(sel.genes[i])
    end(gr[length(gr)]) <- end(sel.genes[i])
    return(gr)
  }))
  
  names(sel.region.tiled) <- sel.genes$gene_id
  sel.region.tiled.dt <- as.data.table(sel.region.tiled)
  setnames(sel.region.tiled.dt, "group_name", "gene_id")
  
  # Exons
  sel.exons.dt <- as.data.table(sel.exons)
  sel.exons.dt$gene_id <- as.character(sel.exons.dt$gene_id) # otherwise type is AsIs
  sel.exons.dt <- merge(sel.exons.dt, unique(sel.region.tiled.dt[, .(group, gene_id)]), by = "gene_id")
  
  # Plot
  p.annot <- ggplot() +
    geom_segment(data = sel.region.tiled.dt, mapping = aes(x = start, xend = end, y = group, yend = group), arrow = arrow(length = unit(0.1, "cm")), colour = "grey50") +
    geom_rect(data = sel.exons.dt, mapping = aes(xmin = start, xmax = end, ymin = group - 0.25, ymax = group + 0.25, fill = gene_id)) +
    theme_cowplot() + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position = "bottom") +
    scale_fill_tableau() +
    labs(x = "Coordinate",
         y = "",
         fill = "") +
    coord_cartesian(xlim = c(start(region.gr), end(region.gr)))

}

# ==========
# Plot transcript kind of structure
# ==========

if(opt$annotation == "transcript") {
  
  # Get transcripts that overlap region and order for plotting
  sel.tx_genes <- transcriptsByOverlaps(TxDb, region.gr, columns = c("gene_id", "tx_name"))
  sel.tx_genes <- sel.tx_genes[order(width(sel.tx_genes), decreasing = TRUE)]
  tx.order.dt <- data.table(transcript_id = sel.tx_genes$tx_name,
                            gene_id = unlist(sel.tx_genes$gene_id),
                            centre = start(sel.tx_genes) + width(sel.tx_genes)/2)[, group := 1:.N]
  setkey(tx.order.dt, gene_id)
  tx.order.dt <- rosetta.dt[tx.order.dt]
  tx.order.dt[, gene := paste0(gene_name, " | ", gene_id)]
  
  # Region (for arrows)
  region <- region.gr
  region$gene_id <- NULL
  region.tiled <- rep(region, length(sel.tx_genes))
  region.tiled$tx_name <- sel.tx_genes$tx_name
  region.tiled <- tile(region.tiled, width = round(width(region)/15, -1))

  sel.region.tiled <- GRangesList(lapply(seq_along(sel.tx_genes), function(i) {
    gr <- subsetByOverlaps(region.tiled[[i]], sel.tx_genes[i], type = "any")
    start(gr[1]) <- start(sel.tx_genes[i])
    end(gr[length(gr)]) <- end(sel.tx_genes[i])
    return(gr)
  }))
  
  names(sel.region.tiled) <- sel.tx_genes$tx_name
  sel.region.tiled.dt <- as.data.table(sel.region.tiled)
  sel.region.tiled.dt[, group := NULL]
  setnames(sel.region.tiled.dt, "group_name", "transcript_id")
  sel.region.tiled.dt <- merge(sel.region.tiled.dt, tx.order.dt, by = "transcript_id")
  
  # CDS
  cds_tx <- cdsBy(TxDb, by = "tx", use.names = TRUE)
  sel.cds_tx <- cds_tx[names(cds_tx) %in% sel.tx_genes$tx_name]
  sel.cds_tx.dt <- as.data.table(sel.cds_tx)
  sel.cds_tx.dt[, group := NULL]
  setnames(sel.cds_tx.dt, "group_name", "transcript_id")
  sel.cds_tx.dt <- merge(sel.cds_tx.dt, tx.order.dt, by = "transcript_id")
  
  # UTR5
  utr5_tx <- fiveUTRsByTranscript(TxDb, use.names = TRUE)
  sel.utr5_tx <- utr5_tx[names(utr5_tx) %in% sel.tx_genes$tx_name]
  sel.utr5_tx.dt <- as.data.table(sel.utr5_tx)
  sel.utr5_tx.dt[, group := NULL]
  setnames(sel.utr5_tx.dt, "group_name", "transcript_id") 
  sel.utr5_tx.dt <- merge(sel.utr5_tx.dt, tx.order.dt, by = "transcript_id")
    
  # UTR3
  utr3_tx <- threeUTRsByTranscript(TxDb, use.names = TRUE)
  sel.utr3_tx <- utr3_tx[names(utr3_tx) %in% sel.tx_genes$tx_name]
  sel.utr3_tx.dt <- as.data.table(sel.utr3_tx)
  sel.utr3_tx.dt[, group := NULL]
  setnames(sel.utr3_tx.dt, "group_name", "transcript_id")
  sel.utr3_tx.dt <- merge(sel.utr3_tx.dt, tx.order.dt, by = "transcript_id")
  
  # Exons
  exons_tx <- exonsBy(TxDb, by = "tx", use.names = TRUE)
  sel.exons_tx <- exons_tx[names(exons_tx) %in% sel.tx_genes$tx_name[!sel.tx_genes$tx_name %in% names(sel.cds_tx)]] # Don't want genes with CDS, just e.g. ncRNA
  sel.exons_tx.dt <- as.data.table(sel.exons_tx)
  sel.exons_tx.dt[, group := NULL]
  setnames(sel.exons_tx.dt, "group_name", "transcript_id")
  sel.exons_tx.dt <- merge(sel.exons_tx.dt, tx.order.dt, by = "transcript_id")
  
  # Plot
  p.annot <- ggplot() +
    geom_segment(data = sel.region.tiled.dt, mapping = aes(x = start, xend = end, y = group, yend = group), arrow = arrow(length = unit(0.1, "cm")), colour = "grey50") +
    # geom_rect(data = sel.exons_tx.dt, mapping = aes(xmin = start, xmax = end, ymin = group - 0.25, ymax = group + 0.25), fill = "black") +
    # geom_rect(data = sel.cds_tx.dt, mapping = aes(xmin = start, xmax = end, ymin = group - 0.25, ymax = group + 0.25), fill = "black") +
    # geom_rect(data = sel.utr5_tx.dt, mapping = aes(xmin = start, xmax = end, ymin = group - 0.15, ymax = group + 0.15), fill = "black") +
    # geom_rect(data = sel.utr3_tx.dt, mapping = aes(xmin = start, xmax = end, ymin = group - 0.15, ymax = group + 0.15), fill = "black") +
    # geom_text(data = tx.order.dt, aes(label = transcript_id, y = group, x = centre), size = 3, vjust = 0) +
    geom_rect(data = sel.exons_tx.dt, mapping = aes(xmin = start, xmax = end, ymin = group - 0.25, ymax = group + 0.25, fill = gene)) +
    geom_rect(data = sel.cds_tx.dt, mapping = aes(xmin = start, xmax = end, ymin = group - 0.25, ymax = group + 0.25, fill = gene)) +
    geom_rect(data = sel.utr5_tx.dt, mapping = aes(xmin = start, xmax = end, ymin = group - 0.15, ymax = group + 0.15, fill = gene)) +
    geom_rect(data = sel.utr3_tx.dt, mapping = aes(xmin = start, xmax = end, ymin = group - 0.15, ymax = group + 0.15, fill = gene)) +
    scale_fill_tableau() +
    theme_cowplot() + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position = "bottom") +
    labs(x = "Coordinate",
         y = "",
         fill = "") +
    coord_cartesian(xlim = c(start(region), end(region)))
  
  
  if(opt$show_tx_name == TRUE) {
    
    p.annot <- p.annot + geom_text(data = tx.order.dt, aes(label = transcript_id, y = group + 0.3, x = centre), size = 3, vjust = 0)
    # Need to fix ones that are outside of window
    
  }
  
}

# ==========
# Omit annotation plot
# ==========

if(opt$annotation == "none") {

  p.annot <- ggplot() + theme_cowplot() + theme(axis.line = element_blank())

}

# ==========
# Old style
# ==========

if(opt$annotation == "original") {

suppressPackageStartupMessages(library(ggbio))

annot.gr <- biovizBase::crunch(TxDb, which = region.gr)
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

}

# ==========
# Part 3 - combine plots
# ==========

# plot_grid(p.iclip, p.annot, align = "hv", axis = "tlbr", nrow = 2, rel_heights = c(1, 2))

if(opt$annotation == "original") {
  ggsave(plot_grid(p.iclip, p.peaks, p.annot, align = "hv", axis = "tlbr", nrow = 3, rel_heights = c(1, 0.5, 2)), height = opt$size_y, width = opt$size_x, units = "mm", filename = opt$output)
}

if(opt$annotation == "gene") {
  ggsave(plot_grid(p.iclip, p.peaks, p.annot, align = "hv", axis = "tlbr", nrow = 3, rel_heights = c(2, 1, 1)), height = opt$size_y, width = opt$size_x, units = "mm", filename = opt$output)
}

if(opt$annotation == "transcript") {
  ggsave(plot_grid(p.iclip, p.peaks, p.annot, align = "hv", axis = "tlbr", nrow = 3, rel_heights = c(2, 1, 3)), height = opt$size_y, width = opt$size_x, units = "mm", filename = opt$output)
}

if(opt$annotation == "none") {
  ggsave(plot_grid(p.iclip, p.peaks, p.annot, align = "hv", axis = "tlbr", nrow = 3, rel_heights = c(2, 1, 1)), height = opt$size_y, width = opt$size_x, units = "mm", filename = opt$output)
}


message("Completed")
