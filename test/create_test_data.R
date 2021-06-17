
# Create test data from Zarnack et al. (2013)
# A. M. Chakrabarti
# 27th June 2020

library(rtracklayer)

setwd("~/Dropbox (The Francis Crick)/clipplotr/processed_data/zarnack")

# chr1:207513000:207515000:+

region.gr <- GRanges(seqnames = "chr1",
                     ranges = IRanges(start = 207510000, end = 207520000),
                     strand = "+")

# iCLIP

clip.files <- list.files(pattern = "_all_xlink_events.bedgraph.gz$", full.names = TRUE)
clip.bed <- lapply(clip.files, import.bedGraph, which = region.gr)
lapply(seq_along(clip.bed), function(i) export.bedGraph(GRanges(clip.bed[[i]]), paste0("test_", basename(clip.files[i])))) # GRanges() to remove UCSC track data

# Alu track

subset.alu.bed <- import.bed("Alu_rev.bed.gz", which = region.gr)
export.bed(subset.alu.bed, "test_Alu_rev.bed.gz")

# Bigwigs

rna.files <- list.files(pattern = "_plus.bigwig$", full.names = TRUE)
rna.bw <- lapply(rna.files, import.bw, selection = region.gr)
lapply(seq_along(rna.bw), function(i) export.bw(rna.bw[[i]], paste0("test_", basename(rna.files[i]))))

# Annotation

gtf <- import.gff2("gencode.v34lift37.annotation.gtf.gz")
gtf <- gtf[gtf$gene_name == "CD55"]
export.gff2(gtf, "CD55_gencode.v34lift37.annotation.gtf")
system("pigz CD55_gencode.v34lift37.annotation.gtf")
