# Script to plot gene structure given a GTF
# A. M. Chakrabarti
# Last updated: 28th June 2019

opt <- list(xlinks = "iclip.bedgraph.gz",
	gtf = ".gtf.gz",
	region = "1:123:123:+",
	gene = "ENSG",
	smoothing = "gaussian",
	normalisation = "libSize",
	output = "plot.pdf")