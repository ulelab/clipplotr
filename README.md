# CLIPplotR

## About CLIPplotR

CLIPplotR is a self-contained command-line tool written in R to facilitate comparative visualisation of CLIP data. It plots multiple CLIP tracks across a gene or region of interest with a range of normalisation and smoothing options. It can also optionally plot:

1. A peak track showing CLIP peaks. (This can also be used with the appropriate input to show other interesting features, e.g. _Alus_)
2. A coverage track which plots orthogonal data, e.g. coverage for the region of interest from RNAseq or Quantseq experiments.
3. An annotation track at either the gene or transcript level.

## Installation

CLIPplotR requires R to be installed on your system. It has been tested with version 3.6.2.

To install CLIPplotR, either clone the repository with `git clone ulelab/clipplotr` for the latest version, or download from the [releases](https://github.com/ulelab/clipplotr/releases) page, which may be missing some of the latest features.

CLIPplotR uses some R and Bioconductor packages. If they are not already installed, running the helper script `Rscript install_libraries.R` will do so.

The `CLIPplotR.R` file may need to be made executable on your system depending on the installation method. If you have permission to do so, this can be done by `chmod +x CLIPplotR.R`.

## Quickstart

To get all the options with explanations, run:

```
./CLIPplotR --help
```

The minimum parametets required are:

1. A space-separated list of CLIP tracks in either iCount BEDGRAPH or BED formats
2. The annotation GTF file
3. The region of interest
4. An output filename for the plot

All input files can be gzip compressed.

This can be run with a command such as:

```
./CLIPplotR \
--xlinks "clip1.bedgraph clip2.bedgraph" \
--gtf genome.gtf \
--region ENSMUSG00000037400 \
--output plot.pdf
```

## Slowstart

There is a lot of customisation that can be done to make the desired plot. These are divided into:

1. CLIP plot
2. Peak plot (optional)
3. Coverage plot (optional)
4. Annotation plot (optional)
5. General

All plots are strand aware and only plot signal or features on the same strand as the region of interest. Additional optional plots are added dynamically to the CLIP plot.

Where multiple files are specified for a parameter, these should be space-separated and in `" "`. BED, BEDGRAPH and GTF files can be gzip compressed.

### CLIP plot

* `-x` or `--xlinks` is used to supply the CLIP tracks. These are either in iCount bedgraph format (i.e. a 4-column BED file with the a positive score indicating the positive strand and a negative score the negative strand) or a standard 6-column BED file. In either case the score indicates the number of crosslinks at a given position. 

* `-l` or `--labels` can be used to supply the names for the CLIP tracks and the order should match `--xlinks`. If not provided, the first 10 characters of the CLIP filename is used instead.

* `-c` or `--colours` can be used to supply the colours for the CLIP tracks and the order should match `--xlinks`. If not provided, a default set of colours are automatically generated up to a maximum of 10 tracks.

* `--groups` can be used to supply the group to which each CLIP track belongs and can be used to facet the CLIP plot. There should be an entry for each track and the order should match `--xlinks`. If not provided CLIP tracks are all plotted together.

* `-n` or `--normalisation` can be used to specify how the CLIP tracks should be normalised:

    1. `libsize` - by library size and scaled to crosslinks per million (default)
    2. `maxpeak` - by the maximum peak within the region of interest
    2. `none` - no normalisation, just plot raw crosslink counts.

* `-s` or `--smoothing` can be use to specify how the CLIP track should be smoothed:

    1. `rollmean` - by a rolling mean (default)
    2. `gaussian` - by a Gaussian smooth (this may require a lot of resources for a large region)
    2. `none` - no smoothing

* `w` or `--smoothing_window` can be used to specify the size of the smoothing window in nucleotides (default: 100 nt)

* `--highlight` can be used to specify a region in the format `"start_coordinate:end_coordinate"`. This will be highlighted by grey shading. 

### Peak plot

* `-p` or `--peaks` can be used to supply the peak tracks. These should be in a standard 6-column BED format. Optionally, if a 10-column BED file is supplied then the 9th column `itemRgb` will be used to colour the peaks. This does not necessarily need to be peak intervals, but could be any features of interest. The names are the first 10 characters of the filename.

### Coverage plot

* `--coverage` can be used to supply coverage tracks (e.g. RNA-seq or Quantseq). These should be supplied as BIGWIGs (as these files are not strand aware, ensure the BIGWIG for the correct strand as the region of interest is supplied). If multiple tracks are supplied, each one is plotted separately. Colours for the tracks are automatically generated up to a maximum of 8. The names are the first 10 characters of the filename.

### Annotation plot

* `-g` or `--gtf` should be used to supply the reference GTF file. GENCODE files have been tested.

* `-r` or `--region` should be used to specify the region of interest as:

    1. Colon-separated coordinates: `"chromosome:start_coordinate:end_coordinate:strand"`
    2. The GENCODE/ENSEMBL gene id: `ENSMUSG00000037400`
    3. The gene name: `Atp11b`

* `-a` or `annotation` specifies the type of annotation plot:

    1. `transcript` - all transcripts in the region are plotted and coloured by gene (default)
    2. `gene` - collapsed meta-genes (each containing all exons of a gene) in the region are plotted and coloured by gene
    3. `none` - annotation is not plotted
    4. `original` - plots the original CLIPplotR annotation using `ggbio` (will be deprecated due to some bugs)

### General

* `--size_x` can be used to specify the width of the final plot in mm (default: 210 A4)

* `--size_y` can be used to specify the height of the final plot in mm (default: 297 A4)

* `-o` or `--output` should be used to specify the output filename. The extension (e.g. `.pdf` or `.png`) will determine the output file type.