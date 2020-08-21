# _clipplotr_

## Table of contents

1. [About _clipplotr_](#about-clipplotr)
2. [Installation](#installation)
3. [Quickstart](#quickstart)
4. [Slowstart](#slowstart)
5. [Example](#example)
5. [Test data](#test-data)

## About _clipplotr_

_clipplotr_ is a self-contained command-line tool written in R to facilitate comparative visualisation of CLIP data. It plots multiple CLIP tracks across a gene or region of interest with a range of normalisation and smoothing options. It can also optionally plot:

1. A peak track showing CLIP peaks. (This can also be used with the appropriate input to show other interesting features, e.g. _Alus_)
2. A coverage track which plots orthogonal data, e.g. coverage for the region of interest from RNAseq or Quantseq experiments.
3. An annotation track at either the gene or transcript level.

## Installation

_clipplotr_ requires R to be installed on your system. It has been tested with version 3.6.2.

To install _clipplotr_, either clone the repository with
```
git clone ulelab/clipplotr
```
for the latest version, or download from the [releases](https://github.com/ulelab/clipplotr/releases) page, which may be missing some of the latest features.

_clipplotr_ uses some R (`optparse`, `BiocManager`, `ggplot2`, `ggthemes`, `cowplot`, `patchwork`, `smoother`, `zoo`, `data.table`) and Bioconductor packages (`rtracklayer`, `GenomicFeatures`). If they are not already installed, run the helper script:
````
Rscript install_libraries.R
````

The `clipplotr` file may need to be made executable on your system depending on the installation method. If you have permission to do so, this can be done with:
```
chmod +x clipplotr
```

## Quickstart

The minimum parameters required are:

1. A space-separated list of CLIP tracks in either iCount BEDGRAPH or BED formats
2. The annotation GTF file
3. The region of interest
4. An output filename for the plot

All input files can be gzip compressed.

This can be run with a command such as:

```
./clipplotr \
--xlinks 'clip1.bedgraph clip2.bedgraph' \
--gtf genome.gtf \
--region 'chr1:1000:2000:+' \
--output plot.pdf
```

## Slowstart

To get all the parameters with explanations, run:

```
./clipplotr --help
```

There is a lot of customisation that can be done to make the desired plot. These are divided into:

1. CLIP plot
2. Peak plot (optional)
3. Coverage plot (optional)
4. Annotation plot (optional)
5. General

All plots are strand aware and only plot signal or features on the same strand as the region of interest. Additional optional plots are added dynamically to the CLIP plot.

Where multiple files are specified for a parameter, these should be space-separated and in `' '`. BED, BEDGRAPH and GTF files can be gzip compressed.

### 1. CLIP plot

* `-x` or `--xlinks` is used to supply the CLIP tracks. These are either in iCount bedgraph format (i.e. a 4-column BED file with the a positive score indicating the positive strand and a negative score the negative strand) or a standard 6-column BED file. In either case the score indicates the number of crosslinks at a given position.

* `-l` or `--labels` can be used to supply the unique names for the CLIP tracks and the order should match `--xlinks`. If not provided, the first 10 characters of the CLIP filename is used instead.

* `-c` or `--colours` can be used to supply the colours in hexadecimal format for the CLIP tracks and the order should match `--xlinks`. If not provided, a default set of colours are automatically generated up to a maximum of 10 tracks.

* `--groups` can be used to supply the group to which each CLIP track belongs and can be used to facet the CLIP plot. There should be an entry for each track and the order should match `--xlinks`. If not provided CLIP tracks are all plotted together.

* `-n` or `--normalisation` can be used to specify how the CLIP tracks should be normalised:

    1. `libsize` - by experiment library size and scaled to crosslinks per million (default).
    2. `maxpeak` - by the maximum peak (after smoothing) within the region of interest for the experiment.
    3. `none` - no normalisation, just plot raw crosslink counts.
    4. `custom` - by dividing counts with the provided size factors.

* `--size_factors` can be used to specify the size factors for each CLIP track for `custom` normalisation. The score will be divided by the size factor.

* `--scale_y` can be used to scale the y-axis for different groups.

* `-s` or `--smoothing` can be use to specify how the CLIP track should be smoothed:

    1. `rollmean` - by a rolling mean (default).
    2. `gaussian` - by a Gaussian smooth (this may require a lot of resources for a large region).
    3. `none` - no smoothing.

* `w` or `--smoothing_window` can be used to specify the size of the smoothing window in nucleotides (default: 100 nt)

* `--highlight` can be used to specify a region in the format `'start_coordinate:end_coordinate'`. This will be highlighted by grey shading. 

### 2. Peak plot

* `-p` or `--peaks` can be used to supply the peak tracks. These should be in a standard 6-column BED format. Optionally, if a 9-column BED file is supplied then the 9th column `itemRgb` will be used to colour the peaks. These tracks do not necessarily need to be peak intervals, but could be any features of interest. The names are the first 10 characters of the filename.

### 3. Coverage plot

* `--coverage` can be used to supply coverage tracks (e.g. RNA-seq or Quantseq). These should be supplied as BIGWIGs (as these files are not strand aware, ensure the BIGWIG for the correct strand as the region of interest is supplied). If multiple tracks are supplied, each one is plotted separately. Colours for the tracks are automatically generated up to a maximum of 8. The names are the first 10 characters of the filename.

### 4. Annotation plot

* `-g` or `--gtf` should be used to supply the reference GTF file. GENCODE files have been tested. The first time this GTF file is passed to _clipplotr_ it will generate and save an SQL TxDb database in the same location as the GTF file. This will be used for all future runs with the same GTF file.

* `-r` or `--region` should be used to specify the region of interest as:

    1. Colon-separated coordinates: `'chromosome:start_coordinate:end_coordinate:strand'`
    2. The GENCODE/ENSEMBL gene id: `ENSMUSG00000037400`
    3. The gene name: `Atp11b`

* `-a` or `annotation` specifies the type of annotation plot:

    1. `transcript` - all transcripts in the region are plotted and coloured by gene (default)
    2. `gene` - collapsed meta-genes (each containing all exons of a gene) in the region are plotted and coloured by gene
    3. `none` - annotation is not plotted
    4. `original` - plots the original _clipplotr_ annotation using `ggbio` (will be deprecated due to some bugs)

### 5. General

* `--size_x` can be used to specify the width of the final plot in mm (default: 210 A4)

* `--size_y` can be used to specify the height of the final plot in mm (default: 297 A4)

* `-o` or `--output` should be used to specify the output filename. The extension (e.g. `.pdf` or `.png`) will determine the output file type.

## Example

This is an example which shows many of the features of _clipplotr_ in action.

Here, I have reproduced part of Figure 1C from [Zarnack et al. (2013)](https://doi.org/10.1016/j.cell.2012.12.023) largely using publicly available pre-processed data (the BIGWIGs for the RNA-seq had to be generated from the raw files).

```
./clipplotr \
-x 'hnRNPC_iCLIP_rep1_LUjh03_all_xlink_events.bedgraph.gz hnRNPC_iCLIP_rep2_LUjh25_all_xlink_events.bedgraph.gz U2AF65_iCLIP_ctrl_rep1_all_xlink_events.bedgraph.gz U2AF65_iCLIP_ctrl_rep2_all_xlink_events.bedgraph.gz U2AF65_iCLIP_KD1_rep2_all_xlink_events.bedgraph.gz U2AF65_iCLIP_KD2_rep1_all_xlink_events.bedgraph.gz' \
-l 'hnRNPC_1 hnRNPC_2 U2AF65_WT_1 U2AF65_WT_2 U2AF65_KD_1 U2AF65_KD_2' \
-c '#586BA4 #324376 #0AA398 #067E79 #A54D69 #771434' \
--groups 'hnRNPC hnRNPC U2AF65_WT U2AF65_WT U2AF65_KD U2AF65_KD' \
-n libsize \
-s rollmean \
-w 50 \
-p 'Alu_rev.bed.gz' \
--coverage 'CTRL_plus.bigwig KD1_plus.bigwig KD2_plus.bigwig' \
-g gencode.v34lift37.annotation.gtf.gz \
-r 'chr1:207513000:207515000:+' \
--highlight '207513650:207513800' \
-a transcript \
-o CD55_C.pdf
```

produces the figure:

![Figure 1C from Zarnack et al. (2013)](https://github.com/ulelab/clipplotr/blob/dev/figures/CD55_C.png)

The plot is given the title of the region of interest.

In the CLIP plot we can see in the highlighted region a peak of hnRNPC binding (hnRNPC). There is no overlapping U2AF65 binding when hnRNPC is present (U2AF65_WT), but upon knockdown of hnRNPC, U2AF65 is able to bind to this region as strongly as hnRNPC (U2AF65_KD).

In the peak plot we can see that this binding site falls on the 5' end of reverse orientation _Alu_ element (Alu_rev).

In the coverage track we can see from RNA-seq read coverage that there is little expression of this exon when hnRNPC is present (CTRL_plus), but upon two biological replicates of hnRNPC knockdown there is a large increase in expression (KD1_plus and KD2_plus).

In the annotation track we can see this is contained within in the CD55 gene, which has the ENSEMBL ID ENSG000001962352 with the GENCODE suffix 16_8. Although there are many transcripts where this exon is not expressed, there are two annotated ones where it is in the latest GENCODE V34 annotation. There are no other genes in this region of interest.

## Test data

To try out _clipplotr_ a small test dataset has been created based on the example plot above.

```
cd test

../clipplotr \
-x 'test_hnRNPC_iCLIP_rep1_LUjh03_all_xlink_events.bedgraph.gz test_hnRNPC_iCLIP_rep2_LUjh25_all_xlink_events.bedgraph.gz test_U2AF65_iCLIP_ctrl_rep1_all_xlink_events.bedgraph.gz test_U2AF65_iCLIP_ctrl_rep2_all_xlink_events.bedgraph.gz test_U2AF65_iCLIP_KD1_rep2_all_xlink_events.bedgraph.gz test_U2AF65_iCLIP_KD2_rep1_all_xlink_events.bedgraph.gz' \
-l 'hnRNPC_1 hnRNPC_2 U2AF65_WT_1 U2AF65_WT_2 U2AF65_KD_1 U2AF65_KD_2' \
-c '#586BA4 #324376 #0AA398 #067E79 #A54D69 #771434' \
--groups 'hnRNPC hnRNPC U2AF65_WT U2AF65_WT U2AF65_KD U2AF65_KD' \
-n custom \
--size_factors '4.869687 9.488133 1.781117 10.135903 4.384385 8.227587' \
-s rollmean \
-w 50 \
-p 'test_Alu_rev.bed.gz' \
--coverage 'test_CTRL_plus.bigwig test_KD1_plus.bigwig test_KD2_plus.bigwig' \
-g CD55_gencode.v34lift37.annotation.gtf.gz \
-r 'chr1:207513000:207515000:+' \
--highlight '207513650:207513800' \
-a transcript \
-o test.pdf
```

Note as we cannot calculate library sizes from the small subset dataset, here I have used custom normalisation with pre-calculated library size factors from the full CLIP bedgraph files.
