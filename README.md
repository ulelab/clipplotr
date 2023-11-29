# _clipplotr_

#### 

Chakrabarti, A. M., Capitanchik, C., Ule, J., & Luscombe, N. M. (2023) [clipplotr - a comparative visualisation and analysis tool for CLIP data](https://doi.org/10.1261/rna.079326.122) 
*RNA* 29 (6): 715-723

doi: 10.1261/rna.079326.122

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

To install _clipplotr_, either clone the repository with
```
git clone https://github.com/ulelab/clipplotr.git
```
for the latest version, or download from the [releases](https://github.com/ulelab/clipplotr/releases) page, which may be missing some of the latest features.

There are two options for installing the dependencies.

### 1. Conda option

If you have Conda on your system you can create a virtual environment which installs R and all the dependencies using the provided YAML. First move into the directory into which you clone _clipplotr_ and then run:
```
cd clipplotr
conda env create -f environment.yml
```

You can then activate the enviroment using:
```
conda activate clipplotr
```

### 2. R option

_clipplotr_ requires R to be installed on your system and uses some R (`optparse`, `BiocManager`, `ggplot2`, `ggthemes`, `cowplot`, `patchwork`, `smoother`, `zoo`, `data.table`) and Bioconductor packages (`rtracklayer`, `GenomicFeatures`). If you have R already installed, you can use the helper script to install the packages if needed. First move into the directory into which you clone _clipplotr_ and then run

````
cd clipplotr
Rscript install_libraries.R
````

With either option the `clipplotr` file may need to be made executable on your system depending on the installation method. If you have permission to do so, this can be done with:
```
chmod u+x clipplotr
```

## Quickstart

The minimum parameters required are:

1. A comma-separated list of CLIP tracks in either iCount BEDGRAPH or BED formats
2. The annotation GTF file
3. The region of interest
4. An output filename for the plot

All input files can be gzip compressed.

This can be run with a command such as:

```
./clipplotr \
--xlinks 'clip1.bedgraph,clip2.bedgraph' \
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
2. Auxiliary plot (ie. any annotation in bed format, for example peaks; optional)
3. Coverage plot (optional)
4. Annotation plot (optional)
5. General

All plots are strand aware and only plot signal or features on the same strand as the region of interest. Additional optional plots are added dynamically to the CLIP plot.

Where multiple files are specified for a parameter, these should be comma-separated and in quotation marks `' '`. BED, BEDGRAPH and GTF files can be gzip compressed.

### 1. CLIP plot

* `-x` or `--xlinks` is used to supply the CLIP tracks. These are either in iCount bedgraph format (i.e. a 4-column BED file with the a positive score indicating the positive strand and a negative score the negative strand) or a standard 6-column BED file. In either case the score indicates the number of crosslinks at a given position.

* `-l` or `--labels` can be used to supply the unique names for the CLIP tracks and the order should match `--xlinks`. If not provided, the first 10 characters of the CLIP filename is used instead.

* `-c` or `--colours` can be used to supply the colours in hexadecimal format for the CLIP tracks and the order should match `--xlinks`. If not provided, a default set of colours are automatically generated up to a maximum of 10 tracks.

* `--groups` can be used to supply the group to which each CLIP track belongs and can be used to facet the CLIP plot. There should be an entry for each track and the order should match `--xlinks`. If not provided CLIP tracks are all plotted together.

* `-n` or `--normalisation` can be used to specify how the CLIP tracks should be normalised:

    1. `libsize` - by experiment library size and scaled to crosslinks per million (default).
    2. `maxpeak` - by the maximum peak (after smoothing) within the region of interest for the experiment. If multiple CLIP tracks are plotted, the maxmimum position in the group is taken.
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

* `--tidy_y_labels` can be use to tidy the y-axis labels and keep alternate ones if there are 5 or more.

### 2. Auxiliary plot

* `-y` or `--auxiliary` can be used to supply the auxiliary tracks. These should be in a standard 6-column BED format. Optionally, if a 9-column BED file is supplied then the 9th column `itemRgb` will be used to colour the intervals. These tracks could be peak intervals, but could be any features of interest. You can label the tracks with `--auxiliary_labels`, otherwise the first 10 characters of the filename are used.

### 3. Coverage plot

* `--coverage` can be used to supply coverage tracks (e.g. RNA-seq or Quantseq). These should be supplied as BIGWIGs (as these files are not strand aware, ensure the BIGWIG for the correct strand as the region of interest is supplied). If multiple tracks are supplied, each one is plotted separately. Colours for the tracks are automatically generated up to a maximum of 8. The names are the first 10 characters of the filename.

* `--coverage_labels` can be used to supply the unique names for the coverage tracks and the order should match `--coverage`.

* `--coverage_colours` can be used to supply the colours in hexadecimal format for the coverage tracks and the order should match `--coverage`. If not provided, a default set of colours are automatically generated up to a maximum of 10 tracks.

* `--coverage_groups` can be used to supply the group to which each coverage track belongs and can be used to facet the coverage plot. There should be an entry for each track and the order should match `--xlinks`. If not provided coverage tracks are all plotted together.

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

* `--ratios` allows you to specify the relative sizing of the combined plots. Specify plot ratios in order: xlink track, auxiliary tracks, coverage track, annotation track (comma separated). Put 0 if any of these track types are absent. (default: 2 for xlinks, 0.25 for 1 auxiliary track 0.5 for >1, 2 for coverage, 3 for annotation)

## Example

This is an example which shows many of the features of _clipplotr_ in action.

Here, I have reproduced part of Figure 1C from [Zarnack et al. (2013)](https://doi.org/10.1016/j.cell.2012.12.023) largely using publicly available pre-processed data (the BIGWIGs for the RNA-seq had to be generated from the raw files).

```
./clipplotr \
-x 'hnRNPC_iCLIP_rep1_LUjh03_all_xlink_events.bedgraph.gz,hnRNPC_iCLIP_rep2_LUjh25_all_xlink_events.bedgraph.gz,U2AF65_iCLIP_ctrl_rep1_all_xlink_events.bedgraph.gz,U2AF65_iCLIP_ctrl_rep2_all_xlink_events.bedgraph.gz,U2AF65_iCLIP_KD1_rep2_all_xlink_events.bedgraph.gz,U2AF65_iCLIP_KD2_rep1_all_xlink_events.bedgraph.gz' \
-l 'hnRNPC 1,hnRNPC 2,U2AF2 WT 1,U2AF2 WT 2,U2AF2 KD 1,U2AF2 KD 2' \
-c '#586BA4,#324376,#0AA398,#067E79,#A54D69,#771434' \
--groups 'hnRNPC,hnRNPC,U2AF2 WT,U2AF2 WT,U2AF2 KD,U2AF2 KD' \
-n libsize \
-s rollmean \
-w 50 \
-y 'Alu_rev.bed.gz' \
--auxiliary_labels 'reverse Alu' \
--coverage 'ERR127306_plus.bigwig,ERR127307_plus.bigwig,ERR127308_plus.bigwig,ERR127309_plus.bigwig,ERR127302_plus.bigwig,ERR127303_plus.bigwig,ERR127304_plus.bigwig,ERR127305_plus.bigwig' \
--coverage_labels 'CTRL1 1,CTRL1 2,CTRL2 1,CTRL2 2,KD1 1,KD1 2,KD2 1,KD2 2' \
--coverage_colours '#A1D99B,#74C476,#31A354,#006D2C,#FDAE6B,#E6550D,#FC9272,#DE2D26' \
--coverage_groups 'CTRL,CTRL,CTRL,CTRL,KD,KD,KD,KD' \
-g gencode.v34lift37.annotation.gtf.gz \
-r 'chr1:207513000:207515000:+' \
--highlight '207513650:207513800' \
-a transcript \
-o 'CD55.pdf'
```

produces the figure:

![Figure 1C from Zarnack et al. (2013)](https://github.com/ulelab/clipplotr/blob/dev/figures/CD55.png)

The plot is given the title of the region of interest.

In the CLIP plot we can see in the highlighted region a peak of hnRNPC binding (hnRNPC). There is no overlapping U2AF2 binding when hnRNPC is present (U2AF2 WT), but upon knockdown of hnRNPC, U2AF2 is able to bind to this region as strongly as hnRNPC (U2AF2 KD).

In the peak plot we can see that this binding site falls on the 5' end of reverse orientation _Alu_ element (reverse Alu).

In the coverage track we can see from RNA-seq read coverage that there is little expression of this exon when hnRNPC is present (CTRL), but upon two biological replicates of hnRNPC knockdown there is a large increase in expression (KD1 and KD2).

In the annotation track we can see this is contained within in the CD55 gene, which has the ENSEMBL ID ENSG000001962352 with the GENCODE suffix 16_8. Although there are many transcripts where this exon is not expressed, there are two annotated ones where it is in the latest GENCODE V34 annotation. There are no other genes in this region of interest.

## Test data

To try out _clipplotr_ a small test dataset has been created based on the example plot above.

```
cd test

../clipplotr \
-x 'test_hnRNPC_iCLIP_rep1_LUjh03_all_xlink_events.bedgraph.gz,test_hnRNPC_iCLIP_rep2_LUjh25_all_xlink_events.bedgraph.gz,test_U2AF65_iCLIP_ctrl_rep1_all_xlink_events.bedgraph.gz,test_U2AF65_iCLIP_ctrl_rep2_all_xlink_events.bedgraph.gz,test_U2AF65_iCLIP_KD1_rep2_all_xlink_events.bedgraph.gz,test_U2AF65_iCLIP_KD2_rep1_all_xlink_events.bedgraph.gz' \
-l 'hnRNPC 1,hnRNPC 2,U2AF2 WT 1,U2AF2 WT 2,U2AF2 KD 1,U2AF2 KD 2' \
-c '#586BA4,#324376,#0AA398,#067E79,#A54D69,#771434' \
--groups 'hnRNPC,hnRNPC,U2AF2 WT,U2AF2 WT,U2AF2 KD,U2AF2 KD' \
-n custom \
--size_factors '4.869687,9.488133,1.781117,10.135903,4.384385,8.227587' \
-s rollmean \
-w 50 \
-y 'test_Alu_rev.bed.gz' \
--auxiliary_labels 'reverse Alu' \
--coverage 'test_ERR127306_plus.bigwig,test_ERR127307_plus.bigwig,test_ERR127308_plus.bigwig,test_ERR127309_plus.bigwig,test_ERR127302_plus.bigwig,test_ERR127303_plus.bigwig,test_ERR127304_plus.bigwig,test_ERR127305_plus.bigwig' \
--coverage_labels 'CTRL1 1,CTRL1 2,CTRL2 1,CTRL2 2,KD1 1,KD1 2,KD2 1,KD2 2' \
--coverage_colours '#A1D99B,#74C476,#31A354,#006D2C,#FDAE6B,#E6550D,#FC9272,#DE2D26' \
--coverage_groups 'CTRL,CTRL,CTRL,CTRL,KD,KD,KD,KD' \
-g CD55_gencode.v34lift37.annotation.gtf.gz \
-r 'chr1:207513000:207515000:+' \
--highlight '207513650:207513800' \
-a transcript \
-o test.pdf
```

Note as we cannot calculate library sizes from the small subset dataset, here I have used custom normalisation with pre-calculated library size factors from the full CLIP bedgraph files.
