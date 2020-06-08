# CLIPplotR
Plotting multiple iCLIP tracks across genes/regions with a range of normalisation and smoothing options with peak and annotation tracks.

Run from command line e.g.:
`Rscript --vanilla CLIPplotR.R --xlinks "clip1.bedgraph clip2.bedgraph" --gtf genome.gtf --region ENSMUSG00000037400 --output finalgraph.pdf`

**To get all the options with explanation, run:**
`Rscript --vanilla CLIPplotR.R --help`

**Minimum input:** crosslinks (bedgraphs), genome annotation (gtf), region of interest (as chr3:35754106:35856276:+ or gene as ENSMUSG00000037400 or Atp11b), output file name.

**Important parameters:** 
--normalisation Normalisation options: none, maxpeak, libsize, default = "libsize"
--smoothing Smoothing options: none, rollmean, gaussian, default = "rollmean"
--smoothing_window, default = 100

**Output:** a graph saved in the file format you give as the “output file name” eg. graph.pdf will save as a pdf.

**'Argh this isn't working' checklist**
- Are all the file paths correct relative to where you are running the script?
- Do you have any rogue extra spaces?
- Are the crosslinks files in speech marks?
