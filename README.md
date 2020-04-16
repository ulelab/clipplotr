# CLIPplotR
Plotting multiple iCLIP tracks across genes/regions with a range of normalisation and smoothing options.
Run from command line eg.:
`Rscript --vanilla CLIPplotR.R --xlinks clip1.bedgraph clip2.bedgraph --gtf genome.gtf --region ENSMUSG00000037400 --output finalgraph.pdf`
**To get all the options with explanation, run:**
`Rscript --vanilla CLIPplotR.R --help`
**Minimum input:** crosslinks (bedgraphs), genome annotation (gtf), region of interest (as chr3:35754106:35856276:+ or gene as ENSMUSG00000037400 or Atp11b), output file name.
**Important parameters:** 
--normalisation Normalisation options: none, maxpeak, libsize, default = "libsize"
--smoothing Smoothing options: none, rollmean, gaussian, default = "rollmean"
--smoothing_window, default = 100
**Output:** a graph saved in the file format you give as the “output file name” eg. graph.pdf will save as a pdf.


**Still to do:**
* Make it so you can input the graph output size, at the moment if you want to change that you need to edit the last line of the script.
* Add Rupert’s code so that you can show signal on mature transcripts (ie. without introns).


