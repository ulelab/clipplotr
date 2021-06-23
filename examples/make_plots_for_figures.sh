#!/bin/bash

# Script to make plots for figures
# A. M. Chakrabarti
# 13th April 2021

conda activate clipplotr

# Specify local path to clipplotr
CLIPPLOTR=~/github/clipplotr/clipplotr

# ==========
# Zarnack - Figure 1
# ==========

$CLIPPLOTR \
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
-o '../../plots/CD55.pdf'

# ==========
# Zarnack - Figure 2
# ==========

$CLIPPLOTR \
-x 'hnRNPC_iCLIP_rep1_LUjh03_all_xlink_events.bedgraph.gz,hnRNPC_iCLIP_rep2_LUjh25_all_xlink_events.bedgraph.gz,U2AF65_iCLIP_ctrl_rep1_all_xlink_events.bedgraph.gz,U2AF65_iCLIP_ctrl_rep2_all_xlink_events.bedgraph.gz,U2AF65_iCLIP_KD1_rep2_all_xlink_events.bedgraph.gz,U2AF65_iCLIP_KD1_rep1_all_xlink_events.bedgraph.gz' \
-l 'hnRNPC 1,hnRNPC 2,U2AF2 WT 1,U2AF2 WT 2,U2AF2 KD 1,U2AF2 KD 2' \
-c '#586BA4,#324376,#0AA398,#067E79,#A54D69,#771434' \
--groups 'hnRNPC,hnRNPC,U2AF2 WT,U2AF2 WT,U2AF2 KD,U2AF2 KD' \
-n libsize \
-s none \
-g gencode.v34lift37.annotation.gtf.gz \
-r 'chr1:207513000:207515000:+' \
--highlight '207513650:207513800' \
-a none \
--ratios '1,0,0,0' \
--size_x 180 \
--size_y 100 \
-o '../../plots/CD55_unsmoothed.pdf'

for i in 50 100 200; do

    $CLIPPLOTR \
    -x 'hnRNPC_iCLIP_rep1_LUjh03_all_xlink_events.bedgraph.gz,hnRNPC_iCLIP_rep2_LUjh25_all_xlink_events.bedgraph.gz,U2AF65_iCLIP_ctrl_rep1_all_xlink_events.bedgraph.gz,U2AF65_iCLIP_ctrl_rep2_all_xlink_events.bedgraph.gz,U2AF65_iCLIP_KD1_rep2_all_xlink_events.bedgraph.gz,U2AF65_iCLIP_KD1_rep1_all_xlink_events.bedgraph.gz' \
    -l 'hnRNPC 1,hnRNPC 2,U2AF2 WT 1,U2AF2 WT 2,U2AF2 KD 1,U2AF2 KD 2' \
    -c '#586BA4,#324376,#0AA398,#067E79,#A54D69,#771434' \
    --groups 'hnRNPC,hnRNPC,U2AF2 WT,U2AF2 WT,U2AF2 KD,U2AF2 KD' \
    -n libsize \
    -s rollmean \
    -w $i \
    -g gencode.v34lift37.annotation.gtf.gz \
    -r 'chr1:207513000:207515000:+' \
    --highlight '207513650:207513800' \
    -a none \
    --ratios '1,0,0,0' \
    --size_x 180 \
    --size_y 100 \
    -o "../../plots/CD55_rollmean_$i.pdf"

done

# ==========
# Nostrand - Figure 3
# ==========

$CLIPPLOTR \
-x 'ENCFF239CML.xl.bed.gz,ENCFF170YQV.xl.bed.gz,ENCFF515BTB.xl.bed.gz,ENCFF537RYR.xl.bed.gz,ENCFF296GDR.xl.bed.gz,ENCFF212IIR.xl.bed.gz' \
-l 'HepG2 IP1,HepG2 IP2,HepG2 SMI,K562 IP1,K562_IP2,K562 SMI' \
-c '#324376,#586BA4,#B4BCD6,#771434,#A54D69,#D2A6B4' \
--groups 'HepG2,HepG2,HepG2,K562,K562,K562' \
-n libsize \
-y 'HepG2.bed.gz,K562.bed.gz' \
--auxiliary_labels 'HepG2 peaks,K562 peaks' \
-g gencode.v34.annotation.gtf.gz \
-r 'chr17:8458500:8469500:+' \
--highlight '8462500:8467500' \
-a gene \
--ratio '2,0.5,0,0.25' \
--size_x 200 \
--size_y 200 \
-o '../../plots/NDEL1_libsize_gene.pdf'

$CLIPPLOTR \
-x 'ENCFF239CML.xl.bed.gz,ENCFF170YQV.xl.bed.gz,ENCFF515BTB.xl.bed.gz,ENCFF537RYR.xl.bed.gz,ENCFF296GDR.xl.bed.gz,ENCFF212IIR.xl.bed.gz' \
-l 'HepG2 IP1,HepG2 IP2,HepG2 SMI,K562 IP1,K562_IP2,K562 SMI' \
-c '#324376,#586BA4,#B4BCD6,#771434,#A54D69,#D2A6B4' \
--groups 'HepG2,HepG2,HepG2,K562,K562,K562' \
-n libsize \
--scale_y \
-g gencode.v34.annotation.gtf.gz \
-r 'chr17:8462500:8467500:+' \
-a none \
--ratio '2,0.5,0,0.25' \
--size_x 200 \
--size_y 100 \
--tidy_y_labels 4 \
-o '../../plots/NDEL1_libsize_scaled.pdf'

$CLIPPLOTR \
-x 'ENCFF239CML.xl.bed.gz,ENCFF170YQV.xl.bed.gz,ENCFF515BTB.xl.bed.gz,ENCFF537RYR.xl.bed.gz,ENCFF296GDR.xl.bed.gz,ENCFF212IIR.xl.bed.gz' \
-l 'HepG2 IP1,HepG2 IP2,HepG2 SMI,K562 IP1,K562_IP2,K562 SMI' \
-c '#324376,#586BA4,#B4BCD6,#771434,#A54D69,#D2A6B4' \
--groups 'HepG2,HepG2,HepG2,K562,K562,K562' \
-n maxpeak \
-g gencode.v34.annotation.gtf.gz \
-r 'chr17:8462500:8467500:+' \
-a none \
--ratio '2,0.5,0,0.25' \
--size_x 198 \
--size_y 100 \
--tidy_y_labels 3 \
-o '../../plots/NDEL1_maxpeak.pdf'

$CLIPPLOTR \
-x 'ENCFF239CML.xl.bed.gz,ENCFF170YQV.xl.bed.gz,ENCFF515BTB.xl.bed.gz,ENCFF537RYR.xl.bed.gz,ENCFF296GDR.xl.bed.gz,ENCFF212IIR.xl.bed.gz' \
-l 'HepG2 IP1,HepG2 IP2,HepG2 SMI,K562 IP1,K562_IP2,K562 SMI' \
-c '#324376,#586BA4,#B4BCD6,#771434,#A54D69,#D2A6B4' \
--groups 'HepG2,HepG2,HepG2,K562,K562,K562' \
-n libsize_maxpeak \
-g gencode.v34.annotation.gtf.gz \
-r 'chr17:8462500:8467500:+' \
-a none \
--ratio '2,0.5,0,0.25' \
--size_x 207 \
--size_y 100 \
--tidy_y_labels 3 \
-o '../../plots/NDEL1_libsize_maxpeak.pdf'