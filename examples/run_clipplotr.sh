#!/bin/bash

# Nostrand

cd ~/projects/clipplotr/nostrand

~/Github/clipplotr/CLIPplotR.R -x 'ENCFF239CML.xl.bed.gz ENCFF170YQV.xl.bed.gz ENCFF515BTB.xl.bed.gz ENCFF537RYR.xl.bed.gz ENCFF296GDR.xl.bed.gz ENCFF212IIR.xl.bed.gz' \
-l 'HepG2_IP1 HepG2_IP2 HepG2_SMI K562_IP1 K562_IP2 K562_SMI' \
-c '#324376 #586BA4 #B4BCD6 #771434 #A54D69 #D2A6B4' \
--groups 'HepG2 HepG2 HepG2 K562 K562 K562' \
-p 'HepG2.bed.gz K562.bed.gz' \
-g gencode.v34.annotation.gtf.gz \
-r 'chr17:8458500:8469500:+' \
-a transcript \
-o NDEL1_libsize.pdf

~/Github/clipplotr/CLIPplotR.R -x 'ENCFF239CML.xl.bed.gz ENCFF170YQV.xl.bed.gz ENCFF515BTB.xl.bed.gz ENCFF537RYR.xl.bed.gz ENCFF296GDR.xl.bed.gz ENCFF212IIR.xl.bed.gz' \
-l 'HepG2_IP1 HepG2_IP2 HepG2_SMI K562_IP1 K562_IP2 K562_SMI' \
--groups 'HepG2 HepG2 HepG2 K562 K562 K562' \
-p 'HepG2.bed.gz K562.bed.gz' \
-g gencode.v34.annotation.gtf.gz \
-r 'chr17:8458500:8469500:+' \
-a transcript \
-o NDEL1_libsize.png

~/Github/clipplotr/CLIPplotR.R -x 'ENCFF239CML.xl.bed.gz ENCFF170YQV.xl.bed.gz ENCFF515BTB.xl.bed.gz ENCFF537RYR.xl.bed.gz ENCFF296GDR.xl.bed.gz ENCFF212IIR.xl.bed.gz' \
-l 'HepG2_IP1 HepG2_IP2 HepG2_SMI K562_IP1 K562_IP2 K562_SMI' \
-c '#324376 #586BA4 #B4BCD6 #771434 #A54D69 #D2A6B4' \
--groups 'HepG2 HepG2 HepG2 K562 K562 K562' \
-p 'HepG2.bed.gz K562.bed.gz' \
-g gencode.v34.annotation.gtf.gz \
-r 'chr17:8458500:8469500:+' \
-n maxpeak \
-a transcript \
-o NDEL1_maxpeak.pdf

~/Github/clipplotr/CLIPplotR.R -x 'ENCFF239CML.xl.bed.gz ENCFF170YQV.xl.bed.gz ENCFF515BTB.xl.bed.gz ENCFF537RYR.xl.bed.gz ENCFF296GDR.xl.bed.gz ENCFF212IIR.xl.bed.gz' \
-l 'HepG2_IP1 HepG2_IP2 HepG2_SMI K562_IP1 K562_IP2 K562_SMI' \
-c '#324376 #586BA4 #B4BCD6 #771434 #A54D69 #D2A6B4' \
--groups 'HepG2 HepG2 HepG2 K562 K562 K562' \
-p 'HepG2.bed.gz K562.bed.gz' \
-g gencode.v34.annotation.gtf.gz \
-r 'chr17:8458500:8469500:+' \
-n maxpeak \
-a transcript \
-o NDEL1_maxpeak.png

# Zarnack

cd ~/projects/clipplotr/zarnack

~/Github/clipplotr/CLIPplotR.R -x 'U2AF65_iCLIP_KD1_rep2_all_xlink_events.bedgraph.gz U2AF65_iCLIP_KD2_rep1_all_xlink_events.bedgraph.gz U2AF65_iCLIP_ctrl_rep1_all_xlink_events.bedgraph.gz U2AF65_iCLIP_ctrl_rep2_all_xlink_events.bedgraph.gz hnRNPC_iCLIP_rep1_LUjh03_all_xlink_events.bedgraph.gz hnRNPC_iCLIP_rep2_LUjh25_all_xlink_events.bedgraph.gz' \
-l 'U2AF65_KD_2 U2AF65_KD_1 U2AF65_WT_1 U2AF65_WT_2 hnRNPC_1 hnRNPC_2' \
--groups 'U2AF65_KD U2AF65_KD U2AF65_WT U2AF65_WT hnRNPC hnRNPC' \
-p 'hg19.Alu.reversed.bed.gz' \
--coverage 'CTRL_plus.bigwig KD1_plus.bigwig KD2_plus.bigwig' \
-g gencode.v34lift37.annotation.gtf.gz \
-r 'CD55' \
-a transcript \
-o CD55.pdf

~/Github/clipplotr/CLIPplotR.R -x 'U2AF65_iCLIP_KD1_rep2_all_xlink_events.bedgraph.gz U2AF65_iCLIP_KD2_rep1_all_xlink_events.bedgraph.gz U2AF65_iCLIP_ctrl_rep1_all_xlink_events.bedgraph.gz U2AF65_iCLIP_ctrl_rep2_all_xlink_events.bedgraph.gz hnRNPC_iCLIP_rep1_LUjh03_all_xlink_events.bedgraph.gz hnRNPC_iCLIP_rep2_LUjh25_all_xlink_events.bedgraph.gz' \
-l 'U2AF65_KD_2 U2AF65_KD_1 U2AF65_WT_1 U2AF65_WT_2 hnRNPC_1 hnRNPC_2' \
--groups 'U2AF65_KD U2AF65_KD U2AF65_WT U2AF65_WT hnRNPC hnRNPC' \
-p 'hg19.Alu.reversed.bed.gz' \
--coverage 'CTRL_plus.bigwig KD1_plus.bigwig KD2_plus.bigwig' \
-g gencode.v34lift37.annotation.gtf.gz \
-r 'chr1:207494500:207500000:+' \
-a transcript \
-w 50 \
-o CD55_B.pdf

~/Github/clipplotr/CLIPplotR.R \
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

~/Github/clipplotr/CLIPplotR.R \
-x 'hnRNPC_iCLIP_rep1_LUjh03_all_xlink_events.bedgraph.gz hnRNPC_iCLIP_rep2_LUjh25_all_xlink_events.bedgraph.gz U2AF65_iCLIP_ctrl_rep1_all_xlink_events.bedgraph.gz U2AF65_iCLIP_ctrl_rep2_all_xlink_events.bedgraph.gz U2AF65_iCLIP_KD1_rep2_all_xlink_events.bedgraph.gz U2AF65_iCLIP_KD2_rep1_all_xlink_events.bedgraph.gz' \
-l 'hnRNPC_1 hnRNPC_2 U2AF65_WT_1 U2AF65_WT_2 U2AF65_KD_1 U2AF65_KD_2' \
-c '#586BA4 #324376 #0AA398 #067E79 #A54D69 #771434' \
--groups 'hnRNPC hnRNPC U2AF65_WT U2AF65_WT U2AF65_KD U2AF65_KD' \
-n none \
-s none \
-g gencode.v34lift37.annotation.gtf.gz \
-r 'chr1:207513000:207515000:+' \
--highlight '207513650:207513800' \
-a none \
-o CD55_C_raw.pdf

for i in 5 10 15 20 25 30 40 50 100 200 500; do 

~/Github/clipplotr/CLIPplotR.R \
-x 'hnRNPC_iCLIP_rep1_LUjh03_all_xlink_events.bedgraph.gz hnRNPC_iCLIP_rep2_LUjh25_all_xlink_events.bedgraph.gz U2AF65_iCLIP_ctrl_rep1_all_xlink_events.bedgraph.gz U2AF65_iCLIP_ctrl_rep2_all_xlink_events.bedgraph.gz U2AF65_iCLIP_KD1_rep2_all_xlink_events.bedgraph.gz U2AF65_iCLIP_KD2_rep1_all_xlink_events.bedgraph.gz' \
-l 'hnRNPC_1 hnRNPC_2 U2AF65_WT_1 U2AF65_WT_2 U2AF65_KD_1 U2AF65_KD_2' \
-c '#586BA4 #324376 #0AA398 #067E79 #A54D69 #771434' \
--groups 'hnRNPC hnRNPC U2AF65_WT U2AF65_WT U2AF65_KD U2AF65_KD' \
-n none \
-s rollmean \
-w $i \
-g gencode.v34lift37.annotation.gtf.gz \
-r 'chr1:207513000:207515000:+' \
--highlight '207513650:207513800' \
-a none \
-o CD55_C_raw_rollmean_$i.pdf;

done

~/Github/clipplotr/CLIPplotR.R \
-x 'hnRNPC_iCLIP_rep1_LUjh03_all_xlink_events.bedgraph.gz hnRNPC_iCLIP_rep2_LUjh25_all_xlink_events.bedgraph.gz U2AF65_iCLIP_ctrl_rep1_all_xlink_events.bedgraph.gz U2AF65_iCLIP_ctrl_rep2_all_xlink_events.bedgraph.gz U2AF65_iCLIP_KD1_rep2_all_xlink_events.bedgraph.gz U2AF65_iCLIP_KD2_rep1_all_xlink_events.bedgraph.gz' \
-l 'hnRNPC_1 hnRNPC_2 U2AF65_WT_1 U2AF65_WT_2 U2AF65_KD_1 U2AF65_KD_2' \
-c '#586BA4 #324376 #0AA398 #067E79 #A54D69 #771434' \
--groups 'hnRNPC hnRNPC U2AF65_WT U2AF65_WT U2AF65_KD U2AF65_KD' \
-n libsize \
-s rollmean \
-w 50 \
-g gencode.v34lift37.annotation.gtf.gz \
-r 'chr1:207513000:207515000:+' \
--highlight '207513650:207513800' \
-a none \
-o CD55_C_libsize_rollmean_50.pdf

~/Github/clipplotr/CLIPplotR.R \
-x 'hnRNPC_iCLIP_rep1_LUjh03_all_xlink_events.bedgraph.gz hnRNPC_iCLIP_rep2_LUjh25_all_xlink_events.bedgraph.gz U2AF65_iCLIP_ctrl_rep1_all_xlink_events.bedgraph.gz U2AF65_iCLIP_ctrl_rep2_all_xlink_events.bedgraph.gz U2AF65_iCLIP_KD1_rep2_all_xlink_events.bedgraph.gz U2AF65_iCLIP_KD2_rep1_all_xlink_events.bedgraph.gz' \
-l 'hnRNPC_1 hnRNPC_2 U2AF65_WT_1 U2AF65_WT_2 U2AF65_KD_1 U2AF65_KD_2' \
-c '#586BA4 #324376 #0AA398 #067E79 #A54D69 #771434' \
--groups 'hnRNPC hnRNPC U2AF65_WT U2AF65_WT U2AF65_KD U2AF65_KD' \
-n maxpeak \
-s rollmean \
-w 50 \
-g gencode.v34lift37.annotation.gtf.gz \
-r 'chr1:207513000:207515000:+' \
--highlight '207513650:207513800' \
-a none \
-o CD55_C_maxpeak_rollmean_50.pdf








