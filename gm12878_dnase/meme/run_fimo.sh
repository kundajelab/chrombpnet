#!/bin/bash

#make sequence input
#zcat idr.optimal_peak.narrowPeak.gz | bedtools sort -i - | cut -f1,2,3 | uniq | fastaFromBed -fi /mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta -bed - > meme/GM12878.idr.fa
fimo /mnt/data/meme/motif_databases/HUMAN/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme GM12878.idr.fa
