#!/bin/bash
CUDA_VISIBLE_DEVICES=2 python get_ism_foreground_and_background.py --model /srv/scratch/annashch/chrombpnet/gm12878_dnase/with_bias_unplugged/gm12878.dnase.with.bias.unplugged.0.hdf5 \
       --narrowPeak /srv/scratch/annashch/chrombpnet/gm12878_dnase/idr.optimal_peak.narrowPeak.gz \
       --ref_fasta /mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
       --out_prefix ism.gsm12878.8dil.fold0.75000.80000 \
       --flank 1057 \
       --peak_start 75000 \
       --peak_end 80000

