#!/bin/bash
CUDA_VISIBLE_DEVICES=0 python get_ism.py --model /srv/scratch/annashch/chrombpnet/gm12878_dnase/with_bias_unplugged/gm12878.dnase.with.bias.unplugged.0.hdf5 \
       --narrowPeak /srv/scratch/annashch/chrombpnet/gm12878_dnase/idr.optimal_peak.narrowPeak.gz \
       --ref_fasta /mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
       --out_pickle ism.gsm12878.8dil.background.50000.90000.p \
       --flank 1057 \
       --n_sample_for_background 1000 \
       --get_background \
       --peak_start 50000 \
       --peak_end 90000


