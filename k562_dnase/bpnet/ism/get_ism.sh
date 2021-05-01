#!/bin/bash
CUDA_VISIBLE_DEVICES=3 python get_ism.py --model /srv/scratch/annashch/chrombpnet/k562_dnase/bpnet/with_bias_unplugged/k562.dnase.with.bias.unplugged.0.hdf5 \
       --narrowPeak /srv/scratch/annashch/chrombpnet/k562_dnase/bpnet/K562.dnase.idr.optimal_peak.narrowPeak.gz \
       --ref_fasta /mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
       --out_pickle ism.k562.8dil.0.1000.p \
       --flank 1057 \
       --peak_start 0 \
       --peak_end 1000

