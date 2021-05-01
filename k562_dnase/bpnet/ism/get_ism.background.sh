#!/bin/bash
CUDA_VISIBLE_DEVICES=1 python get_ism.py --model /srv/scratch/annashch/chrombpnet/k562_dnase/bpnet/with_bias_unplugged/k562.dnase.with.bias.unplugged.0.hdf5 \
              --narrowPeak /srv/scratch/annashch/chrombpnet/k562_dnase/bpnet/K562.dnase.idr.optimal_peak.narrowPeak.gz \
	      --ref_fasta /mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
	      --out_pickle ism.k562.8dil.background.10000.90000.p \
	      --flank 1057 \
	      --n_sample_for_background 1000 \
	      --get_background \
	      --peak_start 10000 \
	      --peak_end 90000


