#!/bin/bash
zcat data/overlap.bed.gz | head -n1000 > data/overlap.1000.bed
chrombpnet_predict_to_bigwig  -bm outputs/models/bias_model/bias.h5 -cm outputs/models/chrombpnet_model/chrombpnet.h5  -cmb outputs/models/chrombpnet_model/chrombpnet_wo_bias.h5 -r data/overlap.1000.bed -g data/hg38.fa -c data/hg38.chrom.sizes -o outputs/pred.bigwig -b 500
