#!/bin/bash
./predict.idr.peaks.sh 0 0 bpnet.bias.preds.k562 1234 /srv/scratch/annashch/deeplearning/profile/k562_dnase/bpnet/pred_from_bias/idr_preds
./predict.idr.peaks.sh 0 0 tobias.bias.preds.k562 1234 /srv/scratch/annashch/deeplearning/profile/k562_dnase/bpnet/pred_from_bias/idr_preds
./predict.idr.peaks.sh 0 0 6mer.bias.preds.k562 1234 /srv/scratch/annashch/deeplearning/profile/k562_dnase/bpnet/pred_from_bias/idr_preds

