#!/bin/bash 

./score.sh /srv/scratch/annashch/deeplearning/profile/hepg2_atac/pred_from_bias/idr_preds/seed1234 bpnet.bias.preds.hepg2
./score.sh /srv/scratch/annashch/deeplearning/profile/hepg2_atac/pred_from_bias/idr_preds/seed1234 6mer.bias.preds.hepg2
./score.sh /srv/scratch/annashch/deeplearning/profile/hepg2_atac/pred_from_bias/idr_preds/seed1234 tobias.bias.preds.hepg2

./score.sh /srv/scratch/annashch/deeplearning/profile/hepg2_atac/pred_from_bias/idr_preds/seed2345 bpnet.bias.preds.hepg2
./score.sh /srv/scratch/annashch/deeplearning/profile/hepg2_atac/pred_from_bias/idr_preds/seed2345 6mer.bias.preds.hepg2
./score.sh /srv/scratch/annashch/deeplearning/profile/hepg2_atac/pred_from_bias/idr_preds/seed2345 tobias.bias.preds.hepg2

./score.sh /srv/scratch/annashch/deeplearning/profile/hepg2_atac/pred_from_bias/idr_preds/seed3456 bpnet.bias.preds.hepg2
./score.sh /srv/scratch/annashch/deeplearning/profile/hepg2_atac/pred_from_bias/idr_preds/seed3456 6mer.bias.preds.hepg2
./score.sh /srv/scratch/annashch/deeplearning/profile/hepg2_atac/pred_from_bias/idr_preds/seed3456 tobias.bias.preds.hepg2


