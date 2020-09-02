#!/bin/bash
./predict.idr.sh 0 3 bpnet.bias.preds.gm12878 1234 /srv/scratch/annashch/deeplearning/profile/gm12878_atac/pred_from_bias/idr_preds
./predict.idr.sh 0 3 tobias.bias.preds.gm12878 1234 /srv/scratch/annashch/deeplearning/profile/gm12878_atac/pred_from_bias/idr_preds
./predict.idr.sh 0 3 6mer.bias.preds.gm12878 1234 /srv/scratch/annashch/deeplearning/profile/gm12878_atac/pred_from_bias/idr_preds

