#!/bin/bash
./predict.sh 0 0 bpnet.bias.preds.gm12878 1234 /srv/scratch/annashch/deeplearning/profile/gm12878_atac/pred_from_bias
./predict.sh 0 0 profile.bpnet.bias.preds.gm12878 1234 /srv/scratch/annashch/deeplearning/profile/gm12878_atac/pred_from_bias
./predict.sh 0 0 counts.bpnet.bias.preds.gm12878 1234 /srv/scratch/annashch/deeplearning/profile/gm12878_atac/pred_from_bias
