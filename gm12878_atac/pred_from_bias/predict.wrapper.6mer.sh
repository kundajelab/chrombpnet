#!/bin/bash
./predict.sh 0 1 6mer.bias.preds.gm12878 1234 /srv/scratch/annashch/deeplearning/profile/gm12878_atac/pred_from_bias
./predict.sh 0 1 profile.6mer.bias.preds.gm12878 1234 /srv/scratch/annashch/deeplearning/profile/gm12878_atac/pred_from_bias
./predict.sh 0 1 counts.6mer.bias.preds.gm12878 1234 /srv/scratch/annashch/deeplearning/profile/gm12878_atac/pred_from_bias

