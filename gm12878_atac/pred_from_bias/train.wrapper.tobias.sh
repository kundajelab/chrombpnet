#!/bin/bash
for fold in 2 #`seq 0 4`
do
    ./train.sh $fold 0 tobias.bias.preds.gm12878 1234 /srv/scratch/annashch/chrombpnet/gm12878_atac/pred_from_bias params.tobias.$fold.txt
done
