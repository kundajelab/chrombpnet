#!/bin/bash
for fold in `seq 0 4`
do
    ./train.sh $fold 3 tobias.bias.preds.h1esc 1234 /srv/scratch/annashch/chrombpnet/h1esc_atac/pred_from_bias params.tobias.$fold.txt
done
