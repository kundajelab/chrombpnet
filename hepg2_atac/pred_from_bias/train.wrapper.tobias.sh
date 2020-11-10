#!/bin/bash
for fold in `seq 0 4`
do
    ./train.sh $fold 2 tobias.bias.preds.hepg2 1234 /srv/scratch/annashch/chrombpnet/hepg2_atac/pred_from_bias params.tobias.$fold.txt
done
