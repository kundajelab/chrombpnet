#!/bin/bash
for fold in `seq 0 4`
do
    ./train.sh $fold 3 tobias.bias.preds.imr90 1234 /srv/scratch/annashch/chrombpnet/imr90_atac/pred_from_bias params.tobias.$fold.txt
done
