#!/bin/bash
for fold in 0 1 2 3 4
do
    ./train.bpnet.tobias.sh $fold 0 bpnet.tobias.bias.preds.imr90 1234 /srv/scratch/annashch/chrombpnet/imr90_dnase/pred_from_bias params.bpnet.tobias.$fold.txt
done
