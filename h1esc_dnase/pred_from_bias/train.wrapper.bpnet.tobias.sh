#!/bin/bash
for fold in 0 1 2 3 4
do
    ./train.bpnet.tobias.sh $fold 1 bpnet.tobias.bias.preds.h1esc 1234 /srv/scratch/annashch/chrombpnet/h1esc_dnase/pred_from_bias params.bpnet.tobias.$fold.txt
done
