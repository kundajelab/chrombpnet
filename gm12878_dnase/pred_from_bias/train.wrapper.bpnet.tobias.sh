#!/bin/bash
for fold in 0 1 2 3 4
do
    ./train.bpnet.tobias.sh $fold 2 bpnet.tobias.bias.preds.gm12878 1234 /srv/scratch/annashch/chrombpnet/gm12878_dnase/pred_from_bias params.bpnet.tobias.$fold.txt
done
