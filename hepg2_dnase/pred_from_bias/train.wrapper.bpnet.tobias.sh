#!/bin/bash
for fold in 0 1 2 3 4
do
    ./train.bpnet.tobias.sh $fold 0 bpnet.tobias.bias.preds.hepg2 1234 /srv/scratch/annashch/chrombpnet/hepg2_dnase/pred_from_bias/seed1234 params.bpnet.tobias.$fold.txt
done
