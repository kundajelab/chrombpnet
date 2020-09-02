#!/bin/bash
for fold in 0 1 2 3 4
do
    
    ./train.sh $fold 3 tobias.bias.preds.hepg2 1234 /srv/scratch/annashch/deeplearning/profile/hepg2_atac/pred_from_bias/seed1234 params.tobias.$fold.txt
done
