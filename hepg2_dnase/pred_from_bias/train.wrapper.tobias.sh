#!/bin/bash
for fold in 1 2 3 4
do
    ./train.sh $fold 2 tobias.bias.preds.hepg2 1234 /srv/scratch/annashch/deeplearning/profile/hepg2_dnase/pred_from_bias/seed1234 params.tobias.$fold.txt
done

#./train.counts.sh 0 2 counts.tobias.bias.preds.hepg2 2345 /srv/scratch/annashch/deeplearning/profile/hepg2_dnase/pred_from_bias params.tobias.txt
#./train.profile.sh 0 2 profile.tobias.bias.preds.hepg2 2345 /srv/scratch/annashch/deeplearning/profile/hepg2_dnase/pred_from_bias params.tobias.txt
