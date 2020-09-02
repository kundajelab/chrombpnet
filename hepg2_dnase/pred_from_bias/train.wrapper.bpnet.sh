#!/bin/bash
for fold in 1 2 3 4
do
    ./train.sh $fold 0 bpnet.bias.preds.hepg2 1234 /srv/scratch/annashch/deeplearning/profile/hepg2_dnase/pred_from_bias/seed1234 params.bpnet.$fold.txt
done
#./train.counts.sh 0 2 counts.bpnet.bias.preds.hepg2 1234 /srv/scratch/annashch/deeplearning/profile/hepg2_dnase/pred_from_bias params.bpnet.txt
#./train.profile.sh 0 2 profile.bpnet.bias.preds.hepg2 1234 /srv/scratch/annashch/deeplearning/profile/hepg2_dnase/pred_from_bias params.bpnet.txt
