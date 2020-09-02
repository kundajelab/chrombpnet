#!/bin/bash
for fold in 0 1 2 3 4
do
    ./train.bpnet.tobias.sh $fold 3 bpnet.tobias.bias.preds.k562 1234 /srv/scratch/annashch/deeplearning/profile/k562_dnase/bpnet/pred_from_bias params.bpnet.tobias.$fold.txt
done
#./train.sh 0 0 bpnet.bias.preds.k562 1234 /srv/scratch/annashch/deeplearning/profile/k562_dnase/bpnet/pred_from_bias params.bpnet.txt
#./train.sh 1 0 bpnet.bias.preds.k562 1234 /srv/scratch/annashch/deeplearning/profile/k562_dnase/bpnet/pred_from_bias params.bpnet.1.txt
#./train.sh 2 0 bpnet.bias.preds.k562 1234 /srv/scratch/annashch/deeplearning/profile/k562_dnase/bpnet/pred_from_bias params.bpnet.2.txt
#./train.sh 3 0 bpnet.bias.preds.k562 1234 /srv/scratch/annashch/deeplearning/profile/k562_dnase/bpnet/pred_from_bias params.bpnet.3.txt
#./train.sh 4 0 bpnet.bias.preds.k562 1234 /srv/scratch/annashch/deeplearning/profile/k562_dnase/bpnet/pred_from_bias params.bpnet.4.txt
#./train.counts.sh 0 0 counts.bpnet.bias.preds.k562 1234 /srv/scratch/annashch/deeplearning/profile/k562_dnase/bpnet/pred_from_bias params.bpnet.txt
#./train.profile.sh 0 0 profile.bpnet.bias.preds.k562 1234 /srv/scratch/annashch/deeplearning/profile/k562_dnase/bpnet/pred_from_bias params.bpnet.txt
