#!/bin/bash
#./train.sh 0 3 tobias.bias.preds.k562 1234 /srv/scratch/annashch/deeplearning/profile/k562_dnase/bpnet/pred_from_bias params.tobias.txt
./train.sh 1 3 tobias.bias.preds.k562 1234 /srv/scratch/annashch/deeplearning/profile/k562_dnase/bpnet/pred_from_bias params.tobias.1.txt
./train.sh 2 3 tobias.bias.preds.k562 1234 /srv/scratch/annashch/deeplearning/profile/k562_dnase/bpnet/pred_from_bias params.tobias.2.txt
./train.sh 3 3 tobias.bias.preds.k562 1234 /srv/scratch/annashch/deeplearning/profile/k562_dnase/bpnet/pred_from_bias params.tobias.3.txt
./train.sh 4 3 tobias.bias.preds.k562 1234 /srv/scratch/annashch/deeplearning/profile/k562_dnase/bpnet/pred_from_bias params.tobias.4.txt
#./train.counts.sh 0 0 counts.tobias.bias.preds.k562 1234 /srv/scratch/annashch/deeplearning/profile/k562_dnase/bpnet/pred_from_bias params.tobias.txt
#./train.profile.sh 0 0 profile.tobias.bias.preds.k562 1234 /srv/scratch/annashch/deeplearning/profile/k562_dnase/bpnet/pred_from_bias params.tobias.txt
