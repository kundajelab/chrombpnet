#!/bin/bash
#./train.sh 0 3 6mer.bias.preds.k562 1234 /srv/scratch/annashch/deeplearning/profile/k562_dnase/bpnet/pred_from_bias params.6mer.txt
./train.counts.sh 0 3 counts.6mer.bias.preds.k562 1234 /srv/scratch/annashch/deeplearning/profile/k562_dnase/bpnet/pred_from_bias params.6mer.txt
./train.profile.sh 0 3 profile.6mer.bias.preds.k562 1234 /srv/scratch/annashch/deeplearning/profile/k562_dnase/bpnet/pred_from_bias params.6mer.txt
