#!/bin/bash 
./score.sh /srv/scratch/annashch/deeplearning/profile/k562_dnase/bpnet/pred_from_bias bpnet.bias.preds.k562
./score.sh /srv/scratch/annashch/deeplearning/profile/k562_dnase/bpnet/pred_from_bias profile.bpnet.bias.preds.k562
./score.sh /srv/scratch/annashch/deeplearning/profile/k562_dnase/bpnet/pred_from_bias counts.bpnet.bias.preds.k562

./score.sh /srv/scratch/annashch/deeplearning/profile/k562_dnase/bpnet/pred_from_bias 6mer.bias.preds.k562
./score.sh /srv/scratch/annashch/deeplearning/profile/k562_dnase/bpnet/pred_from_bias profile.6mer.bias.preds.k562
./score.sh /srv/scratch/annashch/deeplearning/profile/k562_dnase/bpnet/pred_from_bias counts.6mer.bias.preds.k562


./score.sh /srv/scratch/annashch/deeplearning/profile/k562_dnase/bpnet/pred_from_bias tobias.bias.preds.k562
./score.sh /srv/scratch/annashch/deeplearning/profile/k562_dnase/bpnet/pred_from_bias profile.tobias.bias.preds.k562
./score.sh /srv/scratch/annashch/deeplearning/profile/k562_dnase/bpnet/pred_from_bias counts.tobias.bias.preds.k562


