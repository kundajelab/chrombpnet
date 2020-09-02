#!/bin/bash
./train.sh 0 2 6mer.bias.preds.hepg2 1234 /srv/scratch/annashch/deeplearning/profile/hepg2_dnase/pred_from_bias params.6mer.txt
./train.counts.sh 0 2 counts.6mer.bias.preds.hepg2 1234 /srv/scratch/annashch/deeplearning/profile/hepg2_dnase/pred_from_bias params.6mer.txt
./train.profile.sh 0 2 profile.6mer.bias.preds.hepg2 1234 /srv/scratch/annashch/deeplearning/profile/hepg2_dnase/pred_from_bias params.6mer.txt
