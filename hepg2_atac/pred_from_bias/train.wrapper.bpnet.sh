#!/bin/bash
./train.sh 0 0 bpnet.bias.preds.hepg2 1234 /srv/scratch/annashch/deeplearning/profile/hepg2_atac/pred_from_bias params.bpnet.txt
./train.counts.sh 0 0 counts.bpnet.bias.preds.hepg2 1234 /srv/scratch/annashch/deeplearning/profile/hepg2_atac/pred_from_bias params.bpnet.txt
./train.profile.sh 0 0 profile.bpnet.bias.preds.hepg2 1234 /srv/scratch/annashch/deeplearning/profile/hepg2_atac/pred_from_bias params.bpnet.txt
