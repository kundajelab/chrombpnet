#!/bin/bash
#./train.sh 0 0 bpnet.bias.preds.gm1287 1234 /srv/scratch/annashch/deeplearning/profile/gm12878_atac/pred_from_bias params.bpnet.txt
./train.sh 1 1 bpnet.bias.preds.gm1287 1234 /srv/scratch/annashch/deeplearning/profile/gm12878_atac/pred_from_bias params.bpnet.1.txt
./train.sh 2 1 bpnet.bias.preds.gm1287 1234 /srv/scratch/annashch/deeplearning/profile/gm12878_atac/pred_from_bias params.bpnet.2.txt
./train.sh 3 1 bpnet.bias.preds.gm1287 1234 /srv/scratch/annashch/deeplearning/profile/gm12878_atac/pred_from_bias params.bpnet.3.txt
./train.sh 4 1 bpnet.bias.preds.gm1287 1234 /srv/scratch/annashch/deeplearning/profile/gm12878_atac/pred_from_bias params.bpnet.4.txt
#./train.counts.sh 0 0 counts.bpnet.bias.preds.gm1287 1234 /srv/scratch/annashch/deeplearning/profile/gm12878_atac/pred_from_bias params.bpnet.txt
#./train.profile.sh 0 0 profile.bpnet.bias.preds.gm1287 1234 /srv/scratch/annashch/deeplearning/profile/gm12878_atac/pred_from_bias params.bpnet.txt
