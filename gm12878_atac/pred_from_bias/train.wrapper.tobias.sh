#!/bin/bash
#./train.sh 0 3 tobias.bias.preds.gm12878 1234 /srv/scratch/annashch/deeplearning/profile/gm12878_atac/pred_from_bias params.tobias.txt
./train.sh 1 2 tobias.bias.preds.gm12878 1234 /srv/scratch/annashch/deeplearning/profile/gm12878_atac/pred_from_bias params.tobias.1.txt
./train.sh 2 2 tobias.bias.preds.gm12878 1234 /srv/scratch/annashch/deeplearning/profile/gm12878_atac/pred_from_bias params.tobias.2.txt
./train.sh 3 2 tobias.bias.preds.gm12878 1234 /srv/scratch/annashch/deeplearning/profile/gm12878_atac/pred_from_bias params.tobias.3.txt
./train.sh 4 2 tobias.bias.preds.gm12878 1234 /srv/scratch/annashch/deeplearning/profile/gm12878_atac/pred_from_bias params.tobias.4.txt
#./train.counts.sh 0 3 counts.tobias.bias.preds.gm12878 1234 /srv/scratch/annashch/deeplearning/profile/gm12878_atac/pred_from_bias params.tobias.txt
#./train.profile.sh 0 3 profile.tobias.bias.preds.gm12878 1234 /srv/scratch/annashch/deeplearning/profile/gm12878_atac/pred_from_bias params.tobias.txt
