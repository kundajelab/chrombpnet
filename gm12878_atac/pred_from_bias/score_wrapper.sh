#!/bin/bash 
./score.sh /srv/scratch/annashch/deeplearning/profile/gm12878_atac/pred_from_bias bpnet.bias.preds.gm12878
./score.sh /srv/scratch/annashch/deeplearning/profile/gm12878_atac/pred_from_bias profile.bpnet.bias.preds.gm12878
./score.sh /srv/scratch/annashch/deeplearning/profile/gm12878_atac/pred_from_bias counts.bpnet.bias.preds.gm12878

./score.sh /srv/scratch/annashch/deeplearning/profile/gm12878_atac/pred_from_bias 6mer.bias.preds.gm12878
./score.sh /srv/scratch/annashch/deeplearning/profile/gm12878_atac/pred_from_bias profile.6mer.bias.preds.gm12878
./score.sh /srv/scratch/annashch/deeplearning/profile/gm12878_atac/pred_from_bias counts.6mer.bias.preds.gm12878


./score.sh /srv/scratch/annashch/deeplearning/profile/gm12878_atac/pred_from_bias tobias.bias.preds.gm12878
./score.sh /srv/scratch/annashch/deeplearning/profile/gm12878_atac/pred_from_bias profile.tobias.bias.preds.gm12878
./score.sh /srv/scratch/annashch/deeplearning/profile/gm12878_atac/pred_from_bias counts.tobias.bias.preds.gm12878


