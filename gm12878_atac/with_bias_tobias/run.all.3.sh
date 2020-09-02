#!/bin/bash
outdir=/srv/scratch/annashch/deeplearning/profile/gm12878_atac/with_bias_tobias
model=gm12878.atac.with.tobias.bias
seed=1234
gpu=3

for fold in  3 #1 2 3 4
do 
    #./train.sh $fold $gpu $model $seed $outdir
    ./predict.sh $fold $gpu $model $seed $outdir 
    ./score.sh $outdir $model $fold
done
