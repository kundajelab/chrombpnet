#!/bin/bash
outdir=/srv/scratch/annashch/deeplearning/profile/hepg2_atac/joint_tier1_train
model=hepg2.atac.with.tobias.bias
seed=1234
gpu=3

for fold in  0 1 2 3 4
do 
    ./train.sh $fold $gpu $model $seed $outdir
    ./predict.sh $fold $gpu $model $seed $outdir 
    ./score.sh $outdir $model $fold 
done
