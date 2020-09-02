#!/bin/bash
outdir=/srv/scratch/annashch/deeplearning/profile/k562_dnase/bpnet/joint_tier1_train
model=k562.dnase.with.bpnet.tobias.bias.joint.tier1.train
seed=1234
gpu=3
for fold in 2
do 
    #./train.sh $fold $gpu $model $seed $outdir
    #./predict.sh $fold $gpu $model $seed $outdir
    title="DNASE with bpnet/tobias, K562, fold $fold"
    ./score.sh $outdir $model $fold $title
done
