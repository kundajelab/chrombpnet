#!/bin/bash
outdir=/srv/scratch/annashch/deeplearning/profile/hepg2_dnase/joint_tier1_train
model=hepg2.dnase.with.bpnet.tobias.bias
seed=1234
gpu=2

for fold in 3 #0 1 2 3 4
do 
    #./train.sh $fold $gpu $model $seed $outdir
    ./predict.sh $fold $gpu $model $seed $outdir
    ./score.sh $outdir $model $fold 
done
