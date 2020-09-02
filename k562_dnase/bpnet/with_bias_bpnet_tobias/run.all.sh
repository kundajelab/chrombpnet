#!/bin/bash
outdir=/srv/scratch/annashch/deeplearning/profile/k562_dnase/bpnet/with_bias_bpnet_tobias
model=k562.dnase.with.bpnet.tobias.bias
seed=1234
gpu=1

for fold in 0 1 2 3 4
do 
    #./train.sh $fold $gpu $model $seed $outdir
    #./predict.sh $fold $gpu $model $seed $outdir
    ./score.sh $outdir $model $fold &
done
