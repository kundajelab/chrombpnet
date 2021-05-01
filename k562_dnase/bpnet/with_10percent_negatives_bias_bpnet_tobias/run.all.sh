#!/bin/bash
outdir=/srv/scratch/annashch/chrombpnet/k562_dnase/bpnet/with_10percent_negatives_bias_bpnet_tobias
model=k562.dnase.with.bpnet.tobias.bias
seed=1234
gpu=2

for fold in 0
do 
    #./train.sh $fold $gpu $model $seed $outdir
    ./predict.sh $fold $gpu $model $seed $outdir
    ./score.sh $outdir $model $fold 
done
