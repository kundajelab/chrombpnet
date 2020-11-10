#!/bin/bash
outdir=/srv/scratch/annashch/chrombpnet/h1esc_dnase/with_bias_bpnet_tobias
model=h1esc.dnase.with.bpnet.tobias.bias
seed=1234
gpu=1

for fold in 2 3 4
do 
    ./train.sh $fold $gpu $model $seed $outdir
    ./predict.sh $fold $gpu $model $seed $outdir
    ./score.sh $outdir $model $fold 
done
