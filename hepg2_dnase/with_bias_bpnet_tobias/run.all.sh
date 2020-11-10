#!/bin/bash
outdir=/srv/scratch/annashch/chrombpnet/hepg2_dnase/with_bias_bpnet_tobias
model=hepg2.dnase.with.bpnet.tobias.bias
seed=1234
gpu=2

for fold in 2 3 4
do 
    ./train.sh $fold $gpu $model $seed $outdir
    ./predict.sh $fold $gpu $model $seed $outdir
    ./score.sh $outdir $model $fold &
done
