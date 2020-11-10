#!/bin/bash
outdir=/srv/scratch/annashch/chrombpnet/hepg2_atac/with_bias_tobias
model=hepg2.atac.with.tobias.bias
seed=1234
gpu=3

for fold in  0 1 2 3 4
do 
    ./train.sh $fold $gpu $model $seed $outdir
    ./predict.sh $fold $gpu $model $seed $outdir 
    ./score.sh $outdir $model $fold &
done
