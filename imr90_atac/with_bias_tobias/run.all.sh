#!/bin/bash
outdir=/srv/scratch/annashch/chrombpnet/imr90_atac/with_bias_tobias
model=imr90.atac.with.tobias.bias
seed=1234
gpu=1

for fold in  0 1 2 3 4
do 
    ./train.sh $fold $gpu $model $seed $outdir
    ./predict.sh $fold $gpu $model $seed $outdir 
    ./score.sh $outdir $model $fold &
done
