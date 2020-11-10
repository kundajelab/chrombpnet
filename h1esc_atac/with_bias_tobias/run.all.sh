#!/bin/bash
outdir=/srv/scratch/annashch/chrombpnet/h1esc_atac/with_bias_tobias
model=h1esc.atac.with.tobias.bias
seed=1234
gpu=2

for fold in `seq 0 4`
do 
    ./train.sh $fold $gpu $model $seed $outdir
    ./predict.sh $fold $gpu $model $seed $outdir 
    ./score.sh $outdir $model $fold &
done
