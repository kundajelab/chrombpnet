#!/bin/bash
outdir=/srv/scratch/annashch/chrombpnet/gm12878_atac/with_bias_tobias
model=gm12878.atac.with.tobias.bias
seed=1234
gpu=0

for fold in `seq 0 4`
do 
    ./train.sh $fold $gpu $model $seed $outdir
    ./predict.sh $fold $gpu $model $seed $outdir 
    ./score.sh $outdir $model $fold &
done

