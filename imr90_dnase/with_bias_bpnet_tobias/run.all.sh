#!/bin/bash
outdir=/srv/scratch/annashch/chrombpnet/imr90_dnase/with_bias_bpnet_tobias
model=imr90.dnase.with.bpnet.tobias.bias
seed=1234
gpu=0

for fold in 0 1 2 3 4
do 
    ./train.sh $fold $gpu $model $seed $outdir
    ./predict.sh $fold $gpu $model $seed $outdir
    ./score.sh $outdir $model $fold &
done
