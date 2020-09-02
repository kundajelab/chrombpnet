#!/bin/bash
outdir=/srv/scratch/annashch/deeplearning/profile/hepg2_atac/with_bias_unplugged
model=hepg2.atac.with.bias.unplugged
seed=1234
gpu=3

for fold in 0 1 2 3 4
do
    title="k562 DNASE bias unplugged fold $fold"
    #CUDA_VISIBLE_DEVICES=$gpu python get_model_with_bias_unplugged.py --model_params model_params.$fold.txt --outf $model.$fold.hdf5
    #./predict.sh $fold $gpu $model $seed $outdir
    ./score.sh $outdir $model $fold $title
done
