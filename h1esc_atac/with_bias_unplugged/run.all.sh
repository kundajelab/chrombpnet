#!/bin/bash
outdir=/srv/scratch/annashch/chrombpnet/h1esc_atac/with_bias_unplugged
model=h1esc.atac.with.bias.unplugged
seed=1234
gpu=1

for fold in 0 1 2 3 4
do
    title="h1esc atac bias unplugged fold $fold"
    CUDA_VISIBLE_DEVICES=$gpu python get_model_with_bias_unplugged.py --model_params model_params.$fold.txt --outf $model.$fold.hdf5
    #./predict.sh $fold $gpu $model $seed $outdir
    #./score.sh $outdir $model $fold $title &
done
