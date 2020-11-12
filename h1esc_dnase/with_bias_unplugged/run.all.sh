#!/bin/bash
outdir=/srv/scratch/annashch/chrombpnet/h1esc_dnase/with_bias_unplugged
model=h1esc.dnase.with.bias.unplugged
seed=1234
gpu=0

for fold in 0 1 2 3 4
do
    CUDA_VISIBLE_DEVICES=$gpu python get_model_with_bias_unplugged.py --model_params params.$fold.txt --outf $model.$fold.hdf5
done
