#!/bin/bash
outdir=/srv/scratch/annashch/chrombpnet/gm12878_dnase/with_10percent_negatives_bias_bpnet_tobias
model=gm12878.dnase.with.bias.unplugged
seed=1234
gpu=3

for fold in 0 #1 2 3 4
do
    #CUDA_VISIBLE_DEVICES=$gpu python get_model_with_bias_unplugged.py --model_params params.$fold.txt --outf $model.$fold.hdf5
    ./predict.sh $fold $gpu $model $seed $outdir
    ./score.sh $outdir $model $fold    
    
done
