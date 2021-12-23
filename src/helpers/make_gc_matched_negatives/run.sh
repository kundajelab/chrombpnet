#!/bin/bash

foreground_bed=$1
exclude_bed=$2
inputlen=$3
output_dir=$4
genome=$5
genomewide_gc=$6
fold=$7
chrom_sizes=$8

echo "get gc content of the positive sequences" 
python $PWD/src/helpers/make_gc_matched_negatives/get_gc_content.py \
        --input_bed $foreground_bed \
        --genome $genome \
        --chrom_sizes $chrom_sizes \
        --output_prefix $output_dir/foreground.gc \
        --inputlen $inputlen \

echo "get candidate negative bed" 
bedtools intersect -v -a $genomewide_gc -b $exclude_bed  > $output_dir/candidate.negatives.bed

echo "find regions in candidate negative bed that gc-match wth foreground" 
python $PWD/src/helpers/make_gc_matched_negatives/get_gc_matched_negatives.py \
        --candidate_negatives $output_dir/candidate.negatives.bed \
        --foreground_gc_bed $output_dir/foreground.gc.bed \
        --output_prefix $output_dir/negatives \
        --chr_fold_path $fold \
        --neg_to_pos_ratio_train 2 \

