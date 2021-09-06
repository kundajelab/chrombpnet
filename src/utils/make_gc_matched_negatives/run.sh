#!/bin/bash

task=$1 
idr=$2 
ref_fasta=$3
chromo_sizes=$4
output_prefix=$task/genomewide_gc_hg38
flank_size=$5
stride=$6

echo "starting $task $idr" 

bash $PWD/main_scripts/make_gc_matched_negatives/get_svm_peak_splits.sh $task $idr
echo "got svm peak splits" 

bash $PWD/main_scripts/make_gc_matched_negatives/get_gc_positives.sh $task $ref_fasta $flank_size
echo "got gc content of the positive sequences" 

bash $PWD/main_scripts/make_gc_matched_negatives/get_genomewide_gc_bins.sh $ref_fasta $chromo_sizes $output_prefix $flank_size $stride
echo "got genome wide GC tracks"

bash $PWD/main_scripts/make_gc_matched_negatives/get_all_negatives.sh $task $idr $output_prefix.tsv $chromo_sizes $flank_size
echo "got candidate negative set" 

bash $PWD/main_scripts/make_gc_matched_negatives/get_chrom_gc_region_dict.sh $task
echo "created python pickle for candidate negatives" 

bash $PWD/main_scripts/make_gc_matched_negatives/form_svm_input_fastas.sh $task $ref_fasta
echo "finished creating SVM inputs" 

 
bash $PWD/main_scripts/make_gc_matched_negatives/script_make_bed.sh $task/bpnet.inputs $chromo_sizes
