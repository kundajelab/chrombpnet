#!/bin/bash

model_name=$1
bed_regions=$2
split=$3
tdb_array=$4
chrom_sizes=$5
output_dir=$6
cell_line=$7
gpu=$8
fold=$9
CUDA_VISIBLE_DEVICES=$gpu python /srv/scratch/anusri/pipeline_script/main_scripts/interpret/bpnet_shap_wrapper.py \
                        --model_hdf5 $model_name \
                        --ref_fasta /mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
                        --bed_regions $bed_regions/$split \
                        --bed_regions_center center \
                        --tdb_array $tdb_array \
                        --chrom_sizes $chrom_sizes \
                        --tdb_input_datasets seq \
                        --tdb_output_datasets $cell_line $cell_line \
                        --batch_size 200 \
                        --tdb_output_source_attribute count_bigwig_unstranded_5p count_bigwig_unstranded_5p \
                        --tdb_output_flank 500 500 \
                        --tdb_output_aggregation None sum \
                        --tdb_output_transformation None log \
                        --tdb_input_source_attribute seq \
                        --tdb_input_flank 1057 \
                        --tdb_input_aggregation None \
                        --tdb_input_transformation None \
                        --out_pickle $output_dir/scores.$split.fold$fold.deepSHAP \
                        --num_threads 10 \
                        --num_inputs 1 \
                        --num_outputs 2 \
                        --task_index 0




