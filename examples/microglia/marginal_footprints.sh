#!/bin/bash

cell_line=microglia
data_type="ATAC"
neg_shift=4
filters=128
n_dil_layers=4

date=$(date +'%m.%d.%Y')
setting=4_$neg_shift"_shifted_"$data_type"_"$date"_bias_filters_"$filters
cur_file_name="microglia_atac_bias_filters_"$filters".sh"
setting=4_4_shifted_ATAC_10.11.2021_bias_filters_128

chrom_sizes=/data/refs/hg38.chrom.sizes
ref_fasta=/data/refs/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta

main_dir=$PWD/results/chrombpnet/$data_type/$cell_line
#mkdir -p $main_dir
data_dir=$PWD
output_dir=$PWD/results/chrombpnet/$data_type/$cell_line/$setting
#mkdir -p $output_dir

### MODEL PARAMS

gpu=1
seed=1234 
model_name=microglia 
neg_dir=$PWD/negatives
flank_size=1057

neg_bed_train=$neg_dir"/train.fold0.negatives.bed"
neg_bed_test=$neg_dir"/test.fold0.negatives.bed"



## MAKE FOOTPRINTS


#CUDA_VISIBLE_DEVICES=$gpu python  $PWD/src/evaluation/marginal_footrprints/main_footprints_check.py --ref_fasta $ref_fasta --gc_neg $neg_bed_test --motifs tn5_c --model_dir $output_dir/final_model_step3/unplug/
#CUDA_VISIBLE_DEVICES=$gpu python  $PWD/src/evaluation/marginal_footrprints/main_footprints_check.py --ref_fasta $ref_fasta --gc_neg $neg_bed_test --motifs gm12878_motifs_set1 --model_dir $output_dir/final_model_step3/unplug/
#CUDA_VISIBLE_DEVICES=$gpu python  $PWD/src/evaluation/marginal_footrprints/main_footprints_check.py --ref_fasta $ref_fasta --gc_neg $neg_bed_test --motifs gm12878_motifs_set2 --model_dir $output_dir/final_model_step3/unplug/

output_dir=/srv/scratch/annashch/chrombpnet/microglia/results/chrombpnet/ATAC/microglia/4_4_shifted_ATAC_10.11.2021_bias_filters_128
CUDA_VISIBLE_DEVICES=$gpu python  $PWD/src/evaluation/marginal_footrprints/main_footprints_check.py --ref_fasta $ref_fasta --gc_neg $neg_bed_test.narrowPeak --motifs tn5 --model_dir $output_dir/invivo_bias_model_step1/
CUDA_VISIBLE_DEVICES=$gpu python  $PWD/src/evaluation/marginal_footrprints/main_footprints_check.py --ref_fasta $ref_fasta --gc_neg $neg_bed_test.narrowPeak --motifs tn5_c --model_dir $output_dir/invivo_bias_model_step1/
CUDA_VISIBLE_DEVICES=$gpu python  $PWD/src/evaluation/marginal_footrprints/main_footprints_check.py --ref_fasta $ref_fasta --gc_neg $neg_bed_test.narrowPeak --motifs microglia_motifs_set --model_dir $output_dir/invivo_bias_model_step1/

output_dir=/srv/scratch/annashch/chrombpnet/microglia/results/chrombpnet/ATAC/microglia/4_4_shifted_ATAC_10.11.2021_bias_filters_500
CUDA_VISIBLE_DEVICES=$gpu python  $PWD/src/evaluation/marginal_footrprints/main_footprints_check.py --ref_fasta $ref_fasta --gc_neg $neg_bed_test.narrowPeak --motifs tn5 --model_dir $output_dir/invivo_bias_model_step1/
CUDA_VISIBLE_DEVICES=$gpu python  $PWD/src/evaluation/marginal_footrprints/main_footprints_check.py --ref_fasta $ref_fasta --gc_neg $neg_bed_test.narrowPeak --motifs tn5_c --model_dir $output_dir/invivo_bias_model_step1/
CUDA_VISIBLE_DEVICES=$gpu python  $PWD/src/evaluation/marginal_footrprints/main_footprints_check.py --ref_fasta $ref_fasta --gc_neg $neg_bed_test.narrowPeak --motifs microglia_motifs_set --model_dir $output_dir/invivo_bias_model_step1/
