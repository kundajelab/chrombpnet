#!/bin/bash

cell_line=K562
data_type="ATAC_PE"

date=$(date +'%m.%d.%Y')
cur_file_name="k562_wih_reduced_jitter.sh"
bias_threshold_factor=0.9
setting=$data_type"_"$date"_testing_gen"
setting=ATAC_PE_01.24.2022_testing_gen
### SIGNAL INPUT

in_bam=/oak/stanford/groups/akundaje/projects/chrombpnet/model_inputs/ENCODE_ATAC_downloads/K562/sorted_merged.bam
overlap_peak=/oak/stanford/groups/akundaje/projects/chrombpnet/model_inputs/ENCODE_ATAC_downloads/K562/peaks.bed.gz

blacklist_region=/mnt/data/annotations/blacklist/GRch38/GRch38_unified_blacklist.bed.gz
chrom_sizes=/mnt/data/annotations/by_release/hg38/hg38.chrom.sizes
ref_fasta=/mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
genomewide_gc="/oak/stanford/groups/akundaje/anusri/refs/genomewide_gc_hg38_stride_50_inputlen_2114_no_header.bed"
fold=/oak/stanford/groups/akundaje/projects/chrombpnet/model_inputs/ENCODE_ATAC_downloads/splits/fold_0.json
oak_dir=/oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/$data_type/$cell_line/
main_dir=$PWD/results/chrombpnet/$data_type/$cell_line
data_dir=$main_dir/data
output_dir=$main_dir/$setting
neg_dir=$main_dir/negatives_data
inputlen=2114
gpu=0

function timestamp {
    # Function to get the current time with the new line character
    # removed 
    
    # current time
    date +"%Y-%m-%d_%H-%M-%S" | tr -d '\n'
}


## CREATE DIRS
if [[ -d $main_dir ]] ; then
    echo "main director already exists"
else
    mkdir $main_dir
fi

if [[ -d $output_dir ]] ; then
    echo "output director already exists"
else
    mkdir $output_dir
fi

overlap_peak=$data_dir/peaks_no_blacklist.bed
bias_model="/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/K562/ATAC_PE_12.30.2021/bias_model/bias.h5"

### STEP 2 - TRAIN CHROMBPNET MODEL
CUDA_VISIBLE_DEVICES=$gpu bash testing_generator.sh $ref_fasta $data_dir"/"$cell_line"_unstranded.bw" $overlap_peak $neg_dir/negatives_with_summit.bed $fold $bias_model $output_dir/chrombpnet_model $data_type

