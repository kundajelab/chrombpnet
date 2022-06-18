#!/bin/bash

cell_line=HEPG2
data_type="ATAC_PE"

date=$(date +'%m.%d.%Y')
cur_file_name="hepg2_negatives_create.sh"
### SIGNAL INPUT

in_bam=/oak/stanford/groups/akundaje/projects/chrombpnet/model_inputs/ENCODE_ATAC_downloads/HEPG2/sorted_merged.bam
overlap_peak=/oak/stanford/groups/akundaje/projects/chrombpnet/model_inputs/ENCODE_ATAC_downloads/HEPG2/peaks.bed.gz

blacklist_region=/mnt/data/annotations/blacklist/GRch38/GRch38_unified_blacklist.bed.gz
chrom_sizes=/mnt/data/annotations/by_release/hg38/hg38.chrom.sizes
ref_fasta=/mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
genomewide_gc="/oak/stanford/groups/akundaje/anusri/refs/genomewide_gc_hg38_stride_1000_inputlen_2114.bed"
fold=/oak/stanford/groups/akundaje/projects/chrombpnet/model_inputs/ENCODE_ATAC_downloads/splits/fold_0.json
oak_dir=/oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/$data_type/$cell_line/
main_dir=$PWD/results/chrombpnet/$data_type/$cell_line
data_dir=$main_dir/data
output_dir=$main_dir/$setting
neg_dir=$main_dir/negatives_data
inputlen=2114
gpu=1

function timestamp {
    # Function to get the current time with the new line character
    # removed 
    
    # current time
    date +"%Y-%m-%d_%H-%M-%S" | tr -d '\n'
}


overlap_peak=$data_dir/peaks_no_blacklist.bed

## MAKE NEGATIVES BED FILE
if [[ -d $neg_dir ]] ; then
    echo "negatives director already exists"
else
    mkdir $neg_dir
    bash step3_get_background_regions.sh $ref_fasta $chrom_sizes $blacklist_region $overlap_peak $inputlen $genomewide_gc $neg_dir $fold 
fi
