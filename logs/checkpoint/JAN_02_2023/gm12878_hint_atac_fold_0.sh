#!/bin/bash

cell_line=GM12878
data_type="ATAC_PE"

date=$(date +'%m.%d.%Y')
setting=$data_type"_"$date"_hint_atac"
cur_file_name="gm12878_atac_fold_0.sh"
### SIGNAL INPUT

in_bigwig=results/hint_atac/ATAC_PE/GM12878/data_new/GM12878_unstranded.bw
overlap_peak=$PWD/results/chrombpnet/ATAC_PE/GM12878/data/peaks_no_blacklist.bed
nonpeaks=$PWD/results/chrombpnet/ATAC_PE/GM12878/negatives_data/negatives_with_summit.bed 

blacklist_region=$PWD/reference/GRch38_unified_blacklist.bed.gz
chrom_sizes=$PWD/reference/chrom.sizes
ref_fasta=$PWD/reference/hg38.genome.fa
fold=$PWD/splits/fold_0.json
#chrom_sizes=/mnt/data/annotations/by_release/hg38/hg38.chrom.sizes
#ref_fasta=/mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
#fold=/oak/stanford/groups/akundaje/projects/chrombpnet/model_inputs/ENCODE_ATAC_downloads/splits/fold_0.json

oak_dir=/oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/$data_type/$cell_line/

main_dir=$PWD/results/hint_atac/$data_type/$cell_line
data_dir=$main_dir/data
output_dir=$main_dir/$setting
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


### STEP 1 - TRAIN BIAS MODEL
if [[ -d $output_dir/hint_atac_model ]] ; then
    echo "skipping model training  - directory present "
else
    mkdir $output_dir/hint_atac_model
    CUDA_VISIBLE_DEVICES=$gpu bash hint_atac_model_no_pos.sh $ref_fasta $in_bigwig $overlap_peak $nonpeaks $fold $output_dir/hint_atac_model 
fi

#CUDA_VISIBLE_DEVICES=$gpu bash hint_atac_model_preds.sh $ref_fasta $in_bigwig $overlap_peak $nonpeaks $fold $output_dir/hint_atac_model $data_type

### INTERPRET BIAS MODEL


if [[ -d $output_dir/hint_atac_model/interpret ]] ; then
    echo "skipping bias model interpretation - directory present "
else
    mkdir $output_dir/hint_atac_model/interpret

    logfile=$output_dir/hint_atac_model/interpret/interpret.log
    touch $logfile

    echo $( timestamp ):"python $PWD/src/evaluation/interpret/interpret.py \
        --genome=$ref_fasta \
        --regions=results/chrombpnet/ATAC_PE/GM12878/data/30K.subsample.overlap.bed \
        --output_prefix=$output_dir/hint_atac_model/interpret/$cell_line \
        --model_h5=$output_dir/hint_atac_model/hint_atac.h5" |  tee -a $logfile

    CUDA_VISIBLE_DEVICES=$gpu python $PWD/src/evaluation/interpret/interpret.py \
        --genome=$ref_fasta \
        --regions=results/chrombpnet/ATAC_PE/GM12878/data/30K.subsample.overlap.bed \
        --output_prefix=$output_dir/hint_atac_model/interpret/$cell_line \
        --model_h5=$output_dir/hint_atac_model/hint_atac.h5  | tee -a $logfile
fi

if [[ -d $oak_dir/$setting/ ]]; then
    echo "dir exists"
else
    mkdir $oak_dir/$setting/
    mkdir $oak_dir/$setting/SIGNAL
fi

if [[ -f $oak_dir/$setting/SIGNAL/$cell_line.counts_scores.h5  && -f $oak_dir/$setting/SIGNAL/$cell_line.profile_scores.h5 ]] ; then
    echo "SIGNAL model files already present on oak"
else
    echo "copying SIGNAL model counts and profile scores to oak"
    cp $output_dir/hint_atac_model/interpret/$cell_line.counts_scores.h5 $oak_dir/$setting/SIGNAL/
    cp $output_dir/hint_atac_model/interpret/$cell_line.profile_scores.h5 $oak_dir/$setting/SIGNAL/
fi

