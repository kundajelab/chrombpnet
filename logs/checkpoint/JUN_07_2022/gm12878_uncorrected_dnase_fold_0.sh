#!/bin/bash

cell_line=GM12878
data_type="DNASE_SE"

date=$(date +'%m.%d.%Y')
cur_file_name="gm12878_uncorrected_atac_fold_0.sh"
### SIGNAL INPUT

chrom_sizes=/mnt/data/annotations/by_release/hg38/hg38.chrom.sizes
ref_fasta=/mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta

main_dir=$PWD/results/chrombpnet/$data_type/$cell_line
data_dir=$main_dir/data
inputlen=2114
#output_dir=$PWD/results/chrombpnet/$data_type/$cell_line/uncorrected_model_$date/
output_dir=$PWD/results/chrombpnet/$data_type/$cell_line/uncorrected_model_05.10.2022/


filters=512
dil=8
seed=1234
bias_threshold_factor=1.0
chrom_sizes=reference/chrom.sizes
ref_fasta=reference/hg38.genome.fa
fold=splits/fold_0.json
gpu=0

overlap_peak=results/chrombpnet/DNASE_SE/GM12878/data/peaks_no_blacklist.bed
negatives_dir=results/chrombpnet/DNASE_SE/GM12878/negatives_data/negatives_with_summit.bed

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

## CREATE DIRS
if [[ -d $output_dir ]] ; then
    echo "output director already exists"
else
    mkdir $output_dir
fi


### CREATE BIGWIGS AND REMOVE BLACKLIST FROM PEAK EDGES


### STEP 1 - TRAIN BIAS MODEL


if [[ -d $output_dir/uncorrected_model ]] ; then
    echo "skipping model training  - directory present "
else
    mkdir $output_dir/uncorrected_model
    CUDA_VISIBLE_DEVICES=$gpu bash hint_atac_model_no_pos.sh $ref_fasta $data_dir"/"$cell_line"_unstranded.bw" $overlap_peak $negatives_dir $fold $output_dir/uncorrected_model $filters $dil $data_type
fi

### INTERPRET BIAS MODEL


if [[ -d $output_dir/uncorrected_model/interpret ]] ; then
    echo "skipping bias model interpretation - directory present "
else
    mkdir $output_dir/uncorrected_model/interpret
    logfile=$output_dir/uncorrected_model/interpret/interpret.log
    touch $logfile

    echo $( timestamp ):"python $PWD/src/evaluation/interpret/interpret.py \
        --genome=$ref_fasta \
        --regions=results/chrombpnet/DNASE_SE/GM12878/data/30K.subsample.overlap.bed \
        --output_prefix=$output_dir/uncorrected_model/interpret/$cell_line \
        --model_h5=$output_dir/uncorrected_model/hint_atac.h5" |  tee -a $logfile

    python $PWD/src/evaluation/interpret/interpret.py \
        --genome=$ref_fasta \
        --regions=results/chrombpnet/DNASE_SE/GM12878/data/30K.subsample.overlap.bed \
        --output_prefix=$output_dir/uncorrected_model/interpret/$cell_line \
        --model_h5=$output_dir/uncorrected_model/hint_atac.h5  | tee -a $logfile
fi
