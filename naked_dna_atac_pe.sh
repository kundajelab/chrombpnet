#!/bin/bash

cell_line=NAKED
data_type="ATAC_PE"

date=$(date +'%m.%d.%Y')
cur_file_name="naked_dna_atac_pe.sh"
### SIGNAL INPUT

in_bam=/oak/stanford/groups/akundaje/projects/enzymatic_bias_correction/pipeline_out/v2.1.2/atac/8728c4b2-29de-4628-9f13-3d38baf96991/call-filter/shard-0/glob-3bcbe4e7489c90f75e0523ac6f3a9385/SRR072187.4_1.srt.nodup.no_chrM_MT.bam
chrom_sizes=/mnt/data/annotations/by_release/hg38/hg38.chrom.sizes
ref_fasta=/mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta

main_dir=$PWD/results/chrombpnet/$data_type/$cell_line
data_dir=$main_dir/data
inputlen=2114
output_dir=$PWD/results/chrombpnet/$data_type/$cell_line/bias_model_$date/

bias_filters=128
bias_dil=4
seed=1234
bias_threshold_factor=1.0
chrom_sizes=reference/chrom.sizes
ref_fasta=reference/hg38.genome.fa
fold=splits/fold_0.json
gpu=0

overlap_peak=results/chrombpnet/ATAC_PE/GM12878/data/peaks_no_blacklist.bed
negatives_dir=results/chrombpnet/ATAC_PE/GM12878/negatives_data/negatives_with_summit.bed

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

if [[ -d $data_dir ]] ; then
    echo "bigwig director already exists"
else
    mkdir $data_dir
    bash step2_make_bigwigs_from_bams.sh $in_bam $data_dir"/"$cell_line $data_type $ref_fasta $chrom_sizes     
fi


### STEP 1 - TRAIN BIAS MODEL


if [[ -d $output_dir/bias_model ]] ; then
    echo "skipping model training  - directory present "
else
    mkdir $output_dir/bias_model
    CUDA_VISIBLE_DEVICES=$gpu bash hint_atac_model_no_pos.sh $ref_fasta $data_dir"/"$cell_line"_unstranded.bw" $overlap_peak $negatives_dir $fold $output_dir/bias_model 128 4 $data_type
fi

### INTERPRET BIAS MODEL


if [[ -d $output_dir/bias_model/interpret ]] ; then
    echo "skipping bias model interpretation - directory present "
else
    mkdir $output_dir/bias_model/interpret
    logfile=$output_dir/bias_model/interpret/interpret.log
    touch $logfile

    echo $( timestamp ):"python $PWD/src/evaluation/interpret/interpret.py \
        --genome=$ref_fasta \
        --regions=results/chrombpnet/ATAC_PE/GM12878/data/30K.subsample.overlap.bed \
        --output_prefix=$output_dir/bias_model/interpret/$cell_line \
        --model_h5=$output_dir/bias_model/hint_atac.h5" |  tee -a $logfile

    python $PWD/src/evaluation/interpret/interpret.py \
        --genome=$ref_fasta \
        --regions=results/chrombpnet/ATAC_PE/GM12878/data/30K.subsample.overlap.bed \
        --output_prefix=$output_dir/bias_model/interpret/$cell_line \
        --model_h5=$output_dir/bias_model/hint_atac.h5  | tee -a $logfile
fi
