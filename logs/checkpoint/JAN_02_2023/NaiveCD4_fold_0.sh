#!/bin/bash

cell_line=NaiveCD4
data_type="ATAC_PE"

date=$(date +'%m.%d.%Y')
bias_threshold_factor=0.6
setting=$data_type"_"$date"_withgm12878bias"
cur_file_name="NaiveCD4_fold_0.sh"
### SIGNAL INPUT


in_bam=/oak/stanford/groups/akundaje/anusri/Mineto-data/NaiveCD4/data/1001_Naive_Teffs-U_R1.trim.merged.srt.nodup.no_chrM_MT.bam
overlap_peak=/oak/stanford/groups/akundaje/anusri/Mineto-data/NaiveCD4/data/overlap.optimal_peak.narrowPeak.gz

chrom_sizes=/mnt/data/annotations/by_release/hg38/hg38.chrom.sizes
ref_fasta=/mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
blacklist_region=/mnt/data/annotations/blacklist/GRch38/GRch38_unified_blacklist.bed.gz
genomewide_gc="/oak/stanford/groups/akundaje/anusri/refs/genomewide_gc_hg38_stride_1000_inputlen_2114.bed"

fold=/oak/stanford/groups/akundaje/projects/chrombpnet/model_inputs/ENCODE_ATAC_downloads/splits/fold_0.json
oak_dir=/oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/ATAC_PE/$cell_line/
mkdir $oak_dir

main_dir=/oak/stanford/groups/akundaje/anusri/Mineto-data/NaiveCD4/

output_dir=$main_dir/$setting
data_dir=$main_dir/data



inputlen=2114
gpu=2

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



### CREATE BIGWIGS AND REMOVE BLACKLIST FROM PEAK EDGES

if [[ -d $data_dir ]] ; then
    echo "bigwig director already exists"
    bash step2_make_bigwigs_from_bams.sh $in_bam $data_dir"/"$cell_line $data_type $ref_fasta $chrom_sizes 

    logfile=$data_dir"/"$cell_line"_preprocessing.log"
    flank_size=$(( inputlen/2 ))
    echo $( timestamp ): "bedtools slop -i $blacklist_region -g $chrom_sizes -b $flank_size > $data_dir/temp.txt" | tee -a $logfile
    bedtools slop -i $blacklist_region -g $chrom_sizes -b $flank_size > $data_dir/temp.txt
    echo $( timestamp ): "bedtools intersect -v -a $overlap_peak -b $data_dir/temp.txt > $data_dir/peaks_no_blacklist.bed" | tee -a $logfile
    bedtools intersect -v -a $overlap_peak -b $data_dir/temp.txt  > $data_dir/peaks_no_blacklist.bed
    echo $( timestamp ): "rm  $data_dir/temp.txt" | tee -a $logfile
    rm  $data_dir/temp.txt

    echo $( timestamp ): "shuf -n 30000 $data_dir/peaks_no_blacklist.bed > $data_dir/30K.subsample.overlap.bed" | tee -a $logfile
    shuf -n 30000 $data_dir/peaks_no_blacklist.bed > $data_dir/30K.subsample.overlap.bed

    cp $PWD/$cur_file_name $data_dir

else
    mkdir $data_dir
    bash step2_make_bigwigs_from_bams.sh $in_bam $data_dir"/"$cell_line $data_type $ref_fasta $chrom_sizes 
    
    logfile=$data_dir"/"$cell_line"_preprocessing.log"
    flank_size=$(( inputlen/2 ))
    echo $( timestamp ): "bedtools slop -i $blacklist_region -g $chrom_sizes -b $flank_size > $data_dir/temp.txt" | tee -a $logfile
    bedtools slop -i $blacklist_region -g $chrom_sizes -b $flank_size > $data_dir/temp.txt
    echo $( timestamp ): "bedtools intersect -v -a $overlap_peak -b $data_dir/temp.txt > $data_dir/peaks_no_blacklist.bed" | tee -a $logfile
    bedtools intersect -v -a $overlap_peak -b $data_dir/temp.txt  > $data_dir/peaks_no_blacklist.bed
    echo $( timestamp ): "rm  $data_dir/temp.txt" | tee -a $logfile
    rm  $data_dir/temp.txt

    echo $( timestamp ): "shuf -n 30000 $data_dir/peaks_no_blacklist.bed > $data_dir/30K.subsample.overlap.bed" | tee -a $logfile
    shuf -n 30000 $data_dir/peaks_no_blacklist.bed > $data_dir/30K.subsample.overlap.bed

    cp $PWD/$cur_file_name $data_dir

fi


overlap_peak=$data_dir/peaks_no_blacklist.bed
neg_dir=$main_dir/negatives_data


# MAKE NEGATIVES BED FILE
if [[ -d $neg_dir ]] ; then
    echo "negatives director already exists"
else
    mkdir $neg_dir
    bash step3_get_background_regions.sh $ref_fasta $chrom_sizes $blacklist_region $overlap_peak $inputlen $genomewide_gc $neg_dir $fold 
fi



transfer_bias_model=$PWD/results/chrombpnet/ATAC_PE/GM12878/nautilus_runs/GM12878_03.01.2022_bias_128_4_1234_0.4_fold_0/bias_model/bias.h5


### STEP 2 - TRAIN CHROMBPNET MODEL

if [[ -d $output_dir/chrombpnet_model ]] ; then
    echo "skipping chrombpnet model training  - directory present "
else
    mkdir $output_dir/chrombpnet_model
    CUDA_VISIBLE_DEVICES=$gpu bash step6_train_chrombpnet_model.sh $ref_fasta $data_dir"/"$cell_line"_unstranded.bw" $overlap_peak $neg_dir/negatives_with_summit.bed $fold $transfer_bias_model $output_dir/chrombpnet_model $data_type
fi



### INTERPRET SEQUENCE MODEL

if [[ -d $output_dir/chrombpnet_model/interpret ]] ; then
    echo "skipping sequence model interpretation - directory present "
else
    mkdir $output_dir/chrombpnet_model/interpret

    logfile=$output_dir/chrombpnet_model/interpret/interpret.log
    touch $logfile

    echo $( timestamp ):"python $PWD/src/evaluation/interpret/interpret.py \
        --genome=$ref_fasta \
        --regions=$data_dir/30K.subsample.overlap.bed \
        --output_prefix=$output_dir/chrombpnet_model/interpret/$cell_line \
        --model_h5=$output_dir/chrombpnet_model/chrombpnet_wo_bias.h5" |  tee -a $logfile

    CUDA_VISIBLE_DEVICES=$gpu python $PWD/src/evaluation/interpret/interpret.py \
        --genome=$ref_fasta \
        --regions=$data_dir/30K.subsample.overlap.bed \
        --output_prefix=$output_dir/chrombpnet_model/interpret/$cell_line \
        --model_h5=$output_dir/chrombpnet_model/chrombpnet_wo_bias.h5  | tee -a $logfile
fi


if [[ -d $oak_dir/$setting/ ]]; then
    echo "dir exists"
else
    mkdir $oak_dir/$setting/
    mkdir $oak_dir/$setting/SIGNAL
fi



if [[ -f $oak_dir/$setting/SIGNAL/$cell_line.counts_scores.h5  && -f $oak_dir/$setting/SIGNAL/$cell_line.profile_scores.h5 ]] ; then
    echo "chrombpnet model files already present on oak"
else
    echo "copying counts and profile scores to oak"
    cp $output_dir/chrombpnet_model/interpret/$cell_line.counts_scores.h5 $oak_dir/$setting/SIGNAL/
    cp $output_dir/chrombpnet_model/interpret/$cell_line.profile_scores.h5 $oak_dir/$setting/SIGNAL/
fi

