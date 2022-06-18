#!/bin/bash

cell_line=GM12878
data_type="DNASE_SE"

date=$(date +'%m.%d.%Y')
setting=$data_type"_"$date
cur_file_name="gm12878_dnase_no_shift.sh"
### SIGNAL INPUT

in_bam=/oak/stanford/groups/akundaje/projects/chrombpnet/model_inputs/DNASE/unfiltered_bams/GM12878/GM12878.unfiltered.sorted.bam
overlap_peak=/oak/stanford/groups/akundaje/projects/chrombpnet/model_inputs/DNASE/optimal_overlap_peaks/GM12878.overlap.optimal_peak.narrowPeak.gz

blacklist_region=/mnt/data/annotations/blacklist/GRch38/GRch38_unified_blacklist.bed.gz
chrom_sizes=/mnt/data/annotations/by_release/hg38/hg38.chrom.sizes
ref_fasta=/mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta

oak_dir=/oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/$data_type/$cell_line/
main_dir=$PWD/results/chrombpnet/$data_type/$cell_line
data_dir=$main_dir/data_no_shift
output_dir=$main_dir/$setting
neg_dir=$main_dir/negatives_data
bias_threshold_factor=0.5
setting=$data_type"_"$date"_"$bias_threshold_factor
inputlen=2114
gpu=7

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



### CREATE BIGWIGS AND REMOVE BLACKLIST FROM PEAK EDGES

if [[ -d $data_dir ]] ; then
    echo "bigwig director already exists"
else
    mkdir $data_dir
    bash step2_make_bigwigs_from_bams_incorrect.sh $in_bam $data_dir"/"$cell_line $data_type $ref_fasta $chrom_sizes 
    
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


