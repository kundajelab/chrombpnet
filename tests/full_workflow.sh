#!/bin/bash

# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG

GPU=${1-0}

#create directory structure for workflow 
data_dir=data
mkdir -p $data_dir
output_dir=outputs
mkdir -p $output_dir
output_data=$output_dir/data
mkdir -p $output_data
mkdir -p $output_data/subsample_peaks
output_models=$output_dir/models
mkdir -p $output_models
mkdir -p $output_models/bias_model
mkdir -p $output_models/bias_model/interpret
mkdir -p $output_models/chrombpnet_model/
mkdir -p $output_models/chrombpnet_model/interpret


step1_download_bams_and_peaks.sh $data_dir
echo "completed step 1"


in_bam=$data_dir/merged.bam
bigwig_prefix=$data_dir/merged
data_type=ATAC
reference_fasta=$data_dir/hg38.fa
chrom_sizes=$data_dir/hg38.chrom.sizes
step2_make_bigwigs_from_bams.sh $in_bam $bigwig_prefix $data_type $reference_fasta $chrom_sizes
echo "completed step 2"

mkdir -p $output_data/splits
chrombpnet_make_splits -o $output_data/splits
fold=$output_data/splits/fold_0.json
echo "made folds"

chrombpnet_genomewide_gc -g $reference_fasta -o $data_dir/genomewide_gc_hg38_stride_1000_inputlen_2114
echo "got genomewide gc"

blacklist_region=$data_dir/blacklist.bed.gz
overlap_peak=$data_dir/overlap.bed.gz
inputlen=2114
genomewide_gc=$data_dir/genomewide_gc_hg38_stride_1000_inputlen_2114.bed
step3_get_background_regions.sh $reference_fasta $chrom_sizes $blacklist_region $overlap_peak $inputlen $genomewide_gc $output_data $fold
echo "completed step 3"

bias_threshold_factor=0.5
nonpeaks=$output_data/negatives_with_summit.bed
bigwig_path=$bigwig_prefix\_unstranded.bw
CUDA_VISIBLE_DEVICES=$GPU step4_train_bias_model.sh $reference_fasta $bigwig_path $overlap_peak $nonpeaks $fold $bias_threshold_factor $output_models/bias_model
echo "completed step 4"

flank_size=$(( $inputlen/2 ))
echo $flank_size
bedtools slop -i $blacklist_region -g $chrom_sizes -b $flank_size > $output_data/subsample_peaks/temp.txt
bedtools intersect -v -a $overlap_peak -b $output_data/subsample_peaks/temp.txt | shuf  > $output_data/subsample_peaks/temp_n.txt
shuf -n 30000 $output_data/subsample_peaks/temp_n.txt > $output_data/subsample_peaks/30K.subsample.overlap.bed
rm  $output_data/subsample_peaks/temp.txt
rm $output_data/subsample_peaks/temp_n.txt

CUDA_VISIBLE_DEVICES=$GPU step5_interpret_bias_model.sh $reference_fasta $output_data/subsample_peaks/30K.subsample.overlap.bed $output_models/bias_model/bias.h5 $output_models/bias_model/interpret/
echo "completed step 5"

CUDA_VISIBLE_DEVICES=$GPU step6_train_chrombpnet_model.sh $reference_fasta $bigwig_path $overlap_peak $nonpeaks $fold $output_models/bias_model/bias.h5 $output_models/chrombpnet_model $data_type
echo "completed step 6" 

CUDA_VISIBLE_DEVICES=$GPU step7_interpret_chrombpnet_model.sh $reference_fasta $output_data/subsample_peaks/30K.subsample.overlap.bed $output_models/chrombpnet_model/chrombpnet_wo_bias.h5 $output_models/chrombpnet_model/interpret/
echo "completed step 7" 
