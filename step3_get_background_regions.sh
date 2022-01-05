#!/bin/bash

reference_fasta=$1
chrom_sizes=$2
blacklist_region=$3
overlap_peak=$4
inputlen=$5
genomewide_gc=$6
output_dir=$7
fold=$8
codedir=$9

function timestamp {
    # Function to get the current time with the new line character
    # removed 
    
    # current time
    date +"%Y-%m-%d_%H-%M-%S" | tr -d '\n'
}

logfile=$output_dir/"make_background_regions.log"
touch $logfile

# We also slop both our blacklist and peak set with the inputlen//2 and merge them. We will use this combined set as exlusion bed.
# we slop to make sure that the regions we choose have no intersection (not even 1bp) with the overlap peaks/blacklist regions
flank_size=$(( inputlen/2 ))
echo $( timestamp ): "bedtools slop -i $blacklist_region -g $chrom_sizes -b $flank_size > $output_dir/blacklist_slop1057.bed" | tee -a $logfile
bedtools slop -i $blacklist_region -g $chrom_sizes -b $flank_size | cut -d$'\t' -f 1-3 > $output_dir/blacklist_slop1057.bed 

echo $( timestamp ): "bedtools slop -i $overlap_peak -g $chrom_sizes -b $flank_size | cut -d$'\t' -f 1-3 > $output_dir/peaks_slop1057.bed" | tee -a $logfile
bedtools slop -i $overlap_peak -g $chrom_sizes -b $flank_size | cut -d$'\t' -f 1-3  > $output_dir/peaks_slop1057.bed 

echo $( timestamp ): "cat $output_dir/blacklist_slop1057.bed $output_dir/peaks_slop1057.bed | bedtools sort | bedtools merge -i stdin > $output_dir/exclude.bed" | tee -a $logfile
cat $output_dir/blacklist_slop1057.bed $output_dir/peaks_slop1057.bed | bedtools sort | bedtools merge -i stdin > $output_dir/exclude.bed

echo $( timestamp ): "rm $output_dir/blacklist_slop1057.bed" | tee -a $logfile
rm $output_dir/blacklist_slop1057.bed

echo $( timestamp ): "rm $output_dir/peaks_slop1057.bed" | tee -a $logfile
rm $output_dir/peaks_slop1057.bed
exclude_bed=$output_dir/exclude.bed

# create regions of size inputlen that do not fall in exclude bed and that gc-match with the given overlap peaks
# create as many regions as there are in overlap peaks bed file
bash $codedir/src/helpers/make_gc_matched_negatives/run.sh $overlap_peak $exclude_bed $inputlen $output_dir $reference_fasta $genomewide_gc $fold $chrom_sizes $logfile
# make a dummy summit file - rest of the pipeline uses this format
echo $( timestamp ): "awk -v OFS=\"\t\" '{print \$1, \$2, \$3, \".\",  \".\", \".\", \".\", \".\", \".\", \"1057\"}\' $output_dir/negatives.bed > $output_dir/negatives_with_summit.bed" | tee -a $logfile
awk -v OFS="\t" '{print $1, $2, $3, ".",  ".", ".", ".", ".", ".", "1057"}' $output_dir/negatives.bed > $output_dir/negatives_with_summit.bed

