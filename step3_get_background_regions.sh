#!/bin/bash

echo "WARNING: If upgrading from v1.0 or v1.1 to v1.2. Note that chrombpnet has undergone linting to generate a modular structure for release on pypi.Hard-coded script paths are no longer necessary. Please refer to the updated README (below) to ensure your script calls are compatible with v1.2"


# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG

cleanup() {
    exit_code=$?
    if [ ${exit_code} == 0 ]
    then
	echo "Completed execution"
    else
	echo "\"${last_command}\" failed with exit code ${exit_code}."
    fi
}

# echo an error message before exiting
trap 'cleanup' EXIT INT TERM


reference_fasta=${1?param missing - reference_fasta}
chrom_sizes=${2?param missing - chrom_sizes}
blacklist_region=${3?param missing - blaklist_region}
overlap_peak=${4?param missing - overlap_peak}
inputlen=${5?param missing - inputlen}
genomewide_gc=${6?param missing - genomewide_gc}
output_dir=${7?param missing - output_dir} 
fold=${8?param missing - fold}

if [[ ! -e $output_dir ]]; then
    mkdir $output_dir
fi


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
make_gc_matched_negatives.sh  $overlap_peak $exclude_bed $inputlen $output_dir $reference_fasta $genomewide_gc $fold $chrom_sizes $logfile
# make a dummy summit file - rest of the pipeline uses this format
echo $( timestamp ): "awk -v OFS=\"\t\" '{print \$1, \$2, \$3, \".\",  \".\", \".\", \".\", \".\", \".\", \"1057\"}\' $output_dir/negatives.bed > $output_dir/negatives_with_summit.bed" | tee -a $logfile
awk -v OFS="\t" '{print $1, $2, $3, ".",  ".", ".", ".", ".", ".", "1057"}' $output_dir/negatives.bed > $output_dir/negatives_with_summit.bed

