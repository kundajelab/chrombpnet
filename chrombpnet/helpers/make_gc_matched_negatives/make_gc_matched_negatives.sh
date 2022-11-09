#!/bin/bash

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

foreground_bed=${1?param missing - foreground_bed}
exclude_bed=${2?param missing - exclude_bed}
inputlen=${3?param missing - inputlen}
output_dir=${4?param missing - output_dir}
genome=${5?param missing - genome}
genomewide_gc=${6?param missing - genomewide_gc}
fold=${7?param missing - fold}
chrom_sizes=${8?param missing - chrom_sizes}
logfile=${9}

function timestamp {
    # Function to get the current time with the new line character
    # removed 
    
    # current time
    date +"%Y-%m-%d_%H-%M-%S" | tr -d '\n'
}

# create the log file
if [ -z "$logfile" ]
  then
    echo "No logfile supplied - creating one"
    logfile=$output_prefix"_make_gc_matching_run.log"
    touch $logfile
fi

echo $( timestamp ): "get gc content of the positive sequences" | tee -a $logfile
echo $( timestamp ): "chrombpnet_gc_content_foreground \\
        --input_bed $foreground_bed \\
        --genome $genome \\
        --chrom_sizes $chrom_sizes \\
        --output_prefix $output_dir/foreground.gc \\
        --inputlen $inputlen" | tee -a $logfile

chrombpnet_gc_content_foreground \
       --input_bed $foreground_bed \
       --genome $genome \
       --chrom_sizes $chrom_sizes \
       --output_prefix $output_dir/foreground.gc \
       --inputlen $inputlen | tee -a $logfile

echo $( timestamp ): "get candidate negative bed" | tee -a $logfile
echo $( timestamp ): "bedtools intersect -v -a $genomewide_gc -b $exclude_bed  > $output_dir/candidate.negatives.bed" | tee -a $logfile
bedtools intersect -v -a $genomewide_gc -b $exclude_bed  > $output_dir/candidate.negatives.bed

echo $( timestamp ): "find regions in candidate negative bed that gc-match wth foreground" | tee -a $logfile
echo $( timestamp ): "chrombpnet_gc_matched_negatives \\
        --candidate_negatives $output_dir/candidate.negatives.bed \\
        --foreground_gc_bed $output_dir/foreground.gc.bed \\
        --output_prefix $output_dir/negatives \\
        --chr_fold_path $fold \\
        --neg_to_pos_ratio_train 2" | tee -a $logfile
chrombpnet_gc_matched_negatives \
    --candidate_negatives $output_dir/candidate.negatives.bed \
    --foreground_gc_bed $output_dir/foreground.gc.bed \
    --output_prefix $output_dir/negatives \
    --chr_fold_path $fold \
    --neg_to_pos_ratio_train 2 | tee -a $logfile


