#!/bin/bash

# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command failed with exit code $?."' EXIT

in_bam=${1?param missing - in_bam}
bigwig_prefix=${2?param missing - bigwig_prefix}
data_type=${3?param missing - data_type}
reference_fasta=${4?param missing - reference_fasta}
chrom_sizes=${5?param missing - chrom_sizes}

function timestamp {
    # Function to get the current time with the new line character
    # removed 
    
    # current time
    date +"%Y-%m-%d_%H-%M-%S" | tr -d '\n'
}

logfile=$bigwig_prefix"_preprocessing.log"
touch $logfile

bam_to_bigwig.sh $in_bam $bigwig_prefix $data_type $chrom_sizes $logfile

echo $( timestamp ): "build_pwm_from_bigwig -i $bigwig_prefix_unstranded.bw -g $reference_fasta -o $bigwig_prefix_bias_pwm.png -c chr20 -cz $chrom_sizes" | tee -a $logfile
build_pwm_from_bigwig -i $bigwig_prefix"_unstranded.bw" -g $reference_fasta -o $bigwig_prefix"_bias_pwm.png" -c "chr20" -cz $chrom_sizes 
