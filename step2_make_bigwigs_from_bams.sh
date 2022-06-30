#!/bin/bash

echo "WARNING: If upgrading from v1.0 or v1.1 to v1.2. Note that chrombpnet has undergone linting to generate a modular structure for release on pypi.Hard-coded script paths are no longer necessary. Please refer to the updated README (below) to ensure your script calls are compatible with v1.2"

# exit when any command fails
set -e
set -o pipefail

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

echo $( timestamp ): "chrombpnet_pwm_from_bigwig -i $bigwig_prefix_unstranded.bw -g $reference_fasta -o $bigwig_prefix_bias_pwm -c chr20 -cz $chrom_sizes" | tee -a $logfile
chrombpnet_pwm_from_bigwig -i $bigwig_prefix"_unstranded.bw" -g $reference_fasta -o $bigwig_prefix"_bias_pwm" -c "chr20" -cz $chrom_sizes 
