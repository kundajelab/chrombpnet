#!/bin/bash

echo "WARNING: If upgrading from v1.0 or v1.1 to v1.2. Note that chrombpnet has undergone linting to generate a modular structure for release on pypi.Hard-coded script paths are no longer necessary. Please refer to the updated README (below) to ensure your script calls are compatible with v1.2"

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
regions=${2?param missing - regions}
model_h5=${3?param missing - model_h5}
output_dir=${4?param missing - output_dir}

## deepshap run

chrombpnet_deepshap \
    --genome=$reference_fasta \
    --regions=$regions \
    --output_prefix=$output_dir/bias \
    --model_h5=$model_h5 \

## modisco run

modisco motifs -i $output_dir"/bias.profile_scores.h5" -n 50000 -o "modisco_results_allChroms_profile.hdf5" -w 500
modisco motifs -i $output_dir"/bias.counts_scores.h5" -n 50000 -o "modisco_results_allChroms_counts.hdf5" -w 500

modisco report -i "modisco_results_allChroms_profile.hdf5" -o  $output_dir/modisco_reports_profile/ -m /home/anusri/chrombpnet/data/motifs.meme.txt
modisco report -i "modisco_results_allChroms_counts.hdf5" -o  $output_dir/modisco_reports_counts/ -m /home/anusri/chrombpnet/data/motifs.meme.txt
