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
regions=${2?param missing - regions}
model_h5=${3?param missing - model_h5}
output_dir=${4?param missing - output_dir}
if [[ ! -e $output_dir ]]; then
    mkdir $output_dir
fi

## deepshap run

chrombpnet_deepshap \
    --genome=$reference_fasta \
    --regions=$regions \
    --output_prefix=$output_dir/corrected \
    --model_h5=$model_h5 \

## modisco run

chrombpnet_modisco -s $output_dir/corrected -p "profile" -o $output_dir -m 50000 -c 500
chrombpnet_modisco -s $output_dir/corrected -p "counts" -o $output_dir -m 50000 -c 500
