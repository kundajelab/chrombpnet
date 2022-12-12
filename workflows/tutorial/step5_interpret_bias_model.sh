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

modisco motifs -i $output_dir"/bias.profile_scores.h5" -n 50000 -o $output_dir/"modisco_results_allChroms_profile.hdf5" -w 500
modisco motifs -i $output_dir"/bias.counts_scores.h5" -n 50000 -o $output_dir/"modisco_results_allChroms_counts.hdf5" -w 500

meme_dir=$(chrombpnet_srcdir)"/../data"

modisco report -i $output_dir/"modisco_results_allChroms_profile.hdf5" -o  $output_dir/modisco_reports_profile/ -s "./" -m $meme_dir/motifs.meme.txt
modisco report -i $output_dir/"modisco_results_allChroms_counts.hdf5" -o  $output_dir/modisco_reports_counts/ -s "./" -m $meme_dir/motifs.meme.txt

chrombpnet_convert_html_to_pdf -i $output_dir/modisco_reports_profile/motifs.html -o $output_dir/modisco_reports_profile/motifs.pdf
chrombpnet_convert_html_to_pdf -i $output_dir/modisco_reports_counts/motifs.html -o $output_dir/modisco_reports_counts/motifs.pdf


