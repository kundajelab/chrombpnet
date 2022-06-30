#!/bin/bash

# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT


scores_prefix=${1?param missing - scores_prefix}
output_dir=${2?param missing - output_dir}
score_type=${3?param missing - score_type}
seqlets=${4?param missing - seqlets}
crop=${5?param missing - crop} 
meme_db=${6?param missing - meme_db}
meme_logos=${7?param missing - meme_logos}
vier_logos=${8?param_missing - vier_logos}
vier_html=${9?param_missing - vier_html}
html_link=${10?param_missing - html_link}

chrombpnet_modisco -s $scores_prefix -p $score_type -o $output_dir -m $seqlets -c $crop
chrombpnet_tomtom_hits -m $output_dir/modisco_results_allChroms_counts.hdf5 -o $output_dir/$score_type.tomtom.tsv -d $meme_db -n 10 -th 0.3
chrombpnet_visualize_motif_matches -m $output_dir/modisco_results_allChroms_counts.hdf5 -t $output_dir/$score_type.tomtom.tsv -o $output_dir \
     -vd $vier_logos -th 0.3 -hl $html_link -vhl $vier_html \
      -s $score_type -d $meme_logos

