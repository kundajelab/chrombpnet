# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT


scores_prefix=$1
output_dir=$2
score_type=$3
seqlets=$4
crop=$5
meme_db=$6
meme_logos=$7
vier_logos=$8
vier_html=$9
html_link=${10}

python run_modisco.py -s $scores_prefix -p $score_type -o $output_dir -m $seqlets -c $crop
python fetch_tomtom.py -m $output_dir/modisco_results_allChroms_counts.hdf5 -o $output_dir/$score_type.tomtom.tsv -d $meme_db -n 10 -th 0.3
python visualize_motif_matches.py -m $output_dir/modisco_results_allChroms_counts.hdf5 -t $output_dir/$score_type.tomtom.tsv -o $output_dir \
     -vd $vier_logos -th 0.3 -hl  $html_link -vhl $vier_html \
      -s $score_type -d $meme_logos

