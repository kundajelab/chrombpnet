


scores_prefix=$1
output_dir=$2
score_type=$3
seqlets=$4
crop=$5
#meme_db=$6
#meme_logos=$7
#vier_logos=$8
#vier_html=$9
#html_link=${10}

scores_prefix=/oak/stanford/groups/akundaje/projects/chrombpnet_paper/testing/testing
output_dir=/oak/stanford/groups/akundaje/projects/chrombpnet_paper/testing/
score_type=counts
seqlets=50000
crop=1000
meme_db=/oak/stanford/groups/akundaje/soumyak/motifs/motifs.meme.txt
meme_logos=/oak/stanford/groups/akundaje/soumyak/motifs/pfms/
vier_logos=/oak/stanford/groups/akundaje/projects/chromatin-atlas/vierstra_logos/

vier_html=http://mitra.stanford.edu/kundaje/anusri/chrombpnet_paper/vierstra_logos/
html_link=http://mitra.stanford.edu/kundaje/anusri/chrombpnet_paper/testing/

python run_modisco.py -s $scores_prefix -p $score_type -o $output_dir -m $seqlets -c $crop
#python fetch_tomtom.py -m $output_dir/modisco_results_allChroms_counts.hdf5 -o $output_dir/$score_type.tomtom.tsv -d $meme_db -n 10 -th 0.3
#python visualize_motif_matches.py -m $output_dir/modisco_results_allChroms_counts.hdf5 -t $output_dir/$score_type.tomtom.tsv -o $output_dir \
#     -vd $vier_logos -th 0.3 -hl  $html_link -vhl $vier_html \
#      -s $score_type -d $meme_logos

