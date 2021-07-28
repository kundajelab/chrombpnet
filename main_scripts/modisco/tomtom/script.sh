flank=500
input_dir=/oak/stanford/groups/akundaje/projects/chrombpnet_paper/importance_scores/SIGNAL/
output_dir=/oak/stanford/groups/akundaje/projects/chrombpnet_paper/modisco/flank_len_$flank/SIGNAL/
tomtom_dir=/oak/stanford/groups/akundaje/projects/chrombpnet_paper/tomtom/flank_len_$flank/SIGNAL/

bash run_tomtom.sh GM12878 ATAC_07.22.2021 $input_dir $output_dir $tomtom_dir
bash run_tomtom.sh HEPG2 ATAC_07.22.2021 $input_dir $output_dir $tomtom_dir
bash run_tomtom.sh IMR90 ATAC_07.22.2021 $input_dir $output_dir $tomtom_dir
bash run_tomtom.sh H1 ATAC_07.22.2021 $input_dir $output_dir $tomtom_dir


input_dir=/oak/stanford/groups/akundaje/projects/chrombpnet_paper/importance_scores/BIAS/
output_dir=/oak/stanford/groups/akundaje/projects/chrombpnet_paper/modisco/flank_len_$flank/BIAS/
tomtom_dir=/oak/stanford/groups/akundaje/projects/chrombpnet_paper/tomtom/flank_len_$flank/BIAS/

bash run_tomtom.sh GM12878 ATAC_07.22.2021 $input_dir $output_dir $tomtom_dir
bash run_tomtom.sh HEPG2 ATAC_07.22.2021 $input_dir $output_dir $tomtom_dir
bash run_tomtom.sh IMR90 ATAC_07.22.2021 $input_dir $output_dir $tomtom_dir
bash run_tomtom.sh H1 ATAC_07.22.2021 $input_dir $output_dir $tomtom_dir

