input_dir=/oak/stanford/groups/akundaje/projects/chrombpnet_paper/importance_scores/SIGNAL/
output_dir=/oak/stanford/groups/akundaje/projects/chrombpnet_paper/modisco/SIGNAL/

bash sbatch_modisco.sh GM12878 ATAC_07.22.2021 $input_dir $output_dir
bash sbatch_modisco.sh HEPG2 ATAC_07.22.2021 $input_dir $output_dir
bash sbatch_modisco.sh IMR90 ATAC_07.22.2021 $input_dir $output_dir
bash sbatch_modisco.sh H1 ATAC_07.22.2021 $input_dir $output_dir

