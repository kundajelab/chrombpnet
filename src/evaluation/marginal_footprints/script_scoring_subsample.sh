cell_type=GM12878_250M
model="GM12878_250M/GM12878_250M_07.19.2022_bias_transfer_1234_fold_0_data_type_ATAC_PE"


output="/oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/modisco_jun_30/modisco/ATAC/"$cell_type"/footprints/"
regions=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/GM12878/negatives_data/negatives_with_summit.bed
model_dir="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/"$model

mkdir "/oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/modisco_jun_30/modisco/ATAC/"$cell_type"/
mkdir $output

python footprint_scoring.py \
	-o $output$cell_type \
	-pwm_f gm12878_all_motifs.tsv \
	-tt $cell_type"_ATAC"




