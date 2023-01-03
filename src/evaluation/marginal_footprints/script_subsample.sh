gpu=1

regions=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/GM12878/negatives_data/negatives_with_summit.bed
cell_type=GM12878_250M
model="GM12878_250M/GM12878_250M_07.19.2022_bias_transfer_1234_fold_0_data_type_ATAC_PE"
output="/oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/modisco_jun_30/modisco/ATAC/"$cell_type"/footprints/"
model_dir="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/"$model"/chrombpnet_model"

mkdir "/oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/modisco_jun_30/modisco/ATAC/"$cell_type
mkdir $output

CUDA_VISIBLE_DEVICES=$gpu python marginal_footprinting_with_uncorrected.py -g /mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa \
	-r $regions \
	-fl /mnt/lab_data2/anusri/chrombpnet/splits/fold_0.json \
	-cm $model_dir"/chrombpnet_wo_bias.h5" \
	-um $model_dir"/chrombpnet.h5" \
	-bs 128 \
	-o $output$cell_type \
	-pwm_f gm12878_all_motifs.tsv \
	-tt $cell_type"_ATAC"


cell_type=GM12878_100M
model="GM12878_100M/GM12878_100M_07.19.2022_bias_transfer_1234_fold_0_data_type_ATAC_PE"
output="/oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/modisco_jun_30/modisco/ATAC/"$cell_type"/footprints/"
model_dir="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/"$model"/chrombpnet_model"

mkdir "/oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/modisco_jun_30/modisco/ATAC/"$cell_type

mkdir $output



CUDA_VISIBLE_DEVICES=$gpu python marginal_footprinting_with_uncorrected.py -g /mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa \
	-r $regions \
	-fl /mnt/lab_data2/anusri/chrombpnet/splits/fold_0.json \
	-cm $model_dir"/chrombpnet_wo_bias.h5" \
	-um $model_dir"/chrombpnet.h5" \
	-bs 128 \
	-o $output$cell_type \
	-pwm_f gm12878_all_motifs.tsv \
	-tt $cell_type"_ATAC"

cell_type=GM12878_50M
model="GM12878_50M/GM12878_50M_07.18.2022_bias_transfer_1234_fold_0_data_type_ATAC_PE"
output="/oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/modisco_jun_30/modisco/ATAC/"$cell_type"/footprints/"
model_dir="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/"$model"/chrombpnet_model"

mkdir "/oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/modisco_jun_30/modisco/ATAC/"$cell_type
mkdir $output

CUDA_VISIBLE_DEVICES=$gpu python marginal_footprinting_with_uncorrected.py -g /mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa \
	-r $regions \
	-fl /mnt/lab_data2/anusri/chrombpnet/splits/fold_0.json \
	-cm $model_dir"/chrombpnet_wo_bias.h5" \
	-um $model_dir"/chrombpnet.h5" \
	-bs 128 \
	-o $output$cell_type \
	-pwm_f gm12878_all_motifs.tsv \
	-tt $cell_type"_ATAC"

cell_type=GM12878_25M
model="GM12878_25M/GM12878_25M_07.18.2022_bias_transfer_1234_fold_0_data_type_ATAC_PE"
output="/oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/modisco_jun_30/modisco/ATAC/"$cell_type"/footprints/"
model_dir="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/"$model"/chrombpnet_model"

mkdir "/oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/modisco_jun_30/modisco/ATAC/"$cell_type
mkdir $output


CUDA_VISIBLE_DEVICES=$gpu python marginal_footprinting_with_uncorrected.py -g /mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa \
	-r $regions \
	-fl /mnt/lab_data2/anusri/chrombpnet/splits/fold_0.json \
	-cm $model_dir"/chrombpnet_wo_bias.h5" \
	-um $model_dir"/chrombpnet.h5" \
	-bs 128 \
	-o $output$cell_type \
	-pwm_f gm12878_all_motifs.tsv \
	-tt $cell_type"_ATAC"

cell_type=GM12878_5M
model="GM12878_5M/GM12878_5M_07.18.2022_bias_transfer_1234_fold_0_data_type_ATAC_PE"
output="/oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/modisco_jun_30/modisco/ATAC/"$cell_type"/footprints/"
model_dir="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/"$model"/chrombpnet_model"

mkdir "/oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/modisco_jun_30/modisco/ATAC/"$cell_type
mkdir $output


CUDA_VISIBLE_DEVICES=$gpu python marginal_footprinting_with_uncorrected.py -g /mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa \
	-r $regions \
	-fl /mnt/lab_data2/anusri/chrombpnet/splits/fold_0.json \
	-cm $model_dir"/chrombpnet_wo_bias.h5" \
	-um $model_dir"/chrombpnet.h5" \
	-bs 128 \
	-o $output$cell_type \
	-pwm_f gm12878_all_motifs.tsv \
	-tt $cell_type"_ATAC"

cell_type="GM12878"
model="/nautilus_runs/GM12878_03.01.2022_bias_128_4_1234_0.4_fold_0"
output="/oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/modisco_jun_30/modisco/ATAC/"$cell_type"/footprints/"
model_dir="/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/"$cell_type$model"/chrombpnet_model"
mkdir $output

CUDA_VISIBLE_DEVICES=$gpu python marginal_footprinting_with_uncorrected.py -g /mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa \
	-r $regions \
	-fl /mnt/lab_data2/anusri/chrombpnet/splits/fold_0.json \
	-cm $model_dir"/chrombpnet_wo_bias.h5" \
	-um $model_dir"/chrombpnet.h5" \
	-bs 128 \
	-o $output$cell_type \
	-pwm_f union_motifs.tsv \
	-tt $cell_type"_ATAC"
