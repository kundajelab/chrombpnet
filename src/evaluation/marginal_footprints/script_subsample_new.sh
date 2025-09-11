gpu=2

regions=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/GM12878/negatives_data/negatives_with_summit.bed

reads=2m
cell_type=GM12878_$reads
model_dir="/oak/stanford/groups/akundaje/projects/chrombpnet/model_inputs/ENCODE_ATAC_downloads/GM12878/GM12878_subsampled/GM12878_bams/subsampled_$reads/fold_0/GM12878_$reads/"
output=$model_dir"footprints/"

mkdir $output

CUDA_VISIBLE_DEVICES=$gpu python marginal_footprinting_with_uncorrected.py -g /mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa \
	-r $regions \
	-fl /mnt/lab_data2/anusri/chrombpnet/splits/fold_0.json \
	-cm $model_dir"/chrombpnet_wo_bias.h5" \
	-um $model_dir"/chrombpnet.h5" \
	-bs 128 \
	-o $output"$cell_type" \
	-pwm_f gm12878_all_motifs.tsv \
	-tt $cell_type"_ATAC"


reads=10m
cell_type=GM12878_$reads
model_dir="/oak/stanford/groups/akundaje/projects/chrombpnet/model_inputs/ENCODE_ATAC_downloads/GM12878/GM12878_subsampled/GM12878_bams/subsampled_$reads/fold_0/GM12878_$reads/"
output=$model_dir"footprints/"

mkdir $output

CUDA_VISIBLE_DEVICES=$gpu python marginal_footprinting_with_uncorrected.py -g /mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa \
	-r $regions \
	-fl /mnt/lab_data2/anusri/chrombpnet/splits/fold_0.json \
	-cm $model_dir"/chrombpnet_wo_bias.h5" \
	-um $model_dir"/chrombpnet.h5" \
	-bs 128 \
	-o $output"$cell_type" \
	-pwm_f gm12878_all_motifs.tsv \
	-tt $cell_type"_ATAC"


reads=15m
cell_type=GM12878_$reads
model_dir="/oak/stanford/groups/akundaje/projects/chrombpnet/model_inputs/ENCODE_ATAC_downloads/GM12878/GM12878_subsampled/GM12878_bams/subsampled_$reads/fold_0/GM12878_$reads/"
output=$model_dir"footprints/"

mkdir $output

CUDA_VISIBLE_DEVICES=$gpu python marginal_footprinting_with_uncorrected.py -g /mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa \
	-r $regions \
	-fl /mnt/lab_data2/anusri/chrombpnet/splits/fold_0.json \
	-cm $model_dir"/chrombpnet_wo_bias.h5" \
	-um $model_dir"/chrombpnet.h5" \
	-bs 128 \
	-o $output"$cell_type" \
	-pwm_f gm12878_all_motifs.tsv \
	-tt $cell_type"_ATAC"
