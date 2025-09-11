CUDA_ViSIBLE_DEVICES=0 python $PWD/src/evaluation/interpret/interpret.py \
	 --genome=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/reference/hg38.genome.fa \
	 --regions=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/GM12878_100M/GM12878_100M_07.19.2022_bias_transfer_1234_fold_0_data_type_ATAC_PE/chrombpnet_model/interpret/full_GM12878_100M.interpreted_regions_counts.bed \
	--output_prefix=/oak/stanford/groups/akundaje/projects/chrombpnet/model_inputs/ENCODE_ATAC_downloads/GM12878/GM12878_subsampled/GM12878_bams/subsampled_15m/fold_0/GM12878_15m//interpret/full_GM12878_15m \
	--model_h5=/oak/stanford/groups/akundaje/projects/chrombpnet/model_inputs/ENCODE_ATAC_downloads/GM12878/GM12878_subsampled/GM12878_bams/subsampled_15m/fold_0/GM12878_15m//chrombpnet_wo_bias.h5 \
	--profile_or_counts=counts

