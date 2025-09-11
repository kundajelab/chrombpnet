python src/evaluation/make_bigwigs/make_only_bigwigs.py \
	-cm /oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/GM12878/GM12878_03.01.2022_bias_128_4_1234_0.4_fold_0/chrombpnet_model/preds_atac/GM12878_w_bias_predictions.h5 \
	-cmb /oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/GM12878/GM12878_03.01.2022_bias_128_4_1234_0.4_fold_0/chrombpnet_model/preds_atac/GM12878_wo_bias_predictions.h5 \
	-r /oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/GM12878/peaks_no_blacklist.bed \
	-g reference/hg38.genome.fa \
	-c reference/chrom.sizes \
	-o GM12878/GM12878 \
