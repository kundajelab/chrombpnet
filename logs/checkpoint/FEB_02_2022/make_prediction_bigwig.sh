cell_type=K562
regions=results/chrombpnet/ATAC_PE/$cell_type/data/chr20.peaks.bed
output_dir=results/chrombpnet/ATAC_PE/$cell_type/ATAC_PE_01.24.2022_bigger_jitter_on_bias/chrombpnet_model/

mkdir $output_dir/interpret_chr20/

CUDA_VISIBLE_DEVICES=1 python src/evaluation/make_bigwigs/predict_to_bigwig.py \
	-bm $output_dir/bias_model_scaled.h5 \
	-cm $output_dir/chrombpnet.h5 \
	-cmb $output_dir/chrombpnet_wo_bias.h5 \
	-r $regions \
	-g /mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
	-o $output_dir/interpret_chr20/chr20 \
	-c /mnt/data/annotations/by_release/hg38/hg38.chrom.sizes \
	-t 1 \


