CUDA_VISIBLE_DEVICES=1 python src/evaluation/make_bigwigs/predict_to_bigwig.py \
	-bm results/chrombpnet/ATAC_PE/GM12878/ATAC_PE_12.30.2021/chrombpnet_model/bias_model_scaled.h5 \
	-cm results/chrombpnet/ATAC_PE/GM12878/ATAC_PE_12.30.2021/chrombpnet_model/chrombpnet.h5 \
	-cmb results/chrombpnet/ATAC_PE/GM12878/ATAC_PE_12.30.2021/chrombpnet_model/chrombpnet_wo_bias.h5 \
	-r results/chrombpnet/ATAC_PE/GM12878/ATAC_PE_12.30.2021/chrombpnet_model/testing/ctcf.locus \
	-g /mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
	-o results/chrombpnet/ATAC_PE/GM12878/ATAC_PE_12.30.2021/chrombpnet_model/testing/ctcf \
	-c /mnt/data/annotations/by_release/hg38/hg38.chrom.sizes \
	-t 1 \


