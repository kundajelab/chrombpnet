


CUDA_VISIBLE_DEVICES=1 python $PWD/src/evaluation/interpret/interpret.py \
       --genome=/mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
       --regions=results/chrombpnet/ATAC_PE/GM12878/ATAC_PE_12.30.2021/chrombpnet_model/testing/ctcf.locus \
       --output_prefix=results/chrombpnet/DNASE_SE/GM12878/DNASE_SE_12.30.2021/chrombpnet_model/testing/bias_model_scaled.h5 \
       --model_h5=results/chrombpnet/DNASE_SE/GM12878/DNASE_SE_12.30.2021/chrombpnet_model/bias_model_scaled.h5 \


python src/evaluation/make_bigwigs/importance_hdf5_to_bigwig.py \
	-c /mnt/data/annotations/by_release/hg38/hg38.chrom.sizes \
        -r results/chrombpnet/ATAC_PE/GM12878/ATAC_PE_12.30.2021/chrombpnet_model/testing/ctcf.locus \
	-h5 results/chrombpnet/DNASE_SE/GM12878/DNASE_SE_12.30.2021/chrombpnet_model/testing/bias_model_scaled.h5.profile_scores.h5 \
	-o results/chrombpnet/DNASE_SE/GM12878/DNASE_SE_12.30.2021/chrombpnet_model/testing/bias_model_scaled_profile.bw \
	-s results/chrombpnet/DNASE_SE/GM12878/DNASE_SE_12.30.2021/chrombpnet_model/testing/bias_model_scaled_profile.stats \
	-t 1

python src/evaluation/make_bigwigs/importance_hdf5_to_bigwig.py \
	-c /mnt/data/annotations/by_release/hg38/hg38.chrom.sizes \
        -r results/chrombpnet/ATAC_PE/GM12878/ATAC_PE_12.30.2021/chrombpnet_model/testing/ctcf.locus \
	-h5 results/chrombpnet/DNASE_SE/GM12878/DNASE_SE_12.30.2021/chrombpnet_model/testing/bias_model_scaled.h5.counts_scores.h5 \
	-o results/chrombpnet/DNASE_SE/GM12878/DNASE_SE_12.30.2021/chrombpnet_model/testing/bias_model_scaled_counts.bw \
	-s results/chrombpnet/DNASE_SE/GM12878/DNASE_SE_12.30.2021/chrombpnet_model/testing/bias_model_scaled_counts.stats \
	-t 1



