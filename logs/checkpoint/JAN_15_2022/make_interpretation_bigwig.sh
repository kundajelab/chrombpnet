


CUDA_VISIBLE_DEVICES=1 python $PWD/src/evaluation/interpret/interpret.py \
       --genome=/mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
       --regions=results/chrombpnet/ATAC_PE/GM12878/ATAC_PE_12.30.2021/chrombpnet_model/testing/ctcf.locus \
       --output_prefix=results/chrombpnet/ATAC_PE/GM12878/ATAC_PE_12.30.2021/chrombpnet_model/testing/chrombpnet_wo_bias \
       --model_h5=results/chrombpnet/ATAC_PE/GM12878/ATAC_PE_12.30.2021/chrombpnet_model/chrombpnet_wo_bias.h5 \


CUDA_VISIBLE_DEVICES=1 python $PWD/src/evaluation/interpret/interpret.py \
       --genome=/mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
       --regions=results/chrombpnet/ATAC_PE/GM12878/ATAC_PE_12.30.2021/chrombpnet_model/testing/ctcf.locus \
       --output_prefix=results/chrombpnet/ATAC_PE/GM12878/ATAC_PE_12.30.2021/chrombpnet_model/testing/chrombpnet \
       --model_h5=results/chrombpnet/ATAC_PE/GM12878/ATAC_PE_12.30.2021/chrombpnet_model/chrombpnet.h5 \


python src/evaluation/make_bigwigs/importance_hdf5_to_bigwig.py \
	-c /mnt/data/annotations/by_release/hg38/hg38.chrom.sizes \
        -r results/chrombpnet/ATAC_PE/GM12878/ATAC_PE_12.30.2021/chrombpnet_model/testing/ctcf.locus \
	-h5 results/chrombpnet/ATAC_PE/GM12878/ATAC_PE_12.30.2021/chrombpnet_model/testing/chrombpnet.profile_scores.h5 \
	-o results/chrombpnet/ATAC_PE/GM12878/ATAC_PE_12.30.2021/chrombpnet_model/testing/chrombpnet_profile.bw \
	-s results/chrombpnet/ATAC_PE/GM12878/ATAC_PE_12.30.2021/chrombpnet_model/testing/chrombpnet_profile.stats \
	-t 1

python src/evaluation/make_bigwigs/importance_hdf5_to_bigwig.py \
	-c /mnt/data/annotations/by_release/hg38/hg38.chrom.sizes \
        -r results/chrombpnet/ATAC_PE/GM12878/ATAC_PE_12.30.2021/chrombpnet_model/testing/ctcf.locus \
	-h5 results/chrombpnet/ATAC_PE/GM12878/ATAC_PE_12.30.2021/chrombpnet_model/testing/chrombpnet.counts_scores.h5 \
	-o results/chrombpnet/ATAC_PE/GM12878/ATAC_PE_12.30.2021/chrombpnet_model/testing/chrombpnet_counts.bw \
	-s results/chrombpnet/ATAC_PE/GM12878/ATAC_PE_12.30.2021/chrombpnet_model/testing/chrombpnet_counts.stats \
	-t 1


python src/evaluation/make_bigwigs/importance_hdf5_to_bigwig.py \
	-c /mnt/data/annotations/by_release/hg38/hg38.chrom.sizes \
        -r results/chrombpnet/ATAC_PE/GM12878/ATAC_PE_12.30.2021/chrombpnet_model/testing/ctcf.locus \
	-h5 results/chrombpnet/ATAC_PE/GM12878/ATAC_PE_12.30.2021/chrombpnet_model/testing/chrombpnet_wo_bias.profile_scores.h5 \
	-o results/chrombpnet/ATAC_PE/GM12878/ATAC_PE_12.30.2021/chrombpnet_model/testing/chrombpnet_wo_bias_profile.bw \
	-s results/chrombpnet/ATAC_PE/GM12878/ATAC_PE_12.30.2021/chrombpnet_model/testing/chrombpnet_wo_bias_profile.stats \
	-t 1

python src/evaluation/make_bigwigs/importance_hdf5_to_bigwig.py \
	-c /mnt/data/annotations/by_release/hg38/hg38.chrom.sizes \
        -r results/chrombpnet/ATAC_PE/GM12878/ATAC_PE_12.30.2021/chrombpnet_model/testing/ctcf.locus \
	-h5 results/chrombpnet/ATAC_PE/GM12878/ATAC_PE_12.30.2021/chrombpnet_model/testing/chrombpnet_wo_bias.counts_scores.h5 \
	-o results/chrombpnet/ATAC_PE/GM12878/ATAC_PE_12.30.2021/chrombpnet_model/testing/chrombpnet_wo_bias_counts.bw \
	-s results/chrombpnet/ATAC_PE/GM12878/ATAC_PE_12.30.2021/chrombpnet_model/testing/chrombpnet_wo_bias_counts.stats \
	-t 1



