chrom_sizes=$PWD/reference/chrom.sizes
ref_fasta=$PWD/reference/hg38.genome.fa

reads=10m
modelname=GM12878_$reads
modeldir=/oak/stanford/groups/akundaje/projects/chrombpnet/model_inputs/ENCODE_ATAC_downloads/GM12878/GM12878_subsampled/GM12878_bams/subsampled_$reads/fold_0/GM12878_$reads/

CUDA_VISIBLE_DEVICES=1 python src/evaluation/make_bigwigs/predict_to_bigwig.py \
	-bm $modeldir/bias_model_scaled.h5 \
	-cm $modeldir/chrombpnet.h5 \
	-cmb  $modeldir/chrombpnet_wo_bias.h5 \
	-r /mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/GM12878/nautilus_runs/GM12878_03.01.2022_bias_128_4_1234_0.4_fold_0/chrombpnet_model/filtered.peaks.bed \
	-g $ref_fasta \
	-c $chrom_sizes \
	-o $modeldir/interpret/$modelname


reads=2m
modelname=GM12878_$reads
modeldir=/oak/stanford/groups/akundaje/projects/chrombpnet/model_inputs/ENCODE_ATAC_downloads/GM12878/GM12878_subsampled/GM12878_bams/subsampled_$reads/fold_0/GM12878_$reads/

CUDA_VISIBLE_DEVICES=1 python src/evaluation/make_bigwigs/predict_to_bigwig.py \
	-bm $modeldir/bias_model_scaled.h5 \
	-cm $modeldir/chrombpnet.h5 \
	-cmb  $modeldir/chrombpnet_wo_bias.h5 \
	-r /mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/GM12878/nautilus_runs/GM12878_03.01.2022_bias_128_4_1234_0.4_fold_0/chrombpnet_model/filtered.peaks.bed \
	-g $ref_fasta \
	-c $chrom_sizes \
	-o $modeldir/interpret/$modelname
