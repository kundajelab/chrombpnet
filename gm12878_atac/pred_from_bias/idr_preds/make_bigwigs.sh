kerasAC_bigwigs_from_io --predictions_hdf /srv/scratch/annashch/deeplearning/profile/gm12878_atac/pred_from_bias/idr_preds/tobias.bias.preds.gm12878.0.predictions \
			--out_prefix gm12878_atac_bpnet \
			--chrom_sizes_file hg38.chrom.sizes &
kerasAC_bigwigs_from_io --predictions_hdf /srv/scratch/annashch/deeplearning/profile/gm12878_atac/pred_from_bias/idr_preds/6mer.bias.preds.gm12878.0.predictions \
			--out_prefix gm12878_atac_tobias \
			--chrom_sizes_file hg38.chrom.sizes &
kerasAC_bigwigs_from_io --predictions_hdf /srv/scratch/annashch/deeplearning/profile/gm12878_atac/pred_from_bias/idr_preds/bpnet.bias.preds.gm12878.0.predictions \
			--out_prefix gm12878_atac_6mer \
			--chrom_sizes_file hg38.chrom.sizes &

