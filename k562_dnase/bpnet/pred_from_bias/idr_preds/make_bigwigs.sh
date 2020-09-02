kerasAC_bigwigs_from_io --predictions_hdf /srv/scratch/annashch/deeplearning/profile/k562_dnase/bpnet/pred_from_bias/idr_preds/tobias.bias.preds.k562.0.predictions \
			--out_prefix k562_dnase_bpnet \
			--chrom_sizes_file hg38.chrom.sizes
kerasAC_bigwigs_from_io --predictions_hdf /srv/scratch/annashch/deeplearning/profile/k562_dnase/bpnet/pred_from_bias/idr_preds/6mer.bias.preds.k562.0.predictions \
			--out_prefix k562_dnase_tobias \
			--chrom_sizes_file hg38.chrom.sizes
kerasAC_bigwigs_from_io --predictions_hdf /srv/scratch/annashch/deeplearning/profile/k562_dnase/bpnet/pred_from_bias/idr_preds/bpnet.bias.preds.k562.0.predictions \
			--out_prefix k562_dnase_6mer \
			--chrom_sizes_file hg38.chrom.sizes

